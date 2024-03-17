/**
 * @file   fracture_2d_example.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Fri Feb 5 2021
 *
 * @brief  TODO
 *
 *
 * Copyright (C) 2021 ETH Zurich (David S. Kammer)
 *
 * This file is part of uguca.
 *
 * uguca is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * uguca is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with uguca.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "static_communicator_mpi.hh"
#include "uca_parameter_reader.hh"
#include "material.hh"
#include "uca_simple_mesh.hh"

#include "bimat_interface.hh"
#include "barras_law.hh"

#include <iostream>
#include <cmath>

using namespace uguca;

int main(int argc, char *argv[]) {

  // communicator for parallel simulation
  StaticCommunicatorMPI * comm = StaticCommunicatorMPI::getInstance();
  int world_rank = comm->whoAmI();

  // get input file name
  std::string fname;
  if(argc<2) {
    if (world_rank==0) {
      std::cerr << "Not enough arguments:"
		<< " ./fracture_2d_example <input_file>" << std::endl;
    }
    return 1;
  }
  else {
    fname = argv[1];
  }

  // read input file
  ParameterReader data;
  data.readInputFile(fname);

  // mesh
  double length   = data.get<double>("length");
  int nb_elements = data.getOrUse<int>("nb_elements",512); // getOrUse let's you set default value
  SimpleMesh mesh(length, nb_elements);

  // constitutive interface law
  BarrasLaw law(mesh,
		data.get<double>("tauc", "cohesion"),
		data.get<double>("dc", "cohesion"));

  // materials
  Material top_mat = Material(data.get<double>("E","top"),
			      data.get<double>("nu","top"),
			      data.get<double>("rho","top"));
  top_mat.readPrecomputedKernels();
  Material bot_mat = Material(data.get<double>("E","bot"),
			      data.get<double>("nu","bot"),
			      data.get<double>("rho","bot"));
  bot_mat.readPrecomputedKernels();

  // interface
  BimatInterface interface(mesh, {_x,_y}, top_mat, bot_mat, law);

  // external loading
  interface.getLoad().setAllValuesTo(data.get<double>("shear_load"),0);
  interface.getLoad().setAllValuesTo(data.get<double>("normal_load"),1);

  // time step
  double duration = data.get<double>("duration");
  double time_step = data.get<double>("tsf") * interface.getStableTimeStep();
  int nb_time_steps = duration/time_step;
  interface.setTimeStep(time_step);

  // initialization
  interface.init();

  // heterogeneity for nucleation: decreased strength
  const TwoDVector & coords = mesh.getLocalCoords();
  NodalField & tau_max = law.getTauMax();
  double a0 = data.get<double>("a0");
  for (int i=0;i<mesh.getNbLocalNodes(); ++i)
    if (std::abs(coords(i,0) - length/2.) < a0/2.)
      tau_max(i) = 0.;

  // dumping
  interface.initDump(data.get<std::string>("sim_name"),".");
  interface.registerDumpFields(data.get<std::string>("dump_fields"));
  unsigned int dump_int = std::max(1, nb_time_steps/data.get<int>("nb_dumps"));

  // setup of restart to dump
  Restart restart(data.get<std::string>("sim_name"),".");
  interface.registerToRestart(restart);
  unsigned int restart_int = std::max(1, nb_time_steps/data.get<int>("nb_restarts"));
  
  // time stepping
  int s=1;

  // check if restart from file is required
  if (data.has("restart_step")) {
    s = data.get<int>("restart_step");
    // different from dumper not to overwrite it
    Restart restart_load = restart;
    restart_load.initIO(data.get<std::string>("restart_name"),".");
    restart_load.load(s);
  }
  else { // only dump when not restarted
    interface.dump(0,0);
  }
  
  for (; s<=nb_time_steps; ++s) {
    
    if (world_rank==0) {
      std::cout << "s=" << s << "/" << nb_time_steps << "\r";
      std::cout.flush();
    }

    // time integration
    interface.advanceTimeStep();

    // dump
    if (s % dump_int == 0)
      interface.dump(s,s*time_step);

    // dump restart
    if (s % restart_int == 0)
      restart.dump(s);
  }

  // dump restart
  restart.dump(nb_time_steps);
  
  comm->finalize();
  return 0;
}
