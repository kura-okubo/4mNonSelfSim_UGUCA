/**
 * @file   basic_example.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Fri Feb 5 2021
 *
 * @brief  basic example simulation using the uguca library
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
#include "material.hh"
#include "uca_simple_mesh.hh"
#include "unimat_shear_interface.hh"
#include "linear_shear_cohesive_law.hh"

#include <cmath>

using namespace uguca;

int main() {

  // materials - args: E, nu, rho
  Material top_mat = Material(3e9, 0.25, 2000);
  top_mat.readPrecomputedKernels();

  // mesh
  double length = 1.;
  int nb_elements = 512;
  SimpleMesh mesh(length, nb_elements);

  // constitutive interface law - args: mesh, Gc, tauc
  LinearShearCohesiveLaw law(mesh, 10., 1e6);

  // weak interface - args: mesh, material, interface law
  UnimatShearInterface interface(mesh, top_mat, law);

  // external loading
  interface.getLoad().component(1).setAllValuesTo(-5e6);

  // heterogeneity for nucleation
  double * X = mesh.getLocalCoords()[0];
  NodalFieldComponent & ext_shear = interface.getLoad().component(0);
  double Xnuc = length/2.;
  for (int i=0;i<mesh.getNbLocalNodes(); ++i)
    ext_shear(i) = 1.1e6 - 0.7e6*std::tanh(80*std::abs(X[i] - Xnuc) - 2.);

  // time step
  double duration = 5e-4;
  double time_step = 0.3 * interface.getStableTimeStep();
  int nb_time_steps = duration/time_step;
  interface.setTimeStep(time_step);

  // initialization
  interface.init();
  
  // dumping
  interface.initDump("basic_example","."); // args: name, path
  interface.registerDumpField("cohesion_0");
  interface.dump(0,0); // args: time-step, time
  unsigned int dump_int = std::max(1, nb_time_steps/300);

  // restart
  Restart restart("basic_example", "."); // args: name, path
  interface.registerToRestart(restart);
  
  // time stepping
  for (int s=1; s<=nb_time_steps; ++s) {
    interface.advanceTimeStep();
    if (s % dump_int == 0)
      interface.dump(s,s*time_step); // args: time-step, time
  }

  // write restart files
  restart.dump(nb_time_steps); // args: step
  
  StaticCommunicatorMPI::getInstance()->finalize();
  return 0;
}
