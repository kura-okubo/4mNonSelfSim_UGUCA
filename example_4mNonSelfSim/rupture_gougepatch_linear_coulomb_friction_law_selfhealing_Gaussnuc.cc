/*
Conduct the dynamic rupture modeling on the gouge-mediated seismic event.
This code refers the example of fracture_3d, TPV3 benchmark, and the project of breakdown_energy_scaling by:

Ke, C.-Y., McLaskey, G. C., and Kammer, D. S. Earthquake breakdown energy scaling despite constant fracture energy.
 Nature Communications, 13(1):1005, 2022, doi:10.1038/s41467-022-28647-4.

https://gitlab.com/uguca/projects/breakdown_energy_scaling

- new functions: self-healing friction law and the Gaussian stress distribution in the nucleation area.
*/

/**
 * @file   fracture_3d_example.cc
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
#include "unimat_shear_interface.hh"
// #include "linear_shear_cohesive_law.hh"
#include "linear_coulomb_friction_law_selfhealing.hh"

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
		<< " ./rupture_gougepatch <input_file>" << std::endl;
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
  double x_length   = data.get<double>("x_length");
  double z_length   = data.get<double>("z_length");
  int nb_x_elements = data.get<int>("nb_x_elements");
  int nb_z_elements = data.get<int>("nb_z_elements");
  SimpleMesh mesh(x_length, nb_x_elements,
	    z_length, nb_z_elements);

  // friction parameters
  // double Gc_nuc = data.get<double>("Gc_nuc");
  double dc_nuc = data.get<double>("dc_nuc");
  double ds_nuc = data.get<double>("ds_nuc");
  double fp_nuc = data.get<double>("fp_nuc");
  double fr_nuc = data.get<double>("fr_nuc");
  // double tau_c_nuc = data.get<double>("tau_c_nuc");
  // double tau_r_nuc = data.get<double>("tau_r_nuc");

  // double Gc_patch = data.get<double>("Gc_patch");
  double dc_patch = data.get<double>("dc_patch");
  double ds_patch = data.get<double>("ds_patch");
  double fp_patch = data.get<double>("fp_patch");
  double fr_patch = data.get<double>("fr_patch");
  // double tau_c_patch = data.get<double>("tau_c_patch");
  // double tau_r_patch = data.get<double>("tau_r_patch");

  // double Gc_background = data.get<double>("Gc_background");
  double dc_background = data.get<double>("dc_background");
  double ds_background = data.get<double>("ds_background");
  double fp_background = data.get<double>("fp_background");
  double fr_background = data.get<double>("fr_background");
  // double tau_c_background = data.get<double>("tau_c_background");
  // double tau_r_background = data.get<double>("tau_r_background");

  // initial stress state
  double sn_nuc = data.get<double>("sn_nuc");
  double sn_patch = data.get<double>("sn_patch");
  double sn_background = data.get<double>("sn_background");

  // double tau_nuc = data.get<double>("tau_nuc");
  double tau_patch = data.get<double>("tau_patch");
  double tau_background = data.get<double>("tau_background");

  double R_nuc = data.get<double>("R_nuc");
  double R_patch = data.get<double>("R_patch");
  double R_margin = data.get<double>("R_margin");

  double nuc_x = data.get<double>("nuc_x");
  double nuc_z = data.get<double>("nuc_z");

  double c_nucexcess = data.get<double>("c_nucexcess"); // percentage for the excess of initial stress: (1+c)taup

  // constitutive interface law
  // BarrasLaw law(mesh,
	// 	data.get<double>("tauc"),
	// 	data.get<double>("dc"));

  // Initialize the friction law with the values of patch
  // LinearShearCohesiveLaw law(mesh, Gc_patch, tau_c_patch, tau_r_patch);
  double char_reg_time = data.get<double>("char_reg_time");
  LinearCoulombFrictionLawSelfHealing law(mesh, fp_patch, fr_patch, dc_patch, ds_patch, char_reg_time);
	
  // materials
  Material top_mat = Material(data.get<double>("E_top"),
			      data.get<double>("nu_top"),
			      data.get<double>("rho_top"));
  top_mat.readPrecomputedKernels();

  // interface
  // DefRigInterface interface(mesh, top_mat, law);
  UnimatShearInterface interface(mesh, top_mat, law);

  // Initialize external loading with the values of patch
  interface.getLoad().component(0).setAllValuesTo(tau_patch);
  interface.getLoad().component(1).setAllValuesTo(sn_patch);

  // time step
  double duration = data.get<double>("duration");
  double time_step = data.get<double>("tsf") * interface.getStableTimeStep();
  int nb_time_steps = duration/time_step;
  interface.setTimeStep(time_step);

  // initialization
  interface.init();

  // Set the nucleation zone
  double * X = mesh.getLocalCoords()[0];
  double * Z = mesh.getLocalCoords()[2];
  NodalFieldComponent & load_0 = interface.getLoad().component(0);
  NodalFieldComponent & load_1 = interface.getLoad().component(1);
  NodalFieldComponent & cohesion_0 = interface.getCohesion().component(0);
  NodalFieldComponent & cohesion_1 = interface.getCohesion().component(1);
  // NodalFieldComponent & Gamma_c = law.getGc();
  // NodalFieldComponent & tauc = law.getTauc();
  // NodalFieldComponent & taur = law.getTaur();
  NodalFieldComponent & MuS = law.getMuS();
  NodalFieldComponent & MuK = law.getMuK();
  NodalFieldComponent & Dc = law.getDc();
  NodalFieldComponent & Ds = law.getDs();

  // for (int i=0;i<mesh.getNbLocalNodes(); ++i){
  //   double x0 = std::abs(X[i] - x_length/2. - nuc_x);
  //   double z0 = std::abs(Z[i] - z_length/2. - nuc_z);
  //   double r = std::sqrt(x0 * x0 + z0 * z0);
  //   if (r <= R_nuc) {
  //   load_0(i) = tau_nuc;
  //   load_1(i) = sn_nuc;
  //   cohesion_0(i) = tau_nuc; // initialize the cohesion same as the loading just for the sake of dumping at t=0
  //   cohesion_1(i) = sn_nuc;
  //   // Gamma_c(i) = Gc_nuc;
  //   // tauc(i) = tau_c_nuc;
  //   // taur(i) = tau_r_nuc;
  //   MuS(i) = fp_nuc;
  //   MuK(i) = fr_nuc;
  //   Dc(i) = dc_nuc;
  //   Ds(i) = ds_nuc;
  //   }
  // }

  // 2024.05.13 compute Gaussian distribution in the nucleation area
  double tau_c_nuc = -sn_nuc * fp_nuc; // normal stress is defined as negative
  double sig_nuc = R_nuc * std::sqrt(1/(2*std::log( (1+c_nucexcess)*tau_c_nuc/tau_patch )));
  double k_nuc = (2*M_PI*sig_nuc*sig_nuc)*(1+c_nucexcess)*tau_c_nuc;
  printf("debug M_PI=%f, sig_nuc = %f, k_nuc = %f\n", M_PI, sig_nuc, k_nuc);

  for (int i=0;i<mesh.getNbLocalNodes(); ++i){
    double x0 = std::abs(X[i] - x_length/2. - nuc_x);
    double z0 = std::abs(Z[i] - z_length/2. - nuc_z);
    double r = std::sqrt(x0 * x0 + z0 * z0);
    if (r <= R_nuc) {
    double tau_nuc_Gauss =  (k_nuc/(2*M_PI*sig_nuc*sig_nuc)) * std::exp(-(r*r)/(2*sig_nuc*sig_nuc));
    load_0(i) = tau_nuc_Gauss;
    load_1(i) = sn_nuc;
    cohesion_0(i) = tau_nuc_Gauss; // initialize the cohesion same as the loading just for the sake of dumping at t=0
    cohesion_1(i) = sn_nuc;
    // Gamma_c(i) = Gc_nuc;
    // tauc(i) = tau_c_nuc;
    // taur(i) = tau_r_nuc;
    MuS(i) = fp_nuc; // fp_nuc is constant even though tau_nuc is in Gaussian shape
    MuK(i) = fr_nuc;
    Dc(i) = dc_nuc;
    Ds(i) = ds_nuc;
    }
  }

  // Set the margin 
  for (int i=0;i<mesh.getNbLocalNodes(); ++i){
    double x0 = std::abs(X[i] - x_length/2.);
    double z0 = std::abs(Z[i] - z_length/2.);
    double r = std::sqrt(x0 * x0 + z0 * z0);
    if ( (R_patch < r)  && (r <= R_margin) ) {
    load_0(i) = 0; // mimic the free traction
    load_1(i) = 0;
    cohesion_0(i) = 0;
    cohesion_1(i) = 0;
    // Gamma_c(i) = 0;
    // tauc(i) = 0;
    // taur(i) = 0;
    MuS(i) = fp_patch; // apply the friction law same as patch
    MuK(i) = fr_patch;
    Dc(i) = dc_patch;
    Ds(i) = ds_patch;
    }
  }

  // Set background
  for (int i=0;i<mesh.getNbLocalNodes(); ++i){
    double x0 = std::abs(X[i] - x_length/2.);
    double z0 = std::abs(Z[i] - z_length/2.);
    double r = std::sqrt(x0 * x0 + z0 * z0);
    if (R_margin < r) {
    load_0(i) = tau_background;
    load_1(i) = sn_background;
    cohesion_0(i) = tau_background;
    cohesion_1(i) = sn_background;
    // Gamma_c(i) = Gc_background;
    // tauc(i) = tau_c_background; // set residual level 
    // taur(i) = tau_r_background;
    MuS(i) = fp_background;
    MuK(i) = fr_background;
    Dc(i) = dc_background;
    Ds(i) = ds_background;
    }
  }

  // dumping
  std::string bname = data.get<std::string>("simulation_id");
  // interface.initDump(bname,".");
  interface.initDump(bname,".",Dumper::Format::Binary);
  if (world_rank==0) {
    printf("Size of float: %d\n", sizeof(float));
  }
  
  interface.registerDumpFields(data.get<std::string>("dump_fields"));
  interface.dump(0,0);
  unsigned int dump_int = std::max(1, nb_time_steps/data.get<int>("nb_dumps"));

  // time stepping
  for (int s=1; s<=nb_time_steps; ++s) {
    if (world_rank==0) {
      // std::cout << "s=" << s << "/" << nb_time_steps << "\r";
      // Not work
      // NodalFieldComponent & velo0_top = interface.getTop().getVelo().component(0);
      // double * max_velo0_top = std::max({velo0_top});
      // std::cout << "max vx" << max_velo0_top << "m/s" << "\r";
      // std::cout.flush();
      printf("%d/%d\n", s, nb_time_steps);
    }

    // time integration
    interface.advanceTimeStep();

    // dump
    if (s % dump_int == 0)
      interface.dump(s,s*time_step);
  }

  comm->finalize();
  return 0;
}
