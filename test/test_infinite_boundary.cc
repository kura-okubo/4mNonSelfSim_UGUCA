/**
 * @file   test_inf_bc.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Tue Jun 30 2020
 * @date last modification: Tue Jun 30 2020
 *
 * @brief  TODO
 *
 *
 * Copyright (C) 2020 ETH Zurich (David S. Kammer)
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
#include "uca_custom_mesh.hh"
#include "infinite_boundary.hh"

#include <iostream>
#include <iomanip>      // std::setprecision
#include <cmath>
#include <sstream>
#include <vector>
#include <unistd.h>

using namespace uguca;

int main() {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();

  //----------------------------------------------------------------//
  // initialize infinite boundary

  //material
  double inf_E   = 5e6;
  double inf_nu  = 0.25;
  double inf_rho = 1e3;

  Material mat = Material(inf_E,inf_nu,inf_rho);
  mat.readPrecomputedKernels();

  //dimension
  unsigned int inf_nb_nodes_x = 3;
  double inf_length_x = 3.0;
  unsigned int inf_nb_nodes_y = 4;
  double inf_length_y = 4.0;

  double time_step=0.001;
  if (prank==0)
    std::cout << "time_step = " << time_step << std::endl;
  int side_factor=1;

  std::vector<double> x_coord_tmp;
  std::vector<double> z_coord_tmp;

  if (psize==1) {
    x_coord_tmp = {2.0,  // 9
		   0.0,  // 1
		   0.0,  // 0
		   0.0,  // 2
		   0.0,  // 3
		   1.0,  // 4
		   1.0,  // 5
		   1.0,  // 6
		   1.0,  // 7
		   2.0,  // 8
		   2.0,  // 10
		   2.0}; // 11

    z_coord_tmp = {1.0,  // 9
		   1.0,  // 1
		   0.0,  // 0
		   2.0,  // 2
		   3.0,  // 3
		   0.0,  // 4
		   1.0,  // 5
		   2.0,  // 6
		   3.0,  // 7
		   0.0,  // 8
		   2.0,  // 10
		   3.0}; // 11
  }
  else {
    if (prank==0) {
      x_coord_tmp = {2.0,   // 9
		     0.0,   // 1
		     0.0,   // 0
		     0.0,   // 2
		     0.0};  // 3
      
      z_coord_tmp = {1.0,   // 9
		     1.0,   // 1
		     0.0,   // 0
		     2.0,   // 2
		     3.0};  // 3
    }
    if (prank==1) {
      x_coord_tmp = {1.0,   // 4
		     1.0,   // 5
		     1.0,   // 6
		     1.0,   // 7
		     2.0,   // 8
		     2.0,   // 10
		     2.0};  // 11
      
      z_coord_tmp = {0.0,   // 4
		     1.0,   // 5
		     2.0,   // 6
		     3.0,   // 7
		     0.0,   // 8
		     2.0,   // 10
		     3.0};  // 11
    }
  }
  
  CustomMesh mesh(inf_length_x, inf_nb_nodes_x,
		  inf_length_y, inf_nb_nodes_y,
		  x_coord_tmp, z_coord_tmp);
  int dim = mesh.getDim();
  const TwoDVector & coord_tmp = mesh.getLocalCoords();

  InfiniteBoundary infinite_boundary(mesh,
				     {_x,_y,_z},
				     side_factor,
				     mat);

  // time step
  infinite_boundary.setTimeStep(time_step);
  infinite_boundary.initPredictorCorrector();
  infinite_boundary.initConvolutions();
  

  NodalField & inf_bc_ext = infinite_boundary.getExternal();
  FFTableNodalField & inf_bc_dsp = infinite_boundary.getDisp();
  NodalField & inf_bc_vel = infinite_boundary.getVelo();

  inf_bc_ext.zeros();
  inf_bc_dsp.zeros();
  
  // populate fields
  for (int i=0; i<dim; ++i) {
    for (int n=0; n<inf_bc_vel.getNbNodes(); ++n){
      inf_bc_vel(n,i)=coord_tmp(n,0)*coord_tmp(n,2)+coord_tmp(n,0)+coord_tmp(n,2);
    }
  }

  
  std::cout<<"test advance time step Neumann "<<prank<<std::endl;
  infinite_boundary.advanceTimeStepNeumann();

  double mu = mat.getShearModulus();
  double Cs = mat.getCs();
  double Cp = mat.getCp();
  std::vector<double> eta = {1.0, Cp / Cs, 1.0};
  for (int i=0; i<dim; ++i) {
    for (int n=0; n<inf_bc_ext.getNbNodes(); ++n){
      if(std::abs(inf_bc_ext(n,i)-
		  (- side_factor * mu/Cs*eta[i]*inf_bc_vel(n,i)))>1e-12){
	std::cout<<"error "<<std::endl;
	return 1;
      }
    }
  }
  std::cout<<"advance time step Neumann success! "<<prank<<std::endl;  


  std::cout<<"test advance time step Dirichlet "<<prank<<std::endl;
  infinite_boundary.advanceTimeStepDirichlet();
  {
    NodalField & disp = infinite_boundary.getDisp();
    
    if (false) {
      std::cout<<"solution"<<std::endl
	       << std::setprecision(12)
	       << disp(4,0) << ", "
	       << disp(9,0) << std::endl;
    }
    else {
      if (psize==1) {
	if (std::abs(disp(4,0)- 0.003)>1e-6 ||
	    std::abs(disp(9,0)- 0.002)>1e-6) {
	  std::cout << "failed" << std::endl
		    << disp(4,0) << ", "
		    << disp(9,0) << std::endl;
	  return 1; // failure
	}
      }
      else {
	if ((prank==0 && std::abs(disp(4,0)- 0.003)>1e-6) ||
	    (prank==1 && std::abs(disp(4,0)- 0.002)>1e-6)) {
	  std::cout << "prank=" << prank << ": " << disp(4,0) << std::endl << std::flush;
	  std::cout << "failed" << std::endl;
	  return 1; // failure
	}
      }
    }
  }

  
  std::cout<<"test predict time step Dirichlet "<<prank<<std::endl;  
  infinite_boundary.predictTimeStepDirichlet();
  {
    NodalField & disp = infinite_boundary.getDisp(1);

    if (false) {
      std::cout<<"solution"<<std::endl
	       << std::setprecision(12)
	       << disp(4,0) << ", "
	       << disp(9,0) << std::endl;
    }
    else {
      if (psize==1) {
	if (std::abs(disp(4,0)- (-2.00839640955e-05))>1e-12 ||
	    std::abs(disp(9,0)- (-2.86716308883e-07))>1e-12) {
	  std::cout << "failed" << std::endl
		    << disp(4,0) << ", "
		    << disp(9,0) << std::endl;
	  return 1; // failure
	}
      }
      else {
	if ((prank==0 && std::abs(disp(4,0)- (-2.00839640955e-05))>1e-12) ||
	    (prank==1 && std::abs(disp(4,0)- (-2.86716308883e-07))>1e-12)) {
	  std::cout << "prank=" << prank << ": " << disp(4,0) << std::endl << std::flush;
	  std::cout << "failed" << std::endl;
	  return 1; // failure
	}
      }
    }
  }

  std::cout<<"went til the end rank "<<prank<<"\n";

  std::cout << "all checks passed -> overall success "<<prank << std::endl;
  
  StaticCommunicatorMPI::getInstance()->finalize();
  std::cout << "finalized " << prank << std::endl;
  return 0;
}
