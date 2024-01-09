/**
 * @file   test_half_space_quasi_dynamic.cc
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

#include "uca_simple_mesh.hh"
#include "half_space_quasi_dynamic.hh"
#include "static_communicator_mpi.hh"

#include <iostream>
#include <iomanip>      // std::setprecision
#include <unistd.h>

using namespace uguca;

int main(){

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  //sleep(prank);  
  
  if (prank==0)
    std::cout << "start test: test_half_space_quasi_dynamic" << std::endl;
  // --------------------------------------------------------------
  // init 2D half_space

  {
    double length = 2.0;
    int nb_elements = 16;

    SimpleMesh msh(length,nb_elements);
    Material mat = Material(71e9, 0.33, 2777);
    mat.readPrecomputedKernels();

    HalfSpaceQuasiDynamic hs2(mat,msh,1,{_x,_y});

    NodalField & disp = hs2.getDisp();
    for (int i=0; i<disp.getNbNodes(); ++i){
      disp(i,0)=0.3*cos(i*3)+sin(i*2);
      disp(i,1)=0.5*cos(i*6)+0.4*(sin(i));
      // std::cout << (*u0)(i)<<std::endl;
    }

    double dx = msh.getDelta(0);
    double cs = mat.getCs();
    double dt = dx/cs*0.1;

    hs2.setTimeStep(dt);

    hs2.initPredictorCorrector();
    hs2.initConvolutions();

    if (prank==0) {
      std::cout << "initConvolutions success (prank==0)" << std::endl;
      std::cout << "check computeStressFourierCoeff 2D (prank==0)" << std::endl;
    }
    else {
      std::cout << "initConvolutions success (prank>0)" << std::endl;
      std::cout << "check computeStressFourierCoeff 2D (prank>0)" << std::endl;
    }

    hs2.computeInternal(false,false,false); // predicting, correcting, dynamic

    FFTableNodalField & inter = hs2.getInternal();

    if (prank==0) { // real space computations are on 0 rank process
      if (false) {
	std::cout<<"solution"<<std::endl
		 << std::setprecision(12) << inter.fd(4,0)[0] << ", " << inter.fd(4,0)[1] << std::endl
		 << inter.fd(2,1)[0] << ", " << inter.fd(2,1)[1] << std::endl;
      }
      else {
	if (std::abs(inter.fd(4,0)[0]- (-16538131902.7))>1e0 ||
	    std::abs(inter.fd(4,0)[1]- (196679920694))>1e0) {
	  std::cout << "prank == " << prank << std::endl;
	  std::cout << "failed 4" << std::endl
		    << inter.fd(4,0)[0] << ", " << inter.fd(4,0)[1] << std::endl;
	  return 1; // failure
	}
	if (std::abs(inter.fd(2,1)[0]- (-698803824158))>1e0 ||
	    std::abs(inter.fd(2,1)[1]- (153102174931))>1e0) {
	  std::cout << "prank == " << prank << std::endl;
	  std::cout << std::setprecision(12) << "failed 2" << std::endl
		    << inter.fd(2,1)[0] << ", " << inter.fd(2,1)[1] << std::endl;
	  return 1; // failure
	}
      }
    }
  }

    // --------------------------------------------------------------
{
    if (prank==0)
      std::cout << "check computeStressFourierCoeff 3D (prank==0)" << std::endl;
    else
      std::cout << "check computeStressFourierCoeff 3D (prank>0)" << std::endl;

    // init 3D half_space
    double lx = 2.0;
    double lz = 1.5;
    int nb_x = 16;
    int nb_z = 8;
    SimpleMesh msh3(lx,nb_x,
		    lz,nb_z);
    Material mat = Material(71e9, 0.33, 2777);
    mat.readPrecomputedKernels();

    HalfSpaceQuasiDynamic hs3(mat,msh3,1,{_x,_y,_z});

    double dx = msh3.getDelta(0);
    double cs = mat.getCs();
    double dt = dx/cs*0.1;

    hs3.setTimeStep(dt);

    hs3.initPredictorCorrector();
    hs3.initConvolutions();

    NodalField & disp = hs3.getDisp();

    if (prank==0) {
      for (int i=0; i<nb_x; ++i){
	for (int j=0; j<nb_z; ++j){
	  int ij =i*nb_z+j;
	  disp(ij,0)=0.3*cos(i*3+6)+sin(i*2+1) +2*cos(j*2)+sin(j*6-2);
	  disp(ij,1)=0.5*cos(i*6)+0.4*(sin(i+2))-1*cos(j*7)+sin(j*9-5);
	  disp(ij,2)=0.5*cos(i*6+5)+0.4*(sin(i-2))+0.5*cos(5-j*2)+sin(j);;
	}
      }
    }

    // destroys the fourier space
    hs3.computeInternal(false,false,false); // predicting, correcting, dynamic

    if (prank==0) // complete data is gathered to process 0
    {
      std::cout<<"test computeStressFourierCoeff 3D (prank==0)"<<std::endl;
      
      FFTableNodalField & inter = hs3.getInternal();
      
      if (false) {
	std::cout<<"solution"<<std::endl
		 << std::setprecision(12)
		 << inter(4,0) << std::endl
		 << inter(62,1) << std::endl
		 << inter(47,2) << std::endl;
      }
      else {
	if (std::abs(inter(4,0) - (-493509032930))>1e0) {
	  std::cout << "failed 4" << std::endl
		    << inter(4,0) << std::endl;
	  return 1; // failure
	}
	if (std::abs(inter(62,1) - (627247523903))>1e0) {
	  std::cout << "failed 62" << std::endl
		    << inter(62,1) << std::endl;
	  return 1; // failure
	}
	if (std::abs(inter(47,2) - (66306881961.8))>1e0) {
	  std::cout << "failed 47" << std::endl
		    << inter(47,2) << std::endl;
	  return 1; // failure
	}
      }
    }
  }

  if (prank==0)
    std::cout << "all checks passed -> overall success (prank==0)" << std::endl;
  else
    std::cout << "all checks passed -> overall success (prank>0)" << std::endl;
  
  StaticCommunicatorMPI::getInstance()->finalize();
  return 0;
}
