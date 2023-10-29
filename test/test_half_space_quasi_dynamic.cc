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

    HalfSpaceQuasiDynamic hs2(mat,msh,1);

    NodalFieldComponent & u0 = hs2.getDisp().component(0);
    NodalFieldComponent & u1 = hs2.getDisp().component(1);
    for (int i=0; i<u0.getNbNodes(); ++i){
      u0.set(i)=0.3*cos(i*3)+sin(i*2);
      u1.set(i)=0.5*cos(i*6)+0.4*(sin(i));
      // std::cout << (*u0)(i)<<std::endl;
    }

    double dx = msh.getDeltaX();
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

    FFTableNodalFieldComponent & s0 = hs2.getInternal().component(0);
    FFTableNodalFieldComponent & s1 = hs2.getInternal().component(1);

    if (prank==0) { // real space computations are on 0 rank process
      if (false) {
	std::cout<<"solution"<<std::endl
		 << std::setprecision(12) << s0.fd(4)[0] << ", " << s0.fd(4)[1] << std::endl
		 << s1.fd(2)[0] << ", " << s1.fd(2)[1] << std::endl;
      }
      else {
	if (std::abs(s0.fd(4)[0]- (-16538131902.7))>1e0 ||
	    std::abs(s0.fd(4)[1]- (196679920694))>1e0) {
	  std::cout << "prank == " << prank << std::endl;
	  std::cout << "failed 4" << std::endl
		    << s0.fd(4)[0] << ", " << s0.fd(4)[1] << std::endl;
	  return 1; // failure
	}
	if (std::abs(s1.fd(2)[0]- (-698803824158))>1e0 ||
	    std::abs(s1.fd(2)[1]- (153102174931))>1e0) {
	  std::cout << "prank == " << prank << std::endl;
	  std::cout << std::setprecision(12) << "failed 2" << std::endl
		    << s1.fd(2)[0] << ", " << s1.fd(2)[1] << std::endl;
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

    HalfSpaceQuasiDynamic hs3(mat,msh3,1);

    double dx = msh3.getDeltaX();
    double cs = mat.getCs();
    double dt = dx/cs*0.1;

    hs3.setTimeStep(dt);

    hs3.initPredictorCorrector();
    hs3.initConvolutions();

    NodalFieldComponent & u0 = hs3.getDisp().component(0);
    NodalFieldComponent & u1 = hs3.getDisp().component(1);
    NodalFieldComponent & u2 = hs3.getDisp().component(2);

    if (prank==0) {
      for (int i=0; i<nb_x; ++i){
	for (int j=0; j<nb_z; ++j){
	  int ij =i*nb_z+j;
	  u0.set(ij)=0.3*cos(i*3+6)+sin(i*2+1) +2*cos(j*2)+sin(j*6-2);
	  u1.set(ij)=0.5*cos(i*6)+0.4*(sin(i+2))-1*cos(j*7)+sin(j*9-5);
	  u2.set(ij)=0.5*cos(i*6+5)+0.4*(sin(i-2))+0.5*cos(5-j*2)+sin(j);;
	}
      }
    }

    // destroys the fourier space
    hs3.computeInternal(false,false,false); // predicting, correcting, dynamic

    if (prank==0) // complete data is gathered to process 0
    {
      std::cout<<"test computeStressFourierCoeff 3D (prank==0)"<<std::endl;
      
      FFTableNodalFieldComponent & s0 = hs3.getInternal().component(0);
      FFTableNodalFieldComponent & s1 = hs3.getInternal().component(1);
      FFTableNodalFieldComponent & s2 = hs3.getInternal().component(2);
      
      if (false) {
	std::cout<<"solution"<<std::endl
		 << std::setprecision(12)
		 << s0.at(4) << std::endl
		 << s1.at(62) << std::endl
		 << s2.at(47) << std::endl;
      }
      else {
	if (std::abs(s0.at(4) - (-493509032930))>1e0) {
	  std::cout << "failed 4" << std::endl
		    << s0.at(4) << std::endl;
	  return 1; // failure
	}
	if (std::abs(s1.at(62) - (627247523903))>1e0) {
	  std::cout << "failed 62" << std::endl
		    << s1.at(62) << std::endl;
	  return 1; // failure
	}
	if (std::abs(s2.at(47) - (66306881961.8))>1e0) {
	  std::cout << "failed 47" << std::endl
		    << s2.at(47) << std::endl;
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
