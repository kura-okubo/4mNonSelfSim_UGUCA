/**
 * @file   test_half_space.cc
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
#include "half_space_dynamic.hh"
#include "static_communicator_mpi.hh"

#include <iostream>
#include <iomanip>      // std::setprecision
#include <unistd.h>

using namespace uguca;

int main(){

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  //sleep(prank);  
  
  if (prank==0)
    std::cout << "start test: test_half_space" << std::endl;

  /* no need to test:
     initConvolutions
     forwardFFT
     backwardFFT
     gatherCostumMeshForwardFFT
     backwardFFTscatterCostumMesh
  */

  // --------------------------------------------------------------
  // init 2D half_space

  {
    double length = 2.0;
    int nb_elements = 16;

    SimpleMesh msh(length,nb_elements);
    Material mat = Material(71e9, 0.33, 2777);
    mat.readPrecomputedKernels();

    HalfSpaceDynamic hs2(mat,msh,1);

    // --------------------------------------------------------------
    if (prank==0)
      std::cout << "check getStableTimeStep (prank==0)" << std::endl;
    else
      std::cout << "check getStableTimeStep (prank>0)" << std::endl;
    
    double delta_x = msh.getDeltaX();
    double delta_z = msh.getDeltaZ();

    if (msh.getDim()==2)
      delta_z = 1e100;

    if (std::abs(hs2.getStableTimeStep() -
		 std::min(delta_x,delta_z) / mat.getCs())>1e-12) {
      std::cout << "wrong stable time step (prank=" << prank
		<< ")" << std::endl;
      return 1; // failure
    }


    double dt=hs2.getStableTimeStep();
    bool caught_exception = true;

    try {
      hs2.setTimeStep(2*dt);
      caught_exception = false;
    }
    catch (std::runtime_error &e) {
      if (prank==0)
	std::cout << "caught exception -> success (prank==0)" << std::endl;
      else
	std::cout << "caught exception -> success (prank>0)" << std::endl;
    }
    if (!caught_exception) {
      std::cout << "prank=" << prank << std::endl;
      std::cout << "failed" << std::endl;
      return 1; // failure
    }

    if (prank==0)
      std::cout << "check set/get TimeStep (prank==0)" << std::endl;
    else
      std::cout << "check set/get TimeStep (prank>0)" << std::endl;

    dt*=0.25;
    hs2.setTimeStep(dt);

    if (prank==0)
      std::cout << "check computeDisplacement (prank==0)" << std::endl;
    else
      std::cout << "check computeDisplacement (prank>0)" << std::endl;

    NodalFieldComponent & v0 = hs2.getVelo().component(0);
    NodalFieldComponent & v1 = hs2.getVelo().component(1);
    
    for (int i=0; i<v0.getNbNodes(); ++i){
      v0.set(i)=0.3*cos(i*3)+sin(i*2);
      v1.set(i)=0.5*cos(i*6)+0.4*(sin(i));
      // std::cout << (*v0)(i)<<std::endl;
    }
    
    hs2.computeDisplacement();
    NodalFieldComponent & u0 = hs2.getDisp().component(0);
    NodalFieldComponent & u1 = hs2.getDisp().component(1);

    for (int i=0; i<u0.getNbNodes(); ++i){
      double err = std::abs(u0.at(i) - dt*v0.at(i));
      if (err>1e-12) {
	std::cout << u0.at(i) << " != " << dt*v0.at(i)
		  << " diff=" << err << " (prank=="
		  << prank << ")" << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      err = std::abs(u1.at(i) - dt*v1.at(i));
      if (err>1e-12) {
	std::cout << u1.at(i) << " != " << dt*v1.at(i)
		  << " diff=" << err << " (prank=="
		  << prank << ")" << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    if (prank==0)
      std::cout << "check initPredictorCorrector (prank==0)" << std::endl;
    else
      std::cout << "check initPredictorCorrector (prank>0)" << std::endl;
    
    hs2.initPredictorCorrector();
    NodalFieldComponent & up0 = hs2.getDisp(true).component(0);
    NodalFieldComponent & up1 = hs2.getDisp(true).component(1);
    NodalFieldComponent & vp0 = hs2.getVelo(true).component(0);
    NodalFieldComponent & vp1 = hs2.getVelo(true).component(1);

    if (prank==0) {
      std::cout << "initPredictorCorrector success (prank==0)" << std::endl;
      std::cout << "check updateVelocity (prank==0)" << std::endl;
    }
    else {
      std::cout << "initPredictorCorrector success (prank>0)" << std::endl;
      std::cout << "check updateVelocity (prank>0)" << std::endl;
    }

    hs2.updateVelocity();
    for (int i=0; i<vp0.getNbNodes(); ++i){
      if (std::abs(vp0.at(i) - v0.at(i))>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs(vp1.at(i) - v1.at(i))>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    if (prank==0)
      std::cout << "check computeDisplacement pred (prank==0)" << std::endl;
    else
      std::cout << "check computeDisplacement pred (prank>0)" << std::endl;

    hs2.computeDisplacement(true);

    for (int i=0; i<up0.getNbNodes(); ++i){
      if (std::abs(up0.at(i) - dt*v0.at(i)-u0.at(i))>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs(up1.at(i) - dt*v1.at(i)-u1.at(i))>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    if (prank==0)
      std::cout << "check correct velocity (prank==0)" << std::endl;
    else
      std::cout << "check correct velocity (prank>0)" << std::endl;

    for (int i=0; i<v0.getNbNodes(); ++i){
      v0.set(i)*=3;
      v1.set(i)*=3;
    }
    hs2.correctVelocity(false);
    for (int i=0; i<vp0.getNbNodes(); ++i){
      if (std::abs(vp0.at(i)/2-v0.at(i)/3)>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs(vp1.at(i)/2-v1.at(i)/3)>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    if (prank==0)
      std::cout << "check correct velocity last (prank==0)" << std::endl;
    else
      std::cout << "check correct velocity last (prank>0)" << std::endl;
    
    hs2.correctVelocity(true);
    for (int i=0; i<vp0.getNbNodes(); ++i){
      if (std::abs(vp0.at(i)/2-v0.at(i)/2.5)>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs(vp1.at(i)/2-v1.at(i)/2.5)>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    if (prank==0)
      std::cout << "check initConvolutions (prank==0)" << std::endl;
    else
      std::cout << "check initConvolutions (prank>0)" << std::endl;
    
    hs2.initConvolutions();

    if (prank==0) {
      std::cout << "initConvolutions success (prank==0)" << std::endl;
      std::cout << "check computeStressFourierCoeff 2D (prank==0)" << std::endl;
    }
    else {
      std::cout << "initConvolutions success (prank>0)" << std::endl;
      std::cout << "check computeStressFourierCoeff 2D (prank>0)" << std::endl;
    }

    hs2.computeInternal();

    FFTableNodalFieldComponent & s0 = hs2.getInternal().component(0);
    FFTableNodalFieldComponent & s1 = hs2.getInternal().component(1);

    if (prank==0) { // real space computations are on 0 rank process
      if (false) {
	std::cout<<"solution"<<std::endl
		 << std::setprecision(12) << s0.fd(4)[0] << ", " << s0.fd(4)[1] << std::endl
		 << s1.fd(2)[0] << ", " << s1.fd(2)[1] << std::endl;
      }
      else {
	if (std::abs(s0.fd(4)[0]- (-138070.216931))>1e-6 ||
	    std::abs(s0.fd(4)[1]- (731156.525273))>1e-6) {
	  std::cout << "prank == " << prank << std::endl;
	  std::cout << "failed 4" << std::endl
		    << s0.fd(4)[0] << ", " << s0.fd(4)[1] << std::endl;
	  return 1; // failure
	}
	if (std::abs(s1.fd(2)[0]- (-44634.1170066))>1e-6 ||
	    std::abs(s1.fd(2)[1]- (7529.73118974))>1e-6) {
	  std::cout << "prank == " << prank << std::endl;
	  std::cout << std::setprecision(12) << "failed 2" << std::endl
		    << s1.fd(2)[0] << ", " << s1.fd(2)[1] << std::endl;
	  return 1; // failure
	}
      }
    }

    if (prank==0)
      std::cout << "check computeResidual (prank==0)" << std::endl;
    else
      std::cout << "check computeResidual (prank>0)" << std::endl;

    NodalField & u = hs2.getDisp();
    hs2.computeResidual(u);

    NodalFieldComponent & r0 = hs2.getResidual().component(0);
    NodalFieldComponent & r1 = hs2.getResidual().component(1);

    for (int i=0;i<r0.getNbNodes();++i){
      if (std::abs(r1.at(i)-(s1.at(i)+u1.at(i)))>1e-10) {
	std::cout<<r1.at(i)<<" != "<<(s1.at(i)+u1.at(i))<<std::endl;
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    if (prank==0)
      std::cout << "check computeVelocity (prank==0)" << std::endl;
    else
      std::cout << "check computeVelocity (prank>0)" << std::endl;
    
    hs2.computeVelocity();
    if (prank==0) // real space computations are on 0 rank process
    for (int i=0;i<v0.getNbNodes();++i){
      if (std::abs(v0.at(i) - ( mat.getCs()/mat.getShearModulus()*r0.at(i)) )>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed 0" << std::endl
		  << v0.at(i) << " != " << ( mat.getCs()/mat.getShearModulus()*r0.at(i)) << std::endl;
	return 1; // failure
      }
      if (std::abs(v1.at(i) - ( mat.getCs()/mat.getShearModulus()/(mat.getCp()/mat.getCs())*r1.at(i)) )>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed 1" << std::endl
		  << v1.at(i) << " != " << ( mat.getCs()/mat.getShearModulus()/ (mat.getCp()/mat.getCs()) *r1.at(i) ) << std::endl;
	return 1; // failure
      }
    }



  // --------------------------------------------------------------

  }
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

    HalfSpaceDynamic hs3(mat,msh3,1);

    double dt=hs3.getStableTimeStep()*0.1;

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

    hs3.computeInternal(); // destroys the fourier space

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
	if (std::abs(s0.at(4) - (-99443675301.1))>1e3) {
	  std::cout << "failed 4" << std::endl
		    << s0.at(4) << std::endl;
	  return 1; // failure
	}
	if (std::abs(s1.at(62) - (17019761277.7))>1e-1) {
	  std::cout << "failed 62" << std::endl
		    << s1.at(62) << std::endl;
	  return 1; // failure
	}
	if (std::abs(s2.at(47) - (5687137976.34))>1e-2) {
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
  return 0; // success
}
