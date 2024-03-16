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
    SpatialDirectionSet sds{_x,_y};
    Material mat = Material(71e9, 0.33, 2777);
    mat.readPrecomputedKernels();

    HalfSpaceDynamic hs2(mat,msh,1,sds);

    // --------------------------------------------------------------
    if (prank==0)
      std::cout << "check getStableTimeStep (prank==0)" << std::endl;
    else
      std::cout << "check getStableTimeStep (prank>0)" << std::endl;
    
    double delta_x = msh.getDelta(0);
    double delta_z = msh.getDelta(2);

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

    NodalField & v = hs2.getVelo();
    
    for (int i=0; i<v.getNbNodes(); ++i){
      v(i,0)=0.3*cos(i*3)+sin(i*2);
      v(i,1)=0.5*cos(i*6)+0.4*(sin(i));
      // std::cout << (*v0)(i)<<std::endl;
    }
    
    hs2.computeDisplacement();
    NodalField & u = hs2.getDisp();

    for (int i=0; i<u.getNbNodes(); ++i){
      double err = std::abs(u(i,0) - dt*v(i,0));
      if (err>1e-12) {
	std::cout << u(i,0) << " != " << dt*v(i,0)
		  << " diff=" << err << " (prank=="
		  << prank << ")" << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      err = std::abs(u(i,1) - dt*v(i,1));
      if (err>1e-12) {
	std::cout << u(i,1) << " != " << dt*v(i,1)
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
    NodalField & up = hs2.getDisp(true);
    NodalField & vp = hs2.getVelo(true);

    if (prank==0) {
      std::cout << "initPredictorCorrector success (prank==0)" << std::endl;
      std::cout << "check updateVelocity (prank==0)" << std::endl;
    }
    else {
      std::cout << "initPredictorCorrector success (prank>0)" << std::endl;
      std::cout << "check updateVelocity (prank>0)" << std::endl;
    }

    hs2.updateVelocity();
    for (int i=0; i<vp.getNbNodes(); ++i){
      if (std::abs(vp(i,0) - v(i,0))>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs(vp(i,1) - v(i,1))>1e-12) {
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

    for (int i=0; i<up.getNbNodes(); ++i){
      if (std::abs(up(i,0) - dt*v(i,0)-u(i,0))>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs(up(i,1) - dt*v(i,1)-u(i,1))>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    if (prank==0)
      std::cout << "check correct velocity (prank==0)" << std::endl;
    else
      std::cout << "check correct velocity (prank>0)" << std::endl;

    for (int i=0; i<v.getNbNodes(); ++i){
      v(i,0)*=3;
      v(i,1)*=3;
    }
    hs2.correctVelocity(false);
    for (int i=0; i<vp.getNbNodes(); ++i){
      if (std::abs(vp(i,0)/2-v(i,0)/3)>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs(vp(i,1)/2-v(i,1)/3)>1e-12) {
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
    for (int i=0; i<vp.getNbNodes(); ++i){
      if (std::abs(vp(i,0)/2-v(i,0)/2.5)>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs(vp(i,1)/2-v(i,1)/2.5)>1e-12) {
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

    FFTableNodalField & s = hs2.getInternal();

    if (prank==0) { // real space computations are on 0 rank process
      if (false) {
	std::cout<<"solution"<<std::endl
		 << std::setprecision(12) << s.fd(4,0)[0] << ", " << s.fd(4,0)[1] << std::endl
		 << s.fd(2,1)[0] << ", " << s.fd(2,1)[1] << std::endl;
      }
      else {
	if (std::abs(s.fd(4,0)[0]- (-138070.216931))>1e-6 ||
	    std::abs(s.fd(4,0)[1]- (731156.525273))>1e-6) {
	  std::cout << "prank == " << prank << std::endl;
	  std::cout << "failed 4" << std::endl
		    << s.fd(4,0)[0] << ", " << s.fd(4,0)[1] << std::endl;
	  return 1; // failure
	}
	if (std::abs(s.fd(2,1)[0]- (-44634.1170066))>1e-6 ||
	    std::abs(s.fd(2,1)[1]- (7529.73118974))>1e-6) {
	  std::cout << "prank == " << prank << std::endl;
	  std::cout << std::setprecision(12) << "failed 2" << std::endl
		    << s.fd(2,1)[0] << ", " << s.fd(2,1)[1] << std::endl;
	  return 1; // failure
	}
      }
    }

    if (prank==0)
      std::cout << "check computeResidual (prank==0)" << std::endl;
    else
      std::cout << "check computeResidual (prank>0)" << std::endl;

    NodalField & hs2_u = hs2.getDisp();
    hs2.computeResidual(u);

    NodalField & res = hs2.getResidual();

    for (int i=0;i<res.getNbNodes();++i){
      if (std::abs(res(i,1)-(s(i,1)+hs2_u(i,1)))>1e-10) {
	std::cout<<res(i,1)<<" != "<<(s(i,1)+hs2_u(i,1))<<std::endl;
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
    for (int i=0;i<v.getNbNodes();++i){
      if (std::abs(v(i,0) - ( mat.getCs()/mat.getShearModulus()*res(i,0)) )>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed 0" << std::endl
		  << v(i,0) << " != " << ( mat.getCs()/mat.getShearModulus()*res(i,0)) << std::endl;
	return 1; // failure
      }
      if (std::abs(v(i,1) - ( mat.getCs()/mat.getShearModulus()/(mat.getCp()/mat.getCs())*res(i,1)) )>1e-12) {
	std::cout << "prank=" << prank << std::endl;
	std::cout << "failed 1" << std::endl
		  << v(i,1) << " != " << ( mat.getCs()/mat.getShearModulus()/ (mat.getCp()/mat.getCs()) *res(i,1) ) << std::endl;
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

    HalfSpaceDynamic hs3(mat,msh3,1,{_x,_y,_z});

    double dt=hs3.getStableTimeStep()*0.1;

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

    hs3.computeInternal(); // destroys the fourier space

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
	if (std::abs(inter(4,0) - (-99443675301.1))>1e3) {
	  std::cout << "failed 4" << std::endl
		    << inter(4,0) << std::endl;
	  return 1; // failure
	}
	if (std::abs(inter(62,1) - (17019761277.7))>1e-1) {
	  std::cout << "failed 62" << std::endl
		    << inter(62,1) << std::endl;
	  return 1; // failure
	}
	if (std::abs(inter(47,2) - (5687137976.34))>1e-2) {
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
  return 0; // success
}
