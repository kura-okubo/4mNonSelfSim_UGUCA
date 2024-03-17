/**
 * @file   test_restart.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Sat Sept 25 2021
 * @date last modification: Sat Sept 25 2021
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

#include "uca_restart.hh"
#include "uca_simple_mesh.hh"
#include "material.hh"
#include "linear_shear_cohesive_law.hh"
#include "bimat_interface.hh"
#include "half_space_dynamic.hh"

#include <iostream>
#include <cmath>
#include <cstring>
#include <stdexcept>

using namespace uguca;


int main(int argc, char *argv[]) {

  std::string folder = ".";
  if (argc==2) {
    folder = argv[1];
  }
    
  std::cout << "start test: test_restart" << std::endl;
  
  double Lx = 1;
  int Nx = 4;
  double Lz = 2;
  int Nz = 8;
  SimpleMesh mesh(Lx, Nx, Lz, Nz);
  
  // materials
  Material top_mat = Material(1e9, 0.25, 1000);
  top_mat.readPrecomputedKernels();
  Material bot_mat = Material(2e9, 0.33, 2000);
  bot_mat.readPrecomputedKernels();

  LinearShearCohesiveLaw law(mesh, 1., 2e6);
  
  BimatInterface interface(mesh, {_x,_y,_z}, top_mat, bot_mat, law);
  interface.setTimeStep(0.3*interface.getStableTimeStep());
  interface.init();
  
  // test init of restart
  std::cout << "start: init" << std::endl;
  Restart restart_dump("rs1",folder);
  Restart restart_load("rs1",folder);
  std::cout << "end: init" << std::endl;

  // test dump and read of NodalFieldComponent
  std::cout << "start: dump and reload NodalFieldComponent" << std::endl;
  int rs_number = 1;
  NodalField nf1(mesh);
  nf1.setName("nf1");
  restart_dump.registerIO(nf1);
  restart_load.registerIO(nf1);
  
  double nf1v = 55.5;
  nf1.setAllValuesTo(nf1v);
  restart_dump.dump(rs_number);
  nf1.setAllValuesTo(nf1v*2); // set different value
  restart_load.load(rs_number);

  // check
  for (int i=0; i<nf1.getNbNodes(); ++i) {
    if (std::abs((nf1(i) - nf1v) / nf1v) > 1e-6) {
      std::cerr << "should be " << nf1v << ": " << nf1(i) << std::endl;
      return 1; // failure
    }
  }
  std::cout << "NodalFieldComponent correct -> success" << std::endl;
  
  // test binary format
  std::cout << "start: dump and reload NodalFieldComponent in binary" << std::endl;
  Restart restart_dump_binary("rs_binary",folder,BaseIO::Format::Binary);
  Restart restart_load_binary("rs_binary",folder,BaseIO::Format::Binary);
  restart_dump_binary.registerIO(nf1);
  restart_load_binary.registerIO(nf1);

  rs_number = 2;
  double nf2v = 66.6;
  nf1.setAllValuesTo(nf2v);
  restart_dump_binary.dump(rs_number);
  nf1.setAllValuesTo(nf2v*2); // set different value
  restart_load_binary.load(rs_number);

  // check
  for (int i=0; i<nf1.getNbNodes(); ++i) {
    if (std::abs((nf1(i) - nf2v)/nf2v) > 1e-6) {
      std::cerr << "should be " << nf2v << ": " << nf1(i) << std::endl;
      return 1; // failure
    }
  }
  std::cout << "NodalFieldComponent bindary correct -> success" << std::endl;


  // test dump and read of NodalField
  std::cout << "start: dump and reload NodalField" << std::endl;
  NodalField nf3(mesh);
  nf3.setName("nf3");
  restart_dump.registerIO(nf3);
  restart_load.registerIO(nf3);
  
  rs_number = 3;
  double nf3v = 77.7;
  nf3.setAllValuesTo(nf3v);
  restart_dump.dump(rs_number);
  nf3.setAllValuesTo(nf3v*2); // set different value
  restart_load.load(rs_number);

  // check
  for (const auto& d : nf3.getComponents()) {
    for (int i=0; i<nf3.getNbNodes(); ++i) {
      if (std::abs((nf3(i,d) - nf3v) / nf3v) > 1e-6) {
	std::cerr << "should be " << nf3v << ": " << nf3(i,d) << std::endl;
	return 1; // failure
      }
    }
  }
  std::cout << "NodalField correct -> success" << std::endl;

  // ---------------------------------------------------------------------------------
  //   INTERFACE
  // ---------------------------------------------------------------------------------
  
  // test dump and read of interface
  std::cout << "start: dump and reload Interface" << std::endl;
  rs_number = 4;
  double nf4v = 3.2;
  double nf4v2 = 1.2;
  interface.getCohesion().setAllValuesTo(nf4v);
  HalfSpaceDynamic & top = dynamic_cast<HalfSpaceDynamic&>(interface.getTop());
  HistFFTableNodalField & top_disp = dynamic_cast<HistFFTableNodalField&>(top.getDisp());
  double * Ur00 = const_cast<double*>(top_disp.hist(1,0).real());
  int Ur00_size = top_disp.hist(1,0).getSize();
  std::fill_n(Ur00,Ur00_size,nf4v);
  double * Ui00 = const_cast<double*>(top_disp.hist(1,0).imag());
  int Ui00_size = top_disp.hist(1,0).getSize();
  std::fill_n(Ui00,Ui00_size,nf4v2);
  unsigned int nb_history_correct = top_disp.hist(1,0).getNbHistoryPoints();
  unsigned int index_now_correct = top_disp.hist(1,0).getIndexNow();

  // dump the solution
  interface.registerToRestart(restart_dump);
  interface.registerToRestart(restart_dump_binary);
  restart_dump.dump(rs_number);
  restart_dump_binary.dump(rs_number);

  // fill with other information
  fftw_complex fortyfour = {44., 0};
  const_cast<ModalLimitedHistory&>(top_disp.hist(1,0)).addCurrentValue(fortyfour);
  interface.getCohesion().setAllValuesTo(nf4v*2);
  std::fill_n(Ur00,Ur00_size,2*nf4v);
  std::fill_n(Ui00,Ui00_size,2*nf4v2);
  // reload from dump
  interface.registerToRestart(restart_load);
  restart_load.load(rs_number);

  // check ACII
  if (nb_history_correct != top_disp.hist(1,0).getNbHistoryPoints()) {
    std::cerr << "should be " << nb_history_correct
	      << ": " << top_disp.hist(1,0).getNbHistoryPoints()
	      << std::endl;
    return 1; // failure
  }
  if (index_now_correct != top_disp.hist(1,0).getIndexNow()) {
    std::cerr << "should be " << index_now_correct
	      << ": " << top_disp.hist(1,0).getIndexNow()
	      << std::endl;
    return 1; // failure
  }
      
  NodalField & to_check = interface.getCohesion();
  for (const auto& d : to_check.getComponents()) {
    for (int i=0; i<to_check.getNbNodes(); ++i) {
      if (std::abs((to_check(i,d) - nf4v) / nf4v) > 1e-6) {
	std::cerr << "should be " << nf4v << ": " << to_check(i,d) << std::endl;
	return 1; // failure
      }
    }
  }

  HalfSpaceDynamic & top_to_check = dynamic_cast<HalfSpaceDynamic&>(interface.getTop());
  HistFFTableNodalField & top_disp_to_check = dynamic_cast<HistFFTableNodalField&>(top_to_check.getDisp());
  Ur00 = const_cast<double*>(top_disp_to_check.hist(1,0).real());
  Ur00_size = top_disp_to_check.hist(1,0).getSize();
  Ui00 = const_cast<double*>(top_disp_to_check.hist(1,0).imag());
  Ui00_size = top_disp_to_check.hist(1,0).getSize();
  
  for (int i=0; i<Ur00_size; ++i) {
    if (std::abs((Ur00[i] - nf4v) / nf4v) > 1e-6) {
      std::cerr << "Ur should be " << nf4v << ": " << Ur00[i] << std::endl;
      return 1; // failure
    }
  }
  
  for (int i=0; i<Ui00_size; ++i) {
    if (std::abs((Ui00[i] - nf4v2) / nf4v2) > 1e-6) {
      std::cerr << "Ui should be " << nf4v2 << ": " << Ui00[i] << std::endl;
      return 1; // failure
    }
  }

  // fill with other information
  const_cast<ModalLimitedHistory&>(top_disp.hist(1,0)).addCurrentValue(fortyfour);
  interface.getCohesion().setAllValuesTo(nf4v*2);
  std::fill_n(Ur00,Ur00_size,2*nf4v);
  std::fill_n(Ui00,Ui00_size,2*nf4v2);

  // reload from dump
  interface.registerToRestart(restart_load_binary);
  restart_load_binary.load(rs_number);

  // check binary
  if (nb_history_correct != top_disp.hist(1,0).getNbHistoryPoints()) {
    std::cerr << "nb_history_points should be " << nb_history_correct
	      << ": " << top_disp.hist(1,0).getNbHistoryPoints()
	      << std::endl;
    return 1; // failure
  }
  if (index_now_correct != top_disp.hist(1,0).getIndexNow()) {
    std::cerr << "index_now should be " << index_now_correct
	      << ": " << top_disp.hist(1,0).getIndexNow()
	      << std::endl;
    return 1; // failure
  }

  for (const auto& d : to_check.getComponents()) {
    for (int i=0; i<to_check.getNbNodes(); ++i) {
      if (std::abs((to_check(i,d) - nf4v) / nf4v) > 1e-6) {
	std::cerr << "should be " << nf4v << ": " << to_check(i,d) << std::endl;
	return 1; // failure
      }
    }
  }

  Ur00 = const_cast<double*>(top_disp_to_check.hist(1,0).real());
  Ur00_size = top_disp_to_check.hist(1,0).getSize();
  Ui00 = const_cast<double*>(top_disp_to_check.hist(1,0).imag());
  Ui00_size = top_disp_to_check.hist(1,0).getSize();
  
  for (int i=0; i<Ur00_size; ++i) {
    if (std::abs((Ur00[i] - nf4v) / nf4v) > 1e-6) {
      std::cerr << "Ur should be " << nf4v << ": " << Ur00[i] << std::endl;
      return 1; // failure
    }
  }
  
  for (int i=0; i<Ui00_size; ++i) {
    if (std::abs((Ui00[i] - nf4v2) / nf4v2) > 1e-6) {
      std::cerr << "Ui should be " << nf4v2 << ": " << Ui00[i] << std::endl;
      return 1; // failure
    }
  }
  std::cout << "Interface correct -> success" << std::endl;


  // load into wrong simulation
  SimpleMesh wrong_mesh(Lx, 2*Nx, Lz, 2*Nz);
  Restart wrong_restart_load("rs1",folder);
  Restart wrong_restart_load_binary("rs_binary",folder,BaseIO::Format::Binary);
  
  NodalField wrong_nf1(wrong_mesh);
  wrong_nf1.setName("nf1");
  wrong_restart_load.registerIO(wrong_nf1);
  wrong_restart_load_binary.registerIO(wrong_nf1);
  
  rs_number = 1;
  bool caught_exception = true;
  try {
    wrong_restart_load.load(rs_number);
    caught_exception = false;
  }
  catch (std::runtime_error &e) {
    std::cout << "caught exception -> success" << std::endl;
  }
  if (!caught_exception) {
    std::cerr << "did not find wrong load of NFC" << std::endl;
    return 1; // failure
  }

  // not possible anymore since Limitedhistory is dumped in single file
  /*
  caught_exception = true;
  try {
    wrong_restart_load_binary.load(rs_number);
    caught_exception = false;
  }
  catch (std::runtime_error &e) {
    std::cout << "caught exception -> success" << std::endl;
  }
  if (!caught_exception) {
    std::cerr << "did not find wrong load of NFC (binary)" << std::endl;
    return 1; // failure
  }
  */

  
  LinearShearCohesiveLaw wrong_law(wrong_mesh, 1., 2e6);
  BimatInterface wrong_interface(wrong_mesh, {_x,_y,_z}, top_mat, bot_mat, wrong_law);
  wrong_interface.setTimeStep(0.3*wrong_interface.getStableTimeStep());
  wrong_interface.init();
  wrong_interface.registerToRestart(wrong_restart_load);
  wrong_interface.registerToRestart(wrong_restart_load_binary);
  
  rs_number = 4;
  caught_exception = true;
  try {
    wrong_restart_load.load(rs_number);
    caught_exception = false;
  }
  catch (std::runtime_error &e) {
    std::cout << "caught exception -> success" << std::endl;
  }
  if (!caught_exception) {
    std::cerr << "did not find wrong load of interface" << std::endl;
    return 1; // failure
  }

  caught_exception = true;
  try {
    wrong_restart_load_binary.load(rs_number);
    caught_exception = false;
  }
  catch (std::runtime_error &e) {
    std::cout << "caught exception -> success" << std::endl;
  }
  if (!caught_exception) {
    std::cerr << "did not find wrong load of interface (binary)" << std::endl;
    return 1; // failure
  }

  
  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
