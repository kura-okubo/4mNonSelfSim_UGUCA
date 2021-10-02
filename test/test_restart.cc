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
#include "nodal_field_component.hh"
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
  
  BimatInterface interface(mesh, top_mat, bot_mat, law);
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
  NodalFieldComponent nf1(mesh,"nf1");
  nf1.registerToRestart(restart_dump);
  nf1.registerToRestart(restart_load);
  
  double nf1v = 55.5;
  nf1.setAllValuesTo(nf1v);
  restart_dump.dump(rs_number);
  nf1.setAllValuesTo(nf1v*2); // set different value
  restart_load.load(rs_number);

  // check
  for (int i=0; i<nf1.getNbNodes(); ++i) {
    if (std::abs((nf1.at(i) - nf1v) / nf1v) > 1e-6) {
      std::cerr << "should be " << nf1v << ": " << nf1.at(i) << std::endl;
      return 1; // failure
    }
  }
  std::cout << "NodalFieldComponent correct -> success" << std::endl;
  
  // test binary format
  std::cout << "start: dump and reload NodalFieldComponent in binary" << std::endl;
  Restart restart_dump_binary("rs_binary",folder,BaseIO::Format::Binary);
  Restart restart_load_binary("rs_binary",folder,BaseIO::Format::Binary);
  nf1.registerToRestart(restart_dump_binary);
  nf1.registerToRestart(restart_load_binary);

  rs_number = 2;
  double nf2v = 66.6;
  nf1.setAllValuesTo(nf2v);
  restart_dump_binary.dump(rs_number);
  nf1.setAllValuesTo(nf2v*2); // set different value
  restart_load_binary.load(rs_number);

  // check
  for (int i=0; i<nf1.getNbNodes(); ++i) {
    if (std::abs((nf1.at(i) - nf2v)/nf2v) > 1e-6) {
      std::cerr << "should be " << nf2v << ": " << nf1.at(i) << std::endl;
      return 1; // failure
    }
  }
  std::cout << "NodalFieldComponent bindary correct -> success" << std::endl;


  // test dump and read of NodalField
  std::cout << "start: dump and reload NodalField" << std::endl;
  NodalField nf3(mesh, "nf3");
  nf3.registerToRestart(restart_dump);
  nf3.registerToRestart(restart_load);
  
  rs_number = 3;
  double nf3v = 77.7;
  nf3.setAllValuesTo(nf3v);
  restart_dump.dump(rs_number);
  nf3.setAllValuesTo(nf3v*2); // set different value
  restart_load.load(rs_number);

  // check
  for (int d=0; d<nf3.getDim(); ++d) {
    for (int i=0; i<nf3.getNbNodes(); ++i) {
      if (std::abs((nf3.component(d).at(i) - nf3v) / nf3v) > 1e-6) {
	std::cerr << "should be " << nf3v << ": " << nf3.component(d).at(i) << std::endl;
	return 1; // failure
      }
    }
  }
  std::cout << "NodalField correct -> success" << std::endl;

  // ---------------------------------------------------------------------------------------
  //   INTERFACE
  // ---------------------------------------------------------------------------------------
  
  // test dump and read of interface
  std::cout << "start: dump and reload Interface" << std::endl;
  rs_number = 4;
  double nf4v = 3.2;
  double nf4v2 = 1.2;
  interface.getCohesion().setAllValuesTo(nf4v);
  HalfSpaceDynamic & top = dynamic_cast<HalfSpaceDynamic&>(interface.getTop());
  double * Ur00 = top.getLimitedHistoryReal(0,1).getValues();
  int Ur00_size = top.getLimitedHistoryReal(0,1).getSize();
  std::fill_n(Ur00,Ur00_size,nf4v);
  double * Ui00 = top.getLimitedHistoryImag(0,1).getValues();
  int Ui00_size = top.getLimitedHistoryImag(0,1).getSize();
  std::fill_n(Ui00,Ui00_size,nf4v2);
  unsigned int nb_history_correct = top.getLimitedHistoryReal(0,1).getNbHistoryPoints();
  unsigned int index_now_correct = top.getLimitedHistoryReal(0,1).getIndexNow();

  // dump the solution
  interface.registerToRestart(restart_dump);
  interface.registerToRestart(restart_dump_binary);
  restart_dump.dump(rs_number);
  restart_dump_binary.dump(rs_number);

  // fill with other information
  top.getLimitedHistoryReal(0,1).addCurrentValue(44.);
  interface.getCohesion().setAllValuesTo(nf4v*2);
  std::fill_n(Ur00,Ur00_size,2*nf4v);
  std::fill_n(Ui00,Ui00_size,2*nf4v2);
  
  // reload from dump
  interface.registerToRestart(restart_load);
  restart_load.load(rs_number);

  // check ACII
  if (nb_history_correct != top.getLimitedHistoryReal(0,1).getNbHistoryPoints()) {
    std::cerr << "should be " << nb_history_correct
	      << ": " << top.getLimitedHistoryReal(0,1).getNbHistoryPoints()
	      << std::endl;
    return 1; // failure
  }
  if (index_now_correct != top.getLimitedHistoryReal(0,1).getIndexNow()) {
    std::cerr << "should be " << index_now_correct
	      << ": " << top.getLimitedHistoryReal(0,1).getIndexNow()
	      << std::endl;
    return 1; // failure
  }
      
  NodalField & to_check = interface.getCohesion();
  for (int d=0; d<to_check.getDim(); ++d) {
    for (int i=0; i<to_check.getNbNodes(); ++i) {
      if (std::abs((to_check.component(d).at(i) - nf4v) / nf4v) > 1e-6) {
	std::cerr << "should be " << nf4v << ": " << to_check.component(d).at(i) << std::endl;
	return 1; // failure
      }
    }
  }

  HalfSpaceDynamic & top_to_check = dynamic_cast<HalfSpaceDynamic&>(interface.getTop());
  Ur00 = top_to_check.getLimitedHistoryReal(0,1).getValues();
  Ur00_size = top_to_check.getLimitedHistoryReal(0,1).getSize();
  Ui00 = top_to_check.getLimitedHistoryImag(0,1).getValues();
  Ui00_size = top_to_check.getLimitedHistoryImag(0,1).getSize();
  
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
  top.getLimitedHistoryReal(0,1).addCurrentValue(44.);
  interface.getCohesion().setAllValuesTo(nf4v*2);
  std::fill_n(Ur00,Ur00_size,2*nf4v);
  std::fill_n(Ui00,Ui00_size,2*nf4v2);

  // reload from dump
  interface.registerToRestart(restart_load_binary);
  restart_load_binary.load(rs_number);

  // check binary
  if (nb_history_correct != top.getLimitedHistoryReal(0,1).getNbHistoryPoints()) {
    std::cerr << "nb_history_points should be " << nb_history_correct
	      << ": " << top.getLimitedHistoryReal(0,1).getNbHistoryPoints()
	      << std::endl;
    return 1; // failure
  }
  if (index_now_correct != top.getLimitedHistoryReal(0,1).getIndexNow()) {
    std::cerr << "index_now should be " << index_now_correct
	      << ": " << top.getLimitedHistoryReal(0,1).getIndexNow()
	      << std::endl;
    return 1; // failure
  }

  for (int d=0; d<to_check.getDim(); ++d) {
    for (int i=0; i<to_check.getNbNodes(); ++i) {
      if (std::abs((to_check.component(d).at(i) - nf4v) / nf4v) > 1e-6) {
	std::cerr << "should be " << nf4v << ": " << to_check.component(d).at(i) << std::endl;
	return 1; // failure
      }
    }
  }

  Ur00 = top_to_check.getLimitedHistoryReal(0,1).getValues();
  Ur00_size = top_to_check.getLimitedHistoryReal(0,1).getSize();
  Ui00 = top_to_check.getLimitedHistoryImag(0,1).getValues();
  Ui00_size = top_to_check.getLimitedHistoryImag(0,1).getSize();
  
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
  
  NodalFieldComponent wrong_nf1(wrong_mesh,"nf1");
  wrong_nf1.registerToRestart(wrong_restart_load);
  wrong_nf1.registerToRestart(wrong_restart_load_binary);
  
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

  
  LinearShearCohesiveLaw wrong_law(wrong_mesh, 1., 2e6);
  BimatInterface wrong_interface(wrong_mesh, top_mat, bot_mat, wrong_law);
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
