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
#include "uca_fftable_mesh.hh"
#include "nodal_field_component.hh"

#include <iostream>
/*
#include <stdio.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>
*/

using namespace uguca;


int main(){
  std::cout << "start test: test_restart" << std::endl;

  double Lx = 1;
  int Nx = 4;
  double Lz = 2;
  int Nz = 8;
  FFTableMesh mesh(Lx, Nx, Lz, Nz);

  // test init of restart
  std::cout << "start: init" << std::endl;
  Restart restart_dump("test_restart",".");
  Restart restart_load("test_restart",".");
  std::cout << "end: init" << std::endl;

  // test dump and read of NodalFieldComponent
  std::cout << "start: dump and reload NodalFieldComponent" << std::endl;
  int rs_number = 1;
  NodalFieldComponent nf1(mesh);
  double nf1v = 55.5;
  nf1.setAllValuesTo(nf1v);
  restart_dump.registerIO("nf1", nf1);
  restart_dump.dump(rs_number);
  nf1.setAllValuesTo(nf1v*2); // set different value

  restart_load.registerIO("nf1", nf1);
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
  rs_number = 2;
  double nf2v = 66.6;
  nf1.setAllValuesTo(nf2v);
  Restart restart_dump_binary("test_restart_binary",".", BaseIO::Format::Binary);
  restart_dump_binary.registerIO("nf1", nf1);
  restart_dump_binary.dump(rs_number);
  nf1.setAllValuesTo(nf2v*2); // set different value

  Restart restart_load_binary("test_restart_binary",".",BaseIO::Format::Binary);
  restart_load_binary.registerIO("nf1", nf1);
  restart_load_binary.load(rs_number);

  // check
  for (int i=0; i<nf1.getNbNodes(); ++i) {
    if (std::abs((nf1.at(i) - nf2v)/nf2v) > 1e-6) {
      std::cerr << "should be " << nf2v << ": " << nf1.at(i) << std::endl;
      return 1; // failure
    }
  }
  std::cout << "NodalFieldComponent bindary correct -> success" << std::endl;


  
  
  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
