/**
 * @file   test_nodal_field.cc
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

#include "uca_fftable_mesh.hh"
#include "nodal_field.hh"

#include <iostream>
#include <cmath>

using namespace uguca;

int main(){
  std::cout << "start test: nodal_field" << std::endl;

  int dim = 2;
  int sz = 4;
  FFTableMesh mesh(dim,sz);
  
  std::cout << "check initialization and getSize" << std::endl;

  NodalField nf1(mesh, {_x,_y});
  if (nf1.getNbNodes() != sz) {
    std::cerr << "wrong size" << std::endl;
    return 1; // failure
  }
  if (nf1.getNbComponents() != dim) {
    std::cerr << "wrong dim" << std::endl;
    return 1; // failure
  }
  std::cout << "size correct -> success" << std::endl;

  std::cout << "check zeros" << std::endl;
  nf1.zeros();
  for (const auto& d : nf1.getComponents()) {
    for (int i=0; i<nf1.getNbNodes(); ++i) {
      if (nf1(i,d) != 0) {
	std::cerr << "should be zero: " << nf1(i,d) << std::endl;
	return 1; // failure
      }
    }
  }
  std::cout << "zeros correct -> success" << std::endl;

  std::cout << "check setAllValuesTo" << std::endl;
  double vl=3.5;
  nf1.setAllValuesTo(vl);
  for (const auto& d : nf1.getComponents()) {
    for (int i=0; i<nf1.getNbNodes(); ++i) {
      if (nf1(i,d) != vl) {
	std::cerr << "should be " << vl << ": " << nf1(i,d) << std::endl;
	return 1; // failure
      }
    }
  }

  std::vector<double> vls = {5.6,7.5};
  nf1.setAllValuesTo(vls[0],0);
  nf1.setAllValuesTo(vls[1],1);
  for (const auto& d : nf1.getComponents()) {
    for (int i=0; i<nf1.getNbNodes(); ++i) {
      if (nf1(i,d) != vls[d]) {
	std::cerr << "should be " << vls[d] << ": " << nf1(i,d) << std::endl;
	return 1; // failure
      }
    }
  }
  std::cout << "setAllValuesTo correct -> success" << std::endl;

  // set non-uniform values to check
  double fct = 2.;
  for (const auto& d : nf1.getComponents()) {
    for (int i=0; i<nf1.getNbNodes(); ++i) {
      double vv = i*fct + d*100;
      nf1(i,d) = vv;
      std::cout << "nfc(" << i << "," << d << ")="
		<< vv << std::endl;
    }
  }

  std::cout << "check accessing data" << std::endl;
  for (const auto& d : nf1.getComponents()) {
    for (int i=0; i<nf1.getNbNodes(); ++i) {
      if (nf1(i,d) != i*fct + d*100) {
	std::cerr << "() operator failed" << std::endl;
	return 1; // failure
      }
    }
  }

  for (const auto& d : nf1.getComponents()) {
    for (int i=0; i<nf1.getNbNodes(); ++i) {
      if (nf1.data(d)[i] != i*fct + d*100) {
	std::cerr << "data access failed" << std::endl;
	return 1; // failure
      }
    }
  }
  std::cout << "data accessing correct -> success" << std::endl;


  
  // Check NodalField::computeNorm
  // --------------------------------------------------------------
  std::cout << "check NodalField::computeNorm" << std::endl;
  
  double Lx = 0.4;
  int Nx = 3;
  double Lz = 0.3;
  int Nz = 2;
  FFTableMesh mesh3d(Lx, Nx, Lz, Nz);
  
  NodalField ans(mesh3d);
  NodalField field(mesh3d, {_x,_y,_z});
  field.setAllValuesTo(1,0);
  field.setAllValuesTo(2,1);
  field.setAllValuesTo(3,2);
  bool checks_out = true;
  double tol = 1e-16;
  double ref_ans = std::sqrt(1 * 1 + 2 * 2 + 3 * 3);
  field.computeNorm(ans);
  for (int i = 0; i < mesh3d.getNbLocalNodes(); ++i) {
    if (std::abs(ans(i) - ref_ans) > tol) {
      checks_out = false;
      break;
    }
  }
  ref_ans = std::sqrt(1 * 1 + 3 * 3);
  field.computeNorm(ans, 1);
  for (int i = 0; i < mesh3d.getNbLocalNodes(); ++i) {
    if (std::abs(ans(i) - ref_ans) > tol) {
      checks_out = false;
      break;
    }
  }
  if (checks_out) {
    std::cout << "NodalField::computeNorm correct -> success" << std::endl;
  } else {
    std::cout << "NodalField::computeNorm failed" << std::endl;
    return 1;  // failed
  }
  // --------------------------------------------------------------

  // Check Interface::multiplyFieldByScalar
  // --------------------------------------------------------------
  std::cout << "check NodalField::multiplyByScalar" << std::endl;
  checks_out = true;
  tol = 1e-16;
  ref_ans = 0;
  ans.setAllValuesTo(0);
  field.multiplyByScalarField(ans);
  for (const auto& d : nf1.getComponents()) {
    for (int i = 0; i < mesh3d.getNbLocalNodes(); ++i) {
      if (std::abs(field(i,d) - ref_ans) > tol) {
        checks_out = false;
        break;
      }
    }
  }
  if (checks_out) {
    std::cout << "NodalField::multiplyByScalar correct -> success"
              << std::endl;
  } else {
    std::cout << "NodalField::multiplyByScalar failed" << std::endl;
    return 1;  // failed
  }

  
  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
