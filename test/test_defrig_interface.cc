/**
 * @file   test_defrig_interface.cc
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

#include "defrig_interface.hh"
#include "fftable_nodal_field.hh"
#include "linear_shear_cohesive_law.hh"
#include "material.hh"
#include "nodal_field.hh"
#include "uca_simple_mesh.hh"

#include <iostream>
#include <random>
#include <vector>

using namespace uguca;

int checkClosingNormalGapForce(SimpleMesh &mesh, DefRigInterface &interface,
                               NodalField &fields) {
  // assign random u_1 and load_1
  std::vector<double> u_1;
  u_1.reserve((size_t)mesh.getNbLocalNodes());
  std::vector<double> load_1;
  load_1.reserve((size_t)mesh.getNbLocalNodes());
  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;
  HalfSpace &top = interface.getTop();
  NodalFieldComponent & disp_1_top = top.getDisp(false).component(1);
  for (size_t i = 0; i < (size_t)mesh.getNbLocalNodes(); ++i) {
    double random_value = unif(re);
    u_1.push_back(random_value);
    disp_1_top.set(i) = random_value;
    random_value = unif(re);
    load_1.push_back(random_value);
    (interface.getLoad().component(1)).set(i) = random_value;
  }
  // compute reference closing forces, note that internal is zero.
  std::vector<double> ref_close_force;
  ref_close_force.reserve((size_t)mesh.getNbLocalNodes());
  const Material &mat_top = top.getMaterial();
  double fact_top = interface.getTimeStep() * mat_top.getCs() *
                    mat_top.getCs() / mat_top.getCp() /
                    mat_top.getShearModulus();
  for (size_t i = 0; i < (size_t)mesh.getNbLocalNodes(); ++i) {
    ref_close_force.push_back(u_1[i] / fact_top + load_1[i] + 0.0);
  }
  // check results
  double tol = 1e-10;
  interface.closingNormalGapForce(fields.component(1), false);  // not predicting
  for (size_t i = 0; i < (size_t)mesh.getNbLocalNodes(); ++i) {
    if (std::abs((fields.component(1)).set(i)-ref_close_force[i]) > tol) {
      std::cout << "discrepancy found in closingNormalGapForce" << std::endl;
      return 1;  // failure
    }
  }
  return 0;  // success
}

int checkMaintainShearGapForce(SimpleMesh &mesh, DefRigInterface &interface,
                               NodalField &fields) {
  // assign random load_0 and load_2
  std::vector<double> load_0;
  load_0.reserve((size_t)mesh.getNbLocalNodes());
  std::vector<double> load_2;
  load_2.reserve((size_t)mesh.getNbLocalNodes());
  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;
  for (size_t i = 0; i < (size_t)mesh.getNbLocalNodes(); ++i) {
    double random_value = unif(re);
    load_0.push_back(random_value);
    (interface.getLoad().component(0)).set(i) = random_value;
    random_value = unif(re);
    load_2.push_back(random_value);
    (interface.getLoad().component(2)).set(i) = random_value;
  }
  // compute reference solutions, note that internal is zero.
  std::vector<double> ref_maintain_force_0;
  ref_maintain_force_0.reserve((size_t)mesh.getNbLocalNodes());
  std::vector<double> ref_maintain_force_2;
  ref_maintain_force_2.reserve((size_t)mesh.getNbLocalNodes());
  for (size_t i = 0; i < (size_t)mesh.getNbLocalNodes(); ++i) {
    ref_maintain_force_0.push_back(load_0[i] + 0.0);
    ref_maintain_force_2.push_back(load_2[i] + 0.0);
  }
  // check results
  interface.maintainShearGapForce(fields);  // not predicting
  double tol = 1e-10;
  for (size_t i = 0; i < (size_t)mesh.getNbLocalNodes(); ++i) {
    if (std::abs((fields.component(0)).set(i)-ref_maintain_force_0[i]) > tol) {
      std::cout << "discrepancy found in maintainShearGapForce" << std::endl;
      return 1;  // failure
    }
    if (std::abs((fields.component(2)).set(i)-ref_maintain_force_2[i]) > tol) {
      std::cout << "discrepancy found in maintainShearGapForce" << std::endl;
      return 1;  // failure
    }
  }
  return 0;  // success
}

int checkComputeGaps(SimpleMesh &mesh, DefRigInterface &interface,
                     NodalField &fields) {
  // assign random u
  HalfSpace &top = interface.getTop();
  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;
  std::vector<std::vector<double>> disp(3);
  std::vector<std::vector<double>> velo(3);
  for (size_t j = 0; j < 3; ++j) {
    disp[j].reserve((size_t)mesh.getNbLocalNodes());
    velo[j].reserve((size_t)mesh.getNbLocalNodes());
    for (size_t i = 0; i < (size_t)mesh.getNbLocalNodes(); ++i) {
      double random_value = unif(re);
      disp[j].push_back(random_value);
      top.getDisp().component(j).set(i) = random_value;
      random_value = unif(re);
      velo[j].push_back(random_value);
      top.getVelo().component(j).set(i) = random_value;
    }
  }
  // check results
  interface.computeGap(fields, false);  // not predicting
  double tol = 1e-10;
  for (size_t j = 0; j < 3; ++j) {
    for (size_t i = 0; i < (size_t)mesh.getNbLocalNodes(); ++i) {
      if (std::abs(fields.component(j).set(i)-disp[j][i]) > tol) {
        std::cout << "discrepancy found in computeGap" << std::endl;
        return 1;  // failure
      }
    }
  }
  interface.computeGapVelocity(fields, false);  // not predicting
  for (size_t j = 0; j < 3; ++j) {
    for (size_t i = 0; i < (size_t)mesh.getNbLocalNodes(); ++i) {
      if (std::abs(fields.component(j).set(i)-velo[j][i]) > tol) {
        std::cout << "discrepancy found in computeGapVelocity" << std::endl;
        return 1;  // failure
      }
    }
  }
  return 0;  // success
}

int main(){

  std::cout << "start test: test_defrig_interface" << std::endl;

  // --------------------------------------------------------------
  // INITIALIZE
  // --------------------------------------------------------------
  // mesh
  // --------------------------------------------------------------
  double Lx = 0.4;
  int Nx = 3;
  double Lz = 0.3;
  int Nz = 2;
  SimpleMesh mesh(Lx, Nx, Lz, Nz);
  // --------------------------------------------------------------
  // material
  // --------------------------------------------------------------
  double E = 2.1;
  double nu = 0.25;
  double rho = 10.0;
  Material material(E, nu, rho, false);
  // --------------------------------------------------------------
  // interface law
  // --------------------------------------------------------------
  double Gc = 1.0;
  double tau_c = 3.0;
  double tau_r = 1.5;
  LinearShearCohesiveLaw law(mesh, Gc, tau_c, tau_r);
  // --------------------------------------------------------------
  // interface
  // --------------------------------------------------------------
  double time_step = 1.0e-1;
  DefRigInterface interface(mesh, material, law);
  interface.setTimeStep(time_step);
  // --------------------------------------------------------------
  // container
  // --------------------------------------------------------------
  NodalField fields(mesh);
  // --------------------------------------------------------------

  // Check DefRigInterface::closingNormalGapForce
  // --------------------------------------------------------------
  std::cout << "check DefRigInterface::closingNormalGapForce" << std::endl;
  if (checkClosingNormalGapForce(mesh, interface,fields)) {
    std::cout << "DefRigInterface::closingNormalGapForce failed" << std::endl;
    return 1;  // failed
  } else {
    std::cout << "DefRigInterface::closingNormalGapForce correct -> success" << std::endl;
  }
  // --------------------------------------------------------------

  // Check DefRigInterface::maintainShearGapForce
  // --------------------------------------------------------------
  std::cout << "check DefRigInterface::maintainShearGapForce" << std::endl;
  if (checkMaintainShearGapForce(mesh, interface, fields)) {
    std::cout << "DefRigInterface::maintainShearGapForce failed" << std::endl;
    return 1;  // failed
  } else {
    std::cout << "DefRigInterface::maintainShearGapForce correct -> success" << std::endl;
  }
  // --------------------------------------------------------------

  // Check DefRigInterface::computeGap*
  // --------------------------------------------------------------
  std::cout << "check DefRigInterface::computeGap*" << std::endl;
  if (checkComputeGaps(mesh, interface, fields)) {
    std::cout << "DefRigInterface::computeGap* failed" << std::endl;
    return 1;  // failed
  } else {
    std::cout << "DefRigInterface::computeGap* correct -> success" << std::endl;
  }
  // --------------------------------------------------------------


  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
