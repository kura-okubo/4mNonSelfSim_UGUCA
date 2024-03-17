/**
 * @file   test_barras_law.cc
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

#include "material.hh"
#include "uca_simple_mesh.hh"
#include "barras_law.hh"
#include "defrig_interface.hh"

#include <iostream>
#include <cmath>

using namespace uguca;

int main(){

  std::cout << "start test: test_barras_law" << std::endl;

  // information for checks
  SimpleMesh mesh(1., 2);
  double tauc = 4e6;
  double dc = 0.01;
  Material material(1e9, 0.25, 1000, false);

  std::cout << "check constructor" << std::endl;
  BarrasLaw law(mesh, tauc, dc);
  DefRigInterface interface(mesh, {_x,_y}, material, law);
  interface.setTimeStep(1e-8);
  std::cout << "constructor correct -> success" << std::endl;

  std::cout << "check data" << std::endl;
  NodalField & TauMax_tmp = law.getTauMax();
  if (std::abs(TauMax_tmp(0) - tauc) / tauc > 1e-5){
    std::cout << "wrong tau_max: " << TauMax_tmp(0) << std::endl;
    return 1; // failure
  }
  NodalField & Dc_tmp = law.getDc();
  if (std::abs(Dc_tmp(0) - dc) / dc > 1e-5){
    std::cout << "wrong d_c: " << Dc_tmp(0) << std::endl;
    return 1; // failure
  }
  std::cout << "data correct -> success" << std::endl;

  std::cout << "check computeCohesiveForces" << std::endl;

  // fill empty cohesion vector for testing
  NodalField cohesion(mesh, {_x,_y});

  // access to various properties needed to apply values
  NodalField & load = interface.getLoad();
  NodalField & disp = interface.getTop().getDisp();

  // check shear: tau0 < tauc & u=0
  double tau0v = 0.9*tauc;
  double sig0v = 0;
  load.setAllValuesTo(tau0v,0);
  load.setAllValuesTo(sig0v,1);
  double u0v = 0.;
  double u1v = 0.;
  disp.setAllValuesTo(u0v,0);
  disp.setAllValuesTo(u1v,1);
  double val = tau0v;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(cohesion(0,0) - val) / val > 1e-5) || (cohesion(0,0) * tau0v < 0)) {
    std::cout << "shear failed (" << val << "): " << cohesion(0,0) << std::endl;
    return 1; // failure
  }

  // check normal: 0 < sig0 < tauc & u=0
  tau0v = 0.;
  sig0v = 0.9*tauc;
  load.setAllValuesTo(tau0v,0);
  load.setAllValuesTo(sig0v,1);
  u0v = 0.;
  u1v = 0.;
  disp.setAllValuesTo(u0v,0);
  disp.setAllValuesTo(u1v,1);
  val = sig0v;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(cohesion(0,1) - val) / val > 1e-5) || (cohesion(0,1) * sig0v < 0)) {
    std::cout << "normal failed (" << val << "): " << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  // check contact: sig0 < 0 & u=0
  tau0v = 0.;
  sig0v = -1.5*tauc;
  load.setAllValuesTo(tau0v,0);
  load.setAllValuesTo(sig0v,1);
  u0v = 0.;
  u1v = 0.;
  disp.setAllValuesTo(u0v,0);
  disp.setAllValuesTo(u1v,1);
  val = sig0v;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(cohesion(0,1) - val) / val > 1e-5) || (cohesion(0,1) * sig0v < 0)) {
    std::cout << "contact failed (" << val << "): " << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  // check decohesion: tau0 > tauc & sig0 > tauc & u=0
  tau0v = 1.2*tauc;
  sig0v = 1.5*tauc;
  load.setAllValuesTo(tau0v,0);
  load.setAllValuesTo(sig0v,1);
  u0v = 0.;
  u1v = 0.;
  disp.setAllValuesTo(u0v,0);
  disp.setAllValuesTo(u1v,1);
  val = tauc;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(cohesion(0,0) - val) / val > 1e-5) || (cohesion(0,0) * tau0v < 0) ||
      (std::abs(cohesion(0,1) - val) / val > 1e-5) || (cohesion(0,1) * sig0v < 0)) {
    std::cout << "decohesion failed (" << val << "): "
	      << cohesion(0,0) << "," << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  // check weakening: tau0 > tauc & sig0 > tauc & u < dc
  tau0v = 1.2*tauc;
  sig0v = 1.5*tauc;
  load.setAllValuesTo(tau0v,0);
  load.setAllValuesTo(sig0v,1);
  u0v = 0.1*dc;
  u1v = 0.1*dc;
  disp.setAllValuesTo(u0v,0);
  disp.setAllValuesTo(u1v,1);
  val = tauc * (1 - std::sqrt(u0v*u0v+u1v*u1v)/dc);
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(cohesion(0,0) - val) / val > 1e-5) || (cohesion(0,0) * tau0v < 0) ||
      (std::abs(cohesion(0,1) - val) / val > 1e-5) || (cohesion(0,1) * sig0v < 0)) {
    std::cout << "weakening failed (" << val << "): "
	      << cohesion(0,0) << "," << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  // check open: tau0 > tauc & sig0 > tauc & u > dc
  tau0v = 1.2*tauc;
  sig0v = 1.5*tauc;
  load.setAllValuesTo(tau0v,0);
  load.setAllValuesTo(sig0v,1);
  u0v = 0.9*dc;
  u1v = 0.9*dc;
  disp.setAllValuesTo(u0v,0);
  disp.setAllValuesTo(u1v,1);
  val = 0.;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(cohesion(0,0) - val) / val > 1e-5) || (cohesion(0,0) * tau0v < 0) ||
      (std::abs(cohesion(0,1) - val) / val > 1e-5) || (cohesion(0,1) * sig0v < 0)) {
    std::cout << "open failed (" << val << "): "
	      << cohesion(0,0) << "," << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  std::cout << "computeCohesiveForces -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
