/**
 * @file   test_linear_shear_cohesive_law.cc
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
#include "material.hh"
#include "linear_shear_cohesive_law.hh"
#include "defrig_interface.hh"

#include <iostream>
#include <cmath>

using namespace uguca;

int main(){

  std::cout << "start test: test_liner_shear_cohesive_law" << std::endl;

  // information for checks
  SimpleMesh mesh(1., 2);
  double Gc=2.;
  double tauc=4e6;
  double taur=1e6;
  double dc = 2*Gc/(tauc-taur);
  Material material(1e9, 0.25, 1000, false);

  std::cout << "check constructor" << std::endl;
  LinearShearCohesiveLaw law(mesh, Gc, tauc, taur);
  DefRigInterface interface(mesh, material, law);
  interface.setTimeStep(1e-8);
  std::cout << "constructor correct -> success" << std::endl;

  std::cout << "check data" << std::endl;
  NodalFieldComponent & Gc_tmp = law.getGc();
  if (std::abs(Gc_tmp.at(0) - Gc) / Gc > 1e-5){
    std::cout << "wrong Gc: " << Gc_tmp.at(0) << std::endl;
    return 1; // failure
  }
  NodalFieldComponent & Tauc_tmp = law.getTauc();
  if (std::abs(Tauc_tmp.at(0) - tauc) / tauc > 1e-5){
    std::cout << "wrong tauc: " << Tauc_tmp.at(0) << std::endl;
    return 1; // failure
  }
  NodalFieldComponent & Taur_tmp = law.getTaur();
  if (std::abs(Taur_tmp.at(0) - taur) / taur > 1e-5){
    std::cout << "wrong taur: " << Taur_tmp.at(0) << std::endl;
    return 1; // failure
  }
  std::cout << "data correct -> success" << std::endl;

  
  std::cout << "check computeCohesiveForces" << std::endl;

  // fill empty cohesion vector for testing
  NodalField cohesion(mesh);
  NodalFieldComponent & coh0 = cohesion.component(0);
  NodalFieldComponent & coh1 = cohesion.component(1);

  // access to various properties needed to apply values
  NodalField & load = interface.getLoad();
  NodalFieldComponent & tau0 = load.component(0);
  NodalFieldComponent & sig0 = load.component(1);
  HalfSpace & top = interface.getTop();
  NodalFieldComponent & u0 = top.getDisp().component(0);

  // check stick: tau0 < tauc & u=0
  double tau0v = 0.9*tauc;
  tau0.setAllValuesTo(tau0v);
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - tau0v) / tau0v > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "stick failed (" << tau0v << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check stick-to-slip: tau0 > tauc & u=0
  tau0v = 1.1*tauc;
  tau0.setAllValuesTo(tau0v);
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - tauc) / tauc > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "stick-to-slip failed (" << tauc << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check slip during weakening: tau0 > tauc & 0 < u < dc
  tau0v = 1.1*tauc;
  tau0.setAllValuesTo(tau0v);
  double u0v = 0.7*dc;
  u0.setAllValuesTo(u0v);
  double val = tauc - u0v/dc*(tauc - taur);
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - val) / val > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "weakening failed (" << val << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check residual: tau0 > tauc & dc < u
  tau0v = 1.1*tauc;
  tau0.setAllValuesTo(tau0v);
  u0v = 1.1*dc;
  u0.setAllValuesTo(u0v);
  val = taur;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - val) / val > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "residual failed (" << val << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check negative direction: tau0 > tauc & dc < u
  tau0v = -1.1*tauc;
  tau0.setAllValuesTo(tau0v);
  u0v = -1.1*dc;
  u0.setAllValuesTo(u0v);
  val = -taur;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - val) / val > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "negative failed (" << val << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check no penetration: sig0 < 0
  double sig0v = -2*tauc;
  sig0.setAllValuesTo(sig0v);
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh1.at(0) - sig0v) / sig0v > 1e-5) || (coh1.at(0) * sig0v < 0)) {
    std::cout << "no penetration failed (" << sig0v << "): " << coh1.at(0) << std::endl;
    return 1; // failure
  }

  // check no adhesion: sig0 > 0
  sig0v = 2*tauc;
  sig0.setAllValuesTo(sig0v);
  law.computeCohesiveForces(cohesion, false);
  if (std::abs(coh1.at(0)) > 1e-8) {
    std::cout << "no adhesion failed (" << 0. << "): " << coh1.at(0) << std::endl;
    return 1; // failure
  }

  std::cout << "computeCohesiveForces -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
