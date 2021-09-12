/**
 * @file   test_linear_coulomb_friction_law.cc
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
#include "linear_coulomb_friction_law.hh"
#include "defrig_interface.hh"

#include <iostream>

using namespace uguca;

int main(){

  std::cout << "start test: test_liner_coulomb_friction_law" << std::endl;

  // information for checks
  SimpleMesh mesh(1., 2);
  double mus=0.8;
  double muk=0.3;
  double dc=2e-6;
  double T=1e-7;
  double sig=-1e6;
  double dt=1e-8;
  Material material(1e9, 0.25, 1000, false);

  std::cout << "check constructor" << std::endl;
  LinearCoulombFrictionLaw law(mesh, mus, muk, dc, T);
  DefRigInterface interface(mesh, material, law);
  interface.setTimeStep(dt);
  std::cout << "constructor correct -> success" << std::endl;

  std::cout << "check data" << std::endl;
  NodalFieldComponent & MuS_tmp = law.getMuS();
  if (std::abs(MuS_tmp.at(0) - mus) / mus > 1e-5){
    std::cout << "wrong mus: " << MuS_tmp.at(0) << std::endl;
    return 1; // failure
  }
  NodalFieldComponent & MuK_tmp = law.getMuK();
  if (std::abs(MuK_tmp.at(0) - muk) / muk > 1e-5){
    std::cout << "wrong muk: " << MuK_tmp.at(0) << std::endl;
    return 1; // failure
  }
  NodalFieldComponent & Dc_tmp = law.getDc();
  if (std::abs(Dc_tmp.at(0) - dc) / dc > 1e-5){
    std::cout << "wrong dc: " << Dc_tmp.at(0) << std::endl;
    return 1; // failure
  }
  NodalFieldComponent & T_tmp = law.getCharacteristicTime();
  if (std::abs(T_tmp.at(0) - T) / T > 1e-5){
    std::cout << "wrong T: " << T_tmp.at(0) << std::endl;
    return 1; // failure
  }
  std::cout << "data correct -> success" << std::endl;

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

  sig0.setAllValuesTo(sig);

  std::cout << "check computeRegContactPressure" << std::endl;
  NodalFieldComponent rcp(mesh);
  double prev_rcp = 0.5e6;
  rcp.setAllValuesTo(prev_rcp);
  coh1.setAllValuesTo(std::abs(sig));
  law.computeRegContactPressure(coh1, rcp);

  double val = (prev_rcp + dt/T * std::abs(sig)) / (1 + dt/T);
  if (std::abs(rcp.at(0) - val) / val > 1e-5) {
    std::cout << "regularized contact pressure failed ("
	      << val << "): " << rcp.at(0) << std::endl;
    return 1; // failure
  }
  std::cout << "computeRegContactPressure correct -> success" << std::endl;

  std::cout << "check computeCohesiveForces" << std::endl;

  // check stick: tau0 < mus*sig & u=0
  double tau0v = 0.9*mus*sig;
  tau0.setAllValuesTo(tau0v);
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - tau0v) / tau0v > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "stick failed (" << tau0v << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check stick-to-slip: tau0 > tauc & u=0
  tau0v = 1.1*mus*sig;
  tau0.setAllValuesTo(tau0v);
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - mus*sig) / (mus*sig) > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "stick-to-slip failed (" << (mus*sig) << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check slip during weakening: tau0 > tauc & 0 < u < dc
  tau0v = 1.1*mus*sig;
  tau0.setAllValuesTo(tau0v);
  double u0v = 0.7*dc;
  u0.setAllValuesTo(u0v);
  val = (mus - u0v/dc*(mus - muk)) * sig;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - val) / val > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "weakening failed (" << val << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check residual: tau0 > tauc & dc < u
  tau0v = 1.1*mus*sig;
  tau0.setAllValuesTo(tau0v);
  u0v = 1.1*dc;
  u0.setAllValuesTo(u0v);
  val = muk*sig;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - val) / val > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "residual failed (" << val << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check negative direction: tau0 > tauc & dc < u
  tau0v = -1.1*mus*sig;
  tau0.setAllValuesTo(tau0v);
  u0v = -1.1*dc;
  u0.setAllValuesTo(u0v);
  val = -muk*sig;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0.at(0) - val) / val > 1e-5) || (coh0.at(0) * tau0v < 0)) {
    std::cout << "negative failed (" << val << "): " << coh0.at(0) << std::endl;
    return 1; // failure
  }

  // check no penetration: sig0 < 0
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh1.at(0) - sig) / sig > 1e-5) || (coh1.at(0) * sig < 0)) {
    std::cout << "no penetration failed (" << sig << "): " << coh1.at(0) << std::endl;
    return 1; // failure
  }

  // check no adhesion: sig0 > 0
  double sig0v = -sig;
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
