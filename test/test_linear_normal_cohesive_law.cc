/**
 * @file   test_linear_normal_cohesive_law.cc
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
#include "linear_normal_cohesive_law.hh"
#include "defrig_interface.hh"

#include <iostream>
#include <cmath>

using namespace uguca;

int main(){

  std::cout << "start test: test_liner_normal_cohesive_law" << std::endl;

  // information for checks
  SimpleMesh mesh(1., 2);
  double Gc=2.;
  double sigmac=4e6;
  double dc = 2*Gc/sigmac;
  Material material(1e9, 0.25, 1000, false);

  std::cout << "check constructor" << std::endl;
  LinearNormalCohesiveLaw law(mesh, Gc, sigmac);
  DefRigInterface interface(mesh, {_x,_y}, material, law);
  interface.setTimeStep(1e-8);
  std::cout << "constructor correct -> success" << std::endl;

  std::cout << "check data" << std::endl;
  NodalField & Gc_tmp = law.getGc();
  if (std::abs(Gc_tmp(0) - Gc) / Gc > 1e-5){
    std::cout << "wrong Gc: " << Gc_tmp(0) << std::endl;
    return 1; // failure
  }
  NodalField & Sigmac_tmp = law.getSigmac();
  if (std::abs(Sigmac_tmp(0) - sigmac) / sigmac > 1e-5){
    std::cout << "wrong sigmac: " << Sigmac_tmp(0) << std::endl;
    return 1; // failure
  }
  std::cout << "data correct -> success" << std::endl;

  
  std::cout << "check computeCohesiveForces" << std::endl;

  // fill empty cohesion vector for testing
  //NodalField cohesion(mesh, {_x,_y});

  // access to various properties needed to apply values
  NodalField & load = interface.getLoad();
  HalfSpace & top = interface.getTop();
  NodalField & u = top.getDisp();

  // check: sigma0 < sigmac & u=0
  double sigma0v = 0.9*sigmac;
  load.setAllValuesTo(sigma0v,1);
  law.computeCohesiveForces(false);
  NodalField & cohesion =  law.getCohesion();

  if ((std::abs(cohesion(0,1) - sigma0v) / sigma0v > 1e-5)
      || (cohesion(0,1) * sigma0v < 0)) {
    std::cout << "cohesion failed (" << sigma0v << "): "
	      << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  // check: sigma0 > sigmac & u=0
  sigma0v = 1.1*sigmac;
  load.setAllValuesTo(sigma0v,1);
  law.computeCohesiveForces(false);
  
  if ((std::abs(cohesion(0,1) - sigmac) / sigmac > 1e-5)
      || (cohesion(0,1) * sigma0v < 0)) {
    std::cout << "onset failed (" << sigmac << "): "
	      << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  // check during weakening: sigma0 > sigmac & 0 < u < dc
  sigma0v = 1.1*sigmac;
  load.setAllValuesTo(sigma0v,1);
  double u1v = 0.7*dc;
  u.setAllValuesTo(u1v,1);
  double val = sigmac - u1v/dc*sigmac;
  law.computeCohesiveForces(false);
  if ((std::abs(cohesion(0,1) - val) / val > 1e-5)
      || (cohesion(0,1) * sigma0v < 0)) {
    std::cout << "weakening failed (" << val << "): "
	      << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  // check: sigma0 > sigmac & dc < u
  sigma0v = 1.1*sigmac;
  load.setAllValuesTo(sigma0v,1);
  u1v = 1.1*dc;
  u.setAllValuesTo(u1v,1);
  val = 0;
  law.computeCohesiveForces(false);
  if ((std::abs(cohesion(0,1) - val) / val > 1e-5)
      || (cohesion(0,1) * sigma0v < 0)) {
    std::cout << "residual failed (" << val << "): "
	      << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  // check no penetration: sig0 < 0
  double sig0v = -2*sigmac;
  load.setAllValuesTo(sig0v,1);
  law.computeCohesiveForces(false);
  if ((std::abs(cohesion(0,1) - sig0v) / sig0v > 1e-5)
      || (cohesion(0,1) * sig0v < 0)) {
    std::cout << "no penetration failed (" << sig0v << "): "
	      << cohesion(0,1) << std::endl;
    return 1; // failure
  }

  std::cout << "computeCohesiveForces -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
