/**
 * @file   linear_shear_cohesive_law.cc
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
#include "linear_shear_cohesive_law.hh"
#include "interface.hh"

#include <cmath>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

LinearShearCohesiveLaw::LinearShearCohesiveLaw(BaseMesh & mesh,
					       double Gc_default,
					       double tau_c_default,
					       double tau_r_default,
					       const std::string & name) :
  InterfaceLaw(mesh,name),
  G_c(mesh),
  tau_c(mesh),
  tau_r(mesh)
{
  this->G_c.setAllValuesTo(Gc_default);
  this->G_c.setName(name+"_G_c");
  
  this->tau_c.setAllValuesTo(tau_c_default);
  this->tau_c.setName(name+"_tau_c");
  
  this->tau_r.setAllValuesTo(tau_r_default);
  this->tau_r.setName(name+"_tau_r");
}

/* -------------------------------------------------------------------------- */
void LinearShearCohesiveLaw::computeCohesiveForces(NodalField & cohesion,
                                                   bool predicting) {

  // find forces needed to close normal gap
  this->interface->closingNormalGapForce(cohesion, predicting);

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion);

  // get norm of shear cohesion
  NodalField shear_trac_norm(this->mesh);
  cohesion.computeNorm(shear_trac_norm, 1);

  // find current gap
  NodalField gap(this->mesh, cohesion.getComponents());
  this->interface->computeGap(gap, predicting);

  // compute norm of shear gap
  NodalField shear_gap_norm(this->mesh);
  gap.computeNorm(shear_gap_norm, 1);

  // to be filled
  NodalField alpha(this->mesh);

  // coh1 > 0 is a adhesive force
  // coh1 < 0 is a contact pressure
  for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n) {

    // avoid penetration "at any cost"
    // apply no normal cohesive force
    cohesion(n,1) = std::min(cohesion(n,1), 0.);

    double slope = pow(tau_c(n) - tau_r(n),2) / 2. / G_c(n);
    double strength = std::max(tau_c(n) - shear_gap_norm(n) * slope,
			       tau_r(n));

    // maximal shear cohesive force given by strength.
    // keep orientation of shear force
    alpha(n) = std::min(1.,std::abs(strength / shear_trac_norm(n)));
  }

  // only in shear direction
  for (const auto& d : cohesion.getComponents()) {
    if (d==1) // ignore normal direction
      continue;
    cohesion.multiplyByScalarField(alpha,d);
  }
}

/* -------------------------------------------------------------------------- */
void LinearShearCohesiveLaw::registerDumpField(const std::string & field_name) {

  // G_c
  if (field_name == "G_c") {
    this->interface->registerIO(field_name,
				     this->G_c);
  }

  // tau_c
  else if (field_name == "tau_c") {
    this->interface->registerIO(field_name,
				     this->tau_c);
  }

  // tau_r
  else if (field_name == "tau_r") {
    this->interface->registerIO(field_name,
				     this->tau_r);
  }

  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }

}

/* -------------------------------------------------------------------------- */
void LinearShearCohesiveLaw::registerToRestart(Restart & restart) {

  restart.registerIO(this->G_c);
  restart.registerIO(this->tau_c);
  restart.registerIO(this->tau_r);
  
  InterfaceLaw::registerToRestart(restart);
}

__END_UGUCA__
