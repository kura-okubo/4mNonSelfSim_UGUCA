/**
 * @file   barras_law.cc
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
#include "barras_law.hh"
#include "interface.hh"
#include <iostream>
#include <cmath>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

BarrasLaw::BarrasLaw(BaseMesh & mesh,
		     double tau_max_default,
		     double delta_c_default,
		     const std::string & name) :
  InterfaceLaw(mesh,name),
  tau_max(mesh),
  delta_c(mesh),
  gap_norm(mesh)
{
  this->tau_max.setAllValuesTo(tau_max_default);
  this->tau_max.setName(name+"_tau_max");
  
  this->delta_c.setAllValuesTo(delta_c_default);
  this->delta_c.setName(name+"_delta_c");
  
  this->gap_norm.zeros();
  this->gap_norm.setName(name+"_gap_norm");
}

/* -------------------------------------------------------------------------- */
void BarrasLaw::computeCohesiveForces(NodalField & cohesion,
                                      bool predicting) {

  // find current gap
  NodalField gap(this->mesh, cohesion.getComponents());
  this->interface->computeGap(gap, predicting);
  gap.computeNorm(this->gap_norm);

  // find forces needed to close normal gap
  this->interface->closingNormalGapForce(cohesion, predicting);

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion);

  // get norm of shear cohesion
  NodalField shear_trac_norm(this->mesh);
  cohesion.computeNorm(shear_trac_norm, 1);

  // to be filled
  NodalField alpha(this->mesh);

  // coh1 > 0 is a adhesive force
  // coh1 < 0 is a contact pressure
  for (int n = 0; n < this->mesh.getNbLocalNodes(); ++n) {

    double rel_gap = this->gap_norm(n) / delta_c(n);
    double strength = tau_max(n) * std::max(0., 1 - rel_gap);

    cohesion(n,1) = std::min(cohesion(n,1), strength);

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

/* --------------------------------------------------------------------------*/
void BarrasLaw::registerDumpField(const std::string & field_name) {

  // tau_max
  if (field_name == "tau_max") {
    this->interface->registerIO(field_name, this->tau_max);
  }

  // delta_c
  else if (field_name == "delta_c") {
    this->interface->registerIO(field_name, this->delta_c);
  }

  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }
}


/* -------------------------------------------------------------------------- */
void BarrasLaw::registerToRestart(Restart & restart) {

  restart.registerIO(this->tau_max);
  restart.registerIO(this->delta_c);
  restart.registerIO(this->gap_norm);
  
  InterfaceLaw::registerToRestart(restart);
}

__END_UGUCA__
