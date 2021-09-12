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
		     double delta_c_default) :
  InterfaceLaw(mesh),
  tau_max(mesh),
  delta_c(mesh),
  gap_norm(mesh)
{
  this->tau_max.setAllValuesTo(tau_max_default);
  this->delta_c.setAllValuesTo(delta_c_default);
  this->gap_norm.zeros();
}

/* -------------------------------------------------------------------------- */
void BarrasLaw::computeCohesiveForces(NodalField & cohesion,
                                      bool predicting) {

  // find current gap
  //std::vector<NodalField *> gap = this->interface->getBufferField();
  NodalField gap(this->mesh);
  this->interface->computeGap(gap, predicting);
  gap.computeNorm(this->gap_norm);
  double *gap_norm_p = this->gap_norm.storage();

  // find forces needed to close normal gap
  NodalFieldComponent & coh1 = cohesion.component(1);
  this->interface->closingNormalGapForce(coh1, predicting);
  double * p_coh1 = coh1.storage();

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion);

  // get norm of shear cohesion
  NodalFieldComponent shear_trac_norm(this->mesh);
  cohesion.computeNorm(shear_trac_norm, 1);
  double * tau_shear = shear_trac_norm.storage();

  // interface properties
  double *tauc = this->tau_max.storage();
  double *dc = this->delta_c.storage();

  // to be filled
  NodalFieldComponent alpha_field(this->mesh);
  double * alpha = alpha_field.storage();

  // coh1 > 0 is a adhesive force
  // coh1 < 0 is a contact pressure
  for (int n = 0; n < this->mesh.getNbLocalNodes(); ++n) {

    double rel_gap = gap_norm_p[n] / dc[n];
    double strength = tauc[n] * std::max(0., 1 - rel_gap);

    p_coh1[n] = std::min(p_coh1[n], strength);

    // maximal shear cohesive force given by strength.
    // keep orientation of shear force
    alpha[n] = std::min(1.,std::abs(strength / tau_shear[n]));
  }

  // only in shear direction
  cohesion.multiplyByScalar(alpha_field, 1);
}

/* --------------------------------------------------------------------------*/
void BarrasLaw::registerDumpField(const std::string & field_name) {

  // tau_max
  if (field_name == "tau_max") {
    this->interface->registerForDump(field_name, this->tau_max);
  }

  // delta_c
  else if (field_name == "delta_c") {
    this->interface->registerForDump(field_name, this->delta_c);
  }

  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }
}

__END_UGUCA__
