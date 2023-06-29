/**
 * @file   continuous_shear_weakening_law.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Thu Jun 29 2023
 * @date last modification: Thu Jun 29 2023
 *
 * @brief  TODO
 *
 *
 * Copyright (C) 2023 ETH Zurich (David S. Kammer)
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
#include "continuous_shear_weakening_law.hh"
#include "interface.hh"

#include <cmath>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

ContinuousShearWeakeningLaw::ContinuousShearWeakeningLaw(BaseMesh & mesh,
							 double tau_s_default,
							 double tau_d_default,
							 double d_c_default,
							 double alpha_default,
							 const std::string & name) :
  InterfaceLaw(mesh,name),
  tau_s(mesh,name+"_tau_s"),
  tau_d(mesh,name+"_tau_d"),
  d_c(mesh,name+"_d_c"),
  alpha(mesh,name+"_alpha")
{
  this->tau_s.setAllValuesTo(tau_s_default);
  this->tau_d.setAllValuesTo(tau_d_default);
  this->d_c.setAllValuesTo(d_c_default);
  this->alpha.setAllValuesTo(alpha_default);
}

/* -------------------------------------------------------------------------- */
void ContinuousShearWeakeningLaw::computeCohesiveForces(NodalField & cohesion,
							bool predicting) {

  // find forces needed to close normal gap
  NodalFieldComponent & coh1 = cohesion.component(1);
  this->interface->closingNormalGapForce(coh1, predicting);

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion);

  // get norm of shear cohesion
  NodalFieldComponent shear_trac_norm(this->mesh, "shear_trac_norm");
  cohesion.computeNorm(shear_trac_norm, 1);
  double * tau_shear = shear_trac_norm.storage();

  // find current gap
  NodalField gap(this->mesh, "gap");
  this->interface->computeGap(gap, predicting);

  // compute norm of shear gap
  NodalFieldComponent shear_gap_norm(this->mesh, "shear_gap_norm");
  gap.computeNorm(shear_gap_norm, 1);
  double * shear_gap = shear_gap_norm.storage();

  // interface properties
  double * taus = this->tau_s.storage();
  double * taud = this->tau_d.storage();
  double * dc = this->d_c.storage();
  double * alpha = this->alpha.storage();

  double * p_coh1 = coh1.storage();

  // to be filled
  NodalFieldComponent factor_field(this->mesh, "factor");
  double * factor = factor_field.storage();

  // coh1 > 0 is a adhesive force
  // coh1 < 0 is a contact pressure
  for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n) {

    // avoid penetration "at any cost"
    // apply no normal cohesive force
    p_coh1[n] = std::min(p_coh1[n], 0.);

    double strength = taud[n] + (taus[n] - taud[n]) / pow(1+shear_gap[n]/dc[n],alpha[n]);
    
    // maximal shear cohesive force given by strength.
    // keep orientation of shear force
    factor[n] = std::min(1.,std::abs(strength / tau_shear[n]));
  }

  // only in shear direction
  for (int d=0; d<cohesion.getDim(); ++d) {
    if (d==1) // ignore normal direction
      continue;
    cohesion.multiplyByScalar(d,factor_field);
  }
}

/* -------------------------------------------------------------------------- */
void ContinuousShearWeakeningLaw::registerDumpField(const std::string & field_name) {

  // tau_s
  if (field_name == "tau_s") {
    this->interface->registerIO(field_name,
				this->tau_s);
  }

  // tau_d
  else if (field_name == "tau_d") {
    this->interface->registerIO(field_name,
				this->tau_d);
  }

  // d_c
  else if (field_name == "d_c") {
    this->interface->registerIO(field_name,
				this->d_c);
  }

  // alpha
  else if (field_name == "alpha") {
    this->interface->registerIO(field_name,
				this->alpha);
  }

  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }

}

/* -------------------------------------------------------------------------- */
void ContinuousShearWeakeningLaw::registerToRestart(Restart & restart) {

  restart.registerIO(this->tau_s);
  restart.registerIO(this->tau_d);
  restart.registerIO(this->d_c);
  restart.registerIO(this->alpha);
  
  InterfaceLaw::registerToRestart(restart);
}

__END_UGUCA__
