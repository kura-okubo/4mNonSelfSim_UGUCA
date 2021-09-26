/**
 * @file   linear_coulomb_friction_law.cc
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
#include "linear_coulomb_friction_law.hh"
#include "interface.hh"

#include <cmath>
#include <iostream>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
LinearCoulombFrictionLaw::LinearCoulombFrictionLaw(BaseMesh & mesh,
						   double mu_s_default,
						   double mu_k_default,
						   double d_c_default,
						   double char_reg_time) :
  InterfaceLaw(mesh),
  reg_contact_pressure(mesh),
  mu_s(mesh),
  mu_k(mesh),
  d_c(mesh),
  char_time(mesh),
  reg_cont_pres_tmp(mesh)
{
  if (d_c_default < 1e-12) {
    std::cerr << "d_c cannot be zero, and it is currently: " << d_c_default << std::endl;
    throw;
  }

  this->initialized = false;
  this->reg_contact_pressure.setAllValuesTo(0. );
  this->mu_s.setAllValuesTo(mu_s_default);
  this->mu_k.setAllValuesTo(mu_k_default);
  this->d_c.setAllValuesTo(d_c_default);
  this->char_time.setAllValuesTo(char_reg_time);
}

/* -------------------------------------------------------------------------- */
void LinearCoulombFrictionLaw::computeCohesiveForces(NodalField & cohesion,
						     bool predicting) {

  // find forces needed to close normal gap
  NodalFieldComponent & coh1 = cohesion.component(1);
  this->interface->closingNormalGapForce(coh1, predicting);

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion);

  // get norm of shear cohesion
  NodalFieldComponent shear_trac_norm(this->mesh);
  cohesion.computeNorm(shear_trac_norm, 1);
  double * tau_shear = shear_trac_norm.storage();

  // find current gap
  //NodalField gap = this->interface->getBufferField();
  NodalField gap(this->mesh);
  this->interface->computeGap(gap, predicting);

  // compute norm of shear gap
  NodalFieldComponent shear_gap_norm(this->mesh);
  gap.computeNorm(shear_gap_norm, 1);
  double * shear_gap = shear_gap_norm.storage();

  // interface properties
  double * mus = this->mu_s.storage();
  double * muk = this->mu_k.storage();
  double * dc  = this->d_c.storage();

  double * p_coh1 = coh1.storage();
  
  // initialize regularized contact pressure
  if (!this->initialized) {
    for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n) {
      this->reg_contact_pressure(n) = p_coh1[n];
    }
    this->initialized = true;
  }

  // coh1 > 0 is a adhesive force
  // coh1 < 0 is a contact pressure
  for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n) {
    // avoid penetration "at any cost"
    // apply no normal cohesive force
    p_coh1[n] = std::min(p_coh1[n], 0.);
  }

  // regularized contact pressure
  double * reg_sig = NULL;
  if (predicting) {
    for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n)
      this->reg_cont_pres_tmp(n) = this->reg_contact_pressure(n);
    this->computeRegContactPressure(cohesion.component(1),
				    this->reg_cont_pres_tmp);
    reg_sig = this->reg_cont_pres_tmp.storage();
  }
  else {
    this->computeRegContactPressure(cohesion.component(1),
				    this->reg_contact_pressure);
    reg_sig = this->reg_contact_pressure.storage();
  }

  // to be filled
  NodalFieldComponent alpha_field(this->mesh);
  double * alpha = alpha_field.storage();

  for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n) {

    // compute friction coefficient
    double dmu = (1 - shear_gap[n] / dc[n]) * (mus[n] - muk[n]);
    double mu = std::max(muk[n], muk[n] + dmu);

    // maximal shear cohesive force given by strength.
    // keep orientation of shear force
    double strength = std::abs(mu * reg_sig[n]);

    // maximal shear cohesive force given by strength.
    // keep orientation of shear force
    alpha[n] = std::min(1.,std::abs(strength / tau_shear[n]));
  }

  // only in shear direction
  cohesion.multiplyByScalar(alpha_field, 1);
}

/* -------------------------------------------------------------------------- */
void LinearCoulombFrictionLaw::computeRegContactPressure(NodalFieldComponent & cohesion_1,
							 NodalFieldComponent & reg_cont_pres) {

  double dt = this->interface->getTimeStep();

  double * reg_sig = reg_cont_pres.storage();
  double * coh1 = cohesion_1.storage();
  double * tc  = this->char_time.storage();

  for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {

    // regularized
    if (tc[n] > 0) {
      // interface opened -> no history to preserve
      if (std::abs(coh1[n]) < 1e-12)
	reg_sig[n] = 0.;
      else
	reg_sig[n] = (reg_sig[n] + dt / tc[n] * coh1[n]) / (1 + dt / tc[n]);
    }
    // not regularized
    else {
      reg_sig[n] = coh1[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void LinearCoulombFrictionLaw::registerDumpField(const std::string & field_name) {

  // mu_s
  if (field_name == "mu_s") {
    this->interface->registerIO(field_name,
				     this->mu_s);
  }

  // mu_k
  else if (field_name == "mu_k") {
    this->interface->registerIO(field_name,
				     this->mu_k);
  }

  // d_c
  else if (field_name == "d_c") {
    this->interface->registerIO(field_name,
				     this->d_c);
  }

  // reg_cont_pres
  else if (field_name == "reg_cont_pres") {
    this->interface->registerIO(field_name,
				     this->reg_contact_pressure);
  }

  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }

}

__END_UGUCA__

