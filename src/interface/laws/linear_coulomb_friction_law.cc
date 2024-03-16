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
						   double char_reg_time,
						   const std::string & name) :
  InterfaceLaw(mesh,name),
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

  this->reg_contact_pressure.setAllValuesTo(0.);
  this->reg_contact_pressure.setName(name+"_reg_cont_pres");
  
  this->mu_s.setAllValuesTo(mu_s_default);
  this->mu_s.setName(name+"_mus_s");
  
  this->mu_k.setAllValuesTo(mu_k_default);
  this->mu_k.setName(name+"_mu_k");
  
  this->d_c.setAllValuesTo(d_c_default);
  this->d_c.setName(name+"_d_c");
  
  this->char_time.setAllValuesTo(char_reg_time);
  this->char_time.setName(name+"_chart_time");
}

/* -------------------------------------------------------------------------- */
void LinearCoulombFrictionLaw::computeCohesiveForces(NodalField & cohesion,
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

  // initialize regularized contact pressure
  if (!this->initialized) {
    for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n) {
      this->reg_contact_pressure(n) = cohesion(n,1);
    }
    this->initialized = true;
  }
  
  // coh1 > 0 is a adhesive force
  // coh1 < 0 is a contact pressure
  for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n) {
    // avoid penetration "at any cost"
    // apply no normal cohesive force
    cohesion(n,1) = std::min(cohesion(n,1), 0.);
  }

  // regularized contact pressure
  NodalField * reg_sig = NULL;
  if (predicting) {
    for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n)
      this->reg_cont_pres_tmp(n) = this->reg_contact_pressure(n);
    this->computeRegContactPressure(cohesion, this->reg_cont_pres_tmp);
    reg_sig = &(this->reg_cont_pres_tmp);
  }
  else {
    this->computeRegContactPressure(cohesion, this->reg_contact_pressure);
    reg_sig = &(this->reg_contact_pressure);
  }

  // to be filled
  NodalField alpha(this->mesh);

  for (int n = 0; n<this->mesh.getNbLocalNodes(); ++n) {

    // compute friction coefficient
    double dmu = (1 - shear_gap_norm(n) / d_c(n)) * (mu_s(n) - mu_k(n));
    double mu = std::max(mu_k(n), mu_k(n) + dmu);

    // maximal shear cohesive force given by strength.
    // keep orientation of shear force
    double strength = std::abs(mu * (*reg_sig)(n));

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
void LinearCoulombFrictionLaw::computeRegContactPressure(NodalField & cohesion,
							 NodalField & reg_cont_pres) {

  double dt = this->interface->getTimeStep();

  for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {

    // regularized
    if (this->char_time(n) > 0) {
      // interface opened -> no history to preserve
      if (std::abs(cohesion(n,1)) < 1e-12)
	reg_cont_pres(n) = 0.;
      else
	reg_cont_pres(n) = (reg_cont_pres(n) + dt / this->char_time(n) * cohesion(n,1))
	  / (1 + dt / this->char_time(n));
    }
    // not regularized
    else {
      reg_cont_pres(n) = cohesion(n,1);
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

/* -------------------------------------------------------------------------- */
void LinearCoulombFrictionLaw::registerToRestart(Restart & restart) {

  restart.registerIO(this->reg_contact_pressure);
  restart.registerIO(this->mu_s);
  restart.registerIO(this->mu_k);
  restart.registerIO(this->d_c);
  restart.registerIO(this->char_time);
  restart.registerIO(this->reg_cont_pres_tmp);
  
  InterfaceLaw::registerToRestart(restart);
}

__END_UGUCA__

