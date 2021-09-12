/**
 * @file   rate_and_state_law.cc
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

/**
 * Note: Only computes cohesion in 0 direction.
 */

#include "rate_and_state_law.hh"
#include "interface.hh"

#include <cmath>
#include <iostream>

__BEGIN_UGUCA__

/*
  Rate and state slip law in shear direction only.
  No interpenetration allowed
  but also no opening allowed.
  Thus: should only be used for pure mode II fracture
  Note: Currently only computes cohesion in 0 direction.
        Currently only supports UnimatShearInterface and BimatInterface.
*/

RateAndStateLaw::RateAndStateLaw(
       BaseMesh & mesh,
       double a_default,
       double b_default,
       double Dc_default,
       double V0,
       double f0,
       double theta_default,
       double sigma_default,
       EvolutionLaw evolution_law,
       bool predictor_corrector,
       double plate_velocity):
  InterfaceLaw(mesh),
  theta(mesh),
  theta_pc(),
  sigma(mesh),
  V(mesh),
  iterations(mesh),
  rel_error(mesh),
  V0(V0),
  f0(f0),
  a(mesh),
  b(mesh),
  Dc(mesh),
  predictor_corrector(predictor_corrector),
  evolution_law(evolution_law),
  Vplate(plate_velocity),
  Vw()
{
  this->theta.setAllValuesTo(theta_default);
  this->sigma.setAllValuesTo(sigma_default);
  this->a.setAllValuesTo(a_default);
  this->b.setAllValuesTo(b_default);
  this->Dc.setAllValuesTo(Dc_default);

  if (evolution_law == EvolutionLaw::SlipLawWithStrongRateWeakening) {
    this->Vw.init(mesh);
  }
  if (predictor_corrector) {
    this->theta_pc.init(mesh);
  }
}

/* -------------------------------------------------------------------------- */
RateAndStateLaw::~RateAndStateLaw() {}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::init() {
  // check if external load is not all zero on direction 2
  if (mesh.getDim() == 2)
    return;
  unsigned int nb = this->mesh.getNbLocalNodes();
  const double *external_2 = this->interface->getLoad().storage(2);
  bool checks_out = true;
  for (unsigned i = 0; i < nb; ++i) {
    if (external_2[i] != 0.0) {
      checks_out = false;
      break;
    }
  }
  if (!checks_out) {
    std::string message = "[ERROR] RateAndStateLaw currently only supports friction in 0 direction.";
    std::cerr << message << std::endl;
    throw message;
  }
  // check fw for EvolutionLaw::SlipLawWithStrongRateWeakening
  if (evolution_law == EvolutionLaw::SlipLawWithStrongRateWeakening) {
    if (fw < 0) {
      std::string message =
          "[ERROR] fw needs to be set before init() when using "
          "EvolutionLaw::SlipLawWithStrongRateWeakening";
      std::cerr << message << std::endl;
      throw message;
    }
  }
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::computeCohesiveForces(NodalField & cohesion,
					    bool predicting) {

  // find forces needed to close normal gap
  NodalFieldComponent & coh1 = cohesion.component(1);
  this->interface->closingNormalGapForce(coh1, predicting);
  // double * coh_1_p = cohesion[1]->storage();

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion);
  double * coh_0_p = cohesion.storage(0);

  // interface properties
  const HalfSpace & top = this->interface->getTop();
  const HalfSpace & bot = this->interface->getBot();
  const Material & mat_top = top.getMaterial();
  const Material & mat_bot = bot.getMaterial();
  double fact_top = mat_top.getCs() / mat_top.getShearModulus();
  double fact_bot = mat_bot.getCs() / mat_bot.getShearModulus();
  double fact_both = fact_top + fact_bot;

  const double *a_p = this->a.storage();
  const double *b_p = this->b.storage();
  const double *Dc_p = this->Dc.storage();
  const double *sigma_p = this->sigma.storage();
  const double *int0top_p = top.getInternal().storage(0);
  const double *int0bot_p = bot.getInternal().storage(0);
  const double *ext0_p = this->interface->getLoad().storage(0);

  //  std::vector<NodalField *> gap_velo = this->interface->getBufferField();
  NodalField gap_velo(this->mesh);
  // compute delta_dot: do not use v* here -- predicting = false
  this->interface->computeGapVelocity(gap_velo, false);
  // pass true to compute slip rate only in shear directions vectorially
  gap_velo.computeNorm(this->V, 1);
  const double * V_p = this->V.storage();

  // compute theta
  this->computeTheta(predicting ? this->theta_pc : this->theta, this->V);
  double * theta_p = (predicting ? this->theta_pc : this->theta).storage();

  double * iterations_p = this->iterations.storage();
  double * rel_error_p = this->rel_error.storage();
  // double * abs_error_p = this->abs_error->storage();

  // solve tau_coh using Newton-Raphson
  // V is then solved in Interface::advanceTimeStep()
  for (int i = 0; i < this->mesh.getNbLocalNodes(); ++i) {
    double Z;
    if (evolution_law == EvolutionLaw::SlipLawWithStrongRateWeakening){
      Z = 0.5 / V0 * std::exp(theta_p[i] / a_p[i]);
    } else {
      Z = 0.5 / V0 * std::exp((f0 + b_p[i] * std::log(V0 * theta_p[i] / Dc_p[i])) / a_p[i]);
    }

    // initial guess
    double v_prev = V_p[i];
    double tau, dtau, F, dF, v;
    double rel_change = 1;
    double rel_tol = 1e-8;
    unsigned max_iter = 1000;
    unsigned min_iter = 0;
    unsigned iter = 0;
    unsigned sign_change_count = 0;
    // Newton-Raphson
    do {
      ++iter;
      tau = a_p[i] * sigma_p[i] * std::asinh(Z * (v_prev + Vplate));
      dtau = a_p[i] * sigma_p[i] * Z / std::sqrt(1.0 + Z * Z * (v_prev + Vplate) * (v_prev + Vplate));
      F = fact_both * (ext0_p[i] - tau) + fact_top * int0top_p[i] +
          fact_bot * int0bot_p[i] - v_prev;
      dF = -fact_both * dtau - 1.0;
      v = v_prev - F / dF;
      rel_change = std::abs(F / dF / v_prev);
      // catching infinite loop when v is jumping across 0 more than 4 times
      // set v to 0 (Vp) and continue the iterations
      if ((v_prev + Vplate) * (v + Vplate) < 0) ++sign_change_count;
      if (sign_change_count >= 4) {
        v = -Vplate;
        sign_change_count = 0;
        rel_change = 1.0;
      }
      if (!std::isfinite(v)) {
        v = -Vplate;
        rel_change = 1.0;
      } else if (!std::isfinite(rel_change)) {
        rel_change = 1.0;
      }
      if (!std::isfinite(v)){
        v = -Vplate;
        rel_change = 1;
      }
      if (!std::isfinite(rel_change)) rel_change = 1.0;
      v_prev = v;
    } while ((rel_change > rel_tol && iter < max_iter) || iter < min_iter);
    if (iter == max_iter) {
      throw std::runtime_error("Newton-Raphson not converged in RateAndStateLaw::computeCohesiveForces");
    }
    tau = a_p[i] * sigma_p[i] * std::asinh(Z * (v + Vplate));
    coh_0_p[i] = tau;
    iterations_p[i] = iter;
    rel_error_p[i] = rel_change;
  }
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::computeTheta(NodalFieldComponent & target,
				   NodalFieldComponent & delta_dot) {

  double dt = this->interface->getTimeStep();

  const double *theta_p = this->theta.storage();
  const double *Dc_p = this->Dc.storage();
  const double *a_p = this->a.storage();
  const double *b_p = this->b.storage();

  double *theta_target = target.storage();
  double *v_p = delta_dot.storage();

  switch (evolution_law) {
    case EvolutionLaw::AgingLaw: {
      for (int i = 0; i < this->mesh.getNbLocalNodes(); ++i) {
        double v_i = std::abs(v_p[i] + Vplate);
        v_i = std::isfinite(v_i) ?
          std::max(v_i, this->Vguard) : this->Vguard;
        // assuming constant v within time step:
        double D = exp(-v_i * dt / Dc_p[i]);
        theta_target[i] = theta_p[i] * D + Dc_p[i] / v_i * (1.0 - D);
        // assuming constant d(theta)/dt within time step:
        // theta_target[i] = theta_p[i] + (1 - v_i * theta_p[i] / Dc_p[i]) * dt;
      }
      break;
    }
    case EvolutionLaw::SlipLaw: {
      for (int i = 0; i < this->mesh.getNbLocalNodes(); ++i) {
        double v_i = std::abs(v_p[i] + Vplate);
        v_i = std::isfinite(v_i) ?
          std::max(v_i, this->Vguard) : this->Vguard;
        // assuming constant v within time step:
        double D = exp(-v_i * dt / Dc_p[i]);
        theta_target[i] = Dc_p[i] / v_i * pow(v_i * theta_p[i] / Dc_p[i], D);
        // assuming constant d(theta)/dt within time step:
        // double D = v_i * theta_p[i] / Dc_p[i];
        // theta_target[i] = theta_p[i] - D * std::log(D) * dt;
      }
      break;
    }
    case EvolutionLaw::SlipLawWithStrongRateWeakening: {
      const double *Vw_p = this->Vw.storage();
      for (int i = 0; i < this->mesh.getNbLocalNodes(); ++i) {
        double v_i = std::abs(v_p[i] + Vplate);
        v_i = std::isfinite(v_i) ? std::max(v_i, this->Vguard) : this->Vguard;
        double f_LV_i = f0 - (b_p[i] - a_p[i]) * std::log(v_i / V0);
        double f_ss_i =
            fw + (f_LV_i - fw) / std::pow(1 + std::pow(v_i / Vw_p[i], 8.0), 0.125);
        double theta_ss_i = a_p[i] * std::log(2.0 * V0 / v_i *
                                         std::sinh(f_ss_i / a_p[i]));
        double D = exp(-v_i * dt / Dc_p[i]);
        theta_target[i] = D * theta_p[i] + (1.0 - D) * theta_ss_i;
      }
      break;
    }
    default:
      std::string message = "[ERROR] RateAndStateLaw found unsupported evolution law option.";
      std::cerr << message << std::endl;
      throw message;
  }
}

/* -------------------------------------------------------------------------- */
void RateAndStateLaw::registerDumpField(const std::string &field_name) {
  // theta
  if (field_name == "theta") {
    this->interface->registerForDump(field_name,
				     this->theta);
  }
  // iterations
  else if (field_name == "iterations") {
    this->interface->registerForDump(field_name,
				     this->iterations);
  }
  // rel_error in Newton-Raphson
  else if (field_name == "rel_error") {
    this->interface->registerForDump(field_name,
				     this->rel_error);
  }
  // a
  else if (field_name == "a") {
    this->interface->registerForDump(field_name,
				     this->a);
  }
  // b
  else if (field_name == "b") {
    this->interface->registerForDump(field_name,
				     this->b);
  }
  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }

}

NodalFieldComponent & RateAndStateLaw::getVw() {
  if (evolution_law == EvolutionLaw::SlipLawWithStrongRateWeakening)
    return this->Vw;
  else {
    throw std::runtime_error(
        "Vw is only used in EvolutionLaw::SlipLawWithStrongRateWeakening");
  }
}

void RateAndStateLaw::setFw(double fw) {
  if (evolution_law == EvolutionLaw::SlipLawWithStrongRateWeakening)
    this->fw = fw;
  else {
    throw std::runtime_error(
        "fw is only used in EvolutionLaw::SlipLawWithStrongRateWeakening");
  }
}
__END_UGUCA__
