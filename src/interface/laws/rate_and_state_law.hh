/**
 * @file   rate_and_state_law.hh
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

#ifndef __RATE_AND_STATE_LAW_H__
#define __RATE_AND_STATE_LAW_H__
/* -------------------------------------------------------------------------- */

#include "interface_law.hh"

#include <vector>

/*
  Rate and state slip law in shear direction only.
  No interpenetration allowed
  but also no opening allowed.
  Thus: should only be used for pure mode II fracture
  Note: Only computes cohesion in 0 direction.
        Slip rate MUST NOT be zero.
*/

__BEGIN_UGUCA__

class RateAndStateLaw : public InterfaceLaw {
public:
  enum class EvolutionLaw{AgingLaw, SlipLaw, SlipLawWithStrongRateWeakening};
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  RateAndStateLaw(
       BaseMesh & mesh,
       double a_default,
       double b_default,
       double Dc_default,
       double V0,
       double f0,
       double theta_default,
       EvolutionLaw evolution_law = EvolutionLaw::AgingLaw,
       bool predictor_corrector = true,
       double plate_velocity = 0.0
       );
  virtual ~RateAndStateLaw();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void computeCohesiveForces(NodalField & cohesion,
			     bool predicting = false);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);
  virtual void init();

protected:
  virtual void computeTheta(NodalFieldComponent & target,
			    NodalFieldComponent & delta_dot);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  NodalFieldComponent & getTheta() { return this->theta; };
  NodalFieldComponent & getA() { return this->a; };
  NodalFieldComponent & getB() { return this->b; };
  NodalFieldComponent & getVw();    // for slip law with strong rate weakening
  void setFw(double fw);  // for slip law with strong rate weakening

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
 protected:
  NodalFieldComponent theta;          // state variable
  NodalFieldComponent theta_pc;       // state variable for predictor
  //NodalFieldComponent sigma;          // normal stress (compression is positive)
  NodalFieldComponent V;              // slip rate
  NodalFieldComponent iterations;     // iterations took in Newton-Raphson algorithm
  NodalFieldComponent rel_error;      // for Newton-Raphson algorithm debugging
  // NodalField * abs_error;          // for Newton-Raphson algorithm debugging
  double V0;                          // reference slip rate
  double f0;                          // reference parameter
  NodalFieldComponent a;              // friction parameter
  NodalFieldComponent b;              // friction parameter
  NodalFieldComponent Dc;             // friction parameter
  bool predictor_corrector;           // indicates whether predictor-corrector is activated
  EvolutionLaw evolution_law;         // indicates which state evolution law to use
  double Vguard = 1.0e-20;            // minimum absolute velocity
  double Vplate;                      // plate velocity
  NodalFieldComponent Vw;             // for slip law with strong rate weakening
  double fw = -1.0;                   // for slip law with strong rate weakening
};

__END_UGUCA__

#endif /* __RATE_AND_STATE_LAW_H__ */
