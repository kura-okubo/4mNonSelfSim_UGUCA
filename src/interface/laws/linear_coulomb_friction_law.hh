/**
 * @file   linear_coulomb_friction_law.hh
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
#ifndef __LINEAR_COULOMB_FRICTION_LAW_H__
#define __LINEAR_COULOMB_FRICTION_LAW_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"

/*
   Linear coulomb friction law.
   No interpenetration allowed
   but opening allowed.

   Parameter:
   mu_s static friction coefficient
   mu_k kinetic friction coefficient
   d_c characteristic weakening length
   Tstar regularization time

   Notes:
   if mu_s is equal to mu_k and d_c>0, this is static constant coulomb friction
 */

__BEGIN_UGUCA__

class LinearCoulombFrictionLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  LinearCoulombFrictionLaw(BaseMesh & mesh,
			   double mu_s_default,
			   double mu_k_default,
			   double d_c_default,
			   double char_reg_time = 0.);

  virtual ~LinearCoulombFrictionLaw() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void computeCohesiveForces(NodalField & cohesion,
			     bool predicting = false);

  void computeRegContactPressure(NodalFieldComponent & cohesion_1,
				 NodalFieldComponent & reg_cont_pres);

 // dumper function
 virtual void registerDumpField(const std::string & field_name);

 /* ------------------------------------------------------------------------ */
 /* Accessors                                                                */
 /* ------------------------------------------------------------------------ */
public:
  NodalFieldComponent & getMuS() { return this->mu_s; };
  NodalFieldComponent & getMuK() { return this->mu_k; };
  NodalFieldComponent & getDc() { return this->d_c; };
  NodalFieldComponent & getCharacteristicTime() { return this->char_time; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  bool initialized;

  NodalFieldComponent reg_contact_pressure;
  NodalFieldComponent mu_s;
  NodalFieldComponent mu_k;
  NodalFieldComponent d_c;
  NodalFieldComponent char_time;

  // for predictor-corrector approach
  NodalFieldComponent reg_cont_pres_tmp;
};

__END_UGUCA__

//#include "linear_coulomb_friction_law_impl.cc"

#endif /* __LINEAR_COULOMB_FRICTION_LAW_H__ */
