/**
 * @file   continuous_shear_weakening_law.hh
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
#ifndef __CONTINUOUS_SHEAR_WEAKENING_LAW_H__
#define __CONTINUOUS_SHEAR_WEAKENING_LAW_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"
/*
   Continuous shear weakening law.
   No interpenetration allowed
   but also no opening allowed.
   Thus: should only be used for pure mode II fracture

   Parameter:
   tau_s peak
   tau_d dynamic
   d_c characterstic length
   alpha order
 */

__BEGIN_UGUCA__

class ContinuousShearWeakeningLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:


  ContinuousShearWeakeningLaw(BaseMesh & mesh,
			      double tau_s_default,
			      double tau_d_default,
			      double d_c_default,
			      double alpha,
			      const std::string & name = "cswlaw");

  virtual ~ContinuousShearWeakeningLaw() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
 void computeCohesiveForces(NodalField & cohesion,
                            bool predicting = false);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  // restart
  virtual void registerToRestart(Restart & restart);

 /* ------------------------------------------------------------------------ */
 /* Accessors                                                                */
 /* ------------------------------------------------------------------------ */
public:
  NodalFieldComponent & getTauc()  { return this->tau_s; }
  NodalFieldComponent & getTaur()  { return this->tau_d; }
  NodalFieldComponent & getDc()    { return this->d_c; }
  NodalFieldComponent & getAlpha() { return this->alpha; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  NodalFieldComponent tau_s;
  NodalFieldComponent tau_d;
  NodalFieldComponent d_c;
  NodalFieldComponent alpha;
};

__END_UGUCA__

//#include "continuous_shear_weakening_law_impl.cc"

#endif /* __CONTINUOUS_SHEAR_WEAKENING_LAW_H__ */
