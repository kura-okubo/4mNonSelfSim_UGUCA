/**
 * @file   linear_normal_cohesive_law.hh
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Fri Jun 3 2022
 * @date last modification: Fri Jun 3 2022
 *
 * @brief  linear cohesive law for fracture
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
#ifndef __LINEAR_NORMAL_COHESIVE_LAW_H__
#define __LINEAR_NORMAL_COHESIVE_LAW_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"
/*
   Linear cohesive law in normal direction only.
   No interpenetration allowed

   Parameter:
   Gc - fracture energy
   sigma_c peak strength
 */

__BEGIN_UGUCA__

class LinearNormalCohesiveLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:


  LinearNormalCohesiveLaw(BaseMesh & mesh,
			  double Gc_default,
			  double sigma_c_default,
			  const std::string & name = "lnclaw");

  virtual ~LinearNormalCohesiveLaw() {};

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
  NodalFieldComponent & getGc()     { return this->G_c;   }
  NodalFieldComponent & getSigmac() { return this->sigma_c; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  NodalFieldComponent G_c;
  NodalFieldComponent sigma_c;
};

__END_UGUCA__

//#include "linear_normal_cohesive_law_impl.cc"

#endif /* __LINEAR_NORMAL_COHESIVE_LAW_H__ */
