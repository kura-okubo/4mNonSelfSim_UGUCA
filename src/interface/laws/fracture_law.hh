/**
 * @file   fracture_law.hh
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
#ifndef __FRACTURE_LAW_H__
#define __FRACTURE_LAW_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"

/*
   This is not exactly the same law than used in Barras et al. 2014.
   Within the cohesive zone, this law can close in normal direction,
   whereas Fracture' law will only maintain the normal gap.
   Also, here no friction behind the cohesive zone is applied.

   strength = tau_max ( 1 - |delta| / delta_c )
 */

__BEGIN_UGUCA__

class FractureLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  FractureLaw(BaseMesh & mesh,
	    double tau_max_default, double delta_c_default,
	    const std::string & name = "fraclaw");
  virtual ~FractureLaw() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void computeCohesiveForces(NodalField & cohesion,
			     bool predicting = false);
  
  virtual void registerDumpField(const std::string & field_name);
  
  // restart
  virtual void registerToRestart(Restart & restart);

 /* ------------------------------------------------------------------------ */
 /* Accessors                                                                */
 /* ------------------------------------------------------------------------ */
public:
  NodalFieldComponent & getTauMax() { return this->tau_max; }
  NodalFieldComponent & getDc()     { return this->delta_c; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  NodalFieldComponent tau_max;
  NodalFieldComponent delta_c;
  NodalFieldComponent gap_norm;
  NodalFieldComponent strength;
  
};

__END_UGUCA__

//#include "fracture_law_impl.cc"

#endif /* __FRACTURE_LAW_H__ */
