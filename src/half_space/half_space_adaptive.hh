/**
 * @file   half_space_adaptive.hh
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
#ifndef __HALF_SPACE_ADAPTIVE_H__
#define __HALF_SPACE_ADAPTIVE_H__
/* -------------------------------------------------------------------------- */
#include "half_space_quasi_dynamic.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class HalfSpaceAdaptive : public HalfSpaceQuasiDynamic {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  HalfSpaceAdaptive(Material & material, FFTableMesh & mesh, int side_factor,
			SpatialDirectionSet components,
			const std::string & name = "half_space");

  virtual ~HalfSpaceAdaptive();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // init convolutions
  virtual void initConvolutions();

  // restart
  virtual void registerToRestart(Restart & restart);
  
  // set to steady state
  virtual void setSteadyState(bool predicting = false);

protected:
  virtual void computeStressFourierCoeff(bool predicting = false,
					 bool correcting = false,
					 SolverMethod sm = _quasi_dynamic,
					 unsigned int ts_factor = 1);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // set time step
  virtual void setTimeStep(double time_step);

  // get stable time step
  virtual double getStableTimeStep();


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  // keeps information if previous step was dynamic
  bool previously_dynamic;
};

__END_UGUCA__

//#include "half_space_adaptive_impl.cc"

#endif /* __HALF_SPACE_ADAPTIVE_H__ */
