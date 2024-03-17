/**
 * @file   half_space_dynamic.hh
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
#ifndef __HALF_SPACE_DYNAMIC_H__
#define __HALF_SPACE_DYNAMIC_H__
/* -------------------------------------------------------------------------- */
#include "half_space.hh"
#include "convolutions.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class HalfSpaceDynamic : public HalfSpace {
  
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  HalfSpaceDynamic(Material & material, FFTableMesh & mesh, int side_factor,
		   SpatialDirectionSet components,
		   const std::string & name = "half_space");

  virtual ~HalfSpaceDynamic();

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
  /// preintegrate kernels
  virtual void preintegrateKernels();

  /// compute the stress fourier coefficients (fails if not dynamic)
  virtual void computeStressFourierCoeff(bool predicting = false,
					 bool correcting = false,
					 bool dynamic = true);

  /// compute the stress fourier coefficients for dynamic case
  void computeStressFourierCoeffDynamic(bool predicting,
					bool correcting);

  /// compute F from U in fourier space
  void computeF(FFTableNodalField & F,
		const FFTableNodalField & U,
		const Convolutions & cnvls);
  
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

  // convolutions 
  Convolutions convols;
};

__END_UGUCA__

//#include "half_space_dynamic_impl.cc"

#endif /* __HALF_SPACE_DYNAMIC_H__ */
