/**
 * @file   half_space_quasidynamic.hh
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
#ifndef __HALF_SPACE_QUASIDYNAMIC_H__
#define __HALF_SPACE_QUASIDYNAMIC_H__
/* -------------------------------------------------------------------------- */
#include "half_space.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class HalfSpaceQuasiDynamic : public HalfSpace {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  HalfSpaceQuasiDynamic(FFTableMesh & mesh, int side_factor) :
    HalfSpace(mesh, side_factor) {}

  virtual ~HalfSpaceQuasiDynamic() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  
protected:
  virtual void computeStressFourierCoeff(bool /*predicting = false*/,
					 bool /*correcting = false*/) {
    throw std::runtime_error(
	  "HalfSpaceQuasiDynamic::computeStressFourierCoeff not implemented.");
  }

  void computeStressFourierCoeffQuasiDynamic(bool /*predicting*/,
					     bool /*correcting*/) {
    throw std::runtime_error(
	  "HalfSpaceQuasiDynamic::computeStressFourierCoeffQuasiDynamic not implemented.");
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

};

__END_UGUCA__

//#include "half_space_quasidynamic_impl.cc"

#endif /* __HALF_SPACE_QUASIDYNAMIC_H__ */
