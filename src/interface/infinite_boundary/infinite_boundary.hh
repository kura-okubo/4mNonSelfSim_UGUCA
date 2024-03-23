/**
 * @file   infinite_boundary.hh
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
#ifndef __INFINITE_BOUNDARY_H__
#define __INFINITE_BOUNDARY_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "interface.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class InfiniteBoundary : public Interface {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  // side factor top=1 bot=-1
  InfiniteBoundary(FFTableMesh & mesh,
		   SpatialDirectionSet components,
		   int side_factor,
		   Material & material,
		   const std::string & name = "inf_boundary",
		   const SolverMethod & method = _dynamic);

  virtual ~InfiniteBoundary();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initConvolutions() { this->hs->initConvolutions(); }
  
  void advanceTimeStepDirichlet();
  void predictTimeStepDirichlet();
  void advanceTimeStepNeumann();

  virtual void computeResidual();
  
  // dumper function
  virtual void registerDumpField(const std::string & field_name);

private:
  void computeExternal();

  // due to inheritance from interface
  void closingNormalGapForce(NodalField &, bool, unsigned int) {
    throw std::runtime_error(
	 "InfiniteBoundary::closingNormalGapForce not implemented.");
  }
  void maintainShearGapForce(NodalField & ) {
    throw std::runtime_error(
	 "InfiniteBoundary::maintainShearGapForce not implemented.");
  }
  void computeGap(NodalField &, bool) {
    throw std::runtime_error(
	 "InfiniteBoundary::computeGap not implemented.");
  }
  void computeGapVelocity(NodalField &, bool) {
    throw std::runtime_error(
	 "InfiniteBoundary::computeGapVelocity not implemented.");
  }
  HalfSpace & getTop() {
    throw std::runtime_error(
	 "InfiniteBoundary::getTop not implemented.");
  }
  HalfSpace & getBot() {
    throw std::runtime_error(
	 "InfiniteBoundary::getBot not implemented.");
  }
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  NodalField & getExternal() { return this->external; }

  FFTableNodalField & getDisp(bool predicting = false) {
    return this->hs->getDisp(predicting);
  }
  NodalField & getVelo(bool predicting = false) {
    return this->hs->getVelo(predicting); 
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // half space
  HalfSpace * hs;
  NodalField external;
};

__END_UGUCA__

//#include "infinite_boundary_impl.cc"

#endif /* __INFINITE_BOUNDARY_H__ */
