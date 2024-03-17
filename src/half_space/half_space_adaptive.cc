/**
 * @file   half_space_adaptive.cc
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
#include "half_space_adaptive.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
HalfSpaceAdaptive::HalfSpaceAdaptive(Material & material,
				     FFTableMesh & mesh,
				     int side_factor,
				     SpatialDirectionSet components,
				     const std::string & name) :
  HalfSpaceQuasiDynamic(material, mesh, side_factor, components, name),
  previously_dynamic(true) {}

/* -------------------------------------------------------------------------- */
HalfSpaceAdaptive::~HalfSpaceAdaptive() {}

/* -------------------------------------------------------------------------- */
double HalfSpaceAdaptive::getStableTimeStep() {
  return HalfSpaceDynamic::getStableTimeStep();
}

/* -------------------------------------------------------------------------- */
void HalfSpaceAdaptive::setTimeStep(double time_step) {
  HalfSpaceDynamic::setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceAdaptive::computeStressFourierCoeff(bool predicting,
						  bool correcting,
						  bool dynamic) {
  if (!dynamic) {
    this->computeStressFourierCoeffQuasiDynamic(predicting, correcting);
    this->previously_dynamic = false;
  }
  else {
    if (!this->previously_dynamic)
      this->setSteadyState(predicting);
    this->computeStressFourierCoeffDynamic(predicting, correcting);
    this->previously_dynamic = true;
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpaceAdaptive::initConvolutions() {
  HalfSpaceDynamic::initConvolutions();
}

/* -------------------------------------------------------------------------- */
void HalfSpaceAdaptive::registerToRestart(Restart & restart) {
  //HalfSpaceQuasiDynamic::registerToRestart(restart); // does nothing but cause bug
  HalfSpaceDynamic::registerToRestart(restart);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceAdaptive::setSteadyState(bool predicting) {
  HalfSpaceDynamic::setSteadyState(predicting);
}


__END_UGUCA__
