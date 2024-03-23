/**
 * @file   half_space_quasi_dynamic.cc
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
#include "half_space_quasi_dynamic.hh"
#include "static_communicator_mpi.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
HalfSpaceQuasiDynamic::HalfSpaceQuasiDynamic(Material & material,
					     FFTableMesh & mesh,
					     int side_factor,
					     SpatialDirectionSet components,
					     const std::string & name) :
  HalfSpaceDynamic(material, mesh, side_factor, components, name) {

#ifdef UCA_VERBOSE
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::cout << "HSQD construct (prank="
	    << prank << "): " << this->name
	    << " : " << this->mesh.getNbLocalFFT() << std::endl;
#endif // UCA_VERBOSE
}

/* -------------------------------------------------------------------------- */
HalfSpaceQuasiDynamic::~HalfSpaceQuasiDynamic() {}

/* -------------------------------------------------------------------------- */
double HalfSpaceQuasiDynamic::getStableTimeStep() {
  return HalfSpaceDynamic::getStableTimeStep();
}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::setTimeStep(double time_step) {
  HalfSpace::setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::initConvolutions() {
  HalfSpaceDynamic::initConvolutions();
}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::computeStressFourierCoeff(bool predicting,
						      bool correcting,
						      SolverMethod sm,
						      unsigned int /*ts_factor*/) {
  if (sm == _quasi_dynamic)
    this->computeStressFourierCoeffQuasiDynamic(predicting, correcting);
  else
    throw std::runtime_error("HalfSpaceQuasiDynamic cannot compute stress for dynamic problem!");
}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::computeStressFourierCoeffQuasiDynamic(bool predicting,
								  bool /*correcting*/) {

  // compute convolutions (with current U instead of history of U)
  this->convols.convolveSteadyState();

  // compute F
  FFTableNodalField & _disp = predicting ? this->disp_pc : this->disp;
  this->computeF(this->internal, _disp, this->convols);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::registerToRestart(Restart & restart) {
  HalfSpace::registerToRestart(restart);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::setSteadyState(bool /*predicting*/) {
  // do nothing, it is steady state at all times
}


__END_UGUCA__
