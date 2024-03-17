/**
 * @file   infinite_boundary.cc
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
#include "infinite_boundary.hh"
#include "half_space.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
InfiniteBoundary::InfiniteBoundary(FFTableMesh & mesh,
				   SpatialDirectionSet components,
				   int side_factor,
				   Material & material,
				   const std::string & name,
				   const SolverMethod & method) :
  Interface(mesh,components,name),
  external(mesh,components,"external") {

  this->hs = HalfSpace::newHalfSpace(material, mesh, side_factor, components, this->name+"_hs", method);
  
  this->half_spaces.resize(1);
  this->half_spaces[0] = this->hs;
}

/* -------------------------------------------------------------------------- */
InfiniteBoundary::~InfiniteBoundary() {
  delete this->hs;
}
				   
/* -------------------------------------------------------------------------- */
void InfiniteBoundary::computeResidual() {
  this->hs->computeResidual(this->external);
}

/* -------------------------------------------------------------------------- */
void InfiniteBoundary::predictTimeStepDirichlet() {

  // copy v to scratch memory
  this->updateVelocity();
  
  // Predict
  // u* = u + v * dt
  this->computeDisplacement(true);

  // tau_ext* -> compute residual -> tau_res*
  this->computeResidual();
  // tau_res* -> compute velocity -> v*
  this->computeVelocity(true);

  // Correct
  // v** = (v + v*) / 2 ---> overwrite storage if reached last prediction step
  this->correctVelocity(true);
}

/* -------------------------------------------------------------------------- */
void InfiniteBoundary::advanceTimeStepDirichlet() {

  // copy FEM displacement and internal
  this->computeDisplacement();
  this->computeInternal();
  this->computeResidual();
  this->computeVelocity();
}

// --------------------------------------------------------------------------
// NEUMANN BC on FEM
void InfiniteBoundary::computeExternal() {

  double mu = this->hs->getMaterial().getShearModulus();
  double Cs = this->hs->getMaterial().getCs();
  double Cp = this->hs->getMaterial().getCp();
  std::vector<double> eta = {1.0, Cp / Cs, 1.0};

  int sf = this->hs->getSideFactor();

  for (const auto& d : this->external.getComponents()) {
    double *int_p =  this->hs->getInternal().data(d);
    double *ext_p =  this->external.data(d);
    double *velo_p = this->hs->getVelo().data(d);
    double eta_d = eta[d];

    for (int n = 0; n < this->external.getNbNodes(); ++n) {
      ext_p[n] = - sf * mu / Cs * eta_d * velo_p[n] + int_p[n];
    }
  }
}


/* --------------------------------------------------------------------------
 *
 * external = traction to keep gap closed
 *
 * displacements and velocities are use in FEM as Dirichlet boundary condition
 * for the next time step
 *
 */
void InfiniteBoundary::advanceTimeStepNeumann() {

  // set displacement and velocity from FEM
  this->computeInternal();
  this->computeExternal();
}

/* -------------------------------------------------------------------------- */
void InfiniteBoundary::registerDumpField(const std::string &field_name) {

  bool registered = false;

  registered = this->hs->registerDumpFieldToDumper(field_name,
						   field_name,
						   this);

  if (!registered)
    Interface::registerDumpField(field_name);
}


__END_UGUCA__
