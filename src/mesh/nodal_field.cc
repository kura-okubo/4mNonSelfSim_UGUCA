/**
 * @file   nodal_field.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Fri Jan 5 2024
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
#include "nodal_field.hh"
#include <stdexcept>
#include <cmath>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
void NodalField::init(BaseMesh & mesh, SpatialDirectionSet components) {
  // do not initialize twice
  if (this->initialized)
    throw std::runtime_error("NodalField: do not initialize twice\n");
  
  this->mesh = &mesh;
  this->components = components;

  // Loop over the components to define their start
  this->start = {0};
  for (size_t i = 0; i < _spatial_dir_count; ++i) {
    int new_start = this->start.back();
    if (this->components.test(i)) {
      new_start += mesh.getNbLocalNodesAlloc();
    }
    this->start.push_back(new_start);
  }
  
#ifdef UCA_VERBOSE
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::cout << "NF init (prank="
	    << prank << "): " << this->name
	    << " : " << this->start.back() << std::endl;
#endif // UCA_VERBOSE

  // allocate space in storage and fill with zeros
  this->storage.resize(this->start.back(), 0);

  // done
  this->initialized = true;
}

/* -------------------------------------------------------------------------- */
/*
void NodalField::free() {
  if (this->initialized) {
    for (int d=0; d<this->mesh->getDim(); ++d) {
      delete this->field[d];
    }
  }
  this->field.resize(0);
  this->initialized = false;
}
*/

/* -------------------------------------------------------------------------- */
void NodalField::zeros() {
  this->setAllValuesTo(0.);
}

/* -------------------------------------------------------------------------- */
void NodalField::setAllValuesTo(double value, int d) {
  if (!this->initialized)
    throw std::runtime_error("NodalField: cannot set value\n");

  // check that this component exists
  if (d!=-1 && !this->components.test(d)) 
    throw std::runtime_error("NodalField "+this->name+" has no component "+std::to_string(d)+"\n");
  
  // if it should apply to all components
  if (d==-1) 
    std::fill(this->storage.begin(), this->storage.end(), value);
  else { // only on a given component
    std::fill(this->storage.begin()+this->start[d],
	      this->storage.begin()+this->start[d+1], value);
  }
}

/* -------------------------------------------------------------------------- */
/*
void NodalField::computeNorm(NodalFieldComponent & norm,
			     int ignore_dir) const {
  norm.zeros();
  double * norm_p = norm.storage();

  for (int d=0; d<this->getDim(); ++d) {
    if (d == ignore_dir) continue;

    const double * field_d_p = this->storage(d);

    for (int n=0; n<this->getNbNodes(); ++n) {
      norm_p[n] += field_d_p[n] * field_d_p[n];
    }
  }

  for (int n=0; n<this->getNbNodes(); ++n) {
    norm_p[n] = std::sqrt(norm_p[n]);
  }
}

/* -------------------------------------------------------------------------- */
void NodalField::multiplyByScalarField(const NodalField & scalar, int d) {

  // check that this component exists
  if (!this->components.test(d)) 
    throw std::runtime_error("NodalField "+this->name+" has no component "+std::to_string(d)+"\n");

  // check that both have same length
  if (this->getNbNodes() != scalar.getNbNodes())
    throw std::runtime_error("NodelFields "+this->name
			     +" and "+scalar.getName()
			     +" do not have same number of nodes\n");

  for (int n=0; n<this->getNbNodes(); ++n) {
    (*this)(n,d) *= scalar(n); // do not give d to check that scalar is a scalar
  }
}

/* -------------------------------------------------------------------------- */
void NodalField::multiplyByScalarField(const NodalField & scalar) {
  for (size_t i = 0; i < this->components.size(); ++i)
    this->multiplyByScalarField(scalar,static_cast<int>(i));
}


__END_UGUCA__
