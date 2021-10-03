/**
 * @file   nodal_field.cc
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
#include "nodal_field.hh"
#include <stdexcept>
#include <cmath>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
void NodalField::init(BaseMesh & mesh) {
  // do not initialize twice
  if (this->initialized)
    throw std::runtime_error("NodalField: do not initialize twice\n");
  
  this->mesh = &mesh;
  
  // if you want to use less components, simply use the NodalFieldComponent
  this->field.resize(this->mesh->getDim());
  for (int d=0; d<this->mesh->getDim(); ++d) {
    this->field[d] = new NodalFieldComponent(*(this->mesh),d,
					     this->name+"_"+std::to_string(d));
  }

  // done
  this->initialized = true;
}

/* -------------------------------------------------------------------------- */
void NodalField::free() {
  if (this->initialized) {
    for (int d=0; d<this->mesh->getDim(); ++d) {
      delete this->field[d];
    }
  }
  this->field.resize(0);
  this->initialized = false;
}

/* -------------------------------------------------------------------------- */
void NodalField::zeros() {
  this->setAllValuesTo(0.);
}

/* -------------------------------------------------------------------------- */
void NodalField::setAllValuesTo(double value) {
  if (!this->initialized)
    throw std::runtime_error("NodalField: cannot set value\n");

  for (int d=0; d<this->mesh->getDim(); ++d) {
    std::fill_n(this->field[d]->storage(),
		this->mesh->getNbLocalNodes(),
		value);
  }
}

/* -------------------------------------------------------------------------- */
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
void NodalField::multiplyByScalar(const NodalFieldComponent & scalar,
				       int ignore_dir) {

  const double * scalar_p = scalar.storage();

  for (int d=0; d<this->getDim(); ++d) {
    if (d == ignore_dir) continue;

    double * field_d_p = this->storage(d);

    for (int n=0; n<this->getNbNodes(); ++n) {
      field_d_p[n] *= scalar_p[n];
    }
  }
}


__END_UGUCA__
