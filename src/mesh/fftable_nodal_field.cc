/**
 * @file   fftable_nodal_field.cc
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
#include "fftable_nodal_field.hh"
#include <stdexcept>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
FFTableNodalField::FFTableNodalField(FFTableMesh & mesh,
				     const std::string & name) :
  NodalField(mesh,name) {
  // free from nodal field initialization
  this->free();
  this->init(mesh);
}

/* -------------------------------------------------------------------------- */
void FFTableNodalField::init(FFTableMesh & mesh) {

  // do not initialize twice
  if (this->initialized)
    throw std::runtime_error("FFTableNodalField: do not initialize twice\n");

  this->mesh = &mesh;

  // if you want to use less components, simply use the NodalFieldComponent
  this->field.resize(this->mesh->getDim());
  for (int d=0; d<this->mesh->getDim(); ++d) {
    this->field[d] = new FFTableNodalFieldComponent(mesh,d,this->name);
  }

  // done
  this->initialized = true;
}

/* -------------------------------------------------------------------------- */
void FFTableNodalField::forwardFFT() {
  for (int d=0; d<this->mesh->getDim(); ++d) {
    ((FFTableNodalFieldComponent*)this->field[d])->forwardFFT();
  }
}

/* -------------------------------------------------------------------------- */
void FFTableNodalField::backwardFFT() {
  for (int d=0; d<this->mesh->getDim(); ++d) {
    ((FFTableNodalFieldComponent*)this->field[d])->backwardFFT();
  }
}

__END_UGUCA__
