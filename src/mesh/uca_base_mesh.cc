/**
 * @file   uca_base_mesh.cc
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
#include "uca_base_mesh.hh"
#include <stdexcept>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
BaseMesh::BaseMesh(int dim, int N) :
  root(0),
  //  rs_allocated(false),
  dim(dim),
  nb_nodes_local(N),
  nb_nodes_local_alloc(N),
  coords_local(dim,N) {
  //this->allocateRealSpace();
}

/* -------------------------------------------------------------------------- */
BaseMesh::~BaseMesh() {
  //  this->freeRealSpace();
}

/* -------------------------------------------------------------------------- */
void BaseMesh::resize(int nb_nodes, int alloc) {
  this->nb_nodes_local = nb_nodes;
  this->nb_nodes_local_alloc = (alloc < nb_nodes ? nb_nodes : alloc);
  this->coords_local.resize(this->nb_nodes_local_alloc);
}

/* -------------------------------------------------------------------------- */
/*
void BaseMesh::allocateRealSpace() {

  // do not allocate twice
  if (this->rs_allocated)
    throw std::runtime_error("BaseMesh: do not allocate twice\n");

  this->allocateVector(this->coords_local, this->nb_nodes_local_alloc);
  this->rs_allocated = true;
}

/* -------------------------------------------------------------------------- */
/*void BaseMesh::freeRealSpace() {
  /*
  if (this->rs_allocated) {
    this->freeVector(this->coords_local);
  }
  this->rs_allocated = false;
  
}

/* -------------------------------------------------------------------------- */
/*void BaseMesh::allocateVector(double ** vec, int size) {

  for (int d=0; d<this->dim; ++d) {
    vec[d] = new double [size];
    for (int n=0; n<size; ++n) {
      vec[d][n] = 0.; // initialize to zero
    }
  }
}

/* -------------------------------------------------------------------------- */
/*void BaseMesh::freeVector(double ** vec) {
  for (int d=0; d<this->dim; ++d) {
    delete[] vec[d];
  }
  }*/

__END_UGUCA__
