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

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
BaseMesh::BaseMesh(int dim, int N) :
  root(0),
  dim(dim),
  nb_nodes_local(N),
  nb_nodes_local_alloc(N),
  coords_local(dim,N) {}

/* -------------------------------------------------------------------------- */
BaseMesh::~BaseMesh() {}

/* -------------------------------------------------------------------------- */
void BaseMesh::resize(int nb_nodes, int alloc) {
  this->nb_nodes_local = nb_nodes;
  this->nb_nodes_local_alloc = (alloc < nb_nodes ? nb_nodes : alloc);
  this->coords_local.resize(this->nb_nodes_local_alloc);
}

__END_UGUCA__
