/**
 * @file   uca_simple_mesh.cc
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
#include "uca_simple_mesh.hh"
#include "static_communicator_mpi.hh"

#include <cmath>

__BEGIN_UGUCA__

/* --------------------------------------------------------------------------
 * 2D
 *
 * All operations in spatial domain are performed serially:
 * local==global for spatial domain
 * operations in fourier domain can be performed in parallel by default
 */
SimpleMesh::SimpleMesh(double Lx, int Nx, bool initialize) :
  DistributedFFTableMesh(Lx,Nx,false)
{
  if (initialize)
    this->init();
}

/* --------------------------------------------------------------------------
 * 3D with default coords
 *
 * All operations in spatial domain are performed serially:
 * local==global for spatial domain
 * operations in fourier domain can be performed in parallel by default
 *
 */
SimpleMesh::SimpleMesh(double Lx, int Nx,
		       double Lz, int Nz,
		       bool initialize) :
  DistributedFFTableMesh(Lx,Nx,Lz,Nz,false)
{
  if (initialize)
    this->init();
}

/* -------------------------------------------------------------------------- */
SimpleMesh::~SimpleMesh() {
}

/* -------------------------------------------------------------------------- */
void SimpleMesh::init() {
  this->initSimpleCoords(this->coords_local);
  DistributedFFTableMesh::init();
}

/* --------------------------------------------------------------------------
 * FFTW SERAIL METHODS
 * -------------------------------------------------------------------------- */
/*
 * default coordinate initiation as needed for FFTW serial
 *
 *              real data          complex data
 * array size   N0 x N1            2 x N0 x (N1/2+1)
 * stored as    N0 x N1            2 x N0 x (N1/2+1)
 *
 */
void SimpleMesh::initSimpleCoords(double ** coords) {

  // coords only exist on root
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  
  if (prank == this->root) {
  
    // determine element size
    double dx = this->getDeltaX();
    double dz = this->getDeltaZ();
    
    // fill coords for this mesh
    for (int i=0; i<this->nb_nodes_x_global; ++i) {
      for (int j=0; j<this->nb_nodes_z_global; ++j) {
	int ij = i*this->nb_nodes_z_global +j;
	
	coords[0][ij] = i*dx;
	coords[1][ij] = 0.0; // we are on the xz plane
	
	if (this->dim==3)
	  coords[2][ij] = j*dz;
      }
    }
  }
  else { // prank != root
    this->nb_nodes_local = 0;
    this->nb_nodes_local_alloc = 1;
  }
}

__END_UGUCA__
