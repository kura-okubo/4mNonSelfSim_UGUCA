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
SimpleMesh::SimpleMesh(double Lx, int Nx) : 
  DistributedFFTableMesh(Lx,Nx)
{
  // use const cast: not changing length, only content
  this->initSimpleCoords(const_cast<TwoDVector&>(this->getLocalCoords()));
  //this->init();
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
		       double Lz, int Nz) :
  DistributedFFTableMesh(Lx,Nx,Lz,Nz)
{
  // use const cast: not changing length, only content
  this->initSimpleCoords(const_cast<TwoDVector&>(this->getLocalCoords()));
  //this->init();
}

/* -------------------------------------------------------------------------- */
SimpleMesh::~SimpleMesh() {}

/* -------------------------------------------------------------------------- */
/*void SimpleMesh::init() {
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
void SimpleMesh::initSimpleCoords(TwoDVector & coords) {

  // coords only exist on root
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  
  if (prank == this->root) {
  
    // determine element size
    double dx = this->getDelta(0);
    double dz = this->getDelta(2);
    
    // fill coords for this mesh
    for (int i=0; i<this->nb_nodes_global[0]; ++i) {
      for (int j=0; j<this->nb_nodes_global[2]; ++j) {
	int ij = i*this->nb_nodes_global[2] +j;
	
	coords(ij,0) = i*dx;
	coords(ij,1) = 0.0; // we are on the xz plane
	
	if (this->dim==3)
	  coords(ij,2) = j*dz;
      }
    }
  }
  else { // prank != root
    BaseMesh::resize(0);
  }
}

__END_UGUCA__
