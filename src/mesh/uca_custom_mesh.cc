/**
 * @file   uca_custom_mesh.cc
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
#include "uca_custom_mesh.hh"
#include "static_communicator_mpi.hh"
#include "uca_simple_mesh.hh"
#include "fftable_nodal_field.hh"

#include <iomanip>
#include <cmath>
#include <algorithm>
#include <iostream>

__BEGIN_UGUCA__

/* --------------------------------------------------------------------------
 * 2D
 */
CustomMesh::CustomMesh(double Lx, int Nx,
		       std::vector<double> & x_coords) :
  SimpleMesh(Lx,Nx),
  max_nodes_pp(-1)
{
  std::vector<double> empty(0);
  this->init(x_coords, empty);
}

/* --------------------------------------------------------------------------
 * 3D
 */
CustomMesh::CustomMesh(double Lx, int Nx,
		       double Lz, int Nz,
		       std::vector<double> & x_coords,
		       std::vector<double> & z_coords) :
  SimpleMesh(Lx,Nx,Lz,Nz),
  max_nodes_pp(-1)
{
  this->init(x_coords, z_coords);
}

/* -------------------------------------------------------------------------- */
void CustomMesh::init(std::vector<double> & x_coords,
		      std::vector<double> & z_coords) {
    
  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();

  this->initCustomCoords(x_coords, z_coords);

  // all procs could gather and scatter as root
  this->double_buffer.resize(this->max_nodes_pp*psize);
  
  this->initSortCustomNodesMap();
}

/* -------------------------------------------------------------------------- */
void CustomMesh::forwardFFT(FFTableNodalField & nodal_field) {
  
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  NodalField tmp;
  // copy displacement for after gatherAndSort
  if (prank == this->root) {
    tmp.copyDataFrom(nodal_field);
  }

  // loop over all components of the nodal field
  for (const auto& d : nodal_field.getComponents()) {
    this->gatherAndSortCustomNodes(nodal_field.data(d), this->root);
  }
  SimpleMesh::forwardFFT(nodal_field);

  // restore the displacement
  if (prank == this->root) {
    nodal_field.copyDataFrom(tmp);
  }
}

/* -------------------------------------------------------------------------- */
void CustomMesh::backwardFFT(FFTableNodalField & nodal_field) {

  SimpleMesh::backwardFFT(nodal_field);

  // loop over all components of the nodal field
  for (const auto& d : nodal_field.getComponents()) {
    this->sortAndScatterCustomNodes(nodal_field.data(d), this->root);
  }
}

/* -------------------------------------------------------------------------- */
void CustomMesh::initCustomCoords(std::vector<double> & x_coords,
				  std::vector<double> & z_coords) {

  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  
  if ((z_coords.size() > 0) && (x_coords.size() != z_coords.size()))
    throw std::runtime_error("CustomMesh: provided coords do not match in size\n");

  int nb_nodes_local = x_coords.size();

  // find max nb nodes per process
  StaticCommunicatorMPI::getInstance()->allReduce(&nb_nodes_local,
						  &this->max_nodes_pp,
						  1,
						  MPI_MAX);

  // fftforward/backward will be done by the first procs in parallel
  int nb_nodes_local_alloc = -1;
  if (prank == this->root) //< this->dim)
    nb_nodes_local_alloc = this->max_nodes_pp*psize;
  else
    nb_nodes_local_alloc = this->max_nodes_pp;

  // resize the local coords
  BaseMesh::resize(nb_nodes_local, nb_nodes_local_alloc);
  // use const cast: not changing length, only content
  TwoDVector & coords_local = const_cast<TwoDVector&>(this->getLocalCoords());
  
  // copy coordinates
  for (int n=0; n<nb_nodes_local; ++n){
    coords_local(n,0) = x_coords[n];
    coords_local(n,1) = 0.;  // x-z plane is always at y=0
    if (this->dim > 2) {
      coords_local(n,2) = z_coords[n];
    }
  }

  // fill rest with NaN for initSortCustomNodesMap
  for (int n=nb_nodes_local; n<this->max_nodes_pp; ++n) {
    for (int d=0; d<this->dim; ++d)
      coords_local(n,d) = NAN;
  }
}

/* -------------------------------------------------------------------------- */
struct meshcompare {
  meshcompare(const std::vector<double> & coord_x, std::vector<double> & coord_z) :
    coord_x(coord_x),
    coord_z(coord_z)
  { }
  
  bool operator() (int i, int j ) {
    if (std::isnan(this->coord_x[i]))
      return false;
    
    if (std::isnan(this->coord_x[j]))
      return true;

    if ((float)this->coord_x[i]==(float)this->coord_x[j]) // to avoid problems with numerical error in the mesh
      return this->coord_z[i]<this->coord_z[j];
    else
      return this->coord_x[i]<this->coord_x[j];
  }
  
  const std::vector<double> & coord_x;
  const std::vector<double> & coord_z;
};

/* --------------------------------------------------------------------------
 * provide array with nodal indexes - custom domain decomposition
 * needs reordering for compatibility with FFTW
 */
void CustomMesh::initSortCustomNodesMap() {

  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
 
  // allocate sort nodes now that you know needed size
  this->sort_custom_nodes_map.resize(this->getNbGlobalNodes());
  std::fill(this->sort_custom_nodes_map.begin(), this->sort_custom_nodes_map.end(), 0);

  TwoDVector coords_global_tmp(this->dim,this->getNbLocalNodesAlloc());
  
  // gather all local mesh to root
  for (int d=0; d<this->dim; ++d) {
    StaticCommunicatorMPI::getInstance()->gather(coords_global_tmp.data(d).data(),
						 this->getLocalCoordsData(d),
						 this->max_nodes_pp,
						 this->root);
  }
  
  if (prank == this->root) {
    
    meshcompare mymeshcompare(coords_global_tmp.data(0),
			      coords_global_tmp.data(this->dim-1));
    
    // argsort by increasing x, z coords
    std::vector<int> sort_map_vec(0);
    
    // init vector with indexes from 0 to max_nb_nodespp*psize
    for (int i=0; i<this->max_nodes_pp*psize; ++i)
      sort_map_vec.push_back(i);
    
    // reorder indexes -> argsort
    std::sort(sort_map_vec.begin(), sort_map_vec.end(), mymeshcompare);
    
    // create sorting map
    for (int i=0; i<this->getNbGlobalNodes(); ++i)
      this->sort_custom_nodes_map[i] = sort_map_vec[i];

    // verify coords
    this->checkCustomCoords(coords_global_tmp);
  }
  
  StaticCommunicatorMPI::getInstance()->broadcast(this->sort_custom_nodes_map.data(),
						  this->getNbGlobalNodes(),
						  this->root);
}

/* -------------------------------------------------------------------------- */
void CustomMesh::checkCustomCoords(TwoDVector & coords_global) {

  double dx = this->getDelta(0);
  double tolerance = dx*1e-6;

  TwoDVector coords_global_ref(this->dim, this->getNbGlobalNodes());
  this->initSimpleCoords(coords_global_ref);
  
  // declare arrays and alloc memory
  TwoDVector coords_global_sorted(this->dim, this->getNbGlobalNodes());;
  
  //----------------------------------------------------
  // sort coords using sorting map

  for (int d=0; d<this->dim; ++d)
    this->sortCustomNodes(coords_global.data(d).data(),
			  coords_global_sorted.data(d).data(),
			  this->root);

  // find origin
  double x0 = coords_global_sorted(0,0);
  double z0 = 0.;
  if (this->dim == 3)
    z0 = coords_global_sorted(0,2);
  std::vector<double> origin = {x0,0,z0};

  //----------------------------------------------------
  // compare with coords from simple mesh

  // compare and raise error is difference exceeds tolerance
  bool test_passed = true;

  for (int n=0; n<this->getNbGlobalNodes(); ++n) {
    for (int d=0; d<this->dim; d+=2) {
      double error = std::abs(coords_global_sorted(n,d)
			      - origin[d]
			      - coords_global_ref(n,d));
      if (error > tolerance) {
	test_passed = false;
	if (d>0)
	  std::cerr<<n<<" error in custom mesh : is "
		   <<std::setw(6)<<      coords_global_sorted(n,0) - origin[0]
		   <<std::setw(6)<<", "<<coords_global_sorted(n,d) - origin[d]
		   <<" should "
		   <<std::setw(6)<<      coords_global_ref(n,0)
		   <<std::setw(6)<<", "<<coords_global_ref(n,d)
		   <<std::endl;
      }
    }
  }
  if (!test_passed) {
    throw std::runtime_error("Error; custom mesh is not a regular grid\n");
  }
}


/* --------------------------------------------------------------------------
 * SORT CUSTOM NODES
 * -------------------------------------------------------------------------- */
template <typename T>
void CustomMesh::sortCustomNodes(T* un_sorted, T * sorted, int root_rank) {
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (prank == root_rank) {
    for (int n=0; n<this->getNbGlobalNodes(); ++n) {
      sorted[n] = un_sorted[this->sort_custom_nodes_map[n]];
    }
  }
}
/* -------------------------------------------------------------------------- */
template <typename T>
void CustomMesh::unsortCustomNodes(T * sorted, T * un_sorted, int root_rank) {
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (prank == root_rank) {
    for (int n=0; n<this->getNbGlobalNodes(); ++n) {
      un_sorted[this->sort_custom_nodes_map[n]] = sorted[n];
    }
  }
}

/* --------------------------------------------------------------------------
 * SORT SCATTER GATHER CUSTOM NODES TEMPLATE
 * -------------------------------------------------------------------------- */
template <typename T>
void CustomMesh::sortAndScatterCustomNodes(T * Uglobal, T * Ulocal, T * buffer,
					   int size, int root_rank) {

  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();

  this->unsortCustomNodes(Uglobal, buffer, root_rank);
  
  if (psize > 1) {
    StaticCommunicatorMPI::getInstance()->scatter(buffer, Ulocal,
						  size, root_rank);
  }
  else {
    for (int n=0; n<this->getNbLocalNodes(); ++n)
      Ulocal[n] = buffer[n];
  }
}
/* -------------------------------------------------------------------------- */
template <typename T>
void CustomMesh::gatherAndSortCustomNodes(T * Uglobal, T * Ulocal, T * buffer,
					  int size, int root_rank) {
  
  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();

  if (psize > 1) {
    StaticCommunicatorMPI::getInstance()->gather(buffer, Ulocal,
						 size, root_rank);
  }
  else {
    for (int n=0; n<this->getNbLocalNodes(); ++n)
      buffer[n]=Ulocal[n];
  }
  
  this->sortCustomNodes(buffer,Uglobal,root_rank);
}

/* instead of template specialization there are a lot of functions calling
 * the templated one
 */
/* -------------------------------------------------------------------------- */
void CustomMesh::sortAndScatterCustomNodes(int * Uglobal, int * Ulocal,
					   int root_rank) {

  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  // we are allocating and freeing memory here just
  // because this function is used only during initialization
  int * bufferint = NULL; 

  if (prank == root_rank)
    bufferint = new int[this->max_nodes_pp*psize]();

  this->sortAndScatterCustomNodes(Uglobal,Ulocal,bufferint,
				  this->max_nodes_pp,root_rank);

  if (prank==root_rank) delete[] bufferint;
}

/* -------------------------------------------------------------------------- */
void CustomMesh::gatherAndSortCustomNodes(int * Uglobal, int * Ulocal,
					  int root_rank) {
  
  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  // we are allocating and freeing memory here just
  // because this function is used only during initialization
  int * bufferint = NULL; 

  if (prank == root_rank)
    bufferint = new int[this->max_nodes_pp*psize]();

  this->gatherAndSortCustomNodes(Uglobal,Ulocal,bufferint,
				 this->max_nodes_pp, root_rank);

  if (prank==root_rank) delete[] bufferint;
}

/* -------------------------------------------------------------------------- */
void CustomMesh::sortAndScatterCustomNodes(double * Uglobal, double * Ulocal,
					   int root_rank) {
  this->sortAndScatterCustomNodes(Uglobal,Ulocal,this->double_buffer.data(),
				  this->max_nodes_pp,root_rank);
}

/* -------------------------------------------------------------------------- */
void CustomMesh::gatherAndSortCustomNodes(double * Uglobal, double * Ulocal,
					  int root_rank) {

  this->gatherAndSortCustomNodes(Uglobal, Ulocal, this->double_buffer.data(),
				 this->max_nodes_pp, root_rank);
}

__END_UGUCA__
