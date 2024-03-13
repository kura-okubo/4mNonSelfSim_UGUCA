/**
 * @file   uca_custom_mesh.hh
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Fri Feb 5 2021
 *
 * @brief  Contains only coordinates, but does not know anything about them
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
#ifndef __CUSTOM_MESH_H__
#define __CUSTOM_MESH_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_simple_mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

class CustomMesh : public SimpleMesh {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  CustomMesh(double Lx, int Nx,
	     std::vector<double> & x_coords);
  CustomMesh(double Lx, int Nx,
	     double Lz, int Nz,
	     std::vector<double> & x_coords,
	     std::vector<double> & z_coords);

protected:
  // protected constructor for inheritated construction without init
  CustomMesh(double Lx, int Nx);
  CustomMesh(double Lx, int Nx,
	     double Lz, int Nz);
  
public:
  virtual ~CustomMesh();
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void forwardFFT(FFTableNodalField & nodal_field);
  virtual void backwardFFT(FFTableNodalField & nodal_field);

protected:
  virtual void init(std::vector<double> & x_coords,
		    std::vector<double> & z_coords);

  virtual void initCustomCoords(std::vector<double> & x_coords,
				std::vector<double> & z_coords);

  virtual void initSortCustomNodesMap();
  virtual void checkCustomCoords(double ** coords_global);

  // comm for parallel
  template <typename T> void sortCustomNodes  (T * un_sorted, T * sorted,
					       int root_rank);
  template <typename T> void unsortCustomNodes(T * sorted,    T * un_sorted,
					       int root_rank);

  void sortAndScatterCustomNodes(int * Uglobal, int * Ulocal, int root_rank);
  void sortAndScatterCustomNodes(int * U, int root_rank){sortAndScatterCustomNodes(U,U,root_rank);};

  void sortAndScatterCustomNodes(double * Uglobal, double * Ulocal, int root_rank);//used for distributing Nodes
  void sortAndScatterCustomNodes(double * U, int root_rank){sortAndScatterCustomNodes(U,U,root_rank);};

  void gatherAndSortCustomNodes(int * Uglobal, int * Ulocal, int root_rank);
  void gatherAndSortCustomNodes(int * U, int root_rank){gatherAndSortCustomNodes(U,U,root_rank);};

  void gatherAndSortCustomNodes(double * Uglobal, double * Ulocal, int root_rank);
  void gatherAndSortCustomNodes(double * U, int root_rank){gatherAndSortCustomNodes(U,U,root_rank);};

private:
  template <typename T>
  void sortAndScatterCustomNodes(T * Uglobal, T * Ulocal, T * buffer, int size,
				 int root_rank);

  template <typename T>
  void gatherAndSortCustomNodes (T * Uglobal, T * Ulocal, T * buffer, int size,
				 int root_rank);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  // for sorting nodes to processes
  int * sort_custom_nodes_map;

  // maximum of local nodes per proc
  int max_nodes_pp;
  
  // for sorting nodes
  double * double_buffer;
};

__END_UGUCA__

//#include "uca_custom_mesh_impl.cc"

#endif /* __CUSTOM_MESH_H__ */
