/**
 * @file   uca_base_mesh.hh
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
#ifndef __BASE_MESH_H__
#define __BASE_MESH_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

class TwoDVector {
public:
  TwoDVector(int dim, int N) : storage(dim,std::vector<double>(N)) {}

  void resize(int N) {
    for (auto& vec : this->storage)
      vec.resize(N);
  }
  
  inline double & operator()(int n, int d) { return this->storage[d][n]; }
  inline double operator()(int n, int d) const { return this->storage[d][n]; }
  inline std::vector<double> & data(int d) { return this->storage[d]; }
  inline std::vector<std::vector<double>> & getStorage() {return this->storage; }
private:
  std::vector<std::vector<double>> storage;
private:
  // private copy constructor: avoid copying when forgetting &
  TwoDVector(TwoDVector & to_copy);
};

class BaseMesh {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  BaseMesh(int dim, int N);
  virtual ~BaseMesh();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // resize coordinates
  virtual void resize(int nb_nodes, int alloc=-1);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get root of mesh
  int getRoot() const { return this->root; }

  // spatial dimension of mesh
  int getDim() const { return this->dim; }

  // number of nodes
  virtual int getNbLocalNodes() const { return this->nb_nodes_local; }

  // access to coordinates
  //const TwoDVector & getLocalCoords() const { return this->coords_local; }
  
  TwoDVector & getLocalCoords() { return this->coords_local; }

  // inheritate mesh needs to know if it is parallel or not
  virtual int getNbGlobalNodes() const { return this->getNbLocalNodes(); };

  // allocation length
  int getNbLocalNodesAlloc() const { return this->nb_nodes_local_alloc; }

protected:
  double * getLocalCoordsData(int d) { return this->coords_local.data(d).data(); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // root rank for parallel
  int root;
  
  // spatial dimension: either 2 or 3
  int dim;


  // coordinates of local nodes
  TwoDVector coords_local;

private:
  // discretization: local
  int nb_nodes_local;
  int nb_nodes_local_alloc;

};

__END_UGUCA__

//#include "uca_base_mesh_impl.cc"

#endif /* __BASE_MESH_H__ */
