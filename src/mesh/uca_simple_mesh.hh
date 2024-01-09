/**
 * @file   uca_simple_mesh.hh
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
#ifndef __SIMPLE_MESH_H__
#define __SIMPLE_MESH_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_distributed_fftable_mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

class SimpleMesh : public DistributedFFTableMesh {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  SimpleMesh(double Lx, int Nx);
  SimpleMesh(double Lx, int Nx,
	     double Lz, int Nz);

  virtual ~SimpleMesh();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  //virtual void init();
  
protected:
  virtual void initSimpleCoords(TwoDVector & coords);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:


};

__END_UGUCA__

//#include "uca_simple_mesh_impl.cc"

#endif /* __SIMPLE_MESH_H__ */
