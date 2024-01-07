/**
 * @file   fftable_nodal_field.hh
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
/* -------------------------------------------------------------------------- */
#ifndef __FFTABLE_NODAL_FIELD_H__
#define __FFTABLE_NODAL_FIELD_H__
/* -------------------------------------------------------------------------- */
#include <map>
#include "uca_common.hh"
#include "nodal_field.hh"
#include "uca_fftable_mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__
/* **
 * MPI:
 * use default FFTW MPI datastructure N0x(N1/2+1)*2
 * note integer division rounds down
 *
 * Serial:
 * use default FFTW datastructure N0xN1
 */
class FFTableNodalField : public NodalField {

  friend class FFTableMesh;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  FFTableNodalField(const std::string & name = "unnamed") : NodalField(name) {}
  /*FFTableNodalField(FFTableMesh & mesh,
		    SpatialDirectionSet components = {_x},
		    const std::string & name = "unnamed");*/
  FFTableNodalField(FFTableMesh & mesh,
		    SpatialDirectionSet components = {0},
		    const std::string & name = "unnamed");

  virtual ~FFTableNodalField() {}

private:
  // private copy constructor: NodalField cannot be copied (for now to debug)
  FFTableNodalField(FFTableNodalField & to_copy);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // clears the NodalField and reinitializes it
  virtual void init(BaseMesh & mesh, SpatialDirectionSet components) {
    NodalField::init(mesh, components);
    this->init();
  }
private:
  virtual void init();

public:
  void forwardFFT();
  void backwardFFT();
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get number of nodes
  int getNbFFT() const { return ((FFTableMesh *)this->mesh)->getNbLocalFFT(); }

  // get fftw plan id for component
  inline int getFFTWPlanId(int d) { return this->fftw_plan_ids[d]; }

  /*
  inline FFTableNodalFieldComponent & component(int i) {
    return (FFTableNodalFieldComponent&)(*this->field[i]);
  }
  */

  // get one value of frequency domain in direction d
  inline fftw_complex & fd_p(int f, int d=0);                       // < ----------- call it fd_p for now to find all the ones to switch d and f
  //    return ((FFTableNodalFieldComponent*)(this->field[d]))->fd(f);
  
  // get access directly to frequency domain
  // WARNING: convert it to double (assuming that fftw_complex is double[2])
  inline fftw_complex * fd_data(int d=0);
    //    return ((FFTableNodalFieldComponent*)(this->field[d]))->fd_storage();
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  // start indices for each component
  std::vector<int> fd_start;
  
  // values in frequency domain in complex form
  std::vector<fftw_complex> fd_storage;

  // fftw plan id (given by mesh)
  std::map<int, int> fftw_plan_ids;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline fftw_complex & FFTableNodalField::fd_p(int f, int d) {
  if (!this->components.count(d)) 
    throw std::runtime_error("FFTableNodalField "+this->name+" has no component "+std::to_string(d)+"\n");
  return this->fd_storage[this->fd_start[d]+f];
}

inline fftw_complex * FFTableNodalField::fd_data(int d) {
  if (!this->components.count(d))
    throw std::runtime_error("FFTableNodalField "+this->name+" has no component "+std::to_string(d)+"\n");
  return this->fd_storage.data() + this->fd_start[d];
}

__END_UGUCA__

#endif /* __FFTABLE_NODAL_FIELD_H__ */
