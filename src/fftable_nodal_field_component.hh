/**
 * @file   fftable_nodal_field_component.hh
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
#ifndef __FFTABLE_NODAL_FIELD_COMPONENT_H__
#define __FFTABLE_NODAL_FIELD_COMPONENT_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "nodal_field_component.hh"
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
class FFTableNodalFieldComponent : public NodalFieldComponent {

  friend class FFTableMesh;
  
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  FFTableNodalFieldComponent(FFTableMesh & mesh, int direction=0,
			     const std::string & name = "unnamed");
  virtual ~FFTableNodalFieldComponent();

private:
  // private copy constructor: NodalFieldComponent cannot be copied (for now to debug)
  FFTableNodalFieldComponent(FFTableNodalFieldComponent & to_copy);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void init(FFTableMesh & mesh);
  void free();
  void update();
  
  void forwardFFT();
  void backwardFFT();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get number of nodes
  int getNbFFT() const { return ((FFTableMesh *)this->mesh)->getNbLocalFFT(); }

  // get fftw plan id
  inline int getFFTWPlanId() { return this->fftw_plan_id; }
  
  // get one value of frequency domain
  inline fftw_complex & fd(int f);

  // get access directly to frequency domain
  // WARNING: convert it to double (assuming that fftw_complex is double[2])
  inline fftw_complex * fd_storage() { return this->fd_field; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // values in frequency domain in complex form
  fftw_complex * fd_field;

  // fftw plan id (given by mesh)
  int fftw_plan_id;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline fftw_complex & FFTableNodalFieldComponent::fd(int f) {
  return this->fd_field[f];
}


__END_UGUCA__

#endif /* __FFTABLE_NODAL_FIELD_COMPONENT_H__ */
