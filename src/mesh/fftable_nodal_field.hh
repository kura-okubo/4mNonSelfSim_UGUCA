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
#include "uca_common.hh"
#include "nodal_field.hh"
#include "fftable_nodal_field_component.hh"
#include "uca_fftable_mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

class FFTableNodalField : public NodalField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  FFTableNodalField(const std::string & name = "unnamed") : NodalField(name) {}
  FFTableNodalField(FFTableMesh & mesh,
		    const std::string & name = "unnamed");

  virtual ~FFTableNodalField() {}

private:
  // private copy constructor: NodalField cannot be copied (for now to debug)
  FFTableNodalField(FFTableNodalField & to_copy);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void init(FFTableMesh & mesh);

  void forwardFFT();
  void backwardFFT();

protected:
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  inline FFTableNodalFieldComponent & component(int i) {
    return (FFTableNodalFieldComponent&)(*this->field[i]);
  }

  inline fftw_complex * fd_storage(int d) {
    return ((FFTableNodalFieldComponent*)(this->field[d]))->fd_storage();
  }
  
  // get one value of frequency domain in direction d
  inline fftw_complex & fd(int d, int f) {
    return ((FFTableNodalFieldComponent*)(this->field[d]))->fd(f);
  }
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

};

__END_UGUCA__

#endif /* __FFTABLE_NODAL_FIELD_H__ */
