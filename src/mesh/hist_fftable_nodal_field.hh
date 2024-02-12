/**
 * @file   hist_fftable_nodal_field.hh
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
#ifndef __HIST_FFTABLE_NODAL_FIELD_H__
#define __HIST_FFTABLE_NODAL_FIELD_H__
/* -------------------------------------------------------------------------- */
#include "fftable_nodal_field.hh"

#include "modal_limited_history.hh"
#include "convolutions.hh"

#include <memory>

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__
class HistFFTableNodalField : public FFTableNodalField {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  typedef std::vector<ModalLimitedHistory> LHVector;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  HistFFTableNodalField(const std::string & name = "unnamed") : FFTableNodalField(name) {}

  HistFFTableNodalField(FFTableMesh & mesh,
			SpatialDirectionSet components = {0},
			const std::string & name = "unnamed");

  virtual ~HistFFTableNodalField() {}

private:
  // private copy constructor: NodalField cannot be copied (for now to debug)
  HistFFTableNodalField(HistFFTableNodalField & to_copy);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // adds current value to history (for all modes and dimensions)
  void addCurrentValueToHistory();

  // add current value of a different vector to history
  // (needed when predicting)
  void addCurrentValueToHistory(FFTableNodalField & other);

  // change current value to history (for all modes and dimensions)
  void changeCurrentValueOfHistory();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get ModalLimitedHistory of frequency domain in direction d
  inline ModalLimitedHistory & hist(int f, int d=0);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  // start indices for each component
  std::vector<int> hist_start;
  
  // past values of field in frequency domain
  // each LimitedHistory is for a given dimension d and wave number q
  LHVector hist_storage;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline ModalLimitedHistory & HistFFTableNodalField::hist(int f, int d) {
  if (!this->components.count(d)) 
    throw std::runtime_error("HistFFTableNodalField "
			     +this->name
			     +" has no component "
			     +std::to_string(d)+"\n");
  // needs to be structured as fd_storage
  return this->hist_storage[this->hist_start[d]+f]; 
}

__END_UGUCA__

#endif /* __HIST_FFTABLE_NODAL_FIELD_H__ */
