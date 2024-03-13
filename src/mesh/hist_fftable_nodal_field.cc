/**
 * @file   hist_fftable_nodal_field.cc
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
#include "hist_fftable_nodal_field.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
HistFFTableNodalField::HistFFTableNodalField(FFTableMesh & mesh,
					     SpatialDirectionSet components,
					     const std::string & name) :
  FFTableNodalField(mesh,components,name) {

  // copy start from fftable nodal field
  this->hist_start = this->fd_start;

  // resize the storage
  // needs to be structured as fd_storage
  this->hist_storage.resize(this->hist_start.back());
}

/* -------------------------------------------------------------------------- */
void HistFFTableNodalField::addCurrentValueToHistory() {

  // hist_storage and fd_storage have same size (by construction)
  for (size_t i=0; i<this->hist_storage.size(); ++i) {
    this->hist_storage[i].addCurrentValue(this->fd_storage[i]);
  }
}

/* -------------------------------------------------------------------------- */
void HistFFTableNodalField::addCurrentValueToHistory(FFTableNodalField & other) {

  if ((this->getNbFFT() != other.getNbFFT())
      || (this->getNbComponents() != other.getNbComponents()))
    throw std::runtime_error("HistFFTableNodalField don't match!");
      
  // hist_storage and fd_storage have same size (by construction)
  for (size_t i=0; i<this->hist_storage.size(); ++i) {
    this->hist_storage[i].addCurrentValue(other.fd_data(0)[i]);
  }
}

/* -------------------------------------------------------------------------- */
void HistFFTableNodalField::changeCurrentValueOfHistory() {

  // hist_storage and fd_storage have same size (by construction)
  for (size_t i=0; i<this->hist_storage.size(); ++i) {
    this->hist_storage[i].changeCurrentValue(this->fd_storage[i]);
  }
}

/* -------------------------------------------------------------------------- */
void HistFFTableNodalField::fillHistoryWithCurrentValue() {
  // hist_storage and fd_storage have same size (by construction)
  for (size_t i=0; i<this->hist_storage.size(); ++i) {
    this->hist_storage[i].fillHistory(this->fd_storage[i]);
  }
}

/* -------------------------------------------------------------------------- */
void HistFFTableNodalField::fillHistoryWithCurrentValue(FFTableNodalField & other) {
  if ((this->getNbFFT() != other.getNbFFT())
      || (this->getNbComponents() != other.getNbComponents()))
    throw std::runtime_error("HistFFTableNodalField don't match!");
      
  // hist_storage and fd_storage have same size (by construction)
  for (size_t i=0; i<this->hist_storage.size(); ++i) {
    this->hist_storage[i].fillHistory(other.fd_data(0)[i]);
  }
}

/* -------------------------------------------------------------------------- */
void HistFFTableNodalField::extendHistory(unsigned int size) {
  for (size_t i=0; i<this->hist_storage.size(); ++i) {
    this->hist_storage[i].extend(size);
  }
}

__END_UGUCA__
