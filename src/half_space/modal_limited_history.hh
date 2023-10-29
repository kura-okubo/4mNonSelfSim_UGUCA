/**
 * @file   modal_limited_history.hh
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
#ifndef __MODAL_LIMITED_HISTORY_H__
#define __MODAL_LIMITED_HISTORY_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <complex>

__BEGIN_UGUCA__

class PreintKernel;

class ModalLimitedHistory {

  friend class BaseIO;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ModalLimitedHistory();
  virtual ~ModalLimitedHistory() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // at the current value of the history
  inline void addCurrentValue(std::complex<double> value);
  inline void changeCurrentValue(std::complex<double> value);

  inline void setSteadyState(std::complex<double> value);
  
  // get history value at index with index=0 : now
  inline std::complex<double> at(unsigned int index) const;

  // update to take account for newly added preint kernels
  void resize();
  
  // register preintegrated kernel
  void registerKernel(std::shared_ptr<PreintKernel> pi_kernel);

private:
  void resize(std::vector<double> & vec, bool update_index);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  unsigned int getSize() const { return this->values_real.size(); };
  unsigned int getNbHistoryPoints() const { return std::min(this->nb_history_points,
							    this->values_real.size()); };
  unsigned int getIndexNow() const {return this->index_now; }
  const double * real() const { return this->values_real.data(); }
  const double * imag() const { return this->values_imag.data(); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // number of accumulated history points
  std::vector<double>::size_type nb_history_points;

  // index pointing to the newest entry
  unsigned int index_now;

  // values (keep in separate vectors for BLAS in preint_kernel convolution)
  std::vector<double> values_real;
  std::vector<double> values_imag;

  // preintegrated kernels that use this limited history
  std::vector<std::shared_ptr<PreintKernel>> pi_kernels;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline void ModalLimitedHistory::addCurrentValue(std::complex<double> value) {

  if (this->index_now == 0)
    this->index_now = this->values_real.size();

  this->index_now -= 1;

  this->values_real[this->index_now] = std::real(value);
  this->values_imag[this->index_now] = std::imag(value);
  
  // increase the counter of history points
  this->nb_history_points = std::min(this->nb_history_points + 1,
				     this->values_real.size());
}

/* -------------------------------------------------------------------------- */
inline void ModalLimitedHistory::changeCurrentValue(std::complex<double> value) {
  this->values_real[this->index_now] = std::real(value);
  this->values_imag[this->index_now] = std::imag(value);
}

/* -------------------------------------------------------------------------- */
inline void ModalLimitedHistory::setSteadyState(std::complex<double> value) {
  this->nb_history_points = this->values_real.size();
  std::fill(this->values_real.begin(), this->values_real.end(), std::real(value));
  std::fill(this->values_imag.begin(), this->values_imag.end(), std::imag(value));
}

/* -------------------------------------------------------------------------- */
inline std::complex<double> ModalLimitedHistory::at(unsigned int index) const {
  if (index >= this->values_real.size()) {
    std::cerr << "try to access history value beyond existence" << std::endl;
    throw index;
  }

  unsigned int i = (this->index_now + index) % this->values_real.size();
  return {this->values_real[i], this->values_imag[i]};
}

__END_UGUCA__

#endif /* __MODAL_LIMITED_HISTORY_H__ */
