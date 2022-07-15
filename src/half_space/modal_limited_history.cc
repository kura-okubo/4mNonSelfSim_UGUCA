/**
 * @file   modal_limited_history.cc
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
#include "modal_limited_history.hh"
#include "preint_kernel.hh"
#include <iterator>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
ModalLimitedHistory::ModalLimitedHistory() :
  nb_history_points(0),
  index_now(0),
  values_real(0), values_imag(0) {}

/* -------------------------------------------------------------------------- */
void ModalLimitedHistory::registerKernel(std::shared_ptr<PreintKernel> pi_kernel) {
  this->pi_kernels.push_back(pi_kernel);
  this->resize();
}

/* -------------------------------------------------------------------------- */
void ModalLimitedHistory::resize() {
  this->resize(this->values_real, false);
  this->resize(this->values_imag, true);
}

/* -------------------------------------------------------------------------- */
void ModalLimitedHistory::resize(std::vector<double> & vec,
				 bool update_index) {
  unsigned int new_size = 0;
  for(auto const& pik: this->pi_kernels) {
    new_size = std::max(new_size, pik->getSize());
  }

  // haven't stored any history yet
  if (this->nb_history_points == 0) {
    vec.resize(new_size);
  }
  else {
    // place history to beginning of temporary vector
    std::vector<double> tmp(vec.size());
    std::copy(vec.begin()+this->index_now, vec.end(),
	      tmp.begin());
    std::copy(vec.begin(), vec.begin()+this->index_now,
	      tmp.begin()+(vec.size()-this->index_now));

    // resize
    int copy_size = std::min((unsigned int)(vec.size()), new_size);
    vec.resize(new_size);

    // copy back
    std::copy(tmp.begin(), tmp.begin()+copy_size, vec.begin());
    if (update_index)
      this->index_now = 0;
  }
}


__END_UGUCA__
