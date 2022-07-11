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
  values(0) {
  //std::fill(this->values.begin(), this->values.end(), 0);
}

/* -------------------------------------------------------------------------- */
void ModalLimitedHistory::registerKernel(const PreintKernel * pi_kernel) {
  this->pi_kernels.push_back(pi_kernel);
}

/* -------------------------------------------------------------------------- */
void ModalLimitedHistory::resize() {
  unsigned int new_size = 0;
  for(auto const& pik: this->pi_kernels) {
    new_size = std::max(new_size, pik->getSize());
  }

  // haven't stored any history yet
  if (this->nb_history_points == 0) {
    this->values.resize(new_size);
  }
  else {
    // place history to beginning of temporary vector
    std::vector<double> tmp(this->values.size());
    std::copy(this->values.begin()+this->index_now, this->values.end(),
	      tmp.begin());
    std::copy(this->values.begin(), this->values.begin()+this->index_now,
	      tmp.begin()+(this->values.size()-this->index_now));

    // resize
    int copy_size = std::min((unsigned int)(this->values.size()), new_size);
    this->values.resize(new_size);

    // copy back
    std::copy(tmp.begin(), tmp.begin()+copy_size, this->values.begin());
    this->index_now = 0;
  }
}


__END_UGUCA__
