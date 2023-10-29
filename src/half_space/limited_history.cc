/**
 * @file   limited_history.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Sun Jul 10 2022
 * @date last modification: Sun Jul 10 2022
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
#include "limited_history.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
LimitedHistory::LimitedHistory(FFTableMesh & mesh) :
  dimension(mesh.getDim()), nbfft(mesh.getNbLocalFFT()) {

  this->history.resize(this->dimension*this->nbfft);
  for (unsigned int i=0; i<this->dimension*this->nbfft; ++i) {
    this->history[i] = std::make_shared<ModalLimitedHistory>();
  }
}

/* -------------------------------------------------------------------------- */
void LimitedHistory::registerKernel(Convolutions::PIKernelVector & pi_kernels,
				    unsigned int dim) {

  if (pi_kernels.size() != this->nbfft)
    throw std::runtime_error("incorrect number of pi_kernels provided for register in LimitedHistory");
  
  for (unsigned int j=0; j<this->nbfft; ++j) {
    this->get(dim, j)->registerKernel(pi_kernels[j]);
  }
}

__END_UGUCA__
