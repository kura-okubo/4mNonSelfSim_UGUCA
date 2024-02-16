/**
 * @file   preint_kernel_collection.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Thu Jul 14 2022
 * @date last modification: Thu Jul 14 2022
 *
 * @brief  Contains PreintKernels for each type of kernel and each mode
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
#include "convolutions.hh"

#ifdef UCA_USE_OPENMP
#include <omp.h>
#endif /* UCA_USE_OPENMP */

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
Convolutions::Convolutions(HistFFTableNodalField & field) :
  field(field) {

  // zeros vector
  VecComplex tmp(this->field.getNbFFT());
  for (auto& element : tmp) {
    element = {0.,0.};
  }
  this->complex_zeros.swap(tmp);
}

/* -------------------------------------------------------------------------- */
Convolutions::~Convolutions() {}

/* -------------------------------------------------------------------------- */
void Convolutions::preintegrate(Material & material,
				Krnl kernel,
				double scale_factor,
				double time_step) {

  // preintegrated kernels: create vector and resize it
  this->pi_kernels.insert(std::pair<Krnl,PIKernelVector>(kernel, PIKernelVector()));
  PIKernelVector & pik_vector = this->pi_kernels[kernel];
  pik_vector.resize(this->field.getNbFFT());

  const TwoDVector & wave_numbers = this->field.getMesh().getLocalWaveNumbers();

  // history for q1 is longest q = j*q1
  for (int j=0; j<this->field.getNbFFT(); ++j) { //parallel loop

    pik_vector[j] = std::make_shared<PreintKernel>(material.getKernel(kernel));

    double qq = 0.0;
    for (const auto& d : this->field.getComponents())
      qq += wave_numbers(j,d)*wave_numbers(j,d);

    double qj_cs = std::sqrt(qq) * scale_factor;
    
    pik_vector[j]->preintegrate(qj_cs, time_step);
  }
}

/* -------------------------------------------------------------------------- */
void Convolutions::init(ConvPair conv) {
  
  // make sure that the history of the field has the necessary length
  for (int j=0; j<this->field.getNbFFT(); ++j) {
    for (const auto& d : this->field.getComponents()) {
      this->field.hist(j,d).extend(this->pi_kernels[conv.first][j]->getSize());
    }
  }
  
  // prepare results
  int vec_size = this->field.getNbFFT();
  this->results.insert(std::pair<ConvPair,VecComplex>(conv,VecComplex(vec_size)));
}

/* -------------------------------------------------------------------------- */
void Convolutions::convolve() {

  ConvMap::iterator it;

#ifdef UCA_USE_OPENMP
#pragma omp parallel for
#pragma omp single nowait
#endif
  for(it=this->results.begin(); it!=this->results.end(); ++it) {
#ifdef UCA_USE_OPENMP
#pragma omp task firstprivate(datIt)
#endif
    // history for q1 is longest q = j*q1
    for (int j=0; j<this->field.getNbFFT(); ++j) { //parallel loop

      // kernel 
      Krnl kernel = it->first.first;

      // modal U
      unsigned int U_dim = it->first.second;
      const ModalLimitedHistory & U_j = this->field.hist(j,U_dim);
     
      // std::vector<std::complex<double>> & res = it->second;
      it->second[j] = this->pi_kernels[kernel][j]->convolve(U_j);
    }
  }
}

/* -------------------------------------------------------------------------- */
void Convolutions::convolveSteadyState() {

  // loop over Kernel-dim pairs
  for(ConvMap::iterator it=this->results.begin();
      it!=this->results.end(); ++it) {

    // loop over modes
    for (int j=0; j<this->field.getNbFFT(); ++j) {

      // kernel 
      Krnl kernel = it->first.first;

      // modal U (current value)
      unsigned int U_dim = it->first.second;
      std::complex<double> U_j = this->field.fd_or_zero(j,U_dim);

      // compute convolution with time-constant U_j
      it->second[j] = this->pi_kernels[kernel][j]->getIntegral() * U_j;
    }
  }
}

/* -------------------------------------------------------------------------- */
const VecComplex & Convolutions::getResult(ConvPair pair) const {

  auto it = this->results.find(pair); // attempt to find the result

  if (it != this->results.end()) {
    return it->second; // found result, return it
  }
  else {
    return this->complex_zeros; // not found, return vector full of complex zeros
  }
}

__END_UGUCA__
