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
#include "preint_kernel_collection.hh"
#include "static_communicator_mpi.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
PreintKernelCollection::PreintKernelCollection(FFTableMesh & mesh) :
  mesh(mesh) {
  
}

/* -------------------------------------------------------------------------- */
PreintKernelCollection::~PreintKernelCollection() {

}

/* -------------------------------------------------------------------------- */
void PreintKernelCollection::preintegrate(Material & material,
					  Kernel::Krnl kernel,
					  double scale_factor,
					  double time_step) {

  // create vector and resize it
  this->pi_kernels.insert(std::pair<Kernel::Krnl,PIKernelVector>(kernel, PIKernelVector()));
  PIKernelVector & pik_vector = this->pi_kernels[kernel];
  pik_vector.resize(this->mesh.getNbLocalFFT());
  
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();
  
  double ** wave_numbers = this->mesh.getLocalWaveNumbers();

  // history for q1 is longest q = j*q1
  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { //parallel loop

    // ignore mode 0
    if ((prank == m0_rank) && (j == m0_index)) continue;

    pik_vector[j] = std::make_shared<PreintKernel>(material.getKernel(kernel));

    double qq = 0.0;
    for (int d=0; d<this->mesh.getDim();d+=2)
      qq +=(wave_numbers[d])[j]*(wave_numbers[d])[j];

    double qj_cs = std::sqrt(qq) * scale_factor;
    
    pik_vector[j]->preintegrate(qj_cs, time_step);
  }
}

__END_UGUCA__
