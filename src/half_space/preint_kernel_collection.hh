/**
 * @file   preint_kernel_collection.hh
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
#ifndef __PREINT_KERNEL_COLLECTION_H__
#define __PREINT_KERNEL_COLLECTION_H__
/* -------------------------------------------------------------------------- */

#include "uca_common.hh"
#include "uca_fftable_mesh.hh"
#include "preint_kernel.hh"
#include "material.hh"

#include <map>
#include <memory>

__BEGIN_UGUCA__

class PreintKernelCollection {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  typedef std::vector<std::shared_ptr<PreintKernel>> PIKernelVector;
  typedef std::map<Kernel::Krnl,PIKernelVector> PIKernelMap;

  
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  PreintKernelCollection(FFTableMesh & mesh);

  virtual ~PreintKernelCollection();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // preintegrate kernels
  void preintegrate(Material & material, Kernel::Krnl kernel,
		    double scale_factor, double time_step);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  // hopefully don't need this
  std::shared_ptr<PreintKernel> get(Kernel::Krnl kernel, int mode) { return this->pi_kernels[kernel][mode]; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // reference to mesh
  FFTableMesh & mesh;
  
  // preintegrated kernels [kernel][mode]
  // e.g., pi_kernels[Kernel::Krnl::H00][1]
  PIKernelMap pi_kernels;
  
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

__END_UGUCA__

//#include "preint_kernel_collection_impl.cc"

#endif /* __PREINT_KERNEL_COLLECTION_H__ */
