/**
 * @file   convolutions.hh
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
#ifndef __CONVOLUTIONS_H__
#define __CONVOLUTIONS_H__
/* -------------------------------------------------------------------------- */

#include "uca_common.hh"
#include "material.hh"
#include "uca_fftable_mesh.hh"
#include "preint_kernel.hh"
//#include "limited_history.hh"

#include <map>
#include <memory>

__BEGIN_UGUCA__

class LimitedHistory;

class Convolutions {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef std::vector<std::shared_ptr<PreintKernel>> PIKernelVector;
protected:
  typedef std::map<Kernel::Krnl,PIKernelVector> PIKernelMap;
  typedef std::pair<Kernel::Krnl,unsigned int> ConvPair;
  typedef std::map<ConvPair,std::vector<std::complex<double>>> ConvMap;
  
  
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Convolutions(FFTableMesh & mesh);

  virtual ~Convolutions();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  // preintegrate kernels
  void preintegrate(Material & material, Kernel::Krnl kernel,
		    double scale_factor, double time_step);

  void registerHistory(LimitedHistory & limited_history) {
    this->field = &(limited_history);
  }
  
  // initialize a convolution computation
  void init(ConvPair conv);

  // get convolution results
  void convolve();
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  // hopefully don't need this -> replace this by a
  // getIntegral(kernel, mode)
  std::shared_ptr<PreintKernel> get(Kernel::Krnl kernel, int mode) { return this->pi_kernels[kernel][mode]; }
  
  const ConvMap & getResults() { return this->results; }
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // reference to mesh
  FFTableMesh & mesh;

  // pointer to limited history to convolve kernels with
  LimitedHistory * field;
  
  // preintegrated kernels [kernel][mode]
  // e.g., pi_kernels[Kernel::Krnl::H00][1]
  PIKernelMap pi_kernels;

  /// convolution results
  ConvMap results;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

__END_UGUCA__

//#include "convolutions_impl.cc"

#endif /* __CONVOLUTIONS_H__ */
