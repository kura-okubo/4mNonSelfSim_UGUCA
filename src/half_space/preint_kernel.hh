/**
 * @file   preint_kernel.hh
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
#ifndef __PREINT_KERNEL_H__
#define __PREINT_KERNEL_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "kernel.hh"
#include "modal_limited_history.hh"

__BEGIN_UGUCA__

class PreintKernel {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PreintKernel(const Kernel * kernel);
  virtual ~PreintKernel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // compute pre-integration of kernel (only kernel: no other factors)
  // e.g., integral of H11(q*c_s*t) with time_factor q*c_s
  void preintegrate(double time_factor,
		    double time_step);

  // muliply preintegrated kernel by factor
  void multiplyBy(double factor);

  // compute convolution of kernel with a history
  std::complex<double> convolve(const ModalLimitedHistory & U);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  unsigned int getSize() const { return this->values.size(); }

  // get direct access to values (only used to testing)
  std::vector<double> & getValues() { return this->values; }

  // get entire integral of kernel
  double getIntegral() const { return this->integral; }
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  const Kernel * kernel;

  // integral of trapezoids with a constant time step
  std::vector<double> values;

  // full integral of kernel
  double integral;
};

__END_UGUCA__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */


#endif /* __PREINT_KERNEL_H__ */
