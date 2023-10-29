/**
 * @file   limited_history.hh
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Sun Jul 10 2022
 * @date last modification: Sun Jul 10 2022
 *
 * @brief  Limited history for all dimensions and modes
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
#ifndef __LIMITED_HISTORY_H__
#define __LIMITED_HISTORY_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_fftable_mesh.hh"
#include "modal_limited_history.hh"
#include "convolutions.hh"

#include <memory>

__BEGIN_UGUCA__

class LimitedHistory {
  
  friend class BaseIO;
  
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  typedef std::vector<std::shared_ptr<ModalLimitedHistory>> LHVector;
  
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  LimitedHistory(FFTableMesh & mesh);
  virtual ~LimitedHistory() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // registers Kernel to all modes of a history
  void registerKernel(Convolutions::PIKernelVector & pi_kernels,
		      unsigned int dim);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<ModalLimitedHistory> get(unsigned int dim,
					   unsigned int wave_number) {
    return this->history[dim*this->nbfft+wave_number];
  }
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // copy of spatial dimension for fast use in get accessor
  unsigned int dimension;

  // copy of nbfft for fast useful
  unsigned int nbfft;
  
  // past values of field in frequency domain
  // each LimitedHistory is for a given dimension d and wave number q
  LHVector history;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */


__END_UGUCA__

#endif /* __LIMITED_HISTORY_H__ */
