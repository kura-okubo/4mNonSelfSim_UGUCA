/**
 * @file   half_space_quasidynamic.hh
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
#ifndef __HALF_SPACE_QUASIDYNAMIC_H__
#define __HALF_SPACE_QUASIDYNAMIC_H__
/* -------------------------------------------------------------------------- */
#include "half_space.hh"
#include "convolutions.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class HalfSpaceQuasiDynamic : public HalfSpace {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  HalfSpaceQuasiDynamic(Material & material, FFTableMesh & mesh, int side_factor,
			const std::string & name = "half_space");

  virtual ~HalfSpaceQuasiDynamic();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // init convolutions
  virtual void initConvolutions();

protected:
  virtual void computeStressFourierCoeff(bool predicting = false,
					 bool correcting = false,
					 bool dynamic = false);

  void computeStressFourierCoeffQuasiDynamic(bool predicting,
					     bool correcting);

  void computeF2D(std::vector<std::complex<double>> & F,
		  double q,
		  std::vector<std::complex<double>> & U,
		  std::complex<double> conv_H00_U0_j,
		  std::complex<double> conv_H01_U0_j,
		  std::complex<double> conv_H01_U1_j,
		  std::complex<double> conv_H11_U1_j);
  
  void computeF3D(std::vector<std::complex<double>> & F,
		  double k,
		  double m,
		  std::vector<std::complex<double>> & U,
		  std::complex<double> conv_H00_U0_j,
		  std::complex<double> conv_H00_U2_j,
		  std::complex<double> conv_H01_U0_j,
		  std::complex<double> conv_H01_U2_j,
		  std::complex<double> conv_H01_U1_j,
		  std::complex<double> conv_H11_U1_j,
		  std::complex<double> conv_H22_U0_j,
		  std::complex<double> conv_H22_U2_j);				   

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  
  // convolution kernels
  Convolutions convolutions; 

};

__END_UGUCA__

//#include "half_space_quasidynamic_impl.cc"

#endif /* __HALF_SPACE_QUASIDYNAMIC_H__ */
