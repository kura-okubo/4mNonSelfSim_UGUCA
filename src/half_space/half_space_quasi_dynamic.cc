/**
 * @file   half_space_quasi_dynamic.cc
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
#include "half_space_quasi_dynamic.hh"
#include "static_communicator_mpi.hh"
#include "convolutions.hh"

#ifdef UCA_USE_OPENMP
#include <omp.h>
#endif /* UCA_USE_OPENMP */

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
HalfSpaceQuasiDynamic::HalfSpaceQuasiDynamic(Material & material,
					     FFTableMesh & mesh,
					     int side_factor,
					     SpatialDirectionSet components,
					     const std::string & name) :
  HalfSpace(material, mesh, side_factor, components, name),
  convols(mesh) {

#ifdef UCA_VERBOSE
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::cout << "HSQD construct (prank="
	    << prank << "): " << this->name
	    << " : " << this->mesh.getNbLocalFFT() << std::endl;
#endif // UCA_VERBOSE
}

/* -------------------------------------------------------------------------- */
HalfSpaceQuasiDynamic::~HalfSpaceQuasiDynamic() {}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::initConvolutions() {

  this->convols.preintegrate(this->material, Kernel::Krnl::H00,
				  this->material.getCs(), this->time_step);
  this->convols.preintegrate(this->material, Kernel::Krnl::H01,
				  this->material.getCs(), this->time_step);
  this->convols.preintegrate(this->material, Kernel::Krnl::H11,
				  this->material.getCs(), this->time_step);

  if (this->mesh.getDim()==3)
    this->convols.preintegrate(this->material, Kernel::Krnl::H22,
				    this->material.getCs(), this->time_step);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::computeStressFourierCoeff(bool predicting,
						      bool correcting,
						      bool dynamic) {
  if (!dynamic) {
    this->computeStressFourierCoeffQuasiDynamic(predicting, correcting);
  }
  else {
    throw std::runtime_error("HalfSpaceQuasiDynamic cannot compute stress for dynamic problem!");
  }
}

/* -------------------------------------------------------------------------- */
/*void HalfSpaceQuasiDynamic::computeF2D(std::vector<std::complex<double>> & F,
				       double q,
				       std::vector<std::complex<double>> & U,
				       std::complex<double> conv_H00_U0_j,
				       std::complex<double> conv_H01_U0_j,
				       std::complex<double> conv_H01_U1_j,
				       std::complex<double> conv_H11_U1_j) {
  double mu = this->material.getShearModulus();
  double eta = this->material.getCp() / this->material.getCs();
  
  // imaginary number i
  std::complex<double> imag = {0., 1.};


  // - mu * q * int(H00,U0)
  F[0] = -this->side_factor * mu * q * conv_H00_U0_j;

  // + i * (2 - eta) * mu * q * U1
  F[0] += mu * q * (2 - eta) * (imag * U[1]);

  // + i * mu * q * int(H01, U1)
  F[0] += mu * q * imag * conv_H01_U1_j;

  // - mu * q * int(H11, U1)
  F[1] = -this->side_factor * mu * q * conv_H11_U1_j;

  // - i * (2 - etq) * mu * q * U0
  F[1] -= mu * q * (2 - eta) * (imag * U[0]);

  // - i * mu * q * int(H01, U0)
  F[1] -= mu * q * imag * conv_H01_U0_j;
}
*/

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::computeF3D(std::vector<std::complex<double>> & F,
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
				       std::complex<double> conv_H22_U2_j) {
  double mu = this->material.getShearModulus();
  double eta = this->material.getCp() / this->material.getCs();

  // imaginary number i
  std::complex<double> imag = {0., 1.};

  double q = std::sqrt(k * k + m * m);

  F[0] = imag * mu * (2 - eta) * k * U[1];

  F[0] += imag * mu * k * conv_H01_U1_j;

  F[0] -= this->side_factor * mu *
    ((conv_H00_U0_j * ((k * k) / q) + conv_H00_U2_j * ((k * m) / q)) +
     (conv_H22_U0_j * ((m * m) / q) - conv_H22_U2_j * ((k * m) / q)));

  F[1] = -imag * mu * (2 - eta) * (k * U[0] + m * U[2]);

  F[1] -= mu * imag * (conv_H01_U0_j * k + conv_H01_U2_j * m);

  F[1] -= this->side_factor * mu * q * conv_H11_U1_j;

  F[2] = imag * mu * (2 - eta) * m * U[1];

  F[2] += imag * mu * m * conv_H01_U1_j;

  F[2] -= this->side_factor * mu *
    ((conv_H00_U0_j * ((k * m) / q) + conv_H00_U2_j * ((m * m) / q)) +
     (-conv_H22_U0_j * ((k * m) / q) + conv_H22_U2_j * ((k * k) / q)));
}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::computeF(FFTableNodalField & F,
				     const FFTableNodalField & U,
				     const Convolutions::VecComplex & c_H00_U0,
				     const Convolutions::VecComplex & c_H00_U2,
				     const Convolutions::VecComplex & c_H01_U0,
				     const Convolutions::VecComplex & c_H01_U2,
				     const Convolutions::VecComplex & c_H01_U1,
				     const Convolutions::VecComplex & c_H11_U1,
				     const Convolutions::VecComplex & c_H22_U0,
				     const Convolutions::VecComplex & c_H22_U2) {

  // parallel information
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();

  // material properties
  double mu = this->material.getShearModulus();
  double eta = this->material.getCp() / this->material.getCs();

  // wave numbers 
  const TwoDVector & wave_numbers = this->mesh.getLocalWaveNumbers();
  std::vector<double> wn(_spatial_dir_count,0); // to set those zero that do not exist
  
  // imaginary number i
  std::complex<double> imag = {0., 1.};

  // parallel loop over km modes
  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { 

    // get wave numbers
    for (const auto& d : F.getComponents())
       wn[d] = wave_numbers(j,d);
    double k = wn[0]; double m = wn[2]; // for readibility below
    double q = std::sqrt(k * k + m * m);

    // loop over components
    for (const auto& d : F.getComponents()) {

      std::complex<double> F_tmp = {0,0};

      // mode 0
      if ((prank == m0_rank) && (j == m0_index)) {
	// do nothing: mode 0 should have F = (0,0)
      }

      // X component
      else if (d == 0) {
	F_tmp = imag * mu * (2-eta) * k * U.fd_or_zero(j,1);
	F_tmp += imag * mu * k * c_H01_U1[j];
	F_tmp -= this->side_factor * mu *
	  ((c_H00_U0[j] * ((k * k) / q) + c_H00_U2[j] * ((k * m) / q)) +
	   (c_H22_U0[j] * ((m * m) / q) - c_H22_U2[j] * ((k * m) / q)));
      }
      
      // Y component
      else if (d == 1) {
	F_tmp = -imag * mu * (2-eta) * (k * U.fd_or_zero(j,0) + m * U.fd_or_zero(j,2));
	F_tmp -= mu * imag * (c_H01_U0[j] * k + c_H01_U2[j] * m);
	F_tmp -= this->side_factor * mu * q * c_H11_U1[j];
      }
      
      // Z component
      else if (d == 2) {
	F_tmp = imag * mu * (2-eta) * m * U.fd_or_zero(j,1);
	F_tmp += imag * mu * m * c_H01_U1[j];
	F_tmp -= this->side_factor * mu *
	  ((c_H00_U0[j] * ((k * m) / q) + c_H00_U2[j] * ((m * m) / q)) +
	   (-c_H22_U0[j] * ((k * m) / q) + c_H22_U2[j] * ((k * k) / q)));
      }
      
      else {
	throw std::runtime_error("computeF got component (not implemented)");
      }

      // std::complex<double> -> fftw_complex
      F.fd(j,d)[0] = std::real(F_tmp);
      F.fd(j,d)[1] = std::imag(F_tmp);
    }
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpaceQuasiDynamic::computeStressFourierCoeffQuasiDynamic(bool predicting,
								  bool /*correcting*/) {

  FFTableNodalField & _disp = predicting ? this->disp_pc : this->disp;

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();

  const TwoDVector & wave_numbers = this->mesh.getLocalWaveNumbers();

  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { // parallel loop over km modes

    // ignore mode 0
    if ((prank == m0_rank) && (j == m0_index)) continue;
    
    std::vector<std::complex<double>> U;
    U.resize(3); //this->mesh.getDim()); // <-------------------------------- fix
    for (int d=0; d<3; ++d) U[d] = {0,0}; // <-------------------------------- fix
    for (const auto& d : _disp.getComponents()) {
      U[d] = {_disp.fd(j,d)[0], _disp.fd(j,d)[1]};
    }

    double H00_integrated = this->convols.getKernelIntegral(Kernel::Krnl::H00,j);
    double H01_integrated = this->convols.getKernelIntegral(Kernel::Krnl::H01,j);
    double H11_integrated = this->convols.getKernelIntegral(Kernel::Krnl::H11,j);
    
    std::vector<std::complex<double>> F;
    F.resize(3); //this->mesh.getDim()); // <-------------------------------- fix

    if (this->mesh.getDim() == 2) {
      double q = wave_numbers(j,0);
      this->computeF3D(F,q,0,
		       U,
		       H00_integrated*U[0], // H00-U0
		       {0,0}, // H00-U2
		       H01_integrated*U[0], // H01-U0
		       {0,0}, // H01-U2
		       H01_integrated*U[1], // H01-U1
		       H11_integrated*U[1],  // H11-U1
		       {0,0}, // H22-U0
		       {0,0} // H22-U2
		       );
    } else {
      
      double H22_integrated = this->convols.getKernelIntegral(Kernel::Krnl::H22,j);

      // q = {k,m} wave number in x,y direction
      double k = wave_numbers(j,0);
      double m = wave_numbers(j,2);
      this->computeF3D(F,k,m,
		       U,
		       H00_integrated*U[0],
		       H00_integrated*U[2],
		       H01_integrated*U[0],
		       H01_integrated*U[2],
		       H01_integrated*U[1],
		       H11_integrated*U[1],
		       H22_integrated*U[0],
		       H22_integrated*U[2]);
    }
    
    // set values to internal force
    for (const auto& d : this->internal.getComponents()) {
      this->internal.fd(j,d)[0] = std::real(F[d]);  // real part
      this->internal.fd(j,d)[1] = std::imag(F[d]);  // imag part
    }
  }

  // correct mode 0
  if (prank == m0_rank) {
    for (const auto& d : this->internal.getComponents()) {
      this->internal.fd(m0_index,d)[0] = 0.;  // real part
      this->internal.fd(m0_index,d)[1] = 0.;  // imag part
    }
  }
}


__END_UGUCA__
