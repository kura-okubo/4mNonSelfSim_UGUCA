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

#ifdef UCA_USE_OPENMP
#include <omp.h>
#endif /* UCA_USE_OPENMP */

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
HalfSpaceQuasiDynamic::HalfSpaceQuasiDynamic(FFTableMesh & mesh,
					     int side_factor,
					     const std::string & name) :
  HalfSpace(mesh, side_factor, name),
  pi_kernels(mesh) {

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

  this->pi_kernels.preintegrate(*(this->material), Kernel::Krnl::H00,
				this->material->getCs(), this->time_step);
  this->pi_kernels.preintegrate(*(this->material), Kernel::Krnl::H01,
				this->material->getCs(), this->time_step);
  this->pi_kernels.preintegrate(*(this->material), Kernel::Krnl::H11,
				this->material->getCs(), this->time_step);

  if (this->mesh.getDim()==3)
    this->pi_kernels.preintegrate(*(this->material), Kernel::Krnl::H22,
				  this->material->getCs(), this->time_step);
  
#ifdef UCA_VERBOSE
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  printf("Rank %d has total work %d \n",world_rank,total_work);
#endif /* UCA_VERBOSE */
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
void HalfSpaceQuasiDynamic::computeF2D(std::vector<std::complex<double>> & F,
				       double q,
				       std::vector<std::complex<double>> & U,
				       std::complex<double> conv_H00_U0_j,
				       std::complex<double> conv_H01_U0_j,
				       std::complex<double> conv_H01_U1_j,
				       std::complex<double> conv_H11_U1_j) {
  double mu = this->material->getShearModulus();
  double eta = this->material->getCp() / this->material->getCs();
  
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
  double mu = this->material->getShearModulus();
  double eta = this->material->getCp() / this->material->getCs();

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
void HalfSpaceQuasiDynamic::computeStressFourierCoeffQuasiDynamic(bool predicting,
								  bool /*correcting*/) {

  FFTableNodalField & _disp = predicting ? this->disp_pc : this->disp;

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();

  double ** wave_numbers = this->mesh.getLocalWaveNumbers();

  // access to fourier coefficients of stresses
  fftw_complex * internal_fd[3];
  for (int d = 0; d < this->mesh.getDim(); ++d)
    internal_fd[d] = this->internal.fd_storage(d);

  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { // parallel loop over km modes

    // ignore mode 0
    if ((prank == m0_rank) && (j == m0_index)) continue;
    
    std::vector<std::complex<double>> U;
    U.resize(this->mesh.getDim());

    for (int d = 0; d < this->mesh.getDim(); ++d) {
      U[d] = {_disp.fd(d,j)[0], _disp.fd(d,j)[1]};
    }

    double H00_integrated = this->pi_kernels.get(Kernel::Krnl::H00,j)->getIntegral();
    double H01_integrated = this->pi_kernels.get(Kernel::Krnl::H01,j)->getIntegral();
    double H11_integrated = this->pi_kernels.get(Kernel::Krnl::H11,j)->getIntegral();
    
    std::vector<std::complex<double>> F;
    F.resize(this->mesh.getDim());
    if (this->mesh.getDim() == 2) {
      double q = wave_numbers[0][j];
      this->computeF2D(F,q,
		       U,
		       H00_integrated*U[0],
		       H01_integrated*U[0],
		       H01_integrated*U[1],
		       H11_integrated*U[1]);
    } else {
      
      double H22_integrated = this->pi_kernels.get(Kernel::Krnl::H22,j)->getIntegral();

      // q = {k,m} wave number in x,y direction
      double k = wave_numbers[0][j];
      double m = wave_numbers[2][j];
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
    for (int d = 0; d < this->mesh.getDim(); ++d) {
      internal_fd[d][j][0] = std::real(F[d]);  // real part
      internal_fd[d][j][1] = std::imag(F[d]);  // imag part
    }
  }

  // correct mode 0
  if (prank == m0_rank) {
    for (int d = 0; d < this->mesh.getDim(); ++d) {
      internal_fd[d][m0_index][0] = 0.;  // real part
      internal_fd[d][m0_index][1] = 0.;  // imag part
    }
  }
}


__END_UGUCA__
