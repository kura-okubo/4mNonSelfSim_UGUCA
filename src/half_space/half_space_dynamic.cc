/**
 * @file   half_space_dynamic.cc
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
#include "half_space_dynamic.hh"
#include "static_communicator_mpi.hh"

#ifdef UCA_USE_OPENMP
#include <omp.h>
#endif /* UCA_USE_OPENMP */

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
HalfSpaceDynamic::HalfSpaceDynamic(FFTableMesh & mesh,
				   int side_factor,
				   const std::string & name) :
  HalfSpaceQuasiDynamic(mesh, side_factor, name) {
  
  // allocated pre-integrated kernels
  this->H00_pi.resize(this->mesh.getNbLocalFFT());
  this->H01_pi.resize(this->mesh.getNbLocalFFT());
  this->H11_pi.resize(this->mesh.getNbLocalFFT());
  if (this->mesh.getDim()==3)
    this->H22_pi.resize(this->mesh.getNbLocalFFT());

#ifdef UCA_VERBOSE
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::cout << "HSD construct (prank="
	    << prank << "): " << this->name
	    << " : " << this->mesh.getNbLocalFFT() << std::endl;
#endif // UCA_VERBOSE
  
  // allocate displacement history
  this->U_r.resize(this->mesh.getDim());
  this->U_i.resize(this->mesh.getDim());
  for (int d=0; d<this->mesh.getDim(); ++d) {
    this->U_r[d].resize(this->mesh.getNbLocalFFT());
    this->U_i[d].resize(this->mesh.getNbLocalFFT());
  }
}

/* -------------------------------------------------------------------------- */
HalfSpaceDynamic::~HalfSpaceDynamic() {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();
  
  for (int d=0; d<this->mesh.getDim(); ++d) {
    for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) {

      // ignore mode 0
      if ((prank == m0_rank) && (j == m0_index)) continue;

      delete this->U_r[d][j];
      delete this->U_i[d][j];
    }
  }

  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) {

    // ignore mode 0
    if ((prank == m0_rank) && (j == m0_index)) continue;

    delete this->H00_pi[j];
    delete this->H01_pi[j];
    delete this->H11_pi[j];
    if (this->mesh.getDim()==3)
      delete this->H22_pi[j];
  }
}

/* -------------------------------------------------------------------------- */
double HalfSpaceDynamic::getStableTimeStep() {
  double delta_x = this->mesh.getDeltaX();
  double delta_z = this->mesh.getDeltaZ();

  if (this->mesh.getDim()==2)
    delta_z = 1e100;
  return std::min(delta_x,delta_z) / this->material->getCs();
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::setTimeStep(double time_step) {
  if (((time_step/this->getStableTimeStep()>0.35) && (this->mesh.getDim()==3)) ||
      ((time_step/this->getStableTimeStep()>0.4) && (this->mesh.getDim()==2)))
    throw std::runtime_error("Error: time_step_factor is too large: (<0.4 for 2D and <0.35 for 3D)\n");

  HalfSpace::setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::initConvolutions() {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();
  
  double ** wave_numbers = this->mesh.getLocalWaveNumbers();

  // history for q1 is longest q = j*q1
  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { //parallel loop

    // ignore mode 0
    if ((prank == m0_rank) && (j == m0_index)) continue;
    
    this->H00_pi[j] = new PreintKernel(this->material->getH00());
    this->H01_pi[j] = new PreintKernel(this->material->getH01());
    this->H11_pi[j] = new PreintKernel(this->material->getH11());
    if (this->mesh.getDim()==3)
      this->H22_pi[j] = new PreintKernel(this->material->getH22());

    double qq = 0.0;
    for (int d=0; d<this->mesh.getDim();d+=2)
      qq +=(wave_numbers[d])[j]*(wave_numbers[d])[j];

    double qj_cs = std::sqrt(qq) * this->material->getCs();
    
    this->H00_pi[j]->preintegrate(qj_cs, this->time_step);
    this->H01_pi[j]->preintegrate(qj_cs, this->time_step);
    this->H11_pi[j]->preintegrate(qj_cs, this->time_step);
    if (this->mesh.getDim()==3)
      this->H22_pi[j]->preintegrate(qj_cs, this->time_step);

    // register ModalLimitedHistory
    for (int d=0; d<this->mesh.getDim(); ++d) {
      this->U_r[d][j] = new ModalLimitedHistory();
      this->U_i[d][j] = new ModalLimitedHistory();
    }

    // to be corrected when put into a collection
    this->U_r[0][j]->registerKernel(this->H00_pi[j]);
    this->U_r[0][j]->registerKernel(this->H01_pi[j]);
    this->U_i[0][j]->registerKernel(this->H00_pi[j]);
    this->U_i[0][j]->registerKernel(this->H01_pi[j]);

    this->U_r[1][j]->registerKernel(this->H11_pi[j]);
    this->U_r[1][j]->registerKernel(this->H01_pi[j]);
    this->U_i[1][j]->registerKernel(this->H11_pi[j]);
    this->U_i[1][j]->registerKernel(this->H01_pi[j]);

    if (this->mesh.getDim()==3) {
      this->U_r[0][j]->registerKernel(this->H22_pi[j]);
      this->U_i[0][j]->registerKernel(this->H22_pi[j]);

      this->U_r[2][j]->registerKernel(this->H00_pi[j]);
      this->U_r[2][j]->registerKernel(this->H01_pi[j]);
      this->U_r[2][j]->registerKernel(this->H22_pi[j]);
      this->U_i[2][j]->registerKernel(this->H00_pi[j]);
      this->U_i[2][j]->registerKernel(this->H01_pi[j]);
      this->U_i[2][j]->registerKernel(this->H22_pi[j]);
    }

    for (int d=0; d<this->mesh.getDim(); ++d) {
      this->U_r[d][j]->resize();
      this->U_i[d][j]->resize();
    }
  }

#ifdef UCA_VERBOSE
  int total_work=0;
  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { //parallel loop
    // ignore mode 0
    if ((prank == m0_rank) && (j == m0_index)) continue;
    for (int d=0; d<this->mesh.getDim(); ++d)
      total_work += this->U_r[d][j]->getSize();
  }
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  printf("Rank %d has total work %d \n",world_rank,total_work);
#endif /* UCA_VERBOSE */
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::computeStressFourierCoeff(bool predicting,
						 bool correcting) {
  this->computeStressFourierCoeffDynamic(predicting, correcting);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::computeStressFourierCoeffDynamic(bool predicting,
							bool correcting) {

  FFTableNodalField & _disp = predicting ? this->disp_pc : this->disp;

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();

  double mu = this->material->getShearModulus();
  double eta = this->material->getCp() / this->material->getCs();

  double ** wave_numbers = this->mesh.getLocalWaveNumbers();

  // imaginary number i
  std::complex<double> imag = {0., 1.};

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

      // store current displacement in history
      if (correcting) {
        this->U_r[d][j]->changeCurrentValue(std::real(U[d]));
        this->U_i[d][j]->changeCurrentValue(std::imag(U[d]));
      } else {
        this->U_r[d][j]->addCurrentValue(std::real(U[d]));
        this->U_i[d][j]->addCurrentValue(std::imag(U[d]));
      }
    }

    int nb_conv = (int)std::pow(2.0, this->mesh.getDim());

    std::vector<std::complex<double>> conv;
    conv.resize(nb_conv);

    std::vector<PreintKernel *> krnl = {
      this->H00_pi[j],
      this->H01_pi[j],
      this->H01_pi[j],
      this->H11_pi[j]};

    if (this->mesh.getDim() == 3) {  // add 3d part
      std::vector<PreintKernel *> krnl3d = {
	this->H00_pi[j],
	this->H01_pi[j],
	this->H22_pi[j],
	this->H22_pi[j]};

      krnl.insert(krnl.end(), krnl3d.begin(), krnl3d.end());
    }

    std::vector<ModalLimitedHistory *> U_real = {
      this->U_r[0][j],
      this->U_r[0][j],
      this->U_r[1][j],
      this->U_r[1][j]};

    if (this->mesh.getDim() == 3) {  // add 3d part
      std::vector<ModalLimitedHistory *> U_real3d = {
	this->U_r[2][j],
	this->U_r[2][j],
	this->U_r[0][j],
	this->U_r[2][j]};

      U_real.insert(U_real.end(), U_real3d.begin(), U_real3d.end());
    }

    std::vector<ModalLimitedHistory *> U_imag = {
      this->U_i[0][j],
      this->U_i[0][j],
      this->U_i[1][j],
      this->U_i[1][j]};

    if (this->mesh.getDim() == 3) {  // add 3d part
      std::vector<ModalLimitedHistory *> U_imag3d = {
	this->U_i[2][j],
	this->U_i[2][j],
	this->U_i[0][j],
	this->U_i[2][j]};

      U_imag.insert(U_imag.end(), U_imag3d.begin(), U_imag3d.end());
    }

#ifdef UCA_USE_OPENMP
#pragma omp parallel for
#endif
    for (int n = 0; n < nb_conv; ++n) {
      conv[n] = krnl[n]->convolve(U_real[n], U_imag[n]);
    }
    // convolutions for both 2d and 3d
    std::complex<double> conv_H00_U0_j = conv[0];
    std::complex<double> conv_H01_U0_j = conv[1];
    std::complex<double> conv_H01_U1_j = conv[2];
    std::complex<double> conv_H11_U1_j = conv[3];

    std::vector<std::complex<double>> F;
    F.resize(this->mesh.getDim());
    if (this->mesh.getDim() == 2) {
      double q = wave_numbers[0][j];

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
    } else {
      // convolutions for 3d only
      std::complex<double> conv_H00_U2_j = conv[4];
      std::complex<double> conv_H01_U2_j = conv[5];
      std::complex<double> conv_H22_U0_j = conv[6];
      std::complex<double> conv_H22_U2_j = conv[7];
      // q = {k,m} wave number in x,y direction

      double k = wave_numbers[0][j];
      double m = wave_numbers[2][j];

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

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::registerToRestart(Restart & restart) {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();
  
  for (int d=0; d<this->mesh.getDim(); ++d) {
    for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) {

      // ignore mode 0
      if ((prank == m0_rank) && (j == m0_index)) continue;

      restart.registerIO(this->name+"_Ur_"+std::to_string(d)+"_"+std::to_string(j),
			 *(this->U_r[d][j]));
      restart.registerIO(this->name+"_Ui_"+std::to_string(d)+"_"+std::to_string(j),
			 *(this->U_i[d][j]));
    }
  }

  HalfSpace::registerToRestart(restart);
}


__END_UGUCA__
