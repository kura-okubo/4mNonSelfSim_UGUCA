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
HalfSpaceDynamic::HalfSpaceDynamic(Material & material,
				   FFTableMesh & mesh,
				   int side_factor,
				   SpatialDirectionSet components,
				   const std::string & name) :
  HalfSpace(material, mesh, side_factor, components, name),
  convols(this->disp) {
  
#ifdef UCA_VERBOSE
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::cout << "HSD construct (prank="
	    << prank << "): " << this->name
	    << " : " << this->mesh.getNbLocalFFT() << std::endl;
#endif // UCA_VERBOSE
}

/* -------------------------------------------------------------------------- */
HalfSpaceDynamic::~HalfSpaceDynamic() {}

/* -------------------------------------------------------------------------- */
double HalfSpaceDynamic::getStableTimeStep() {
  double delta_x = this->mesh.getDelta(0);
  double delta_z = this->mesh.getDelta(2);

  if (this->mesh.getDim()==2)
    delta_z = delta_x;
  return std::min(delta_x,delta_z) / this->material.getCs();
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::setTimeStep(double time_step) {
  if (((time_step/this->getStableTimeStep()>0.35) && (this->mesh.getDim()==3)) ||
      ((time_step/this->getStableTimeStep()>0.4) && (this->mesh.getDim()==2)))
    throw std::runtime_error("Error: time_step_factor is too large: (<0.4 for 2D and <0.35 for 3D)\n");
  
  HalfSpace::setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::preintegrateKernels() {

  // for x and y components
  if (this->disp.getComponents().count(_x) &&
      this->disp.getComponents().count(_y)) {
    this->convols.preintegrate(this->material, Krnl::H00,
			       this->material.getCs(), this->time_step);
    this->convols.preintegrate(this->material, Krnl::H01,
			       this->material.getCs(), this->time_step);
    this->convols.preintegrate(this->material, Krnl::H11,
			       this->material.getCs(), this->time_step);
  }

  // for z components
  if (this->disp.getComponents().count(_z)) {
    this->convols.preintegrate(this->material, Krnl::H22,
			       this->material.getCs(), this->time_step);
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::initConvolutions() {

  this->preintegrateKernels();

  // in-plane problem
  if (this->disp.getComponents().count(_x) &&
      this->disp.getComponents().count(_y)) {
    this->convols.init(std::make_pair(Krnl::H00,0)); // H00-U0
    this->convols.init(std::make_pair(Krnl::H01,0)); // H01-U0
    this->convols.init(std::make_pair(Krnl::H11,1)); // H11-U1
    this->convols.init(std::make_pair(Krnl::H01,1)); // H01-U1
  }

  // out-of-plane problem
  if (this->disp.getComponents().count(_z)) {
    this->convols.init(std::make_pair(Krnl::H22,2)); // H22-U2
  }

  // 3D problem
  if (this->disp.getComponents().count(_x) &&
      this->disp.getComponents().count(_y) &&
      this->disp.getComponents().count(_z)) {
    this->convols.init(std::make_pair(Krnl::H22,0)); // H22-U0
    this->convols.init(std::make_pair(Krnl::H00,2)); // H00-U2
    this->convols.init(std::make_pair(Krnl::H01,2)); // H01-U2
  }    

  // history of mode 0 needs to be of length zero
  // otherwise convolution fails
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();
  if (prank == m0_rank) {
    for (const auto& d : this->disp.getComponents())
      this->disp.hist(m0_index,d).resize(0);
  }
  
#ifdef UCA_VERBOSE
  int total_work=0;
  for (int j=0; j<this->disp.getNbFFT(); ++j) { //parallel loop
    for (const auto& d : this->disp.getComponents())
      total_work += this->disp.hist(j,d).getSize();
  }
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  printf("Rank %d has total work %d \n",world_rank,total_work);
#endif /* UCA_VERBOSE */
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::computeStressFourierCoeff(bool predicting,
						 bool correcting,
						 bool dynamic) {
  if (dynamic)
    this->computeStressFourierCoeffDynamic(predicting, correcting);
  else
    throw std::runtime_error("HalfSpaceDynamic cannot compute stress for quasi dynamic problem!");
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::computeStressFourierCoeffDynamic(bool predicting,
							bool correcting) {

  // update history
  if (predicting) {
    this->disp.addCurrentValueToHistory(this->disp_pc);
  }
  else {
    this->disp.addCurrentValueToHistory();
  }
  if (correcting) {
    this->disp.changeCurrentValueOfHistory();
  }
  
  // compute convolutions
  this->convols.convolve();
    
  // compute F
  FFTableNodalField & _disp = predicting ? this->disp_pc : this->disp;
  this->computeF(this->internal, _disp, this->convols);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::computeF(FFTableNodalField & F,
				const FFTableNodalField & U,
				const Convolutions & cnvls) {

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
  for (int j=0; j<F.getNbFFT(); ++j) { 

    // get wave numbers
    for (int d=0; d<this->mesh.getDim(); ++d)
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
	const VecComplex & c_H00_U0 = cnvls.getResult(std::make_pair(Krnl::H00,0));
	const VecComplex & c_H00_U2 = cnvls.getResult(std::make_pair(Krnl::H00,2));
	const VecComplex & c_H01_U1 = cnvls.getResult(std::make_pair(Krnl::H01,1));
	const VecComplex & c_H22_U0 = cnvls.getResult(std::make_pair(Krnl::H22,0));
	const VecComplex & c_H22_U2 = cnvls.getResult(std::make_pair(Krnl::H22,2));
	  
	F_tmp = imag * mu * (2-eta) * k * U.fd_or_zero(j,1);
	F_tmp += imag * mu * k * c_H01_U1[j];
	F_tmp -= this->side_factor * mu *
	  ((c_H00_U0[j] * ((k * k) / q) + c_H00_U2[j] * ((k * m) / q)) +
	   (c_H22_U0[j] * ((m * m) / q) - c_H22_U2[j] * ((k * m) / q)));
      }
      
      // Y component
      else if (d == 1) {
	const VecComplex & c_H01_U0 = cnvls.getResult(std::make_pair(Krnl::H01,0));
	const VecComplex & c_H01_U2 = cnvls.getResult(std::make_pair(Krnl::H01,2));
	const VecComplex & c_H11_U1 = cnvls.getResult(std::make_pair(Krnl::H11,1));
	
	F_tmp = -imag * mu * (2-eta) * (k * U.fd_or_zero(j,0) + m * U.fd_or_zero(j,2));
	F_tmp -= mu * imag * (c_H01_U0[j] * k + c_H01_U2[j] * m);
	F_tmp -= this->side_factor * mu * q * c_H11_U1[j];
      }
      
      // Z component
      else if (d == 2) {
	const VecComplex & c_H00_U0 = cnvls.getResult(std::make_pair(Krnl::H00,0));
	const VecComplex & c_H00_U2 = cnvls.getResult(std::make_pair(Krnl::H00,2));
	const VecComplex & c_H01_U1 = cnvls.getResult(std::make_pair(Krnl::H01,1));
	const VecComplex & c_H22_U0 = cnvls.getResult(std::make_pair(Krnl::H22,0));
	const VecComplex & c_H22_U2 = cnvls.getResult(std::make_pair(Krnl::H22,2));
	
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
void HalfSpaceDynamic::registerToRestart(Restart & restart) {
  restart.registerIO(this->name+"_U",this->disp);
  HalfSpace::registerToRestart(restart);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::setSteadyState(bool predicting) {
  if (predicting)
    this->disp.fillHistoryWithCurrentValue(this->disp_pc);
  else
    this->disp.fillHistoryWithCurrentValue();
}

__END_UGUCA__
