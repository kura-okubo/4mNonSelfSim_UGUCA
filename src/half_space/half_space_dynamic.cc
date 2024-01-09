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

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
HalfSpaceDynamic::HalfSpaceDynamic(Material & material,
				   FFTableMesh & mesh,
				   int side_factor,
				   SpatialDirectionSet components,
				   const std::string & name) :
  HalfSpaceQuasiDynamic(material, mesh, side_factor, components, name),
  U_history(mesh),
  previously_dynamic(true) {
  
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
  
  HalfSpaceQuasiDynamic::setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::initConvolutions() {

  HalfSpaceQuasiDynamic::initConvolutions();

  this->convolutions.registerHistory(this->U_history);
  
  this->convolutions.init(std::make_pair(Kernel::Krnl::H00,0)); // H00-U0
  this->convolutions.init(std::make_pair(Kernel::Krnl::H01,0)); // H01-U0
  this->convolutions.init(std::make_pair(Kernel::Krnl::H11,1)); // H11-U1
  this->convolutions.init(std::make_pair(Kernel::Krnl::H01,1)); // H01-U1

  if (this->mesh.getDim()==3) {
    this->convolutions.init(std::make_pair(Kernel::Krnl::H22,0)); // H22-U0
    this->convolutions.init(std::make_pair(Kernel::Krnl::H00,2)); // H00-U2
    this->convolutions.init(std::make_pair(Kernel::Krnl::H01,2)); // H01-U2
    this->convolutions.init(std::make_pair(Kernel::Krnl::H22,2)); // H22-U2
  }    
  
#ifdef UCA_VERBOSE
  int total_work=0;
  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { //parallel loop
    for (int d=0; d<this->mesh.getDim(); ++d)
      total_work += this->U_history.get(d,j)->getSize();
  }
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  printf("Rank %d has total work %d \n",world_rank,total_work);
#endif /* UCA_VERBOSE */
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::computeStressFourierCoeff(bool predicting,
						 bool correcting,
						 bool dynamic) {
  if (!dynamic) {
    this->computeStressFourierCoeffQuasiDynamic(predicting, correcting);
    this->previously_dynamic = false;
  }
  else {
    if (!this->previously_dynamic)
      this->setSteadyState(predicting);
    this->computeStressFourierCoeffDynamic(predicting, correcting);
    this->previously_dynamic = true;
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::computeStressFourierCoeffDynamic(bool predicting,
							bool correcting) {

  FFTableNodalField & _disp = predicting ? this->disp_pc : this->disp;

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();

  const TwoDVector & wave_numbers = this->mesh.getLocalWaveNumbers();

  // --------------
  // UPDATE HISTORY
  // --------------
  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { // parallel loop over km modes

    // ignore mode 0
    if ((prank == m0_rank) && (j == m0_index)) continue;
    
    for (int d = 0; d < this->mesh.getDim(); ++d) {
      
      std::complex<double> Udj = {_disp.fd(j,d)[0], _disp.fd(j,d)[1]};
      
      if (correcting) {
        this->U_history.get(d,j)->changeCurrentValue(Udj);
      } else {
        this->U_history.get(d,j)->addCurrentValue(Udj);
      }
    }
  }
  
  // --------------
  // CONVOLUTION
  // --------------
  this->convolutions.convolve();
  auto conv = this->convolutions.getResults();
    
  // --------------
  // COMPUTE F
  // --------------

  // access to fourier coefficients of stresses
  fftw_complex * internal_fd[3];
  for (int d = 0; d < this->mesh.getDim(); ++d)
    internal_fd[d] = this->internal.fd_data(d);
 
  
  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { // parallel loop over km modes
    
    // ignore mode 0
    if ((prank == m0_rank) && (j == m0_index)) continue;
    
    // current U
    std::vector<std::complex<double>> U;
    U.resize(this->mesh.getDim());
    for (int d=0; d<this->mesh.getDim(); ++d)
      U[d] = {_disp.fd(j,d)[0], _disp.fd(j,d)[1]};
    
    // to be computed
    std::vector<std::complex<double>> F;
    F.resize(this->mesh.getDim());
    
    if (this->mesh.getDim() == 2) {
      double q = wave_numbers(j,0);
      this->computeF2D(F,
		       q,
		       U,
		       conv[std::make_pair(Kernel::Krnl::H00,0)][j],  // H00-U0
		       conv[std::make_pair(Kernel::Krnl::H01,0)][j],  // H01-U0
		       conv[std::make_pair(Kernel::Krnl::H01,1)][j],  // H01-U1
		       conv[std::make_pair(Kernel::Krnl::H11,1)][j]); // H11-U1
    } else {
      // q = {k,m} wave number in x,y direction
      double k = wave_numbers(j,0);
      double m = wave_numbers(j,2);

      this->computeF3D(F,k,m,
		       U,
		       conv[std::make_pair(Kernel::Krnl::H00,0)][j],  // H00-U0
		       conv[std::make_pair(Kernel::Krnl::H00,2)][j],  // H00-U2
		       conv[std::make_pair(Kernel::Krnl::H01,0)][j],  // H01-U0
		       conv[std::make_pair(Kernel::Krnl::H01,2)][j],  // H01-U2
		       conv[std::make_pair(Kernel::Krnl::H01,1)][j],  // H01-U1
		       conv[std::make_pair(Kernel::Krnl::H11,1)][j],  // H11-U1
		       conv[std::make_pair(Kernel::Krnl::H22,0)][j],  // H22-U0
		       conv[std::make_pair(Kernel::Krnl::H22,2)][j]); // H22-U2
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
  restart.registerIO(this->name+"_U",this->U_history);
  HalfSpace::registerToRestart(restart);
}

/* -------------------------------------------------------------------------- */
void HalfSpaceDynamic::setSteadyState(bool predicting) {

  FFTableNodalField & _disp = predicting ? this->disp_pc : this->disp;

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int m0_rank = this->mesh.getMode0Rank();
  int m0_index = this->mesh.getMode0Index();

  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { // parallel loop over km modes
    
    // ignore mode 0
    if ((prank == m0_rank) && (j == m0_index)) continue;

    for (int d = 0; d < this->mesh.getDim(); ++d) {

      this->U_history.get(d,j)->setSteadyState({_disp.fd(j,d)[0],
	                                        _disp.fd(j,d)[1]});
    }
  }
}

__END_UGUCA__
