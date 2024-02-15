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

  this->convols.registerHistory(this->disp);
  
  this->convols.init(std::make_pair(Kernel::Krnl::H00,0)); // H00-U0
  this->convols.init(std::make_pair(Kernel::Krnl::H01,0)); // H01-U0
  this->convols.init(std::make_pair(Kernel::Krnl::H11,1)); // H11-U1
  this->convols.init(std::make_pair(Kernel::Krnl::H01,1)); // H01-U1

  if (this->mesh.getDim()==3) {
    this->convols.init(std::make_pair(Kernel::Krnl::H22,0)); // H22-U0
    this->convols.init(std::make_pair(Kernel::Krnl::H00,2)); // H00-U2
    this->convols.init(std::make_pair(Kernel::Krnl::H01,2)); // H01-U2
    this->convols.init(std::make_pair(Kernel::Krnl::H22,2)); // H22-U2
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
  for (int j=0; j<this->mesh.getNbLocalFFT(); ++j) { //parallel loop
    for (int d=0; d<this->mesh.getDim(); ++d)
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

  // --------------
  // UPDATE HISTORY
  // --------------
  if (predicting) {
    this->disp.addCurrentValueToHistory(this->disp_pc);
  }
  else {
    this->disp.addCurrentValueToHistory();
  }
  if (correcting) {
    this->disp.changeCurrentValueOfHistory();
  }
  
  // --------------
  // CONVOLUTION
  // --------------
  this->convols.convolve();
    
  // --------------
  // COMPUTE F
  // --------------
  FFTableNodalField & _disp = predicting ? this->disp_pc : this->disp;
  this->computeF(this->internal,_disp,
		 this->convols.getResult(std::make_pair(Kernel::Krnl::H00,0)),  // H00-U0
		 this->convols.getResult(std::make_pair(Kernel::Krnl::H00,2)),  // H00-U2
		 this->convols.getResult(std::make_pair(Kernel::Krnl::H01,0)),  // H01-U0
		 this->convols.getResult(std::make_pair(Kernel::Krnl::H01,2)),  // H01-U2
		 this->convols.getResult(std::make_pair(Kernel::Krnl::H01,1)),  // H01-U1
		 this->convols.getResult(std::make_pair(Kernel::Krnl::H11,1)),  // H11-U1
		 this->convols.getResult(std::make_pair(Kernel::Krnl::H22,0)),  // H22-U0
		 this->convols.getResult(std::make_pair(Kernel::Krnl::H22,2))); // H22-U2
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
