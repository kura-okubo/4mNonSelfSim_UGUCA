/**
 * @file   uca_fftable_mesh.cc
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
#include "uca_fftable_mesh.hh"
#include "fftable_nodal_field.hh"
#include "static_communicator_mpi.hh"

#include <cmath>
#include <stdexcept>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
FFTableMesh::FFTableMesh(double Lx, int Nx, bool initialize) :
  BaseMesh(2, Nx),
  fs_allocated(false),
  length_x(Lx),
  length_z(0.),
  nb_nodes_x_global(Nx),
  nb_nodes_z_global(1),
  mode_zero_rank(0),
  mode_zero_index(0) {
  if (initialize)
    this->init();
}

/* -------------------------------------------------------------------------- */
FFTableMesh::FFTableMesh(double Lx, int Nx,
			 double Lz, int Nz,
			 bool initialize) :
  BaseMesh(3, Nx*Nz),
  fs_allocated(false),
  length_x(Lx),
  length_z(Lz),
  nb_nodes_x_global(Nx),
  nb_nodes_z_global(Nz),
  mode_zero_rank(0),
  mode_zero_index(0) {
  if (initialize)
    this->init();
}

/* -------------------------------------------------------------------------- */
FFTableMesh::~FFTableMesh() {
  this->freeSpectralSpace();

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (prank == this->root) {
    
    for (unsigned int i=0; i<this->forward_plans.size(); ++i) {
      fftw_destroy_plan(this->forward_plans[i]);
      fftw_destroy_plan(this->backward_plans[i]);
      // do not call fftw_cleanup() because we don't know if all plans are destroyed
    }
  }
  //this->forward_plans.resize(0);
  //this->backward_plans.resize(0);
}

/* -------------------------------------------------------------------------- */
void FFTableMesh::init() {
  this->initSpectralSpace();
  this->allocateSpectralSpace();
  this->initWaveNumbersGlobal(this->wave_numbers_local);
}

/* -------------------------------------------------------------------------- */
void FFTableMesh::initSpectralSpace() {

  if (this->dim==2) {
    this->nb_fft_x_global = this->nb_nodes_x_global / 2 + 1;
    this->nb_fft_z_global = 1;
  }
  else if (this->dim==3) {
    this->nb_fft_x_global = this->nb_nodes_x_global;
    this->nb_fft_z_global = this->nb_nodes_z_global / 2 + 1;
  }

  this->nb_fft_local = this->getNbGlobalFFT();
  this->nb_fft_local_alloc = this->nb_fft_local;
}

/* -------------------------------------------------------------------------- */
void FFTableMesh::allocateSpectralSpace() {

  // do not initialize twice
  if (this->fs_allocated)
    throw std::runtime_error("FFTableMesh: do not allocate twice\n");

  this->allocateVector(this->wave_numbers_local, this->nb_fft_local_alloc);
  this->fs_allocated = true;
}

/* -------------------------------------------------------------------------- */
void FFTableMesh::freeSpectralSpace() {

  if (this->fs_allocated) {
    this->freeVector(this->wave_numbers_local);
  }
  this->fs_allocated = false;
}

/* -------------------------------------------------------------------------- */
void FFTableMesh::registerForFFT(FFTableNodalField & nodal_field) {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (prank == this->root) {

    // loop over all components of the nodal field
    for (const auto& d : nodal_field.getComponents()) {

      // access relevant storage
      double * nfc = nodal_field.data(d);
      fftw_complex * fd_nfc = nodal_field.fd_data(d);

      // set plans to NULL
      fftw_plan forward_plan = NULL;
      fftw_plan backward_plan = NULL;

      // create fftw plans for relevant dimension
      if (this->dim==2) {
	forward_plan = fftw_plan_dft_r2c_1d(this->nb_nodes_x_global,
					    nfc, fd_nfc, FFTW_MEASURE);
	backward_plan = fftw_plan_dft_c2r_1d(this->nb_nodes_x_global,
					     fd_nfc, nfc, FFTW_MEASURE);
      }
      else if (this->dim==3) {
	forward_plan = fftw_plan_dft_r2c_2d(this->nb_nodes_x_global,
					    this->nb_nodes_z_global,
					    nfc, fd_nfc, FFTW_MEASURE);
	backward_plan = fftw_plan_dft_c2r_2d(this->nb_nodes_x_global,
					     this->nb_nodes_z_global,
					     fd_nfc, nfc, FFTW_MEASURE);
      }
      this->forward_plans.push_back(forward_plan);
      this->backward_plans.push_back(backward_plan);
    
      // set the fftw_plan_id of the nodal field component (can do because friend)
      nodal_field.fftw_plan_ids[d] = this->forward_plans.size()-1;
    }
  }
}

/* -------------------------------------------------------------------------- */
void FFTableMesh::unregisterForFFT(FFTableNodalField & nodal_field) {
  
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (prank == this->root) {
    // loop over all components of the nodal field
    for (const auto& d : nodal_field.getComponents()) {
      fftw_destroy_plan(this->forward_plans[nodal_field.getFFTWPlanId(d)]);
      fftw_destroy_plan(this->backward_plans[nodal_field.getFFTWPlanId(d)]);
    }
  }
}

/* -------------------------------------------------------------------------- */
void FFTableMesh::forwardFFT(FFTableNodalField & nodal_field) {
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (prank == this->root) {
    // loop over all components of the nodal field
    for (const auto& d : nodal_field.getComponents()) {
      fftw_execute(this->forward_plans[nodal_field.getFFTWPlanId(d)]);
    }
  }
}
  
/* -------------------------------------------------------------------------- */
void FFTableMesh::backwardFFT(FFTableNodalField & nodal_field) {
  
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (prank == this->root) {
    // loop over all components of the nodal field
    for (const auto& d : nodal_field.getComponents()) {

      fftw_execute(this->backward_plans[nodal_field.getFFTWPlanId(d)]);
    
      double nb_nodes_global = this->getNbGlobalNodes();
      double * p_field = nodal_field.data(d);
      for (int i=0; i<this->nb_nodes_local_alloc; ++i) { // for fftw_mpi it includes padding
	p_field[i] /= nb_nodes_global;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
/*
 * Serial version includes mode 0 (not used in computeStressFourierCoeff)
 * to avoid sorting before adding the zero mode back in
 * computeStressFourierCoeff()
 *
 */
void FFTableMesh::initWaveNumbersGlobal(double ** wave_numbers) {

  // fundamental mode
  double k1 = 2*M_PI / this->length_x;
  double m1 = 0.0;
  if (this->dim == 3)
    m1 = 2*M_PI / this->length_z;

  // compute nyquest frequency
  int f_ny_x = this->nb_fft_x_global; // 2D
  if (this->dim == 3)
    f_ny_x = this->nb_fft_x_global/2+1;

  // init wave numbers global
  for (int i=0; i<this->nb_fft_x_global; ++i) {
    for (int j=0; j<this->nb_fft_z_global; ++j) {
      int ij =  i*this->nb_fft_z_global + j;
      
      if (this-> dim == 2) {
	wave_numbers[0][ij] = k1 * i;
      }
      else {
	// after nyquest frequncy modes are negative
	// e.g. 0, 1, 2, ... ny, -ny, -ny+1, ... -2, -1
	wave_numbers[0][ij] = k1 * (i - (i/f_ny_x)*this->nb_fft_x_global);
	wave_numbers[2][ij] = m1 * j;
      }

      wave_numbers[1][ij] = 0.0; // we re on the xz plane
    }
  }
}

__END_UGUCA__
