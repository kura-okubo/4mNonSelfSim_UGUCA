/**
 * @file   uca_distributed_fftable_mesh.cc
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
#include "uca_distributed_fftable_mesh.hh"
#include "fftable_nodal_field.hh"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <numeric>      // std::iota
#include <algorithm>
#include <unistd.h> //sleep

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
DistributedFFTableMesh::DistributedFFTableMesh(double Lx,
					       int Nx) :
  FFTableMesh(Lx, Nx),
  max_fft_pp(0) {
  this->init();
}

/* -------------------------------------------------------------------------- */
DistributedFFTableMesh::DistributedFFTableMesh(double Lx, int Nx,
					       double Lz, int Nz) : 
  FFTableMesh(Lx, Nx, Lz, Nz),
  max_fft_pp(0) {
  this->init();
}

/* -------------------------------------------------------------------------- */
void DistributedFFTableMesh::init() {

  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  
  // if serial: just use FFTableMesh, hence local = global
  if (psize == 1) {} // serial
  
  else { // parallel

    if (psize%2 != 0)
      throw std::runtime_error("ERROR: number of mpi process needs to be even\n");
    
    // global wave numbers
    TwoDVector wave_numbers_global(this->dim,this->getNbGlobalFFT());  // local {k,-,m}
    this->initWaveNumbersGlobal(wave_numbers_global);

    // assign the modes to procs
    // fills the this->sort_fft_modes_map
    this->assignFFTModes(wave_numbers_global);

    // actually sort modes and scatter them
    for (int d=0; d<this->dim; ++d){
      this->sortAndScatterFFTModes(wave_numbers_global.data(d).data(),
				   this->getWaveNumbersLocalData(d),
				   this->root);
    }

    // for sorting (all ranks should be able to do it)
    //this->fftw_complex_buffer.resize(this->max_fft_pp*psize); // does not work because vector of array
    std::vector<fftw_complex> tmp(this->max_fft_pp*psize);
    this->fftw_complex_buffer.swap(tmp);
  }
}

/* -------------------------------------------------------------------------- */
void DistributedFFTableMesh::forwardFFT(FFTableNodalField & nodal_field) {
  
  FFTableMesh::forwardFFT(nodal_field);
  // loop over all components of the nodal field
  for (const auto& d : nodal_field.getComponents()) {
    this->sortAndScatterFFTModes(nodal_field.fd_data(d), this->root);
  }
}
  
/* -------------------------------------------------------------------------- */
void DistributedFFTableMesh::backwardFFT(FFTableNodalField & nodal_field) {

  // loop over all components of the nodal field
  for (const auto& d : nodal_field.getComponents()) {
    this->gatherAndSortFFTModes(nodal_field.fd_data(d), this->root);
  }
  FFTableMesh::backwardFFT(nodal_field);
}

/* --------------------------------------------------------------------------
 * FFT LOAD_BALANCE
 * for small mesh sizes or very big world_size not optimal
 * because 1st mode is very big
 */
void DistributedFFTableMesh::assignFFTModes(TwoDVector & wave_numbers_global) {

  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  if (psize == 1)
    throw std::runtime_error("cannot assign modes in serial\n");
  
  //split modes_km between processes

  // adapt for 3D case order modes with respect to their length
  int nb_fft_global = this->getNbGlobalFFT();
  std::vector<int> mode_assigned(nb_fft_global, this->root); //mode_assinged[job_id] = rank
  std::vector<int> nb_modes_p_rank(psize, 0);

  if (prank == this->root) {
    std::vector<double> totwork(psize, 0);
    std::vector<double> work_p_mode(nb_fft_global,0.0);
    this->computeWorkPerMode(wave_numbers_global,work_p_mode);
    
    // argsort by decreasing work_per_mode
    std::vector<int> sorted_modes_idx(nb_fft_global,0);
    std::iota(sorted_modes_idx.begin(), sorted_modes_idx.end(),0);
    std::sort(sorted_modes_idx.begin(), sorted_modes_idx.end(),
	      [&work_p_mode](int i1, int i2){return work_p_mode[i1]>work_p_mode[i2];});
    
    this->printMaximumNbProc(work_p_mode);

    int nb_jobs = nb_fft_global;//-1; // take away zero frequency
    for (int i=0; i<nb_jobs; ++i) {
      // loop over modes by decreasing work load
      int job_id = sorted_modes_idx[i];
      
      // assign job to process rank with minimum work
      int min_work_rank = std::distance(totwork.begin(),
					std::min_element(totwork.begin(),
							 totwork.end()));

      mode_assigned[job_id] = min_work_rank;
      nb_modes_p_rank[min_work_rank]+=1;
      totwork[min_work_rank] += work_p_mode[job_id];
    }
    
    // find max nb jobs per proc and set nb modes per proc
    this->max_fft_pp = *std::max_element(nb_modes_p_rank.begin(),
					 nb_modes_p_rank.end());
  } // prank == root
  
  // scatter individual nb_fft_local to each process
  int nb_fft_local = -1;
  StaticCommunicatorMPI::getInstance()->scatter(nb_modes_p_rank.data(),
						&nb_fft_local,1);
  // broadcast the max of nb_fft_local
  StaticCommunicatorMPI::getInstance()->broadcast(&this->max_fft_pp,1);

  // define how much fft needs to be allocated
  int nb_fft_local_alloc = -1;
  if (prank == this->root)
    nb_fft_local_alloc = this->max_fft_pp*psize;
  else
    nb_fft_local_alloc = this->max_fft_pp;

  // resize the local wave numbers
  FFTableMesh::resize(nb_fft_local, nb_fft_local_alloc);

  // resize the fft modes map
  this->sort_fft_modes_map.resize(nb_fft_global);

  
#ifdef UCA_VERBOSE
  if (prank == this->root) {
    std::cout<<"mode assigned [mode -> rank]:"<<std::endl;
    for (unsigned int i=0; i<mode_assigned.size(); ++i) {
      std::cout<< i << " -> " << mode_assigned[i]<<std::endl;
    }
  }
#endif

  
  // build sorting array from mode_assigned
  // allocate memory now that we know the size
  if (prank == this->root) {

    //this->sort_fft_modes_map[0]=nb_fft_global-1;
    for (int rank=0; rank<psize; ++rank) {
      int j = 0;
      for (int i=0; i<nb_fft_global; ++i) {
	int job_id = i;
	if (rank == mode_assigned[job_id]) {
	  this->sort_fft_modes_map[job_id] = j + this->max_fft_pp*rank;

	  // information about location of mode zero
	  if (job_id == 0) {
	    this->mode_zero_rank = rank;
	    this->mode_zero_index = j;
	  }
	  
	  ++j;
	}
      }
    }
  } // prank == root

  StaticCommunicatorMPI::getInstance()->broadcast(&this->mode_zero_rank,1);
  StaticCommunicatorMPI::getInstance()->broadcast(&this->mode_zero_index,1);

  StaticCommunicatorMPI::getInstance()->broadcast(this->sort_fft_modes_map.data(),
						  nb_fft_global); 

#ifdef UCA_VERBOSE
  if (prank == this->root) {
    std::cout<<"sort_fft_modes_map"<<std::endl;
    for (int i=0; i<nb_fft_global; ++i) {
      std::cout<<this->sort_fft_modes_map[i]<<" ";
    }
    std::cout << std::endl;
  }
#endif
}

/* -------------------------------------------------------------------------- */
void DistributedFFTableMesh::computeWorkPerMode(
						TwoDVector & wave_numbers_global,
						std::vector<double> & work_per_mode) {

  int nb_fft_global = this->getNbGlobalFFT();
  
  // estimate work_per_mode = nb_intrgration_int
  // k,m are components of mode vector q
  work_per_mode[0] = 0.0; // no convolution for 0 mode
  for (int i=1; i<nb_fft_global; ++i) {

    double k = wave_numbers_global(i,0);
    double m = 0.0;
    if (this->dim==3)
      m = wave_numbers_global(i,2);

    double q = std::sqrt(k*k + m*m);
    work_per_mode[i] = 1.0 / q;
  }

#ifdef UCA_VERBOSE
  std::cout<<"work per mode"<<std::endl;
  for (int i=0; i<nb_fft_global; ++i) {
    std::cout<< work_per_mode[i]<<std::endl;
  }
#endif

}

/* -------------------------------------------------------------------------- */
int DistributedFFTableMesh::getMaximumNbProc(std::vector<double> & work_per_mode) {
  /* For 2D:
   * N0 is the number of nodes
   * N0/2+1 is the number of FFT modes
   *
   * The total work during convolution is given by the nth partial sum of the
   * harmonic series,
   * H_n = sum_{k=1}^{n} 1/k
   * where n=N0/2+1
   *
   * For 3D:
   * N0 is the number of nodes along x
   * N1 is the number of nodes along z
   * the number of FFT Modes is N0x(N1/2+1)
   * the total work is give by
   *
   * H = sum_i sum_j 1/q_ij
   * where q_ij = sqrt(k_i^2+m_j^2)
   *
   * for 3D the optimal number of process quickly exceeds the available processes.
   * This is not an issue -> In 3D the are practically no load balance issues.
   * All processes will be assigned the same amount of work.
   *
   * Since our current parallel scheme does not allow us to subdivide the 1st
   * mode between different processes we can prove that the theoretical maximum
   * speedup due to parallelization (neglecting communication overhead) is:
   *
   * Speedup = SetialTime/ParallelTime <= TotalWork/Work1stMode = MaxSpeedup
   *
   * and it is achieved when one process only has the first mode. Therefore, the
   * maximum number of processes is also equivalent to TotalWork/Work1stMode.
   *
   */

  double tot_work = std::accumulate(work_per_mode.begin(),work_per_mode.end(),0.0);
  double max_work_per_mode = *std::max_element(work_per_mode.begin(),work_per_mode.end());

#ifdef UCA_VERBOSE
  std::cout<<"tot work : "<<tot_work<<", max w :"<<max_work_per_mode<<std::endl;
#endif /* UCA_VERBOSE */

  return int(tot_work/max_work_per_mode) + 1;// +1 coz intger division rounds down
}

/* -------------------------------------------------------------------------- */
void DistributedFFTableMesh::printMaximumNbProc(std::vector<double> & work_per_mode) {
  if (StaticCommunicatorMPI::getInstance()->whoAmI()==0) {
    std::cout << "---------------------------------------------------------"
	      << std::endl;
    printf("Maximal number of MPI processes for optimal load balance: \n");
    printf("MaxNbProc = %d \n",
	   this->getMaximumNbProc(work_per_mode));
    std::cout << "---------------------------------------------------------"
	      << std::endl;
  }
}


/* --------------------------------------------------------------------------
 * SORT FFT MODES
 * -------------------------------------------------------------------------- */
template <typename T>
void DistributedFFTableMesh::sortFFT(T* un_sorted, T * sorted, int root_rank) {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int nb_fft_global = this->getNbGlobalFFT();
  
  if (prank==root_rank) {
    for (int n=0; n<nb_fft_global; ++n) {
      sorted[n] = un_sorted[this->sort_fft_modes_map[n]];
    }
  }
}
 
/* -------------------------------------------------------------------------- */
template <typename T>
void DistributedFFTableMesh::unsortFFT(T * sorted, T * un_sorted, int root_rank) {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int nb_fft_global = this->getNbGlobalFFT();

  if (prank==root_rank) {
    for (int n=0; n<nb_fft_global; ++n) {
      un_sorted[this->sort_fft_modes_map[n]] = sorted[n];
    }
  }
}

/* --------------------------------------------------------------------------
 * template specialization for fftw complex
 * -------------------------------------------------------------------------- */
template <>
void DistributedFFTableMesh::sortFFT(fftw_complex * un_sorted,
				     fftw_complex  * sorted,
				     int root_rank) {
  
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int nb_fft_global = this->getNbGlobalFFT();

  if (prank==root_rank) {
    for (int n=0; n<nb_fft_global; ++n) {
      for (int j=0;j<2; ++j) {
	sorted[n][j] = un_sorted[this->sort_fft_modes_map[n]][j];
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
void DistributedFFTableMesh::unsortFFT(fftw_complex * sorted,
				       fftw_complex * un_sorted,
				       int root_rank) {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int nb_fft_global = this->getNbGlobalFFT();

  if (prank==root_rank) {
    for (int n=0; n<nb_fft_global; ++n) {
      for (int j=0;j<2; ++j) {
	un_sorted[this->sort_fft_modes_map[n]][j] = sorted[n][j];
      }
    }
  }
}

/* --------------------------------------------------------------------------
 * SORT SCATTER GATHER FFT MODE TEMPLATE
 * -------------------------------------------------------------------------- */
template <typename T>
void DistributedFFTableMesh::sortAndScatterFFTModes(T * Uglobal,
						    T * Ulocal,
						    T * buffer,
						    int size,
						    int root_rank) {
  
  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();

  if (psize > 1) {
    this->unsortFFT(Uglobal, buffer, root_rank);

    StaticCommunicatorMPI::getInstance()->scatter(
			   buffer, Ulocal, size, root_rank);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void DistributedFFTableMesh::gatherAndSortFFTModes(T * Uglobal,
						   T * Ulocal,
						   T * buffer,
						   int size,
						   int root_rank) {
  
  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();

  if (psize > 1) {
    StaticCommunicatorMPI::getInstance()->gather(
			   buffer, Ulocal, size, root_rank);

    this->sortFFT(buffer,Uglobal,root_rank);
  }
}

/* --------------------------------------------------------------------------
 * SORT AND SCATTER FFTW MODES
 * -------------------------------------------------------------------------- */
/*
void DistributedFFTableMesh::sortAndScatterFFTModes(int * Uglobal,
						    int * Ulocal) {

  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  // we are allocating and freeing memory here just
  // because this function is used only during initialization
  int * bufferint = NULL; 

  if (prank==this->root)
    bufferint = new int[(this->nb_fft_local_alloc)*(psize)]();

  this->sortAndScatterFFTModes(Uglobal,Ulocal,bufferint,this->max_fft_pp);
  if (prank==this->root) delete[] bufferint;
}
*/
/* -------------------------------------------------------------------------- */
void DistributedFFTableMesh::sortAndScatterFFTModes(double * Uglobal,
						    double * Ulocal,
						    int root_rank) {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  // we are allocating and freeing memory here just
  // because this function is used only during initialization
  double * buffer = NULL;

  if (prank==root_rank)
    buffer = new double[this->getNbLocalFFTAlloc()]();

  this->sortAndScatterFFTModes(Uglobal, Ulocal, buffer,
			       this->max_fft_pp, root_rank);
  if (prank==this->root) delete[] buffer;
}

/* -------------------------------------------------------------------------- */
void DistributedFFTableMesh::sortAndScatterFFTModes(fftw_complex * Uglobal,
						    fftw_complex * Ulocal,
						    int root_rank) {
  // 2* accounts for real and imag parts -> fftw_complex aka double[2]
  this->sortAndScatterFFTModes(Uglobal,
			       Ulocal,
			       this->fftw_complex_buffer.data(),
			       (this->max_fft_pp)*2,
			       root_rank);
}

/* --------------------------------------------------------------------------
 * GATHER AND SORT FFTW MODES
 * -------------------------------------------------------------------------- */
void DistributedFFTableMesh::gatherAndSortFFTModes(fftw_complex * Uglobal,
						   fftw_complex * Ulocal,
						   int root_rank) {
  // 2* accounts for real and imag parts -> fftw_complex aka double[2]
  this->gatherAndSortFFTModes(Uglobal,
			      Ulocal,
			      this->fftw_complex_buffer.data(),
			      (this->max_fft_pp)*2,
			      root_rank);
}

__END_UGUCA__
