/**
 * @file   uca_distributed_fftable_mesh.hh
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Fri Feb 5 2021
 *
 * @brief  Contains only coordinates, but does not know anything about them
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
#ifndef __DISTRIBUTED_FFTABLE_MESH_H__
#define __DISTRIBUTED_FFTABLE_MESH_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "static_communicator_mpi.hh"
#include "uca_fftable_mesh.hh"
#include "nodal_field.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

class DistributedFFTableMesh : public FFTableMesh {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  DistributedFFTableMesh(double Lx, int Nx,
			 bool initialize = true);

  DistributedFFTableMesh(double Lx, int Nx,
			 double Lz, int Nz,
			 bool initialize = true);

  virtual ~DistributedFFTableMesh();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void init();

  //virtual void registerForFFT(FFTableNodalFieldComponent & nodal_field_comp);
  virtual void forwardFFT(FFTableNodalField & nodal_field);
  virtual void backwardFFT(FFTableNodalField & nodal_field);
  
protected:
  // for parallel implementation of computeStressFourierCoeff()
  void assignFFTModes(double ** wave_numbers_global);
  void computeWorkPerMode(double ** wave_numbers_global,
			  std::vector<double> & work_per_mode);

  // useful information
  void printMaximumNbProc(std::vector<double> & work_per_mode);
  int getMaximumNbProc(std::vector<double> & work_per_mode);

  // scatter and gather modes
  //void sortAndScatterFFTModes(int * Uglobal, int * Ulocal, int root_rank);
  //void sortAndScatterFFTModes(int * U, int root_rank){sortAndScatterFFTModes(U,U, root_rank);}

  void sortAndScatterFFTModes(double * Uglobal, double * Ulocal, int root_rank); //used for distributing Modes
  //void sortAndScatterFFTModes(double * U, int root_rank){sortAndScatterFFTModes(U,U,root_rank);}

  void sortAndScatterFFTModes(fftw_complex * Uglobal, fftw_complex * Ulocal, int root_rank);
  void sortAndScatterFFTModes(fftw_complex * U, int root_rank){sortAndScatterFFTModes(U,U,root_rank);}

  void gatherAndSortFFTModes (fftw_complex * Uglobal, fftw_complex * Ulocal, int root_rank);
  void gatherAndSortFFTModes (fftw_complex * U, int root_rank){gatherAndSortFFTModes(U,U,root_rank);}
  
  // for load balance in parallel convolution
  template <typename T>
  void sortFFT(T * un_sorted, T * sorted, int root_rank);

  template <typename T>
  void unsortFFT(T * sorted, T * un_sorted, int root_rank);

  template <typename T>
  void sortAndScatterFFTModes(T * Uglobal, T * Ulocal, T * buffer, int size, int root_rank);

  template <typename T>
  void gatherAndSortFFTModes(T * Uglobal, T * Ulocal, T * buffer, int size, int root_rank);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // global wave numbers (only used temporarily)
  double * wave_numbers_global[3];  // local {k,-,m}
  
  // map for distribution of fft modes
  int * sort_fft_modes_map;

  // maximum of local fft per proc
  int max_fft_pp;
  
  // to sort modes
  fftw_complex * fftw_complex_buffer; 
};

__END_UGUCA__

//#include "uca_distributed_fftable_mesh_impl.cc"

#endif /* __DISTRIBUTED_FFTABLE_MESH_H__ */
