/**
 * @file   uca_fftable_mesh.hh
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
#ifndef __FFTABLE_MESH_H__
#define __FFTABLE_MESH_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_base_mesh.hh"

#include <fftw3.h>

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

class FFTableNodalFieldComponent;

class FFTableMesh : public BaseMesh {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  FFTableMesh(double Lx, int Nx,
	      bool initialize=true);

  FFTableMesh(double Lx, int Nx,
	      double Lz, int Nz,
	      bool initialize=true);

  virtual ~FFTableMesh();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void init();

  virtual void registerForFFT(FFTableNodalFieldComponent & nodal_field_comp);
  virtual void unregisterForFFT(FFTableNodalFieldComponent & nodal_field_comp);

  virtual void forwardFFT(FFTableNodalFieldComponent & nodal_field_comp);
  virtual void backwardFFT(FFTableNodalFieldComponent & nodal_field_comp);
  
protected:
  virtual void initWaveNumbersGlobal(double ** wave_numbers);

  virtual void initSpectralSpace();
  virtual void allocateSpectralSpace();
  virtual void freeSpectralSpace();
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // access to the wave numbers
  double ** getLocalWaveNumbers()
  { return this->wave_numbers_local; }

  // geometry of interface
  double getLengthX() const { return this->length_x; }
  double getLengthZ() const { return this->length_z; }

  // element sizes
  double getDeltaX() { return this->length_x / this->nb_nodes_x_global; }
  double getDeltaZ() { return this->length_z / this->nb_nodes_z_global; }
  
  // global nodes that may not be all here
  virtual int getNbGlobalNodes() const {
    return this->nb_nodes_x_global * this->nb_nodes_z_global; }
  int getNbGlobalNodesX() const { return this->nb_nodes_x_global; }
  int getNbGlobalNodesZ() const { return this->nb_nodes_z_global; }

  // global nodes in fourier space
  int getNbGlobalFFTX() const { return this->nb_fft_x_global; }
  int getNbGlobalFFTZ() const { return this->nb_fft_z_global; }
  int getNbGlobalFFT() { return this->nb_fft_x_global * this->nb_fft_z_global; }
  int getNbLocalFFT() { return this->nb_fft_local; }
  int getNbLocalFFTAlloc() { return this->nb_fft_local_alloc; }

  int getMode0Rank() const { return this->mode_zero_rank; }
  int getMode0Index() const { return this->mode_zero_index; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // fourier space allocated
  bool fs_allocated;
  
  // length of domain / replication length
  double length_x;
  double length_z;

  // number of nodes
  int nb_nodes_x_global; // global
  int nb_nodes_z_global; // global

  // Fourier modes
  int nb_fft_x_global; // global
  int nb_fft_z_global; // global
  
  // wave numbers in fourier space: local
  int nb_fft_local;
  int nb_fft_local_alloc;
  double * wave_numbers_local[3];  // local {k,-,m}

  // information on where the zero mode is (it needs to be ignored)
  int mode_zero_rank;
  int mode_zero_index;
  
  // fft plans
  std::vector<fftw_plan> forward_plans;
  std::vector<fftw_plan> backward_plans;
};

__END_UGUCA__

//#include "uca_fftable_mesh_impl.cc"

#endif /* __FFTABLE_MESH_H__ */
