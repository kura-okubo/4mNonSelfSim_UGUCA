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

class FFTableNodalField;

class FFTableMesh : public BaseMesh {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  FFTableMesh(double Lx, int Nx);

  FFTableMesh(double Lx, int Nx,
	      double Lz, int Nz);

  virtual ~FFTableMesh();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void resize(int nb_fft, int alloc=-1);
  
  virtual void registerForFFT(FFTableNodalField & nodal_field);
  virtual void unregisterForFFT(FFTableNodalField & nodal_field);

  virtual void forwardFFT(FFTableNodalField & nodal_field);
  virtual void backwardFFT(FFTableNodalField & nodal_field);
  
protected:
  virtual void initWaveNumbersGlobal(TwoDVector & wave_numbers);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // access to the wave numbers
  const TwoDVector & getLocalWaveNumbers() const { return this->wave_numbers_local; }

  // geometry of interface
  double getLength(int d) const { return this->lengths.at(d); }
  
  // element sizes
  double getDelta(int d) const { return this->lengths.at(d) / this->nb_nodes_global.at(d); }
  
  // global nodes that may not be all here
  virtual int getNbGlobalNodes() const {
    int nb = 1; for (const auto & n : this->nb_nodes_global) nb*=n; return nb; }
  int getNbGlobalNodes(int d) const { return this->nb_nodes_global.at(d); }

  // global nodes in fourier space
  int getNbGlobalFFT() const {
    int nb = 1; for (const auto& e : this->nb_fft_global) nb*=e; return nb; }
  int getNbGlobalFFT(int d) const { return this->nb_fft_global.at(d); }
  
  int getNbLocalFFT() { return this->nb_fft_local; }
  int getNbLocalFFTAlloc() { return this->nb_fft_local_alloc; }

  int getMode0Rank() const { return this->mode_zero_rank; }
  int getMode0Index() const { return this->mode_zero_index; }

protected:
  double * getWaveNumbersLocalData(int d) { return this->wave_numbers_local.data(d).data(); }
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // length of domain / replication length
  std::vector<double> lengths;
  
  // number of nodes
  std::vector<int> nb_nodes_global;
  
  // Fourier modes
  std::vector<int> nb_fft_global;
  
private:
  // wave numbers in fourier space: local
  int nb_fft_local;
  int nb_fft_local_alloc;
  TwoDVector wave_numbers_local;   // local {k,-,m}

protected:
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
