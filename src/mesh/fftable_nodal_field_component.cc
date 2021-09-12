/**
 * @file   fftable_nodal_field.cc
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
#include "fftable_nodal_field_component.hh"
#include "static_communicator_mpi.hh"
#include "uca_distributed_fftable_mesh.hh"
#include "uca_simple_mesh.hh"
#include "uca_custom_mesh.hh"
#include <cstring>
#include <typeinfo>
#include <stdexcept>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
/* use default FFTW datastructure N0xN1                                       */
/* -------------------------------------------------------------------------- */
FFTableNodalFieldComponent::FFTableNodalFieldComponent(FFTableMesh & mesh,
						       int direction,
						       const std::string & name) :
  NodalFieldComponent(mesh,direction,name),
  fd_field(NULL),
  fftw_plan_id(-1)
{
  this->init(mesh);
}

/* -------------------------------------------------------------------------- */
FFTableNodalFieldComponent::~FFTableNodalFieldComponent() {
  delete[] this->fd_field;
}

/* -------------------------------------------------------------------------- */
void FFTableNodalFieldComponent::init(FFTableMesh & mesh) {
  this->fd_field = new fftw_complex[mesh.getNbLocalFFTAlloc()];
  
  // init freq dom field with zeros
  memset(this->fd_field,0.,mesh.getNbLocalFFTAlloc()*sizeof(fftw_complex));

  ((FFTableMesh *)this->mesh)->registerForFFT(*this);
}

/* -------------------------------------------------------------------------- */
void FFTableNodalFieldComponent::free() {
  ((FFTableMesh *)this->mesh)->unregisterForFFT(*this);
  delete[] this->fd_field;
}

/* -------------------------------------------------------------------------- */
void FFTableNodalFieldComponent::update() {
  this->free();
  this->init(*(FFTableMesh *)this->mesh);
}

/* -------------------------------------------------------------------------- */
void FFTableNodalFieldComponent::forwardFFT() {

  if ((typeid(*(this->mesh)) == typeid(SimpleMesh)) ||
      (typeid(*(this->mesh)) == typeid(DistributedFFTableMesh))) {
    ((DistributedFFTableMesh *)this->mesh)->forwardFFT(*this);
  }
  else if (typeid(*(this->mesh)) == typeid(FFTableMesh)) {
    ((FFTableMesh *)this->mesh)->forwardFFT(*this);
  }
  else if (typeid(*(this->mesh)) == typeid(CustomMesh)) {
    ((CustomMesh *)this->mesh)->forwardFFT(*this);
  }
  else {
    throw std::runtime_error("Do not know mesh type\n");
  }
}

/* -------------------------------------------------------------------------- */
void FFTableNodalFieldComponent::backwardFFT() {
  if ((typeid(*(this->mesh)) == typeid(SimpleMesh)) ||
      (typeid(*(this->mesh)) == typeid(DistributedFFTableMesh))) {
    ((DistributedFFTableMesh *)this->mesh)->backwardFFT(*this);
  }
  else if (typeid(*(this->mesh)) == typeid(FFTableMesh)) {
    ((FFTableMesh *)this->mesh)->backwardFFT(*this);
  }
  else if (typeid(*(this->mesh)) == typeid(CustomMesh)) {
    ((CustomMesh *)this->mesh)->backwardFFT(*this);
  }
  else {
    throw std::runtime_error("Do not know mesh type\n");
  }
}

__END_UGUCA__
