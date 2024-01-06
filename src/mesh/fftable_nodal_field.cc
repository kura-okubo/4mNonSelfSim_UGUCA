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
#include <cstring>
#include <typeinfo>
#include <stdexcept>
#include <iostream>
#include "fftable_nodal_field.hh"
#include "uca_distributed_fftable_mesh.hh"
#include "uca_simple_mesh.hh"
#include "uca_custom_mesh.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
/*FFTableNodalField::FFTableNodalField(FFTableMesh & mesh,
				     SpatialDirectionSet components,
				     const std::string & name) :
  NodalField(mesh,components,name) {
  this->init();
}

/* -------------------------------------------------------------------------- */
FFTableNodalField::FFTableNodalField(FFTableMesh & mesh,
				     SpatialDirectionSet components,
				     const std::string & name) :
  NodalField(mesh,components,name) {
  this->init();
}

/* -------------------------------------------------------------------------- */
void FFTableNodalField::init() {
  
  /*
  // do not initialize twice
  if (this->initialized)
    throw std::runtime_error("FFTableNodalField: do not initialize twice\n");

  this->mesh = &mesh;

  // if you want to use less components, simply use the NodalFieldComponent
  this->field.resize(this->mesh->getDim());
  for (int d=0; d<this->mesh->getDim(); ++d) {
    this->field[d] = new FFTableNodalFieldComponent(mesh,d,
						    this->name+"_"+std::to_string(d));
  }

  // done
  this->initialized = true;
  */

  // Loop over the components to define their start
  this->fd_start = {0};
  for (size_t i = 0; i < _spatial_dir_count; ++i) {
    int new_start = this->fd_start.back();
    if (this->components.count(i)) {
      new_start += ((FFTableMesh *)this->mesh)->getNbLocalFFTAlloc();
    }
    this->fd_start.push_back(new_start);
  }
  
#ifdef UCA_VERBOSE
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::cout << "FFTableNF init fd (prank="
	    << prank << "): " << this->name
	    << " : " << this->fd_start.back() << std::endl;
#endif // UCA_VERBOSE

  // allocate space in storage and fill with zeros
  //this->fd_storage.resize(this->fd_start.back()),{0,0}); // does not work because vector of array
  std::vector<fftw_complex> tmp(this->fd_start.back());
  //std::fill(tmp.begin(), tmp.end(), fftw_complex{0.0, 0.0});
  for (auto& element : tmp) {
    element[0] = 0.0;
    element[1] = 0.0;
  }
  this->fd_storage.swap(tmp);

  ((FFTableMesh *)this->mesh)->registerForFFT(*this);
}

/* -------------------------------------------------------------------------- */
void FFTableNodalField::forwardFFT() {
  /*  for (int d=0; d<this->mesh->getDim(); ++d) {
    ((FFTableNodalFieldComponent*)this->field[d])->forwardFFT();
    }*/

  
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
void FFTableNodalField::backwardFFT() {
  /*
  for (int d=0; d<this->mesh->getDim(); ++d) {
    ((FFTableNodalFieldComponent*)this->field[d])->backwardFFT();
    }*/
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
