/**
 * @file   nodal_field.cc
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
#include "nodal_field_component.hh"
#include "static_communicator_mpi.hh"

#include <stdexcept>
#include <iostream>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
NodalFieldComponent::~NodalFieldComponent() {
  if (this->initialized) {
#ifdef UCA_VERBOSE
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::cout << "NFC del (prank="
	    << prank << "): " << this->name << std::endl;
#endif // UCA_VERBOSE

    delete[] this->field;
  }
}

/* -------------------------------------------------------------------------- */
void NodalFieldComponent::init(BaseMesh & mesh) {
  // do not initialize twice
  if (this->initialized)
    throw std::runtime_error("NodalFieldComponent: do not initialize twice\n");

  this->mesh = &mesh;

#ifdef UCA_VERBOSE
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::cout << "NFC init (prank="
	    << prank << "): " << this->name
	    << " : " << mesh.getNbLocalNodesAlloc() << std::endl;
#endif // UCA_VERBOSE
  
  // if nb_nodes==0 -> use 1
  //unsigned int size = std::max(this->mesh->getNbLocalNodesAlloc(),1);
  this->field = new double [this->mesh->getNbLocalNodesAlloc()];
  
  this->initialized = true;
  this->zeros();
}

/* -------------------------------------------------------------------------- */
void NodalFieldComponent::zeros() {
  this->setAllValuesTo(0.);
}

/* -------------------------------------------------------------------------- */
void NodalFieldComponent::multiply(const NodalFieldComponent & other) {

  double * this_p = this->storage();
  const double * other_p = other.storage();
  
  for (int n=0; n<this->getNbNodes(); ++n) {
    this_p[n] *= other_p[n];
  }
}

/* -------------------------------------------------------------------------- */
void NodalFieldComponent::setAllValuesTo(double value) {
  // do not initialize twice
  if (!this->initialized)
    throw std::runtime_error("NodalFieldComponent: not initialize yet\n");

  std::fill_n(this->field,this->mesh->getNbLocalNodes(),value);
}

__END_UGUCA__
