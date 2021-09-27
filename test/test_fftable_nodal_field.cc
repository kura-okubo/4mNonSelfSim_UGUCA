/**
 * @file   test_fftable_nodal_field.cc
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

#include <iostream>
#include <cmath>

using namespace uguca;

int main() {

  std::cout << "start test: fftable_nodal_field" << std::endl;

  int nb_points_x = 64*2;
  int nb_points_z = 64;

  double length_x = 3.0;
  double length_z = 1.0;

  FFTableMesh mesh(length_x,nb_points_x,
		  length_z,nb_points_z);
  
  FFTableNodalField Ftf(mesh);
  NodalField Field(mesh);

  double ** coords = mesh.getLocalCoords();

  // initiate field and ftf
  for (int d=0; d<mesh.getDim(); ++d) {
    FFTableNodalFieldComponent & ftf = Ftf.component(d);
    NodalFieldComponent & field = Field.component(d);
    for (int i=0; i<mesh.getNbGlobalNodesX(); ++i) {
      for (int j=0; j<mesh.getNbGlobalNodesZ(); ++j) {
	int ij = i*mesh.getNbGlobalNodesZ()+j;
	double x = coords[0][ij];
	double z = coords[2][ij];
	field.set(ij) = 1.0+cos(x*4*M_PI)*sin(z*8*M_PI)+cos(x*2*M_PI);
	ftf.set(ij) = field.at(ij);
      }
    }
  }

  std::cout << "fftable field initialized" << std::endl
	    << "check forward and backward fft" << std::endl;

  Ftf.forwardFFT();
  Ftf.setAllValuesTo(0.0);
  Ftf.backwardFFT();

  // check
  for (int d=0; d<mesh.getDim(); ++d) {
    FFTableNodalFieldComponent & ftf = Ftf.component(d);
    NodalFieldComponent & field = Field.component(d);
    for (int i=0; i<mesh.getNbGlobalNodesX(); ++i) {
      for (int j=0; j<mesh.getNbGlobalNodesZ(); ++j) {
	int ij = i*mesh.getNbGlobalNodesZ()+j;
	if (fabs(field.at(ij) - ftf.at(ij))>1e-12) {
	  std::cout <<"("<<i<<","<<j<<") "
		    << field.at(ij) << " != " << ftf.at(ij) << " -> "
		    << fabs(field.at(ij) - ftf.at(ij)) << "\n";
	  std::cout <<"2D FFT failed"<< std::endl;
	  return 1;
	}
      }
    }
  }

  std::cout << "forward and backward FFT correct"<< std::endl
	    << "all checks passed -> overall success" << std::endl;

  return 0; // success
}
