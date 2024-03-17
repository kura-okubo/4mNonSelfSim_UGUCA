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
#include "static_communicator_mpi.hh"

#include <iostream>
#include <cmath>

using namespace uguca;

/* This is serial
 * however it can be executed in parallel and only the 0 rank would perform the operations
 */

int main() {

  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  
  if (prank==0)
    std::cout << "start test: fftable_nodal_field" << std::endl;

  int nb_points_x = 64*2;
  int nb_points_z = 64;

  double length_x = 3.0;
  double length_z = 1.0;

  FFTableMesh mesh(length_x,nb_points_x,
		   length_z,nb_points_z);
  
  FFTableNodalField ftf(mesh, {_x,_y,_z});
  if (prank==0)
    std::cout<<"fftable field contructed"<< std::endl;

  NodalField solution(mesh, {_x,_y,_z});

  const TwoDVector & coords = mesh.getLocalCoords();

  // initiate field and ftf
  for (const auto& d : solution.getComponents()) {
    for (int i=0; i<mesh.getNbGlobalNodes(0); ++i) {
      for (int j=0; j<mesh.getNbGlobalNodes(2); ++j) {
	int ij = i*mesh.getNbGlobalNodes(2)+j;
	double x = coords(ij,0);
	double z = coords(ij,2);
	if (prank!=0) continue;
	solution(ij,d) = 1.0+cos(x*4*M_PI)*sin(z*8*M_PI)+cos(x*2*M_PI);
	ftf(ij,d) = solution(ij,d);
      }
    }
  }

  if (prank==0)
    std::cout << "fftable field initialized" << std::endl
	      << "check forward and backward fft" << std::endl;

  ftf.forwardFFT();
  ftf.setAllValuesTo(0.0);
  ftf.backwardFFT();

  // check
  for (const auto& d : solution.getComponents()) {
    for (int i=0; i<mesh.getNbGlobalNodes(0); ++i) {
      for (int j=0; j<mesh.getNbGlobalNodes(2); ++j) {
	int ij = i*mesh.getNbGlobalNodes(2)+j;
	if (prank!=0) continue;
	if (fabs(solution(ij,d) - ftf(ij,d))>1e-12) {
	  std::cout <<"("<<i<<","<<j<<") "
		    << solution(ij,d) << " != " << ftf(ij,d) << " -> "
		    << fabs(solution(ij,d) - ftf(ij,d)) << "\n";
	  std::cout <<"2D FFT failed"<< std::endl;
	  return 1;
	}
      }
    }
  }

  if (prank==0)
  std::cout << "forward and backward FFT correct"<< std::endl
	    << "all checks passed -> overall success" << std::endl;

  StaticCommunicatorMPI::getInstance()->finalize();
  return 0; // success
}
