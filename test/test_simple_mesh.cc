/**
 * @file   test_simple_mesh.cc
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
#include <iostream>
#include <cmath>

#include "uca_simple_mesh.hh"

using namespace uguca;

class TestMesh : public SimpleMesh {
public:
  TestMesh(double Lx, int Nx) : SimpleMesh(Lx,Nx) {}
  TestMesh(double Lx, int Nx,
	   double Lz, int Nz) : SimpleMesh(Lx,Nx,Lz,Nz) {}
};

int main(){

  //int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  unsigned int nb_nodes_x = 10;
  double length_x = 3;

  unsigned int nb_nodes_z = 4;
  double length_z = 4.0;


  // construct 2D mesh
  if (prank==0)
    printf("==============\n2D mesh\n==============\n");
  TestMesh mesh2d(length_x,nb_nodes_x);

  if (prank==0)
    std::cout << "test 2d coords" << std::endl;

  double ** coords2d = mesh2d.getLocalCoords();

  if (prank==0) {

    int nx = mesh2d.getNbLocalNodes();
    for (int i=0; i<nx;++i){
      if (std::abs(coords2d[0][i] - mesh2d.getLengthX()*i/nx)>1e-6){
	std::cerr << "error0\n";
	return 1; // failure
      }
      if (std::abs(coords2d[1][i] - 0.0)>1e-6){
	std::cerr << "error1\n";
	return 1; // failure
      }
    }
  }

  else if (prank>0) {

    int n = mesh2d.getNbLocalNodes();
    if (n != 0) {
      std::cerr << "rank > 0 should have zero nodes instead of " << n << std::endl;
      return 1; // failure
    }
  }

  
  // construct 3D mesh
  if (prank==0)
    printf("==============\n3D mesh\n==============\n");
  TestMesh mesh3d(length_x,nb_nodes_x,
		  length_z,nb_nodes_z);


  double ** coords3d = mesh3d.getLocalCoords();

  if (prank==0) {

    int nx = mesh3d.getNbGlobalNodesX();
    int nz = mesh3d.getNbGlobalNodesZ();
    for (int i=0; i<nx;++i) {
      for (int j=0; j<nz; ++j) {
	int ij=i*nz+j;
	if (std::abs(coords3d[0][ij] - mesh3d.getLengthX()*i/nx) > 1e-6) {
	  std::cerr << "error0\n";
	  return 1; // failure
	}
	if (std::abs(coords3d[1][ij] - 0.0) > 1e-6) {
	  std::cerr << "error1\n";
	  return 1; // failure
	}
	if (std::abs(coords3d[2][ij] - mesh3d.getLengthZ()*j/nz) > 1e-6) {
	  std::cout<<"error2\n";
	  return 1; // failure
	}
      }
    }
  }

  else if (prank>0) {

    int n = mesh3d.getNbLocalNodes();
    if (n != 0) {
      std::cerr << "rank > 0 should have zero nodes instead of " << n << std::endl;
      return 1; // failure
    }
  }
  
  
  if (prank==0)
    std::cout << "all checks passed -> overall success" << std::endl;

  StaticCommunicatorMPI::getInstance()->finalize();
  return 0;
}

