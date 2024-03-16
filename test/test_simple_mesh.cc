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

#include "uca_simple_mesh.hh"

#include <iostream>
#include <cmath>

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

  const TwoDVector & coords2d = mesh2d.getLocalCoords();

  if (prank==0) {

    int nx = mesh2d.getNbLocalNodes();
    for (int i=0; i<nx;++i){
      if (std::abs(coords2d(i,0) - mesh2d.getLength(0)*i/nx)>1e-6){
	std::cerr << "error0\n";
	return 1; // failure
      }
      if (std::abs(coords2d(i,1) - 0.0)>1e-6){
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


  const TwoDVector & coords3d = mesh3d.getLocalCoords();

  if (prank==0) {

    int nx = mesh3d.getNbGlobalNodes(0);
    int nz = mesh3d.getNbGlobalNodes(2);
    for (int i=0; i<nx;++i) {
      for (int j=0; j<nz; ++j) {
	int ij=i*nz+j;
	if (std::abs(coords3d(ij,0) - mesh3d.getLength(0)*i/nx) > 1e-6) {
	  std::cerr << "error0\n";
	  return 1; // failure
	}
	if (std::abs(coords3d(ij,1) - 0.0) > 1e-6) {
	  std::cerr << "error1\n";
	  return 1; // failure
	}
	if (std::abs(coords3d(ij,2) - mesh3d.getLength(2)*j/nz) > 1e-6) {
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

