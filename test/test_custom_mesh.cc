/**
 * @file   test_custom_mesh.cc
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
#include "uca_custom_mesh.hh"
#include "static_communicator_mpi.hh"
#include "fftable_nodal_field.hh"

#include <random>
#include <iostream>
#include <unistd.h>

using namespace uguca;

void print(double ** fld, int size){
  int dim = 3;//fld.size();
  for (int n=0; n<size;n++) {
    for (int d=0;d<dim;d++) {
      std::cout<<fld[d][n]<<", ";
    }
    std::cout<<std::endl;
  }
  std::cout<<"----"<<std::endl;
}

class TestMesh : public CustomMesh {
public:
  TestMesh(double Lx, int Nx, std::vector<double> & cx) :
    CustomMesh(Lx,Nx,cx) {}
  TestMesh(double Lx, int Nx,
	   double Lz, int Nz,
	   std::vector<double> & cx,
	   std::vector<double> & cz) :
    CustomMesh(Lx,Nx,Lz,Nz,cx,cz) {}
  void gatherAndSortCustomNodes(double * u1,double * u2,int r) {
    CustomMesh::gatherAndSortCustomNodes(u1,u2,r); }
  void sortAndScatterCustomNodes(double * u1,double * u2,int r) {
    CustomMesh::sortAndScatterCustomNodes(u1,u2,r); }
};


int main() {

  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();

  //sleep(prank);
    
  std::cout<<"prank "<<prank<<std::endl;

  unsigned int nb_nodes_x = 3;
  double length_x = 3;

  unsigned int nb_nodes_z = 4;
  double length_z = 4.0;

  if (prank == 0)
    printf("==============\n3D mesh custom\n==============\n");
  if (psize == 1){
    std::cout<<"test check custom mesh "<<std::endl;
    std::vector<double> x_coord_tmp = {2.0,  // 9
				       0.0,  // 1
				       0.0,  // 0
				       0.0,  // 2
				       0.0,  // 3
				       1.0,  // 4
				       1.0,  // 5
				       1.0,  // 6
				       1.0,  // 7
				       2.0,  // 8
				       2.0,  // 10
				       2.0}; // 11

    std::vector<double> z_coord_tmp = {1.1,  // 9
				       1.0,  // 1
				       0.0,  // 0
				       2.0,  // 2
				       3.0,  // 3
				       0.0,  // 4
				       1.0,  // 5
				       2.0,  // 6
				       3.0,  // 7
				       0.0,  // 8
				       2.0,  // 10
				       3.0}; // 11

    bool caught_exception = true;
    try {
      // construct mesh
      CustomMesh mesh3dC(length_x,nb_nodes_x, 
			 length_z,nb_nodes_z,
			 x_coord_tmp, z_coord_tmp);
      std::cout<<"mesh done"<<prank<<std::endl;
      caught_exception = false;
    }
    catch (std::runtime_error &e) {
      std::cout << "caught exception -> success" << std::endl;
    }
    if (!caught_exception) {
      std::cout << "failed" << std::endl;
      return 1; // failure
    }
  } // end test check custom mesh
  
  std::vector<double> x_coord_tmp;
  std::vector<double> z_coord_tmp;

  int max_nb_nodes_pp = -1;

  int root = 0;
  if (psize>1)
    root = 1;
  
  if (psize==1) {
    max_nb_nodes_pp = 12;
    x_coord_tmp = {2.0,  // 9
		   0.0,  // 1
		   0.0,  // 0
		   0.0,  // 2
		   0.0,  // 3
		   1.0,  // 4
		   1.0,  // 5
		   1.0,  // 6
		   1.0,  // 7
		   2.0,  // 8
		   2.0,  // 10
		   2.0}; // 11
    
    z_coord_tmp = {1.0,  // 9
		   1.0,  // 1
		   0.0,  // 0
		   2.0,  // 2
		   3.0,  // 3
		   0.0,  // 4
		   1.0,  // 5
		   2.0,  // 6
		   3.0,  // 7
		   0.0,  // 8
		   2.0,  // 10
		   3.0}; // 11
  }
  else {
    max_nb_nodes_pp = 7;
    if (prank==0) {
      x_coord_tmp = {2.0,   // 9
		     0.0,   // 1
		     0.0,   // 0
		     0.0,   // 2
		     0.0};  // 3

      z_coord_tmp = {1.0,   // 9
		     1.0,   // 1
		     0.0,   // 0
		     2.0,   // 2
		     3.0};  // 3
    }
    if (prank==1) {
      x_coord_tmp = {1.0,   // 4
		     1.0,   // 7
		     2.0,   // 10
		     1.0,   // 5
		     1.0,   // 6
		     2.0,   // 8
		     2.0};  // 11
      
      z_coord_tmp = {0.0,   // 4
		     3.0,   // 7
		     2.0,   // 10
		     1.0,   // 5
		     2.0,   // 6
		     0.0,   // 8
		     3.0};  // 11
    }
  }

  // construct mesh
  TestMesh mesh3dC(length_x,nb_nodes_x,
		   length_z,nb_nodes_z,
		   x_coord_tmp,z_coord_tmp);

 
  // for checking afterwards
  double * coords_local[3];
  for (int d=0; d<mesh3dC.getDim();++d) {
    int size = max_nb_nodes_pp;
    coords_local[d] = new double [size];
    for (int n=mesh3dC.getNbLocalNodes(); n<max_nb_nodes_pp; ++n)
      coords_local[d][n] = NAN;
  }
  for (int n=0; n<mesh3dC.getNbLocalNodes(); ++n) {
      coords_local[0][n] = x_coord_tmp[n];
      coords_local[1][n] = 0.;
      coords_local[2][n] = z_coord_tmp[n];
  }
  
  double ** coords3dC = mesh3dC.getLocalCoords();
  print(coords3dC, mesh3dC.getNbLocalNodes());

  double * coords_global[3];
  for (int d=0; d<mesh3dC.getDim();++d) {
    int size = max_nb_nodes_pp*psize;
    coords_global[d] = new double [size];
    for (int n=0; n<size; ++n)
      coords_global[d][n] = 0.;
  }
  
  mesh3dC.gatherAndSortCustomNodes(coords_global[0],
				   coords_local[0],root);
  mesh3dC.gatherAndSortCustomNodes(coords_global[2],
				   coords_local[2],root);

  std::cout<<std::flush;
  if (prank==root) {
    std::cout<<"coords global sorted"<<std::endl;
    print(coords_global, mesh3dC.getNbGlobalNodes());
  }
  
  // test
  if (prank==root) {
    bool test_pass=true;
    for (int n=0; n<mesh3dC.getNbGlobalNodes();++n) {
      int idx=n/nb_nodes_z;
      int idz=n%nb_nodes_z;
      if (coords_global[0][n] != length_x*idx/nb_nodes_x ||
	  coords_global[2][n] != length_z*idz/nb_nodes_z )
	test_pass=false;
    }
    if (test_pass)
      std::cout<<"Success 3D mesh custom gatherAndSortNodes \n";
    else
      return 1;
  }

  double * coords_local_tmp[3];
  for (int d=0; d<mesh3dC.getDim();d++) {
    int size = max_nb_nodes_pp;
    coords_local_tmp[d] = new double[size];
    for (int n=0; n<size; ++n)
      coords_local_tmp[d][n] = 0.;
  }

  mesh3dC.sortAndScatterCustomNodes(coords_global[0],
				    coords_local_tmp[0],root);
  mesh3dC.sortAndScatterCustomNodes(coords_global[2],
				    coords_local_tmp[2],root);

  //sleep(prank);
  printf("local coords scatter\n");
  print(coords_local_tmp, mesh3dC.getNbLocalNodes());
  
  bool test_pass=true;
  for (int n=0; n<mesh3dC.getNbLocalNodes();++n) {
    if (x_coord_tmp[n] != coords_local_tmp[0][n] ||
	z_coord_tmp[n] != coords_local_tmp[2][n] ) {
      printf("%g != %g or %g != %g \n",
	     x_coord_tmp[n], coords_local_tmp[0][n],
	     z_coord_tmp[n], coords_local_tmp[2][n]);
      test_pass=false;
    }
  }
  if (test_pass)
    std::cout<<"Success 3D mesh custom sortAndScatterNodes \n";
  else {
    std::cout<<"error"<<std::endl;
    return 1;
  }

  for (int d=0; d<mesh3dC.getDim();++d) {
    delete[] coords_global[d];
    delete[] coords_local_tmp[d];
  }

  std::cout<<"went til the end rank "<<prank<<"\n";
  std::cout << "all checks passed -> overall success" << std::endl;
  StaticCommunicatorMPI::getInstance()->finalize();

  return 0;
}

