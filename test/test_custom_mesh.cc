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
#include <stdexcept>
#include <unistd.h>

using namespace uguca;

void print(const TwoDVector & fld, int size){
  int dim = 3;//fld.size();
  for (int n=0; n<size;n++) {
    for (int d=0;d<dim;d++) {
      std::cout<<fld(n,d)<<", ";
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
  TwoDVector coords_local(3,max_nb_nodes_pp);
  for (int d=0; d<mesh3dC.getDim();++d) {
    for (int n=mesh3dC.getNbLocalNodes(); n<max_nb_nodes_pp; ++n)
      coords_local(n,d) = NAN;
  }
  for (int n=0; n<mesh3dC.getNbLocalNodes(); ++n) {
    coords_local(n,0) = x_coord_tmp[n];
    coords_local(n,1) = 0.;
    coords_local(n,2) = z_coord_tmp[n];
  }
  
  const TwoDVector & coords3dC = mesh3dC.getLocalCoords();
  print(coords3dC, mesh3dC.getNbLocalNodes());

  TwoDVector  coords_global(3,max_nb_nodes_pp*psize);
  mesh3dC.gatherAndSortCustomNodes(coords_global.data(0).data(),
				   coords_local.data(0).data(),root);
  mesh3dC.gatherAndSortCustomNodes(coords_global.data(2).data(),
				   coords_local.data(2).data(),root);

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
      if (coords_global(n,0) != length_x*idx/nb_nodes_x ||
	  coords_global(n,2) != length_z*idz/nb_nodes_z )
	test_pass=false;
    }
    if (test_pass)
      std::cout<<"Success 3D mesh custom gatherAndSortNodes \n";
    else
      return 1;
  }

  TwoDVector coords_local_tmp(3,max_nb_nodes_pp);
  for (int d=0; d<mesh3dC.getDim();d++) {
    for (int n=0; n<max_nb_nodes_pp; ++n)
      coords_local_tmp(n,d) = 0.;
  }

  mesh3dC.sortAndScatterCustomNodes(coords_global.data(0).data(),
				    coords_local_tmp.data(0).data(),root);
  mesh3dC.sortAndScatterCustomNodes(coords_global.data(2).data(),
				    coords_local_tmp.data(2).data(),root);

  //sleep(prank);
  printf("local coords scatter\n");
  print(coords_local_tmp, mesh3dC.getNbLocalNodes());
  
  bool test_pass=true;
  for (int n=0; n<mesh3dC.getNbLocalNodes();++n) {
    if (x_coord_tmp[n] != coords_local_tmp(n,0) ||
	z_coord_tmp[n] != coords_local_tmp(n,2) ) {
      printf("%g != %g or %g != %g \n",
	     x_coord_tmp[n], coords_local_tmp(n,0),
	     z_coord_tmp[n], coords_local_tmp(n,2));
      test_pass=false;
    }
  }
  if (test_pass)
    std::cout<<"Success 3D mesh custom sortAndScatterNodes \n";
  else {
    std::cout<<"error"<<std::endl;
    return 1;
  }

  std::cout<<"went til the end rank "<<prank<<"\n";
  std::cout << "all checks passed -> overall success" << std::endl;
  StaticCommunicatorMPI::getInstance()->finalize();

  return 0;
}

