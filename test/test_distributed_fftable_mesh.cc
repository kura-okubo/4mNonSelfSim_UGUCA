/**
 * @file   test_distributed_fftable_mesh.cc
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

#include "uca_distributed_fftable_mesh.hh"
#include "fftable_nodal_field.hh"

#include <iostream>
#include <cmath>
#include <unistd.h>

using namespace uguca;

class TestMesh : public DistributedFFTableMesh {
public:
  TestMesh(double Lx, int Nx) : DistributedFFTableMesh(Lx,Nx) {}
  TestMesh(double Lx, int Nx,
	   double Lz, int Nz) : DistributedFFTableMesh(Lx,Nx,Lz,Nz) {}
  void sortAndScatterFFTModes(fftw_complex * U, int root_rank) {
    DistributedFFTableMesh::sortAndScatterFFTModes(U,root_rank);
  }
  void gatherAndSortFFTModes(fftw_complex * U, int root_rank) {
    DistributedFFTableMesh::gatherAndSortFFTModes(U,root_rank);
  }
  void printSortMap() {
    std::cout << "global=" << this->getNbGlobalFFT()
	      << " local_alloc=" << this->getNbLocalFFTAlloc() << std::endl;
    for (int i=0; i<this->getNbGlobalFFT(); ++i)
      std::cout << this->sort_fft_modes_map[i] << " ";
    std::cout << std::endl;
  }
};

int main(){

  int psize = StaticCommunicatorMPI::getInstance()->getNbProc();
  int prank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int root = 0;
  
  //sleep(prank);  
    
  unsigned int nb_nodes_x = 10;
  double length_x = 3;

  unsigned int nb_nodes_z = 4;
  double length_z = 4.0;

  int itest = -1;
  int isol = -1;
  
  // construct 2D mesh
  if (prank==0)
    printf("==============\n2D mesh\n==============\n");
  TestMesh mesh2d(length_x,nb_nodes_x);

  if (prank==0)
    std::cout << "test 2d fft data" << std::endl;

  itest = mesh2d.getNbGlobalFFT(0);
  isol  = nb_nodes_x / 2 + 1;
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb global fft x failed" << std::endl;
    return 1; // failure
  }
  
  itest = mesh2d.getNbGlobalFFT(2);
  isol  = 1;
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb global fft z failed" << std::endl;
    return 1; // failure
  }

  itest = mesh2d.getNbGlobalFFT();
  isol  = nb_nodes_x / 2 + 1;
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb global fft failed" << std::endl;
    return 1; // failure
  }

  itest = mesh2d.getNbLocalFFT();
  if (psize==1)
    isol = nb_nodes_x / 2 + 1;
  else if (prank==0)
    isol = 2;
  else if (prank==1)
    isol = 4;
  else
    throw;
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb local fft failed" << std::endl;
    return 1; // failure
  }

  itest = mesh2d.getNbLocalFFTAlloc();
  if (psize==1)
    isol = nb_nodes_x / 2 + 1;
  else
    isol = 4;
  if ((prank==0) && (itest != psize*isol)) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb local alloc fft failed (prank==0)" << std::endl;
    return 1; // failure
  }
  if ((prank>0) && (itest != isol)) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb local alloc fft failed (prank>0)" << std::endl;
    return 1; // failure
  }
  
  if (prank==0)
    std::cout << "2d fft data success!" << std::endl;


  FFTableNodalField fnfc(mesh2d);
  if (prank==0) {
    for (int i=0; i<mesh2d.getNbGlobalFFT(); ++i) {
      fnfc.fd(i)[0] = i;
      fnfc.fd(i)[1] = 100+i;
    }
  }
  fftw_complex * p_fnfc = fnfc.fd_data();

  /*
  if (prank==0)
    mesh2d.printSortMap();
  */
  
  if (false) {
    std::cout << "start" << std::endl;
    for (int i=0; i<mesh2d.getNbGlobalFFT(); ++i) {
      std::cout << prank << ": " << p_fnfc[i][0] << "," << p_fnfc[i][1] << std::endl;
    }
    std::cout << "-----------" << std::flush << std::endl;
  }
  
  mesh2d.sortAndScatterFFTModes(fnfc.fd_data(),root);

  if (false) {
    std::cout << "sorted and scattered" << std::endl;
    for (int i=0; i<mesh2d.getNbLocalFFT(); ++i) {
      std::cout << prank << ": " << p_fnfc[i][0] << "," << p_fnfc[i][1] << std::endl;
    }
    std::cout << "-----------" << std::flush << std::endl;
  }
  
  mesh2d.gatherAndSortFFTModes(fnfc.fd_data(),root);

  if (false) {
    std::cout << "gathered and sorted" << std::endl;
    for (int i=0; i<mesh2d.getNbGlobalFFT(); ++i) {
      std::cout << prank << ": " << p_fnfc[i][0] << "," << p_fnfc[i][1] << std::endl;
    }
    std::cout << "-----------" << std::endl;
  }
  
  if (prank==0) {
    for (int i=0; i<mesh2d.getNbGlobalFFT(); ++i) {
      if ((fnfc.fd(i)[0] != i) || (fnfc.fd(i)[1] != 100+i)) {
	std::cerr << "sort scatter gather sort failed" << std::endl;
	return 1; // failure
      }
    }
  }


  // construct 3D mesh
  if (prank==0)
    printf("==============\n3D mesh\n==============\n");
  TestMesh mesh3d(length_x,nb_nodes_x,
		  length_z,nb_nodes_z);
    
  if (prank==0)
    std::cout << "test 3d fft data" << std::endl;

  itest = mesh3d.getNbGlobalFFT(0);
  isol  = nb_nodes_x;
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb global fft x failed" << std::endl;
    return 1; // failure
  }
  
  itest = mesh3d.getNbGlobalFFT(2);
  isol  = nb_nodes_z / 2 + 1;
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb global fft z failed" << std::endl;
    return 1; // failure
  }

  itest = mesh3d.getNbGlobalFFT();
  isol  = nb_nodes_x * (nb_nodes_z / 2 + 1);
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb global fft failed" << std::endl;
    return 1; // failure
  }

  itest = mesh3d.getNbLocalFFT();
  if (psize==1)
    isol = nb_nodes_x * (nb_nodes_z / 2 + 1);
  else if (prank==0)
    isol = 15;
  else if (prank==1)
    isol = 15;
  else
    throw;
  if (itest != isol) {
    std::cerr << isol << " " << itest << "(prank==" << prank << ")" << std::endl;
    std::cerr << "nb local fft failed" << std::endl;
    return 1; // failure
  }

  itest = mesh3d.getNbLocalFFTAlloc();
  if (psize==1)
    isol = nb_nodes_x * (nb_nodes_z / 2 + 1);
  else if (prank==0)
    isol = 30;
  else if (prank==1)
    isol = 15;
  if (itest != isol) {
    std::cerr << isol << " " << itest << "(prank==" << prank << ")" << std::endl;
    std::cerr << "nb local alloc fft failed" << std::endl;
    return 1; // failure
  }

  if (prank==0)
    std::cout << "3d fft data success!" << std::endl;
    
  if (prank==0)
    std::cout << "all checks passed -> overall success" << std::endl;

  StaticCommunicatorMPI::getInstance()->finalize();
  return 0;
}

