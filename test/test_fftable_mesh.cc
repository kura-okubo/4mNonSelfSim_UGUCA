/**
 * @file   test_fftable_mesh.cc
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

#include "uca_fftable_mesh.hh"

using namespace uguca;

class TestMesh : public FFTableMesh {
public:
  TestMesh(double Lx, int Nx, bool init) : FFTableMesh(Lx,Nx,init) {}
  TestMesh(double Lx, int Nx,
	   double Lz, int Nz, bool init) : FFTableMesh(Lx,Nx,Lz,Nz,init) {}
  void initSpectralSpace() { FFTableMesh::initSpectralSpace(); }
  void allocateSpectralSpace() { FFTableMesh::allocateSpectralSpace(); }
  void freeSpectralSpace() { FFTableMesh::freeSpectralSpace(); }
  void initWaveNumbersGlobal(double ** wn) {
    FFTableMesh::initWaveNumbersGlobal(wn);
  }
};

int main(){
  
  unsigned int nb_nodes_x = 10;
  double length_x = 3;

  unsigned int nb_nodes_z = 4;
  double length_z = 4.0;

  int itest = -1;
  int isol = -1;
  
  // construct 2D mesh
  printf("==============\n2D mesh\n==============\n");
  TestMesh mesh2d(length_x,nb_nodes_x,false); // do not initialize

  std::cout << "test 2d fft data" << std::endl;

  // instead of calling .init() we do it step-wise
  mesh2d.initSpectralSpace();

  itest = mesh2d.getNbGlobalFFTX();
  isol  = nb_nodes_x / 2 + 1;
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb global fft x failed" << std::endl;
    return 1; // failure
  }
  
  itest = mesh2d.getNbGlobalFFTZ();
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
  isol  = nb_nodes_x / 2 + 1;
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb fft x failed" << std::endl;
    return 1; // failure
  }

  std::cout << "2d fft data success!" << std::endl;
  
  std::cout << "test 2d allocate" << std::endl;

  mesh2d.allocateSpectralSpace();
  
  bool caught_exception = true;
  try {
    // should not allocate a second time
    mesh2d.allocateSpectralSpace();
    caught_exception = false;
  }
  catch (std::runtime_error &e) {
    std::cout << "caught exception -> success" << std::endl;
  }
  if (!caught_exception) {
    std::cerr << "double allocate failed" << std::endl;
    return 1; // failure
  }
  
  std::cout << "test 2d allocate success!" << std::endl;

  std::cout << "test 2d check wave numbers global" << std::endl;

  // vector to fill with wave numbers
  int dim = mesh2d.getDim();
  int size = mesh2d.getNbLocalFFT();
  double * vec[3];
  for (int d=0; d<dim; ++d) {
    vec[d] = new double [size];
    for (int n=0; n<size; ++n) {
      vec[d][n] = -1.; // initialize to negative value
    }
  }

  // solution
  double * sol2D[2];
  for (int d=0; d<mesh2d.getDim(); ++d) {
    sol2D[d] = new double [mesh2d.getNbLocalFFT()];
  }

  double k1 = 2*M_PI / length_x;
  for (int i=0;  i<mesh2d.getNbLocalFFT(); ++i) {
    sol2D[0][i] = i*k1;
    sol2D[1][i] = 0.;
  }
  
  // fill vector with wave numbers (global=local)
  mesh2d.initWaveNumbersGlobal(vec);

  // test
  for (int d=0; d<mesh2d.getDim(); ++d) {
    for (int i=0; i<mesh2d.getNbLocalFFT(); ++i) {
      std::cout << vec[d][i] << " ";
      if (std::abs(vec[d][i] - sol2D[d][i]) > 1e-4) {
	std::cerr << std::endl << std::abs(vec[d][i] - sol2D[d][i]) << std::endl;
	std::cerr << std::endl << "(" << d << "," << i << ") = "
		  << vec[d][i] << " (should be = "
		  << sol2D[d][i] << ")" << std::endl;
	return 1; // failure
      }
    }
    std::cout << std::endl;
  }
  
  std::cout << "test 2d check wave numbers global success!" << std::endl;
  
  std::cout << "test 2d free" << std::endl;
  
  mesh2d.freeSpectralSpace();
  mesh2d.freeSpectralSpace(); // should do just nothing

  std::cout << "test 2d free success!" << std::endl;

  // forwardFFT and backwardFFT
  // are tested by test_fftable_nodal_field_component.cc

  for (int d=0; d<dim; ++d) {
    delete[] vec[d];
  }

  // construct 3D mesh
  printf("==============\n3D mesh\n==============\n");
  TestMesh mesh3d(length_x,nb_nodes_x,
		  length_z,nb_nodes_z, false); // do not initialize

  std::cout << "test 3d fft data" << std::endl;
  
  // instead of calling .init() we do it step-wise
  mesh3d.initSpectralSpace();

  itest = mesh3d.getNbGlobalFFTX();
  isol  = nb_nodes_x;
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb global fft x failed" << std::endl;
    return 1; // failure
  }
  
  itest = mesh3d.getNbGlobalFFTZ();
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
  isol  = nb_nodes_x * (nb_nodes_z / 2 + 1);
  if (itest != isol) {
    std::cerr << isol << " " << itest << std::endl;
    std::cerr << "nb fft x failed" << std::endl;
    return 1; // failure
  }

  std::cout << "3d fft data success!" << std::endl;
  
  std::cout << "test 3d allocate" << std::endl;

  mesh3d.allocateSpectralSpace();
  
  caught_exception = true;
  try {
    // should not allocate a second time
    mesh3d.allocateSpectralSpace();
    caught_exception = false;
  }
  catch (std::runtime_error &e) {
    std::cout << "caught exception -> success" << std::endl;
  }
  if (!caught_exception) {
    std::cerr << "double allocate failed" << std::endl;
    return 1; // failure
  }
  
  std::cout << "test 3d allocate success!" << std::endl;

  std::cout << "test 3d check wave numbers global" << std::endl;

  // vector to fill with wave numbers
  dim = mesh3d.getDim();
  size = mesh3d.getNbLocalFFT();
  for (int d=0; d<dim; ++d) {
    vec[d] = new double [size];
    for (int n=0; n<size; ++n) {
      vec[d][n] = -1.; // initialize to negative value
    }
  }

  // solution
  std::vector<std::vector<double> > sol3D(3);
  for (int d=0; d<mesh3d.getDim(); ++d) {
    sol3D[d].resize(mesh3d.getNbLocalFFT());
  }

  double m1 = 2*M_PI / length_z;
  int f_ny_x = mesh3d.getNbGlobalFFTX()/2+1;
  for (int i=0;  i<mesh3d.getNbGlobalFFTX(); ++i) {
    for (int j=0;  j<mesh3d.getNbGlobalFFTZ(); ++j) {
      int ij = i*mesh3d.getNbGlobalFFTZ() + j;
      sol3D[0][ij] = k1*(i - (i/f_ny_x)*mesh3d.getNbGlobalFFTX());
      sol3D[1][ij] = 0.;
      sol3D[2][ij] = m1*j;
    }
  }
  
  // fill vector with wave numbers (global=local)
  mesh3d.initWaveNumbersGlobal(vec);

  // test
  for (int d=0; d<mesh3d.getDim(); ++d) {
    for (int i=0; i<mesh3d.getNbLocalFFT(); ++i) {
      std::cout << vec[d][i] << " ";
      if (std::abs(vec[d][i] - sol3D[d][i]) > 1e-4) {
	std::cerr << std::endl << std::abs(vec[d][i] - sol3D[d][i]) << std::endl;
	std::cerr << std::endl << "(" << d << "," << i << ") = "
		  << vec[d][i] << " (should be = "
		  << sol3D[d][i] << ")" << std::endl;
	return 1; // failure
      }
    }
    std::cout << std::endl;
  }
  
  std::cout << "test 3d check wave numbers global success!" << std::endl;
  
  std::cout << "test 3d free" << std::endl;
  
  mesh3d.freeSpectralSpace();
  mesh3d.freeSpectralSpace(); // should do just nothing

  std::cout << "test 3d free success!" << std::endl;

  // forwardFFT and backwardFFT
  // are tested by test_fftable_nodal_field_component.cc

  for (int d=0; d<dim; ++d) {
    delete[] vec[d];
  }
    
  std::cout << "all checks passed -> overall success" << std::endl;

  return 0;
}

