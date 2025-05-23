/**
 * @file   test_interface.cc
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

#include "half_space_dynamic.hh"
#include "interface.hh"
#include "interface_law.hh"
#include "linear_shear_cohesive_law.hh"
#include "uca_simple_mesh.hh"
#include "nodal_field.hh"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using namespace uguca;

class TestInterfaceLaw : public InterfaceLaw {
public:
  TestInterfaceLaw(SimpleMesh & mesh): InterfaceLaw(mesh) {}
  void computeCohesiveForces(NodalField &,
			     bool = false) {
    computeCohesiveForcesCalled = true;
  }
  void registerDumpField(const std::string &) {
    registerDumpFieldCalled = true;
  }
  bool computeCohesiveForcesCalled = false;
  bool registerDumpFieldCalled = false;
};

class TestHalfSpace : public HalfSpaceDynamic {
public:
  TestHalfSpace(SimpleMesh & mesh, int side_factor) : HalfSpaceDynamic(mesh, side_factor) {}
  void computeDisplacement(bool = false) {
    computeDisplacementCalled = true;
  }
  void computeStressFourierCoeff(bool = false,
                                bool = false) {
    computeStressFourierCoeffCalled = true;
  }
  void computeResidual(NodalField &) {
    computeResidualCalled = true;
  }
  void computeVelocity(bool = false) {
    computeVelocityCalled = true;
  }
  void forwardFFT(bool = false) {
    forwardFFTCalled = true;
  }
  void backwardFFT() {
    backwardFFTCalled = true;
  }
  void gatherCostumMeshForwardFFT(NodalField &,
                                  bool = false) {
    gatherCostumMeshForwardFFTCalled = true;
  }
  void backwardFFTscatterCostumMesh() {
    backwardFFTscatterCostumMeshCalled = true;
  }
  void updateVelocity() {
    updateVelocityCalled = true;
  }
  void correctVelocity(bool) {
    correctVelocityCalled = true;
  }
  bool getPredictorCorrector() { return predictor_corrector; }
  bool computeDisplacementCalled = false;
  bool computeStressFourierCoeffCalled = false;
  bool computeResidualCalled = false;
  bool computeVelocityCalled = false;
  bool forwardFFTCalled = false;
  bool backwardFFTCalled = false;
  bool gatherCostumMeshForwardFFTCalled = false;
  bool backwardFFTscatterCostumMeshCalled = false;
  bool updateVelocityCalled = false;
  bool correctVelocityCalled = false;
};

class TestInterface : public Interface {
public:
 TestInterface(SimpleMesh &mesh, InterfaceLaw &law, Material &material)
     : Interface(mesh, law), top(mesh, 1) {
   half_spaces.push_back(&top);
   top.setMaterial(&material);
   constructed = true;
 }
  // pure virtual functions
  void closingNormalGapForce(NodalFieldComponent &, bool = false) {
    throw "TestInterface::closingNormalGapForce() is not implemented.";
  }
  void maintainShearGapForce(NodalField &) {
    throw "TestInterface::maintainShearGapForce() is not implemented.";
  }
  void computeGap(NodalField &, bool = false) {
    throw "TestInterface::computeGap() is not implemented.";
  }
  void computeGapVelocity(NodalField &,
                          bool = false) {
    throw "TestInterface::computeGapVelocity() is not implemented.";
  }
  void initDump(const std::string &, const std::string &,
                const Dumper::Format) {
    initDumpCalled = true;
  }
  HalfSpace &getTop() { throw "TestInterface::getTop() is not implemented."; }
  HalfSpace &getBot() { throw "TestInterface::getBot() is not implemented."; }
  TestHalfSpace &getTestTop() { return this->top; }
  bool constructed = false;
  bool initDumpCalled = false;

 protected:
  TestHalfSpace top;
};

int main(){

  std::cout << "start test: test_interface" << std::endl;

  // --------------------------------------------------------------
  // INITIALIZE
  // --------------------------------------------------------------
  // mesh
  // --------------------------------------------------------------
  double Lx = 0.4;
  int Nx = 3;
  double Lz = 0.3;
  int Nz = 2;
  SimpleMesh mesh(Lx, Nx, Lz, Nz);
  // --------------------------------------------------------------
  // interface law
  // --------------------------------------------------------------
  TestInterfaceLaw law(mesh);
  // --------------------------------------------------------------
  // material
  // --------------------------------------------------------------
  double E = 2.1;
  double nu = 0.25;
  double rho = 10.0;
  Material material(E, nu, rho, false);
  material.readPrecomputedKernels();

  // Check Interface::Interface
  // --------------------------------------------------------------
  std::cout << "check Interface::Interface" << std::endl;
  TestInterface interface(mesh, law, material);
  if (!interface.constructed) {
    std::cout << "Interface::Interface failed" << std::endl;
    return 1;  // failed
  } else {
    std::cout << "Interface::Interface correct -> success" << std::endl;
  }
  // --------------------------------------------------------------

  // Check Interface::setTimeStep
  // --------------------------------------------------------------
  std::cout << "check Interface::setTimeStep & Interface::getTimeStep"
            << std::endl;
  double time_step = 1.0e-1;
  interface.setTimeStep(time_step);
  TestHalfSpace & top = interface.getTestTop();
  if (interface.getTimeStep() == time_step) {
    std::cout
        << "Interface::setTimeStep or Interface::getTimeStep correct -> success"
        << std::endl;
  } else {
    std::cout << "Interface::setTimeStep or Interface::getTimeStep failed"
              << std::endl;
    return 1;  // failed
  }
  // --------------------------------------------------------------

  // Check Interface::initPredictorCorrector
  // --------------------------------------------------------------
  std::cout << "check Interface::initPredictorCorrector" << std::endl;
  interface.initPredictorCorrector(1);
  if (top.getPredictorCorrector()) {
    std::cout << "Interface::initPredictorCorrector correct -> success"
              << std::endl;
  } else {
    std::cout << "Interface::initPredictorCorrector failed" << std::endl;
    return 1;  // failed
  }
  // --------------------------------------------------------------

  // Check Interface::init
  // --------------------------------------------------------------
  std::cout << "check Interface::init" << std::endl;
  interface.init();
  if (top.updateVelocityCalled &&
      top.computeDisplacementCalled &&
      top.computeResidualCalled &&
      top.computeVelocityCalled &&
      top.forwardFFTCalled &&
      top.computeStressFourierCoeffCalled &&
      top.backwardFFTCalled &&
      top.computeResidualCalled &&
      law.computeCohesiveForcesCalled) {
    std::cout << "check Interface::init correct -> success" << std::endl;
    top.updateVelocityCalled = false;
    top.computeDisplacementCalled = false;
    top.computeResidualCalled = false;
    top.computeVelocityCalled = false;
    top.forwardFFTCalled = true;
    top.computeStressFourierCoeffCalled = true;
    top.backwardFFTCalled = true;
    top.computeResidualCalled = true;
    law.computeCohesiveForcesCalled = false;
  } else {
    std::cout << "check Interface::init failed" << std::endl;
    return 1;  // failed
  }
  // --------------------------------------------------------------

  // Check Interface::registerDumpField
  // --------------------------------------------------------------
  std::cout << "check Interface::registerDumpField" << std::endl;
  interface.registerDumpField("nonexisting_field");
  if (law.registerDumpFieldCalled) {
    std::cout << "Interface::registerDumpField correct -> success" << std::endl;
  } else {
    std::cout << "Interface::registerDumpField failed" << std::endl;
    return 1;  // failed
  }
  // --------------------------------------------------------------

  // Check Interface::initDump
  // --------------------------------------------------------------
  std::cout << "check Interface::initDump" << std::endl;
  interface.initDump("test", "test", Dumper::Format::ASCII);
  if (!interface.initDumpCalled) {
    std::cout << "Interface::initDump failed" << std::endl;
    return 1;  // failed
  } else {
    std::cout << "Interface::initDump correct -> success" << std::endl;
  }
  // --------------------------------------------------------------

  // Check Interface::getStableTimeStep
  // --------------------------------------------------------------
  std::cout << "check Interface::getStableTimeStep" << std::endl;
  if (std::abs(interface.getStableTimeStep() - 0.46004370622823615) < 1.0e-17) {
    std::cout << "Interface::getStableTimeStep correct -> success" << std::endl;
  } else {
    std::cout << "Interface::getStableTimeStep failed" << std::endl;
    return 1;  // failed
  }
  // --------------------------------------------------------------

  // Check Interface::advanceTimeStep
  // --------------------------------------------------------------
  std::cout << "check Interface::advanceTimeStep" << std::endl;
  interface.init();
  if (top.updateVelocityCalled &&
      top.computeDisplacementCalled &&
      top.computeResidualCalled &&
      top.computeVelocityCalled &&
      top.forwardFFTCalled &&
      top.computeStressFourierCoeffCalled &&
      top.backwardFFTCalled &&
      top.computeResidualCalled &&
      law.computeCohesiveForcesCalled) {
    std::cout << "check Interface::advanceTimeStep correct -> success"
              << std::endl;
    top.updateVelocityCalled = false;
    top.computeDisplacementCalled = false;
    top.computeResidualCalled = false;
    top.computeVelocityCalled = false;
    top.forwardFFTCalled = true;
    top.computeStressFourierCoeffCalled = true;
    top.backwardFFTCalled = true;
    top.computeResidualCalled = true;
    law.computeCohesiveForcesCalled = false;
  } else {
    std::cout << "check Interface::advanceTimeStep failed" << std::endl;
    return 1;  // failed
  }
  // --------------------------------------------------------------

  // Check Interface::combineLoadAndCohesion
  // --------------------------------------------------------------
  NodalField field(mesh);
  field.component(0).setAllValuesTo(1);
  field.component(1).setAllValuesTo(2);
  field.component(2).setAllValuesTo(3);

  std::cout << "check Interface::combineLoadAndCohesion" << std::endl;
  for (int i = 0; i < 3; ++i) {
    interface.getLoad().component(i).setAllValuesTo(3);
    interface.getCohesion().component(i).setAllValuesTo(42);
  }
  interface.combineLoadAndCohesion(field);
  bool checks_out = true;
  double tol = 1e-16;
  double ref_ans = 3 - 42;
  for (size_t j = 0; j < 3; ++j) {
    for (int i = 0; i < mesh.getNbLocalNodes(); ++i) {
      if (std::abs(field.component(j).at(i) - ref_ans) > tol) {
        checks_out = false;
        break;
      }
    }
  }
  if (checks_out) {
    std::cout << "Interface::combineLoadAndCohesion correct -> success"
              << std::endl;
  } else {
    std::cout << "Interface::combineLoadAndCohesion failed" << std::endl;
    return 1;  // failed
  }
  // --------------------------------------------------------------

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
