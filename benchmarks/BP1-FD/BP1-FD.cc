/**
 * @file   BP1-FD.cc
 *
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Wed Mar 6 2024
 * @date last modification: Wed Mar 6 2024
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

#include <cmath>
#include <sys/time.h>
#include <unistd.h>

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "material.hh"
#include "uca_simple_mesh.hh"
#include "rate_and_state_law.hh"
#include "static_communicator_mpi.hh"
#include "unimat_shear_interface.hh"

using namespace uguca;

int main(int argc, char *argv[]) {

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  // ---------------------------------------------------------------------------
  // default parameters

  double length_x_rpt = 50e3;

  double domain_factor = 2.0;

  double duration = 47304000000; // 1500 years
  double dump_int = 1e6;

  unsigned nb_nodes_x = 2000;  // 50e3 / 25 = 2000
  double time_step_factor = 0.35;

  unsigned n_pc = 1;

  // ---------------------------------------------------------------------------
  // argument processing

  int c;
  extern char* optarg;
  while ((c = getopt(argc, argv, "hN:T:t:s:f:p:")) != -1) {
    switch (c) {
    case 'h':
      fprintf(stderr,
              "%s\n"
              "\t-h: print this message\n"
              "\t-N: number of elements power of 2 (%d)\n"
              "\t-T: duration (%f)\n"
              "\t-t: dump interval (seconds) (%f)\n"
              "\t-s: factor of domain size (%f)\n"
              "\t-f: time step factor (%f)\n"
              "\t-p: predictor-corrector iterations (%d)\n",
              argv[0], nb_nodes_x, duration, dump_int, domain_factor,
              time_step_factor, n_pc);
      return -1;
    case 'N': nb_nodes_x       = atoi(optarg); break;
    case 'T': duration         = atof(optarg); break;
    case 't': dump_int         = atof(optarg); break;
    case 's': domain_factor    = atof(optarg); break;
    case 'f': time_step_factor = atof(optarg); break;
    case 'p': n_pc             = atoi(optarg); break;

    default:
      fprintf(stderr, "Unknown option (-%c)\n", c);
      return -1;
    }
  }

  double length_x = length_x_rpt * domain_factor;

  // ---------------------------------------------------------------------------
  // problem parameters

  // material
  double Cp  = 6000.0;
  double Cs  = 3464.0;
  double rho = 2670.0;

  // rate and state
  double a0 = 0.010;
  double a_max = 0.025;
  double b0 = 0.015;
  double Dc = 0.008;
  double V0 = 1.0e-6;
  double f0 = 0.6;

  // initial conditions
  double sigma_n = 50.0e6;
  double V_init = 1.0e-9;
  double V_p = 1.0e-9;
  double theta_init = 1.606238999213454e9;
  double H = 15e3;
  double h = 3e3;
  double Wf = 40e3;

  // ---------------------------------------------------------------------------
  // mesh
  SimpleMesh mesh(length_x, nb_nodes_x);

  // constitutive interface law
  SpatialDirection slip_dir = _z;
  RateAndStateLaw law(mesh, a_max, b0, Dc, V0, f0, theta_init,
                      RateAndStateLaw::EvolutionLaw::AgingLaw, n_pc > 0, 0.0, slip_dir);
  NodalField & theta = law.getTheta();
  NodalField & a = law.getA();
  NodalField & b = law.getB();

  double mu = Cs * Cs * rho;
  double lambda = Cp * Cp * rho - 2.0 * mu;
  double nu = 0.5 * (lambda / (lambda + mu));
  double E = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
  if (world_rank == 0)  printf("E=%g\nnu=%g\n", E, nu);

  Material mat = Material(E,nu,rho);
  mat.readPrecomputedKernels();
  // Material bot_mat = Material(E,nu,rho);
  if ((std::abs((mat.getCp()-Cp))>1e-15*Cp) ||
      (std::abs((mat.getCs()-Cs))>1e-15*Cs))
    return -1;

  // ---------------------------------------------------------------------------
  // weak interface

  UnimatShearInterface interface(mesh, {_z}, mat, law);

  // ---------------------------------------------------------------------------
  // initial conditions

  // init external load

  NodalField & external = interface.getLoad();
  double tau0 = sigma_n * a_max * std::asinh(V_init / V0 / 2.0 * std::exp((f0 + b0 * std::log(V0 / V_init))/a_max));
  external.setAllValuesTo(tau0, slip_dir);

  // impose a constant contact pressure
  NodalField & sigma = law.getConstantPressure();
  sigma.setAllValuesTo(sigma_n,1);

  // init velocity
  HalfSpace& top = interface.getTop();
  NodalField& velo_top = top.getVelo();
  velo_top.setAllValuesTo(V_init / 2, 2);

  const TwoDVector &  coords = mesh.getLocalCoords();

  // init a
  for (int  i = 0; i < mesh.getNbLocalNodes(); ++i) {
    double x = std::abs(coords(i,0) - length_x / 2);
    if (x < H) {
      a(i) = a0;
    } else if (x < H + h) {
      a(i) = a0 + (a_max - a0) * (x - H) / h;
    }
  }

  // init theta
  for (int  i = 0; i < mesh.getNbLocalNodes(); ++i) {
    theta(i) = Dc / V0 * std::exp(a(i) / b(i) * std::log(2.0 * V0 / V_init * std::sinh(tau0 / a(i) / sigma_n)) - f0 / b(i));
  }

  // time step
  double time_step = time_step_factor * interface.getStableTimeStep();
  interface.setTimeStep(time_step);
  unsigned nb_time_steps = std::ceil(duration / time_step);

  // init interface
  interface.initPredictorCorrector(n_pc);
  law.init();
  interface.init(true);

  // ---------------------------------------------------------------------------
  // dumping
  if (world_rank == 0) std::cout << "dump int = " << dump_int << std::endl;

  std::ostringstream bname_out;
  bname_out << std::fixed << std::setprecision(2)
            << "BP1-FD_Nx" << nb_nodes_x
            << "_Nx" << nb_nodes_x
            << "_s" << domain_factor
            << "_tf" << time_step_factor
            << "_npc" << n_pc;
  std::string bname = bname_out.str();

  if (world_rank == 0) std::cout << bname << std::endl;

  interface.initDump(bname, ".", Dumper::Format::Binary);

  interface.registerDumpField("cohesion");
  interface.registerDumpField("top_disp");
  interface.registerDumpField("top_velo");
  interface.registerDumpField("theta");

  // interface.registerDumpField("iterations");
  // interface.registerDumpField("rel_error");
  // interface.registerDumpField("a");
  // interface.registerDumpField("b");

  interface.dump(0, 0);
  unsigned s_dump = dump_int / time_step + 1;

  NodalField &u_top = top.getDisp();

  if (world_rank == 0) std::cout << "simulation start..." << std::endl;

  // time stepping
  for (unsigned s = 1; s <= nb_time_steps; ++s) {
    if (world_rank == 0) {
      std::cout << "s=" << s << "/" << nb_time_steps << "\r";
      std::cout.flush();
    }

    if (world_rank == mesh.getRoot()) {
      // only works with SimpleMesh
      int nb_nodes_x = mesh.getNbGlobalNodes(0);
      for (int i = 1; i < nb_nodes_x / 2; ++i) {
        // plate rate
        double x = std::abs(coords(i, 0) - length_x / domain_factor);
        if (x > Wf) {
          velo_top(i,2) = V_p;
        }
        // free surface
        u_top(nb_nodes_x - i, 2) = u_top(i, 2);
        velo_top(nb_nodes_x - i, 2) = velo_top(i, 2);
      }
    }


    // time integration
    interface.advanceTimeStep();

    // dump
    if (world_rank == 0 && s % s_dump == 0) interface.dump(s, s * time_step);
  }

  StaticCommunicatorMPI::getInstance()->finalize();

  if (world_rank == 0)
    std::cout << "uguca simulation completed." << std::endl;

  return 0;
}
