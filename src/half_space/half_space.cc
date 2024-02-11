/**
 * @file   half_space.cc
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
#include "half_space.hh"
#include "half_space_dynamic.hh"

#include <limits>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
HalfSpace::HalfSpace(Material & material,
		     FFTableMesh & mesh,
		     int side_factor,
		     SpatialDirectionSet components,
		     const std::string & name) :
  name(name),
  material(material),
  mesh(mesh),
  time_step(0.),
  side_factor(side_factor),
  disp(mesh,components,name+"_disp"),
  velo(mesh,components,name+"_velo"),
  internal(mesh,components,name+"_internal"),
  residual(mesh,components,name+"_residual"),
  disp_pc(name+"_pcdisp"), // not allocated
  velo_pc(name+"_pcvelo") // not allocated
{}

/* -------------------------------------------------------------------------- */
HalfSpace::~HalfSpace() {}

/* -------------------------------------------------------------------------- */
HalfSpace * HalfSpace::newHalfSpace(Material & material,
				    FFTableMesh & mesh,
				    int side_factor,
				    SpatialDirectionSet components,
				    const std::string & name,
				    const SolverMethod & method) {
  HalfSpace * hs = NULL;
  if (method == _dynamic) {
    hs = new HalfSpaceDynamic(material, mesh, side_factor, components, name);
  }
  else {
    throw std::runtime_error("HalfSpace: solver method not implemented");
  }
  return hs;
}

/* -------------------------------------------------------------------------- */
void HalfSpace::initPredictorCorrector() {
  this->predictor_corrector = true;
  this->disp_pc.resize(this->mesh,this->disp.getComponents());
  this->velo_pc.resize(this->mesh,this->velo.getComponents());
}

/* -------------------------------------------------------------------------- */
double HalfSpace::getStableTimeStep() {
#ifdef UCA_VERBOSE
  std::cout << "getStableTimeStep is not implemented for this HalfSpace" << std::endl;
#endif
  return std::numeric_limits<double>::max();
}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeDisplacement(bool predicting) {
  this->computeDisplacement(this->disp,
			    predicting ? this->velo_pc : this->velo,
			    predicting ? this->disp_pc : this->disp);
}

/* -------------------------------------------------------------------------- */
// u_i+1 = u_i + dt * v_i
void HalfSpace::computeDisplacement(NodalField & disp,
				    NodalField & velo,
				    NodalField & target) {

  for (const auto& d : target.getComponents()) {

    double * disp_p   = disp.data(d);
    double * velo_p   = velo.data(d);
    double * target_p = target.data(d);

    for (int n=0; n<target.getNbNodes(); ++n) {
      target_p[n] = disp_p[n] + velo_p[n] * this->time_step;
    }
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeInternal(bool predicting, bool correcting, bool dynamic) {
  this->forwardFFT(predicting);
  this->computeStressFourierCoeff(predicting, correcting, dynamic);
  this->backwardFFT();
}

/* -------------------------------------------------------------------------- */
void HalfSpace::forwardFFT(bool predicting) {
  (predicting ? this->disp_pc : this->disp).forwardFFT();
}

/* -------------------------------------------------------------------------- */
void HalfSpace::backwardFFT() {
  this->internal.backwardFFT();
}

/* -------------------------------------------------------------------------- */
// residual = (internal + external) * side_factor
void HalfSpace::computeResidual(NodalField & external) {

  for (const auto& d : this->residual.getComponents()) {

    double *int_p = this->internal.data(d);
    double *ext_p = external.data(d);
    double *res_p = this->residual.data(d);

    for (int n = 0; n < this->residual.getNbNodes(); ++n) {
      res_p[n] = this->side_factor * (int_p[n] + ext_p[n]);
    }
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeVelocity(bool predicting) {
  this->computeVelocity(predicting ? this->velo_pc : this->velo);
}

/* -------------------------------------------------------------------------- */
void HalfSpace::updateVelocity() {

  for (const auto& d : this->velo.getComponents()) {

    double *velo_p = this->velo.data(d);
    double *velo_pc_p = this->velo_pc.data(d);

    for (int n = 0; n < this->velo_pc.getNbNodes(); ++n) {
      velo_pc_p[n] = velo_p[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::correctVelocity(bool last_step) {
  this->correctVelocity(this->velo, this->velo_pc,
			last_step ? this->velo : this->velo_pc);
}

/* -------------------------------------------------------------------------- */
void HalfSpace::correctVelocity(NodalField & velo_n,
				NodalField & velo_pc,
				NodalField & target) {

  for (const auto& d : target.getComponents()) {

    double * velo_n_p = velo_n.data(d);
    double * velo_pc_p = velo_pc.data(d);
    double * target_p = target.data(d);

    for (int n = 0; n < target.getNbNodes(); ++n) {
      target_p[n] = 0.5 * (velo_n_p[n] + velo_pc_p[n]);
    }
  }
}

/* -------------------------------------------------------------------------- */
// velocity = cs / mu       * residual (for in-plane shear components)
// velocity = cs / mu / eta * residual (for normal component)
// velocity = cs / mu       * residual (for out-of-plane shear components)
void HalfSpace::computeVelocity(NodalField & velo) {
  double mu = this->material.getShearModulus();
  double Cs = this->material.getCs();
  double Cp = this->material.getCp();
  std::vector<double> eta = {1.0, Cp / Cs, 1.0};

  for (const auto& d : velo.getComponents()) {

    double * velo_p = velo.data(d);
    double * res_p = this->residual.data(d);
    double eta_d = eta[d];

    for (int n=0; n<velo.getNbNodes(); ++n)
      velo_p[n] = Cs / mu / eta_d * res_p[n];
  }
}

/* -------------------------------------------------------------------------- */
bool HalfSpace::registerDumpFieldToDumper(const std::string & field_name,
					  const std::string & dump_name,
					  Dumper * const dumper) {
  
  // disp
  if (field_name == "disp") {
    dumper->registerIO(dump_name, this->disp);
    return true;
  }
  // velo
  else if (field_name == "velo") {
    dumper->registerIO(dump_name, this->velo);
    return true;
  }
  // residual
  else if (field_name == "residual") {
    dumper->registerIO(dump_name, this->residual);
    return true;
  }
  // internal
  else if (field_name == "internal") {
    dumper->registerIO(dump_name, this->internal);
    return true;
  }

  return false;
}

/* -------------------------------------------------------------------------- */
void HalfSpace::registerToRestart(Restart & restart) {

  restart.registerIO(this->disp);
  restart.registerIO(this->velo);
  restart.registerIO(this->internal);
  restart.registerIO(this->residual);

  if (this->predictor_corrector) {
    restart.registerIO(this->disp_pc);
    restart.registerIO(this->velo_pc);
  }

}

__END_UGUCA__
