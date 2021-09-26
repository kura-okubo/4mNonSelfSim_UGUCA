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
HalfSpace::HalfSpace(FFTableMesh & mesh,
		     int side_factor,
		     const std::string & name) :
  name(name),
  mesh(mesh),
  time_step(0.),
  side_factor(side_factor),
  disp(mesh,name+"_disp"),
  velo(mesh,name+"_velo"),
  internal(mesh,name+"_internal"),
  residual(mesh,name+"_residual") {}

/* -------------------------------------------------------------------------- */
HalfSpace::~HalfSpace() {}

/* -------------------------------------------------------------------------- */
HalfSpace * HalfSpace::newHalfSpace(FFTableMesh & mesh,
				    int side_factor,
				    const std::string & name,
				    const SolverMethod & method) {
  HalfSpace * hs = NULL;
  if (method == _dynamic) {
    hs = new HalfSpaceDynamic(mesh, side_factor, name);
  }
  else {
    throw std::runtime_error("HalfSpace: solver method not implemented");
  }
  return hs;
}

/* -------------------------------------------------------------------------- */
void HalfSpace::initPredictorCorrector() {
  this->predictor_corrector = true;
  this->disp_pc.init(this->mesh);
  this->velo_pc.init(this->mesh);
}

/* -------------------------------------------------------------------------- */
/*double HalfSpace::getStableTimeStep() {
#ifdef UCA_VERBOSE
  std::cout << "getStableTimeStep is not implemented for this HalfSpace" << std::endl;
#endif
  return std::numeric_limits<double>::max();
}*/

/* -------------------------------------------------------------------------- */
void HalfSpace::computeDisplacement(bool predicting) {
  this->computeDisplacement(this->disp,
			    predicting ? this->velo_pc : this->velo,
			    predicting ? this->disp_pc : this->disp);
}

/* -------------------------------------------------------------------------- */
// u_i+1 = u_i + dt * v_i
void HalfSpace::computeDisplacement(NodalField & _disp,
				    NodalField & _velo,
				    NodalField & target) {

  for (int d = 0; d < this->mesh.getDim(); ++d) {

    double * disp_p = _disp.storage(d);
    double * velo_p = _velo.storage(d);
    double * target_p = target.storage(d);

    for (int n=0; n<target.getNbNodes(); ++n) {
      target_p[n] = disp_p[n] + velo_p[n] * this->time_step;
    }
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeInternal(bool predicting, bool correcting) {
  this->forwardFFT(predicting);
  this->computeStressFourierCoeff(predicting, correcting);
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
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    double *int_p = this->internal.storage(d);
    double *ext_p = external.storage(d);
    double *res_p = this->residual.storage(d);

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
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    double *velo_p = this->velo.storage(d);
    double *velo_pc_p = this->velo_pc.storage(d);
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

  for (int d = 0; d < this->mesh.getDim(); ++d) {
    double * velo_n_p = velo_n.storage(d);
    double * velo_pc_p = velo_pc.storage(d);
    double * target_p = target.storage(d);

    for (int n = 0; n < target.getNbNodes(); ++n) {
      target_p[n] = 0.5 * (velo_n_p[n] + velo_pc_p[n]);
    }
  }
}

/* -------------------------------------------------------------------------- */
// velocity = cs / mu       * residual (for in-plane shear components)
// velocity = cs / mu / eta * residual (for normal component)
// velocity = cs / mu       * residual (for out-of-plane shear components)
void HalfSpace::computeVelocity(NodalField & _velo) {
  double mu = this->material->getShearModulus();
  double Cs = this->material->getCs();
  double Cp = this->material->getCp();
  std::vector<double> eta = {1.0, Cp / Cs, 1.0};

  for (int d=0; d < this->mesh.getDim(); ++d) {
    double * velo_p = _velo.storage(d);
    double * res_p = this->residual.storage(d);
    double eta_d = eta[d];

    for (int n=0; n<_velo.getNbNodes(); ++n)
      velo_p[n] = Cs / mu / eta_d * res_p[n];
  }
}

/* -------------------------------------------------------------------------- */
bool HalfSpace::registerDumpFieldToDumper(const std::string & field_name,
					  const std::string & dump_name,
					  Dumper * const dumper) {

  int d = std::atoi(&field_name[field_name.length() - 1]);

  if (d >= this->mesh.getDim())
    throw std::runtime_error("Field "+field_name
			     +" cannot be dumped, too high dimension");
  
  // disp
  if (field_name == "disp_" + std::to_string(d)) {
    dumper->registerIO(dump_name, this->disp.component(d));
    return true;
  }
  // velo
  else if (field_name == "velo_" + std::to_string(d)) {
    dumper->registerIO(dump_name, this->velo.component(d));
    return true;
  }
  // residual
  else if (field_name == "residual_" + std::to_string(d)) {
    dumper->registerIO(dump_name, this->residual.component(d));
    return true;
  }
  // internal
  else if (field_name == "internal_" + std::to_string(d)) {
    dumper->registerIO(dump_name, this->internal.component(d));
    return true;
  }

  return false;
}

/* -------------------------------------------------------------------------- */
void HalfSpace::registerToRestart(Restart & restart) {

  this->disp.registerToRestart(restart);
  this->velo.registerToRestart(restart);
  this->internal.registerToRestart(restart);
  this->residual.registerToRestart(restart);

  if (this->predictor_corrector) {
    this->disp_pc.registerToRestart(restart);
    this->velo_pc.registerToRestart(restart);
  }

}

__END_UGUCA__
