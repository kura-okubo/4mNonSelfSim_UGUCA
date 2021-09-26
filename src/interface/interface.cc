/**
 * @file   interface.cc
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
#include "interface.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
Interface::Interface(FFTableMesh & mesh,
		     InterfaceLaw & law,
		     const std::string & name) :
  Dumper(mesh),
  name(name),
  mesh(mesh),
  time_step(0.),
  load(mesh,name+"_load"),
  cohesion(mesh,name+"_cohesion"),
  scratch_field(mesh,name+"_scratch"),
  law(&law) {
  
  this->law->setInterface(this);
}

/* -------------------------------------------------------------------------- */
Interface::Interface(FFTableMesh & mesh,
		     const std::string & name) :
  Dumper(mesh),
  name(name),
  mesh(mesh),
  time_step(0.),
  load(mesh,name+"_load"),
  cohesion(mesh,name+"_cohesion"),
  scratch_field(mesh,name+"_scratch"),
  law(NULL) {}

/* -------------------------------------------------------------------------- */
void Interface::init(bool velocity_initial_conditions) {

  // initiate the convolutions in the half spaces
  for (unsigned int i=0;i<this->half_spaces.size();++i)
    this->half_spaces[i]->initConvolutions();

  for (int i = 0; i < this->nb_pc; ++i) {
    if (i == 0)
      this->updateVelocity();
    this->computeDisplacement(true);
    this->computeCohesion(true);
    this->computeResidual();
    if (!velocity_initial_conditions) {
      this->computeVelocity(true);
      this->correctVelocity(i == this->nb_pc - 1);
    }
  }

  this->computeInternal();
  this->computeCohesion();
  this->computeResidual();
  if (!velocity_initial_conditions)
    this->computeVelocity();
}

/* -------------------------------------------------------------------------- */
void Interface::initPredictorCorrector(int iterations) {
  this->nb_pc = iterations;

  if (iterations > 0) {
    for (unsigned int i=0;i<this->half_spaces.size();++i) {
      this->half_spaces[i]->initPredictorCorrector();
    }
  }
}

/* -------------------------------------------------------------------------- */
void Interface::setTimeStep(double time_step) {
  this->time_step = time_step;
  for (unsigned int i=0;i<this->half_spaces.size();++i)
    this->half_spaces[i]->setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
double Interface::getStableTimeStep() {

  if (this->half_spaces.size()==2)
    return std::min(this->half_spaces[0]->getStableTimeStep(),
		    this->half_spaces[1]->getStableTimeStep());
  else
    return this->half_spaces[0]->getStableTimeStep();
}

/* -------------------------------------------------------------------------- */
/*void Interface::setDynamic(bool fully_dynamic) {
  for (unsigned int i = 0; i < this->half_space.size(); ++i)
    this->half_space[i]->setDynamic(fully_dynamic);
    }*/

/* -------------------------------------------------------------------------- */
void Interface::computeDisplacement(bool predicting) {
  for (unsigned int i=0;i<this->half_spaces.size();++i)
    this->half_spaces[i]->computeDisplacement(predicting);
}

/* -------------------------------------------------------------------------- */
void Interface::computeInternal(bool predicting, bool correcting) {
  for (unsigned int i=0;i<this->half_spaces.size();++i)
    this->half_spaces[i]->computeInternal(predicting, correcting);
}

/* -------------------------------------------------------------------------- */
void Interface::computeResidual() {
  this->combineLoadAndCohesion(this->scratch_field);
  for (unsigned int i=0;i<this->half_spaces.size();++i)
    this->half_spaces[i]->computeResidual(this->scratch_field);
}
/* -------------------------------------------------------------------------- */
void Interface::computeVelocity(bool predicting) {
  for (unsigned int i=0;i<this->half_spaces.size();++i)
    this->half_spaces[i]->computeVelocity(predicting);
}

/* -------------------------------------------------------------------------- */
void Interface::updateVelocity() {
  for (unsigned int i=0;i<this->half_spaces.size();++i)
    this->half_spaces[i]->updateVelocity();
}

/* -------------------------------------------------------------------------- */
void Interface::correctVelocity(bool last_step) {
  for (unsigned int i=0;i<this->half_spaces.size();++i)
    this->half_spaces[i]->correctVelocity(last_step);
}

/* -------------------------------------------------------------------------- */
void Interface::advanceTimeStep() {

  // predictor-corrector
  for (int i = 0; i < this->nb_pc; ++i) {
    if (i == 0) {
      // copy v to scratch memory
      this->updateVelocity();
    }
    // Predict
    // u* = u + v * dt
    this->computeDisplacement(true);
    
    // f* -> compute cohesion -> tau_coh*
    this->computeCohesion(true);
    // tau_coh* -> compute residual -> tau_res*
    this->computeResidual();
    // tau_res* -> compute velocity -> v*
    this->computeVelocity(true);
    
    // Correct
    // v** = (v + v*) / 2 ---> overwrite storage if reached last step
    this->correctVelocity(i == this->nb_pc - 1);
  }

  // compute displacement
  this->computeDisplacement();
  this->computeInternal(false,false);
  this->computeCohesion();
  this->computeResidual();
  this->computeVelocity();
}

/* -------------------------------------------------------------------------- */
void Interface::computeCohesion(bool predicting) {
  this->law->computeCohesiveForces(this->cohesion, predicting);
}

/* -------------------------------------------------------------------------- */
void Interface::combineLoadAndCohesion(NodalField & load_and_cohesion) {
  
  for (int d=0; d<this->mesh.getDim(); ++d) {

    // tau_0 - tau_coh
    double * load_and_cohesion_p = load_and_cohesion.storage(d);
    double * coh_p = this->cohesion.storage(d);
    double * load_p = this->load.storage(d);

    for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {
      load_and_cohesion_p[n] = load_p[n] - coh_p[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void Interface::registerDumpField(const std::string & field_name) {

  int d = std::atoi(&field_name[field_name.length()-1]);

  if (d >= this->mesh.getDim())
    throw std::runtime_error("Field "+field_name
			     +" cannot be dumped, too high dimension");

  if (field_name == "load_"+std::to_string(d))
    this->registerIO(field_name,
			  this->load.component(d));

  else if (field_name == "cohesion_"+std::to_string(d))
    this->registerIO(field_name,
			  this->cohesion.component(d));

  else
    this->law->registerDumpField(field_name);
}

/* -------------------------------------------------------------------------- */
void Interface::registerToRestart(Restart & restart) {

  this->cohesion.registerToRestart(restart);
  this->load.registerToRestart(restart);
  this->scratch_field.registerToRestart(restart);

  for (unsigned int i=0;i<this->half_spaces.size();++i)
    this->half_spaces[i]->registerToRestart(restart);
  
  this->law->registerToRestart(restart);
}

__END_UGUCA__
