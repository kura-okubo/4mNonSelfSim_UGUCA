/**
 * @file   interface.hh
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
#ifndef __INTERFACE_H__
#define __INTERFACE_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "interface_law.hh"
#include "uca_dumper.hh"
#include "half_space.hh"
#include "uca_fftable_mesh.hh"
#include "uca_restart.hh"

__BEGIN_UGUCA__

class Interface : public Dumper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Interface(FFTableMesh & mesh,
	    SpatialDirectionSet components,
	    InterfaceLaw & law,
	    const std::string & name = "interface");

  virtual ~Interface() {}

protected:
  // for inheritate object: infinite boundary
  Interface(FFTableMesh & mesh,
	    SpatialDirectionSet components,
	    const std::string & name = "interface");
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // get (limiting) stable time step
  virtual double getStableTimeStep();

  // initiate interface model and half spaces
  virtual void init(bool velocity_initial_conditions = false);

  // initiate predictor-corrector stratch memory
  virtual void initPredictorCorrector(int iterations = 1);

  // iteration of advancing one time step
  virtual void advanceTimeStep(bool dynamic=true);

  // functions used during time stepping for each half-space
  virtual void computeDisplacement(bool predicting = false);
  virtual void computeInternal(bool predicting = false,
			       bool correcting = false,
			       bool dynamic = true);
  virtual void computeCohesion(bool predicting = false);
  virtual void computeResidual();
  virtual void computeVelocity(bool predicting = false);

  // functions used during predictor-corrector time stepping
  virtual void updateVelocity();
  virtual void correctVelocity(bool last_step);

  // function that combine load and cohesionw with correct signs
  void combineLoadAndCohesion(NodalField & load_and_cohesion);

  // compute force needed to close normal gap
  virtual void closingNormalGapForce(NodalField & close_force,
                                     bool predicting = false) = 0;

  // compute force needed to maintain current shear gap
  virtual void maintainShearGapForce(NodalField & maintain_force) = 0;

  // compute gap in displacement
  virtual void computeGap(NodalField & gap,
                          bool predicting = false) = 0;

  // compute gap relative velocity
  virtual void computeGapVelocity(NodalField & gap_velo,
                                  bool predicting = false) = 0;

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  // restart
  virtual void registerToRestart(Restart & restart);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  NodalField & getLoad()        { return this->load; }
  NodalField & getCohesion()    { return this->cohesion; }

  virtual HalfSpace & getTop() = 0;
  virtual HalfSpace & getBot() = 0;

  virtual void setTimeStep(double time_step);
  double getTimeStep() const { return this->time_step; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // name of interface (in case there are more)
  std::string name;
  
  // mesh
  FFTableMesh & mesh;

  // time step
  double time_step;

  // {top, bot}
  std::vector<HalfSpace *> half_spaces;

  // external loading
  NodalField load;

  // interface forces (e.g., cohesion)
  NodalField cohesion;

  // used to add up stuff
  NodalField scratch_field;

  // predictor-corrector
  int nb_pc = 0;

  // interface law: cohesive law and contact/friction law
  InterfaceLaw * law;
};

__END_UGUCA__

//#include "interface_impl.cc"

#endif /* __INTERFACE_H__ */
