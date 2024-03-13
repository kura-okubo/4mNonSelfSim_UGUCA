/**
 * @file   half_space.hh
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
#ifndef __HALF_SPACE_H__
#define __HALF_SPACE_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "material.hh"
#include "nodal_field.hh"
#include "fftable_nodal_field.hh"
#include "hist_fftable_nodal_field.hh"
#include "uca_fftable_mesh.hh"
#include "uca_dumper.hh"
#include "uca_restart.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class HalfSpace {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  // side factor top=1 bot=-1
  HalfSpace(Material & material, FFTableMesh & mesh, int side_factor,
	    SpatialDirectionSet components,
	    const std::string & name = "half_space");

  virtual ~HalfSpace();

  // allocate HalfSpace of type
  static HalfSpace * newHalfSpace(Material & material,
				  FFTableMesh & mesh,
				  int side_factor,
				  SpatialDirectionSet components,
				  const std::string & name,
				  const SolverMethod & method);
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void computeDisplacement(bool predicting = false);
  virtual void computeInternal(bool predicting = false,
			       bool correcting = false,
			       bool dynamic = true);
  virtual void computeResidual(NodalField & external);
  virtual void computeVelocity(bool predicting = false);

  // init convolutions
  virtual void initConvolutions() {}

  // for predictor-corrector implmentation
  void initPredictorCorrector();
  virtual void updateVelocity();
  virtual void correctVelocity(bool last_step);

  // dumper function: returns true if successful
  bool registerDumpFieldToDumper(const std::string & field_name, // what to dump
				 const std::string & dump_name, // how to name it
				 Dumper * const dumper);

  // restart
  virtual void registerToRestart(Restart & restart);

  // set half space to steady state (used for transition from
  // quasi-dynamic to dynamic)
  virtual void setSteadyState(bool /*predicting*/ = false) {};
  
protected:
  // apply fft forward on displacement
  virtual void forwardFFT(bool predicting = false);
  // apply fft backward on internal
  virtual void backwardFFT();

  virtual void computeStressFourierCoeff(bool predicting = false,
					 bool correcting = false,
					 bool dynamic = true) = 0;

private:
  void computeDisplacement(NodalField & disp,
			   NodalField & velo,
			   NodalField & target);
  void computeVelocity(NodalField & _velo);
  void correctVelocity(NodalField & velo_n,
		       NodalField & velo_pc,
		       NodalField & target);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get material of half space
  const Material& getMaterial() const { return this->material; }

  // set time step
  virtual void setTimeStep(double time_step) { this->time_step = time_step; }

  // get side factor
  int getSideFactor() const { return this->side_factor; }
  
  // accessors
  FFTableNodalField & getDisp(bool predicting = false) {
    return predicting ? this->disp_pc : this->disp;
  }
  NodalField & getVelo(bool predicting = false) {
    return predicting ? this->velo_pc : this->velo;
  }

  FFTableNodalField & getInternal() { return this->internal; };

  NodalField & getResidual() { return this->residual; };

  // const accessors
  const FFTableNodalField & getDisp(bool predicting = false) const {
    return predicting ? this->disp_pc : this->disp;
  }
  const NodalField & getVelo(bool predicting = false) const {
    return predicting ? this->velo_pc : this->velo;
  }
  const FFTableNodalField & getInternal() const { return this->internal; }

  virtual double getStableTimeStep();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // name
  std::string name;

  // material properties
  Material & material;

  // mesh
  FFTableMesh & mesh;

  // time step
  double time_step;

  // used to know in which directions the tractions pull
  int side_factor;

  // displacement 0 x "in-plane shear" ; 1 y "normal"; 2 z "out-of-plane shear"
  // with limited history
  HistFFTableNodalField disp;

  // velocity
  NodalField velo;
  
  // tractions due to deformation
  FFTableNodalField internal;

  // all acting forces
  NodalField residual;

  // for predictor-corrector implmentation
  bool predictor_corrector = false;
  FFTableNodalField disp_pc;
  NodalField velo_pc;

};

__END_UGUCA__

//#include "half_space_impl.cc"

#endif /* __HALF_SPACE_H__ */
