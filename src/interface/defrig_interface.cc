/**
 * @file   defrig_interface.cc
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
#include "defrig_interface.hh"
#include "half_space.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
DefRigInterface::DefRigInterface(FFTableMesh & mesh,
				 SpatialDirectionSet components,
				 Material & top_material,
				 InterfaceLaw & law,
				 const SolverMethod & method) :
  Interface(mesh, components, law)
{
  this->top = HalfSpace::newHalfSpace(top_material, mesh, 1, components, this->name+"_top", method);
  
  this->half_spaces.resize(1);
  this->half_spaces[0] = this->top;

  // top material information
  const Material & mat_t = this->top->getMaterial();
  double cp_t = mat_t.getCp();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double eta_t = cp_t / cs_t;
  this->fact_t_2 = cs_t / mu_t / eta_t;
}

/* -------------------------------------------------------------------------- */
DefRigInterface::~DefRigInterface() {
  delete this->top;
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::closingNormalGapForce(NodalField & close_force,
					    bool predicting,
					    unsigned int ts_factor) {
  
  // C factor of notes
  double fact_t = ts_factor * this->time_step * this->fact_t_2;

  // accessors
  double * f_1_t = this->top->getInternal().data(1);
  double * t0_1 = this->load.data(1);
  double * cf = close_force.data(1);

  NodalField * gap = (&this->scratch_field);
  this->computeGap(*gap, predicting);
  double * gap_1_p = gap->data(1);

  for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {
    double u_1_gap = gap_1_p[n] / fact_t;
    double du_1_t = t0_1[n] + f_1_t[n];
    cf[n] = u_1_gap + du_1_t;
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::maintainShearGapForce(NodalField & maintain_force) {

  // only in the shear and anti-shear direction
  SpatialDirectionSet comps = maintain_force.getComponents();
  comps.erase(_y);
  for (const auto& d : comps) {

    // accessors
    double * f_t = this->top->getInternal().data(d);
    double * t0 = this->load.data(d);
    double * mf = maintain_force.data(d);

    for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {
      mf[n] = t0[n] + f_t[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeGap(NodalField & gap,
                                 bool predicting) {

  for (const auto& d : gap.getComponents()) {

    double * top_disp = this->top->getDisp(predicting).data(d);
    double * gap_p = gap.data(d);

    for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {
      gap_p[n] = top_disp[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeGapVelocity(NodalField & gap_velo,
                                         bool predicting) {

  for (const auto& d : gap_velo.getComponents()) {
    
    double *top_velo = this->top->getVelo(predicting).data(d);
    double *gap_velo_p = gap_velo.data(d);

    for (int n = 0; n < this->mesh.getNbLocalNodes(); ++n) {
      gap_velo_p[n] = top_velo[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::registerDumpField(const std::string & field_name) {

  bool registered = false;
  // field_name starts with "top"
  if (field_name.rfind("top", 0) == 0) {
    // cut away "top_" from field_name and give interface as dumper
    registered = this->top->registerDumpFieldToDumper(field_name.substr(4),
						      field_name,
						      this);
  }

  if (!registered)
    Interface::registerDumpField(field_name);
}

__END_UGUCA__
