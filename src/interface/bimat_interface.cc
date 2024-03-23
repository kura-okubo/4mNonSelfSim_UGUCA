/**
 * @file   bimat_interface.cc
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
#include "bimat_interface.hh"
#include "half_space.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

BimatInterface::BimatInterface(FFTableMesh & mesh,
			       SpatialDirectionSet components,
			       Material & top_material,
			       Material & bot_material,
			       InterfaceLaw & law,
			       const SolverMethod & method) :
  Interface(mesh, components, law)
{
  this->top = HalfSpace::newHalfSpace(top_material, mesh,  1, components, this->name+"_top", method);
  this->bot = HalfSpace::newHalfSpace(bot_material, mesh, -1, components, this->name+"_bot", method);
  
  this->half_spaces.resize(2);
  this->half_spaces[0] = this->top;
  this->half_spaces[1] = this->bot;
}

/* -------------------------------------------------------------------------- */
BimatInterface::~BimatInterface() {
  delete this->top;
  delete this->bot;
}

/* -------------------------------------------------------------------------- */
void BimatInterface::closingNormalGapForce(NodalField & close_force,
                                           bool predicting,
					   unsigned int ts_factor) {

  // check if the normal direction exists (it doesn't for antiplane)
  if (!close_force.getComponents().count(_y)) {
    return;
  }

  // top material information
  const Material & mat_t = this->top->getMaterial();
  double cp_t = mat_t.getCp();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double eta_t = cp_t / cs_t;
  double fact_t = ts_factor * this->time_step * cs_t / mu_t / eta_t;

  // bot material information
  const Material & mat_b = this->bot->getMaterial();
  double cp_b = mat_b.getCp();
  double cs_b = mat_b.getCs();
  double mu_b = mat_b.getShearModulus();
  double eta_b = cp_b / cs_b;
  double fact_b = ts_factor * this->time_step * cs_b / mu_b / eta_b;

  // accessors
  double * u_1_t = this->top->getDisp(predicting).data(1);
  double * u_1_b = this->bot->getDisp(predicting).data(1);
  double * f_1_t = this->top->getInternal().data(1);
  double * f_1_b = this->bot->getInternal().data(1);
  double * t0_1 = this->load.data(1);
  double * cf = close_force.data(1);

  double fact_cf = fact_t + fact_b;

  fact_t /= fact_cf;
  fact_b /= fact_cf;

  for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {
    double u_1_gap = u_1_t[n] - u_1_b[n];
    double du_1_t = fact_t * (t0_1[n] + f_1_t[n]);
    double du_1_b = fact_b * (t0_1[n] + f_1_b[n]);
    cf[n] = (u_1_gap/fact_cf + du_1_t + du_1_b);
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::maintainShearGapForce(NodalField & maintain_force) {
  // top material information
  const Material & mat_t = this->top->getMaterial();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double fact_t = cs_t / mu_t;

  // bot material information
  const Material & mat_b = this->bot->getMaterial();
  double cs_b = mat_b.getCs();
  double mu_b = mat_b.getShearModulus();
  double fact_b = cs_b / mu_b;

  double fact_mf = fact_t + fact_b;

  fact_t/=fact_mf;
  fact_b/=fact_mf;

  // only in the shear and anti-shear direction
  SpatialDirectionSet comps = maintain_force.getComponents();
  comps.erase(_y);
  for (const auto& d : comps) {
    // accessors
    double * f_t = this->top->getInternal().data(d);
    double * f_b = this->bot->getInternal().data(d);
    double * t0 = this->load.data(d);
    double * mf = maintain_force.data(d);

    for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {
      //better for numerical error
      mf[n] = t0[n] + (fact_t * f_t[n] + fact_b * f_b[n]);
    }
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeGap(NodalField & gap,
                                bool predicting) {

  for (const auto& d : gap.getComponents()) {
    
    double * top_disp = this->top->getDisp(predicting).data(d);
    double * bot_disp = this->bot->getDisp(predicting).data(d);
    double * gap_p = gap.data(d);

    for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {
      gap_p[n] = top_disp[n] - bot_disp[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeGapVelocity(NodalField & gap_velo,
                                        bool predicting) {

  for (const auto& d : gap_velo.getComponents()) {

    double * top_velo = this->top->getVelo(predicting).data(d);
    double * bot_velo = this->bot->getVelo(predicting).data(d);
    double * gap_velo_p = gap_velo.data(d);

    for (int n = 0; n < this->mesh.getNbLocalNodes(); ++n) {
      gap_velo_p[n] = top_velo[n] - bot_velo[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::registerDumpField(const std::string &field_name) {

  bool registered = false;
  // field_name starts with "top"
  if (field_name.rfind("top", 0) == 0) {
    // cut away "top_" from field_name and give interface as dumper
    registered = this->top->registerDumpFieldToDumper(field_name.substr(4),
						     field_name,
						     this);
  }
  // field_name starts with "bot"
  else if (field_name.rfind("bot", 0) == 0) {
    // cut away "bot_" from field_name and give interface as dumper
    registered = this->bot->registerDumpFieldToDumper(field_name.substr(4),
						      field_name,
						      this);
  }

  if (!registered)
    Interface::registerDumpField(field_name);
}

__END_UGUCA__
