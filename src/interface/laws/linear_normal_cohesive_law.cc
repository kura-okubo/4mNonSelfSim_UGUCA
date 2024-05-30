/**
 * @file   linear_normal_cohesive_law.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Fri Jun 3 2022
 * @date last modification: Fri Jun 3 2022
 *
 * @brief  linear cohesive law for fracture
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
#include "linear_normal_cohesive_law.hh"
#include "interface.hh"

#include <cmath>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

LinearNormalCohesiveLaw::LinearNormalCohesiveLaw(BaseMesh &mesh,
                                                 double Gc_default,
                                                 double sigma_c_default,
                                                 const std::string &name) : InterfaceLaw(mesh, name),
                                                                            G_c(mesh),
                                                                            sigma_c(mesh)
{
  this->G_c.setAllValuesTo(Gc_default);
  this->G_c.setName(name + "_G_c");

  this->sigma_c.setAllValuesTo(sigma_c_default);
  this->sigma_c.setName(name + "_sigma_c");
}

/* -------------------------------------------------------------------------- */
void LinearNormalCohesiveLaw::computeCohesiveForces(bool predicting, unsigned int ts_factor)
{
  NodalField &cohesion = this->interface->getCohesion();
  // find forces needed to close normal gap
  this->interface->closingNormalGapForce(cohesion, predicting, ts_factor);

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion);

  // get norm of normal cohesion
  NodalField normal_trac_norm(this->mesh);
  cohesion.computeNorm(normal_trac_norm, 0);

  // find current gap
  NodalField gap(this->mesh, cohesion.getComponents());
  this->interface->computeGap(gap, predicting);

  // compute norm of normal gap
  NodalField normal_gap_norm(this->mesh);
  gap.computeNorm(normal_gap_norm, 0);

  // to be filled
  NodalField alpha(this->mesh);

  // coh1 > 0 is a adhesive force
  // coh1 < 0 is a contact pressure
  for (int n = 0; n < this->mesh.getNbLocalNodes(); ++n)
  {

    if ((G_c(0) < 1e-12) || (sigma_c(n) < 1e-12))
    {
      alpha(n) = 0.;
    }
    else
    {

      double slope = pow(sigma_c(n), 2) / 2. / G_c(n);
      double strength = std::max(sigma_c(n) - normal_gap_norm(n) * slope, 0.);

      // maximal normal cohesive force given by strength.
      // keep orientation of normal force
      alpha(n) = std::min(1., std::abs(strength / normal_trac_norm(n)));
    }
  }

  // only in normal direction
  cohesion.multiplyByScalarField(alpha, 1);
}

/* -------------------------------------------------------------------------- */
void LinearNormalCohesiveLaw::registerDumpField(const std::string &field_name)
{

  // G_c
  if (field_name == "G_c")
  {
    this->interface->registerIO(field_name,
                                this->G_c);
  }

  // sigma_c
  else if (field_name == "sigma_c")
  {
    this->interface->registerIO(field_name,
                                this->sigma_c);
  }

  // do not know this field
  else
  {
    InterfaceLaw::registerDumpField(field_name);
  }
}

/* -------------------------------------------------------------------------- */
void LinearNormalCohesiveLaw::registerToRestart(Restart &restart)
{

  restart.registerIO(this->G_c);
  restart.registerIO(this->sigma_c);

  InterfaceLaw::registerToRestart(restart);
}

__END_UGUCA__
