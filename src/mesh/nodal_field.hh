/**
 * @file   nodal_field.hh
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
/* -------------------------------------------------------------------------- */
#ifndef __NODAL_FIELD_H__
#define __NODAL_FIELD_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "nodal_field_component.hh"
#include <vector>
/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

class NodalField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NodalField(const std::string & name = "unnamed") : name(name),
						     initialized(false),
						     mesh(NULL) {}
  NodalField(BaseMesh & mesh,
	     const std::string & name = "unnamed") : name(name),
						     initialized(false)
  { this->init(mesh); }
  
  virtual ~NodalField() { this->free(); }

private:
  // private copy constructor: NodalField cannot be copied (for now to debug)
  NodalField(NodalField & to_copy);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void init(BaseMesh & mesh);
  virtual void free();

public:
  void zeros();
  void setAllValuesTo(double value);

  // compute element-wise norm of field
  void computeNorm(NodalFieldComponent & norm,
		   int ignore_dir = -1) const;

  // multiply fields element-wise with scalar
  void multiplyByScalar(const NodalFieldComponent & scalar,
			int ignore_dir = -1);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get name of nodal field
  std::string getName() const { return this->name; }
  
  // get dimension
  int getDim() const { return this->mesh->getDim(); }
  
  // get number of nodes
  int getNbNodes() const { return this->mesh->getNbLocalNodes(); };

  inline NodalFieldComponent & component(int i) { return (*this->field[i]); }
  
  // access to storage
  inline double * storage(int i) { return this->field[i]->storage(); }
  inline const double * storage(int i) const { return this->field[i]->storage(); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // name
  std::string name;
  
  // initialized
  bool initialized;
  
  // associated mesh
  BaseMesh * mesh;

  // nodal field
  std::vector<NodalFieldComponent * > field;
};

__END_UGUCA__

#endif /* __NODAL_FIELD_H__ */
