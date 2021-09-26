/**
 * @file   nodal_field_component.hh
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
#ifndef __NODAL_FIELD_COMPONENT_H__
#define __NODAL_FIELD_COMPONENT_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_base_mesh.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

class Restart;

class NodalFieldComponent {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NodalFieldComponent(const std::string & name = "unnamed") :
    name(name),
    initialized(false),
    mesh(NULL),
    direction(0),
    field(NULL) {}

  NodalFieldComponent(BaseMesh & mesh,
		      int direction=0,
		      const std::string & name = "unnamed") :
    name(name),
    initialized(false),
    mesh(NULL),
    direction(direction),
    field(NULL)
  { this->init(mesh); }

  NodalFieldComponent(BaseMesh & mesh,
		      const std::string & name,
		      int direction=0) :
    name(name),
    initialized(false),
    mesh(NULL),
    direction(direction),
    field(NULL)
  { this->init(mesh); }

  
  virtual ~NodalFieldComponent();

private:
  // private copy constructor: cannot be copied (for compilation error)
  NodalFieldComponent(NodalFieldComponent & to_copy);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void init(BaseMesh & mesh);

  void zeros();
  void setAllValuesTo(double value);

  // restart
  virtual void registerToRestart(Restart & restart);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get number of nodes
  int getNbNodes() const { return this->mesh->getNbLocalNodes(); }

  // get direction
  inline int getDirection() const { return this->direction; }
  
  // access the value of node n (reading and writing)
  inline double & operator()(int node);

  // access the value of node n (only reading for const nodal field)
  inline double & set(int node);

  // access the value of node n (only reading for const nodal field)
  inline double at(int node) const ;

  // access to storage
  inline double * storage() { return this->field; }
  inline const double * storage() const { return this->field; }

  // get name
  const std::string & getName() const { return this->name; }
  
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

  // direction (axis)
  int direction;
  
  // nodal field
  double * field;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline double & NodalFieldComponent::operator()(int node) {
  return this->field[node];
}

inline double & NodalFieldComponent::set(int node) {
  return this->field[node];
}

inline double NodalFieldComponent::at(int node) const {
  return this->field[node];
}

__END_UGUCA__

#endif /* __NODAL_FIELD_COMPONENT_H__ */
