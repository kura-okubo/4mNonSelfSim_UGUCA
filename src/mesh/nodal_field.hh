/**
 * @file   nodal_field.hh
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Fri Jan 5 2024
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
#include <stdexcept>
#include <set>
#include "uca_common.hh"
#include "uca_base_mesh.hh"
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
	     SpatialDirectionSet components = {_x},
	     const std::string & name = "unnamed");
  
  virtual ~NodalField() {}

private:
  // private copy constructor: NodalField cannot be copied (for now to debug)
  NodalField(NodalField & to_copy);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // clears the NodalField and reinitializes it
  virtual void resize(BaseMesh & mesh, SpatialDirectionSet components);

  // set all entries to zero
  void zeros();

  // if d provided, only do it on that component
  void setAllValuesTo(double value, int d = -1);

  // compute element-wise norm of field
  void computeNorm(NodalField & norm,
		   int ignore_dir = -1) const;

  // multiply component of field element-wise with scalar
  void multiplyByScalarField(const NodalField & scalar, int d);
  // multiply all component of field element-wise with scalar
  void multiplyByScalarField(const NodalField & scalar);

  // copy data
  void copyDataFrom(const NodalField & other);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get name of nodal field
  std::string getName() const { return this->name; }

  // set name of nodal field
  void setName(const std::string & name) { this->name = name; }
  
  // get number of nodes
  int getNbNodes() const { return this->mesh->getNbLocalNodes(); }

  // get number of components
  int getNbComponents() const { return this->components.size(); }

  // get the components, i.e. a set of integers
  SpatialDirectionSet getComponents() const { return this->components; }
  
  // access the value of node n (reading and writing)
  inline double & operator()(int n, int d=0);
  inline const double & operator()(int n, int d=0) const;
  
  // access to storage
  inline double * data(int i = 0);
  inline const double * data(int = 0) const;

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

  // components: integers indicating spatial directions
  SpatialDirectionSet components;

  // start indices for each component
  std::vector<int> start;
  
  // storage: x1,x2,...,xn,y1,y2,...,yn,z1,z2,...,zn
  std::vector<double> storage;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline double & NodalField::operator()(int n, int d) {
  if (!this->components.count(d)) 
    throw std::runtime_error("NodalField "+this->name+" has no component "+std::to_string(d)+"\n");
  return this->storage[this->start[d]+n];
}

inline const double & NodalField::operator()(int n, int d) const {
  if (!this->components.count(d)) 
    throw std::runtime_error("NodalField "+this->name+" has no component "+std::to_string(d)+"\n");
  return this->storage[this->start[d]+n];
}

inline double * NodalField::data(int d) {
  if (!this->components.count(d))
    throw std::runtime_error("NodalField "+this->name+" has no component "+std::to_string(d)+"\n");
  return this->storage.data() + this->start[d];
}

inline const double * NodalField::data(int d) const {
  if (!this->components.count(d))
    throw std::runtime_error("NodalField "+this->name+" has no component "+std::to_string(d)+"\n");
  return this->storage.data() + this->start[d];
}

__END_UGUCA__

#endif /* __NODAL_FIELD_H__ */
