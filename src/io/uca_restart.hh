/**
 * @file   uca_restart.hh
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Sat Sept 25 2021
 * @date last modification: Sat Sept 25 2021
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
#ifndef __RESTART_H__
#define __RESTART_H__

/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_base_io.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class Restart : public BaseIO {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Restart(const std::string &bname,
	  const std::string &path,
	  const Format format = Format::ASCII);
  virtual ~Restart();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void dump(unsigned int step);
  virtual void load(unsigned int step);
  
protected:
  void setBaseName(const std::string & bname);
  std::string getFilePath(const std::string & name, unsigned int number);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
 /* ------------------------------------------------------------------------ */
 /* Class Members                                                            */
 /* ------------------------------------------------------------------------ */
protected:
};

__END_UGUCA__

//#include "uca_restart_impl.cc"

#endif /* __RESTART_H__ */
