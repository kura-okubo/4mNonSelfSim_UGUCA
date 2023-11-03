/**
 * @file   uca_input_section.hh
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Sun Jun 26 2022
 * @date last modification: Sun Jun 26 2022
 *
 * @brief  Contains information from a section in the input file
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
#ifndef __INPUT_SECTION_HH__
#define __INPUT_SECTION_HH__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"

// std
#include <map>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class InputSection {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  InputSection(std::string type = "section") : type(type) {};
  virtual ~InputSection() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

private:
  /// get value returns string and checks
  std::string get(std::string) const;
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// insert a new pair
  void insert(std::string key, std::string value);

  /// returns value for key at type T
  template<typename T>
  T get(std::string key) const;

  /// returns value for key at type T or a provided alternative value
  template<typename T>
  T get(std::string key, T alternative_value) const;
  
  /// check if key is in data
  bool has(std::string key) const;

  const std::string getType() const { return this->type; }
  
  /// access to data
  const std::map<std::string,std::string> & getData() const { return this->data; }
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// information on type of input section
  std::string type;
  
  /// data stored as strings
  std::map<std::string,std::string> data;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

__END_UGUCA__

//#include "input_section_inline_impl.cc"

#endif /* __INPUT_SECTION_HH__ */
