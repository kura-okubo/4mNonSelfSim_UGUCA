/**
 * @file   uca_input_section.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 *
 * @date creation: Sun Jun 26 2022
 * @date last modification: Sun Jun 26 2022
 *
 * @brief  contains information from a section in the input file
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

#include "uca_input_section.hh"

// std
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <algorithm>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
void InputSection::insert(std::string key, std::string value) {
  if (this->has(key)) {
    throw std::runtime_error("Section contains already key = "+key);
  }
  this->data.insert(std::pair<std::string,std::string>(key, value));
}

/* -------------------------------------------------------------------------- */
std::string InputSection::get(std::string key) const {
  std::map<std::string,std::string>::const_iterator it;
  it = this->data.find(key);

  // if not in map
  if (it == this->data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. "
	      << "You need the following line in your input file: ";
    std::cerr << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }

  else
    return it->second;
}

/* -------------------------------------------------------------------------- */
// not stricly needed, but so that user always uses template (avoid confusion)
template<>
std::string InputSection::get<std::string>(std::string key) const {
  return this->get(key);
}

/* -------------------------------------------------------------------------- */
template<>
double InputSection::get<double>(std::string key) const {
	return std::stod(this->get(key));
}

/* -------------------------------------------------------------------------- */
template<>
int InputSection::get<int>(std::string key) const {
  std::string value = this->get(key);
  if (std::abs((std::stoi(value)-std::stod(value))/std::stod(value)) > 1e-10) {
    std::cerr << "WARNING: value " << value
	      << " is returned as 'int' with " << std::stoi(value) << std::endl;
  }
  return std::stoi(value);
}

/* -------------------------------------------------------------------------- */
template<>
unsigned int InputSection::get<unsigned int>(std::string key) const {
  int value = this->get<int>(key);
  if (value < 0) {
    throw std::runtime_error("negative value turned into unsigned int"+value);
  }
  return (unsigned int)(value);
}

/* -------------------------------------------------------------------------- */
template<>
bool InputSection::get<bool>(std::string key) const {
  std::string value = this->get(key);
  std::transform(value.begin(),value.end(),value.begin(),::tolower);
  bool b;
  if (value.compare("true") == 0)
    b = true;
  else if (value.compare("false") == 0)
    b = false;
  else {
    throw std::runtime_error("boolean cannot be "+value);
  }
  return b;
}

/* -------------------------------------------------------------------------- */
bool InputSection::has(std::string key) const {
  std::map<std::string,std::string>::const_iterator it;
  it = this->data.find(key);
  return (it != this->data.end());
}

__END_UGUCA__
