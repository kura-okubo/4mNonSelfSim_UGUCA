/**
 * @file   uca_parameter_reader.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Sun Jun 26 2022
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

#include "uca_parameter_reader.hh"

// std
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <algorithm>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
ParameterReader::ParameterReader() {
  this->sections.insert(std::pair<std::string,std::shared_ptr<InputSection>>(ParameterReader::general,
									     std::make_shared<InputSection>()));
}

/* -------------------------------------------------------------------------- */
void ParameterReader::readInputFile(std::string file_name) {

  char comment_char = '#';
  char equal_char = '=';

  // open a file called file name
  std::ifstream infile;
  infile.open(file_name.c_str());

  if(!infile.good()) {
    std::cerr << "Cannot open file " << file_name << "!!!" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::string clean_line;
  while(infile.good()) {
    getline(infile, line);
    clean_line = line;

    // take out comments
    size_t found_comment;
    found_comment = line.find_first_of(comment_char);
    if (found_comment != std::string::npos)
      clean_line = line.substr(0,found_comment);
    if (clean_line.empty())
      continue;

    std::stringstream sstr(clean_line);

    std::string keyword;
    std::string equal;
    std::string value;

    // get keyword
    sstr >> keyword;
    size_t equal_p = keyword.find_first_of(equal_char);
    if (equal_p != std::string::npos) {
      equal = keyword.substr(equal_p,std::string::npos);
      keyword = keyword.substr(0,equal_p);
    }

    // get equal
    if (equal.empty())
      sstr >> equal;
    if (equal.length() != 1) {
      value = equal.substr(1,std::string::npos);
      equal = equal[0];
    }
    if (equal[0] != equal_char) {
      std::cerr << " *** WARNING *** Unrespected convention! Ignore line: ";
      std::cerr << clean_line << std::endl;
      continue;
    }

    // get value
    if (value.empty())
      sstr >> value;

    // no value
    if (value.empty()) {
      std::cerr << " *** WARNING *** No value given! Ignore line: ";
      std::cerr << clean_line << std::endl;
      continue;
    }

    // put value in map
    this->getSection(ParameterReader::general).insert(keyword, value);
  }
}

/* -------------------------------------------------------------------------- */
void ParameterReader::writeInputFile(std::string file_name) const {

  // open file to write input information
  std::ofstream outfile;
  outfile.open(file_name.c_str());

  InputSection & is = this->getSection(ParameterReader::general);
  const std::map<std::string,std::string> & data = is.getData();

  for (std::map<std::string, std::string>::const_iterator it = data.begin();
       it != data.end(); ++it)
    outfile << it->first << " = " << it->second << std::endl;
}

/* -------------------------------------------------------------------------- */
InputSection & ParameterReader::getSection(std::string name) {

  SectionMap::iterator it;
  it = this->sections.find(name);

  // if not in map
  if (it == this->sections.end()) {
    throw std::runtime_error("Don't know input section: "+name);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
InputSection & ParameterReader::getSection(std::string name) const {

  SectionMap::const_iterator it;
  it = this->sections.find(name);

  // if not in map
  if (it == this->sections.end()) {
    throw std::runtime_error("Don't know input section: "+name);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template<>
std::string ParameterReader::get<std::string>(std::string key,
					      std::string section) const {
  return this->getSection(section).get<std::string>(key);
}

/* -------------------------------------------------------------------------- */
template<>
double ParameterReader::get<double>(std::string key,
				    std::string section) const {
  return this->getSection(section).get<double>(key);
}

/* -------------------------------------------------------------------------- */
template<>
int ParameterReader::get<int>(std::string key,
			      std::string section) const {
  return this->getSection(section).get<int>(key);
}

/* -------------------------------------------------------------------------- */
template<>
unsigned int ParameterReader::get<unsigned int>(std::string key,
						std::string section) const {
  return this->getSection(section).get<unsigned int>(key);
}

/* -------------------------------------------------------------------------- */
template<>
bool ParameterReader::get<bool>(std::string key,
				std::string section) const {
  return this->getSection(section).get<bool>(key);
}

/* -------------------------------------------------------------------------- */
bool ParameterReader::has(std::string key, std::string section) const {
  return this->getSection(section).has(key);
}

__END_UGUCA__
