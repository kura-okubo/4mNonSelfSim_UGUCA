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
  // create default section
  this->sections.insert(std::pair<std::string,std::shared_ptr<InputSection>>(ParameterReader::general,
									     std::make_shared<InputSection>()));
}

/* -------------------------------------------------------------------------- */
void ParameterReader::readInputFile(std::string file_name) {

  char comment_char = '#';
  char equal_char = '=';
  char section_start_char = '[';
  char section_end_char = ']';

  // start with general section
  std::string current_section = ParameterReader::general;
  
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

    // this is a line with a parameter
    if (clean_line.find(equal_char) != std::string::npos) {
    
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
      this->getSection(current_section).insert(keyword, value);
    }
    // this is the start of a section
    else if (clean_line.find(section_start_char) != std::string::npos) {
      std::stringstream sstr(clean_line);
      
      std::string sec_type;
      std::string sec_name;
      std::string sec_open;

      sstr >> sec_type;

      size_t open_p = sec_type.find_first_of(section_start_char);
      if (open_p != std::string::npos) { // one word without space to [
	sec_open = sec_type[open_p];
	sec_name = sec_type.substr(0,open_p);
	sec_type = "section";
      }
      else { // two words
	sstr >> sec_name;
	open_p = sec_name.find_first_of(section_start_char);
	if (open_p != std::string::npos) { // no space to [
	  sec_open = sec_name[open_p];
	  sec_name = sec_name.substr(0,open_p);
	}
	else {
	  sstr >> sec_open;
	}
      }

      // check syntax correct
      if ((sec_open.length() != 1) || sec_open[0] != section_start_char) {
	std::cerr << " *** WARNING *** Unrespected convention! Ignore line: ";
	std::cerr << clean_line << std::endl;
	continue;
      }

      // use name to store section
      current_section = sec_name;

      // create new section
      if (this->sections.find(current_section) == this->sections.end()) {
	this->sections.insert(std::pair<std::string,std::shared_ptr<InputSection>>(current_section,
										   std::make_shared<InputSection>(sec_type)));
      }
    }
    // this is the end of a section
    else if (clean_line.find(section_end_char) != std::string::npos) {

      if (clean_line.length() != 1) {
	std::cerr << " *** WARNING *** Unrespected convention! Ignore line: ";
	std::cerr << clean_line << std::endl;
      }
      
      // everything outside a section goes into the general section
      current_section = ParameterReader::general;
    }
    // don't know what this line is
    else {
      std::cerr << " *** WARNING *** Don't know what to do with this line. Ignore it: ";
      std::cerr << clean_line << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
void ParameterReader::writeInputFile(std::string file_name) const {

  // open file to write input information
  std::ofstream outfile;
  outfile.open(file_name.c_str());

  for (auto const& sec : this->sections) {

    outfile << sec.second->getType() << " " << sec.first << " [" << std::endl;

    const std::map<std::string,std::string> & data = sec.second->getData();
    for (auto const& e : data) {
      outfile << "  " << e.first << " = " << e.second << std::endl;
    }

    outfile << "]" << std::endl << std::endl;
  }
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
template<>
std::string ParameterReader::getOrUse<std::string>(std::string key,
						   std::string alt_value,
						   std::string section) const {
  return this->getSection(section).get<std::string>(key, alt_value);
}

/* -------------------------------------------------------------------------- */
template<>
double ParameterReader::getOrUse<double>(std::string key,
					 double alt_value,
					 std::string section) const {
  return this->getSection(section).get<double>(key, alt_value);
}

/* -------------------------------------------------------------------------- */
template<>
int ParameterReader::getOrUse<int>(std::string key,
				   int alt_value,
				   std::string section) const {
  return this->getSection(section).get<int>(key, alt_value);
}

/* -------------------------------------------------------------------------- */
template<>
unsigned int ParameterReader::getOrUse<unsigned int>(std::string key,
						     unsigned int alt_value,
						     std::string section) const {
  return this->getSection(section).get<unsigned int>(key, alt_value);
}

/* -------------------------------------------------------------------------- */
template<>
bool ParameterReader::getOrUse<bool>(std::string key,
				     bool alt_value,
				     std::string section) const {
  return this->getSection(section).get<bool>(key, alt_value);
}


/* -------------------------------------------------------------------------- */
bool ParameterReader::has(std::string key, std::string section) const {
  return this->getSection(section).has(key);
}

__END_UGUCA__
