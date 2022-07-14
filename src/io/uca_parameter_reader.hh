/**
 * @file   uca_parameter_reader.hh
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Sun Jun 26 2022
 *
 * @brief  reads input file and stores input parameters
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
#ifndef __PARAMETER_READER_HH__
#define __PARAMETER_READER_HH__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_input_section.hh"

// std
#include <map>
#include <memory>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class ParameterReader {
  
private:
  inline static const std::string general = "general";

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  typedef std::map<const std::string,std::shared_ptr<InputSection>> SectionMap;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ParameterReader();
  virtual ~ParameterReader() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read input file
  void readInputFile(std::string file_name);

  /// write input file
  void writeInputFile(std::string file_name) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get a section
  InputSection & getSection(std::string name = ParameterReader::general);
  InputSection & getSection(std::string name = ParameterReader::general) const;

  ///
  template<typename T>
  T get(std::string key, std::string section = ParameterReader::general) const;

  bool has(std::string key, std::string section = ParameterReader::general) const;
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// data
  SectionMap sections;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

__END_UGUCA__

//#include "parameter_reader_inline_impl.cc"

#endif /* __PARAMETER_READER_HH__ */
