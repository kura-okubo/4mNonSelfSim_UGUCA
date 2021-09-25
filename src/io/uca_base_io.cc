/**
 * @file   uca_base_io.cc
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
#include "uca_base_io.hh"
#include "static_communicator_mpi.hh"

/*#include "nodal_field_component.hh"
#include "static_communicator_mpi.hh"
#include "uca_custom_mesh.hh"

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <typeinfo>
*/
#include <stdexcept>
	
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
BaseIO::BaseIO() {
}

/* -------------------------------------------------------------------------- */
BaseIO::~BaseIO() {}

/* -------------------------------------------------------------------------- */
void BaseIO::initIO(const std::string &bname,
		    const std::string &path,
		    const Format format) {

  this->setBaseName(bname);
  this->path = path;
  this->dump_format = format;

  switch (this->dump_format) {
  case Format::ASCII: {
    this->separator = " ";
    this->file_extension = ".out";
    break;
  }
  case Format::CSV: {
    this->separator = ",";
    this->file_extension = ".csv";
    break;
  }
  case Format::Binary: {
    this->separator = " ";
    this->file_extension = ".out";
    break;
  }
  default:
    throw std::runtime_error("Unsupported output format.");
  }

  // only one proc creates folder
  int rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  if (rank == 0) {
    // create folder for files (works only on linux)
    // read/write/search permission for owner and group
    // read/search permissions for others
    std::string full_path_to_folder = this->path + BaseIO::directorySeparator() + this->folder_name;
    BaseIO::createDirectory(full_path_to_folder);
  }
}

/* -------------------------------------------------------------------------- */
void BaseIO::setBaseName(const std::string & bname) {
  this->base_name = bname;
  this->folder_name = this->base_name;
}

/* -------------------------------------------------------------------------- */
std::string BaseIO::directorySeparator() {
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  return "\\";
#else
  return "/";
#endif
}

/* -------------------------------------------------------------------------- */
void BaseIO::createDirectory(std::string path) {
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  _mkdir(path.c_str());
#else
  mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}

__END_UGUCA__
