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

#include <iostream>
#include <typeinfo>
*/
#include <stdexcept>
#include <iomanip>
#include <sstream>

#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
BaseIO::BaseIO() :
  initiated(false),
  base_name("uninitialized"),
  folder_name("uninitialized"),
  path("."),
  dump_format(Format::ASCII),
  separator(" "),
  file_extension(".out") {
}

/* -------------------------------------------------------------------------- */
BaseIO::~BaseIO() {
  this->closeFiles(true);
}

/* -------------------------------------------------------------------------- */
void BaseIO::initIO(const std::string &bname,
		    const std::string &path,
		    const Format format) {

  this->initiated = true;

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
void BaseIO::registerIO(const std::string & name,
			NodalFieldComponent & nodal_field) {
  this->registered_fields[name] = (&nodal_field);
}

/* -------------------------------------------------------------------------- */
void BaseIO::setBaseName(const std::string & bname) {
  this->base_name = bname;
  this->folder_name = this->base_name;
}

/* -------------------------------------------------------------------------- */
std::fstream * BaseIO::openFile(const std::string & path_to_file,
				std::fstream::openmode mode) {

  // open file
  std::fstream * new_file = new std::fstream();
  
  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      new_file->open(path_to_file, mode);
      (*new_file) << std::scientific << std::setprecision(this->precision);
      break;
    }
    case Format::Binary: {
      new_file->open(path_to_file, mode | std::ios::binary); // open as binary file
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }

  return new_file;
}

/* -------------------------------------------------------------------------- */
void BaseIO::closeFiles(bool release_memory) {

  FileMap::iterator it = this->open_files.begin();
  FileMap::iterator end = this->open_files.end();
  
  for (; it != end; ++it) {
    it->second->close();
    if (release_memory) delete it->second;
  }
  this->open_files.clear();
}

/* -------------------------------------------------------------------------- */
void BaseIO::dumpField(std::fstream * dump_file,
		       const NodalFieldComponent & nodal_field) {
  if (!this->initiated) return;

  int nb_nodes = nodal_field.getNbNodes();
  
  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      for (int n = 0; n < nb_nodes; ++n) {
        if (n != 0) (*dump_file) << this->separator;
        (*dump_file) << nodal_field.at(n);
      }
      (*dump_file) << std::endl;
      break;
    }
    case Format::Binary: {
      float temp = 0.0;
      for (int n = 0; n < nb_nodes; ++n) {
        temp = (float)(nodal_field.at(n));
        (*dump_file).write((char *)&temp, sizeof(float));
      }
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }
}

/* -------------------------------------------------------------------------- */
void BaseIO::loadField(std::fstream * load_file,
		       NodalFieldComponent & nodal_field) {
  if (!this->initiated) return;

  int nb_nodes = nodal_field.getNbNodes();
  
  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      std::string line;
      std::getline(*load_file,line);
      std::stringstream ss(line);
      for (int n = 0; n < nb_nodes; ++n) {
	ss >> nodal_field.set(n);
      }
      break;
    }
    case Format::Binary: {
      for (int n = 0; n < nb_nodes; ++n) {
	float temp;
	(*load_file).read((char *)&temp, sizeof(float));
	nodal_field.set(n) = (float)(temp);
      }
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }
}

/* -------------------------------------------------------------------------- */
void BaseIO::dump(unsigned int, double) {

  if (!this->initiated) return;

  FieldMap::iterator it = this->registered_fields.begin();
  FieldMap::iterator end = this->registered_fields.end();

  for (; it!=end; ++it) {
    this->dumpField(this->open_files[it->first], *(it->second));
  }
}

/* -------------------------------------------------------------------------- */
void BaseIO::load(unsigned int) {

  if (!this->initiated) return;

  FieldMap::iterator it = this->registered_fields.begin();
  FieldMap::iterator end = this->registered_fields.end();

  for (; it!=end; ++it) {
    this->loadField(this->open_files[it->first], *(it->second));
  }
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
