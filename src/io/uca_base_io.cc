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
  separator(' '),
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
    this->separator = ' ';
    this->file_extension = ".out";
    break;
  }
  case Format::CSV: {
    this->separator = ',';
    this->file_extension = ".csv";
    break;
  }
  case Format::Binary: {
    this->separator = ' ';
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
  if (this->registered_fields.find(name) == this->registered_fields.end())
    this->registered_fields[name] = (&nodal_field);
  else
    throw std::runtime_error("Field already registered: "+name);
}

/* -------------------------------------------------------------------------- */
void BaseIO::registerIO(const std::string & name,
			NodalField & nodal_field) {
  for (int d=0; d<nodal_field.getDim(); ++d) {
    std::string component_name = name + "_" + std::to_string(d);
    this->registerIO(component_name, nodal_field.component(d));
  }
}

/* -------------------------------------------------------------------------- */
void BaseIO::registerIO(const std::string & name,
			ModalLimitedHistory & lim_history) {
  if (this->registered_histories.find(name) == this->registered_histories.end())
    this->registered_histories[name] = (&lim_history);
  else
    throw std::runtime_error("Limited History already registered: "+name);
}

/* -------------------------------------------------------------------------- */
void BaseIO::registerIO(NodalFieldComponent & nodal_field) {
  this->registerIO(nodal_field.getName(), nodal_field);
}

/* -------------------------------------------------------------------------- */
void BaseIO::registerIO(NodalField & nodal_field) {
  for (int d=0; d<nodal_field.getDim(); ++d) {
    this->registerIO(nodal_field.component(d));
  }
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
  const double * nf_data = nodal_field.storage();

  this->write(dump_file, nf_data, nb_nodes);
}

/* -------------------------------------------------------------------------- */
void BaseIO::dumpHistory(std::fstream * dump_file,
			 const ModalLimitedHistory & limited_history) {
  if (!this->initiated) return;

  // write nb_history points and index_now
  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      (*dump_file) << limited_history.getSize() << this->separator;
      (*dump_file) << limited_history.getNbHistoryPoints() << this->separator;
      (*dump_file) << limited_history.getIndexNow() << std::endl;
      break;
    }
    case Format::Binary: {
      float temp = (float)(limited_history.getSize());
      (*dump_file).write((char *)&temp, sizeof(float));
      temp = (float)(limited_history.getNbHistoryPoints());
      (*dump_file).write((char *)&temp, sizeof(float));
      temp = (float)(limited_history.getIndexNow());
      (*dump_file).write((char *)&temp, sizeof(float));
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }

  int size = limited_history.getSize();
  const double * lh_data = limited_history.getValues();
  this->write(dump_file, lh_data, size);
}

/* -------------------------------------------------------------------------- */
void BaseIO::write(std::fstream * dump_file,
		   const double * data, int size) {
  
  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      for (int n = 0; n < size; ++n) {
        if (n != 0) (*dump_file) << this->separator;
        (*dump_file) << data[n];
      }
      (*dump_file) << std::endl;
      break;
    }
    case Format::Binary: {
      float temp = 0.0;
      for (int n = 0; n < size; ++n) {
        temp = (float)(data[n]);
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
  double * nf_data = nodal_field.storage();

  this->read(load_file, nf_data, nb_nodes);
}

/* -------------------------------------------------------------------------- */
void BaseIO::loadHistory(std::fstream * load_file,
			 ModalLimitedHistory & limited_history) {
  if (!this->initiated) return;

  // read nb_history points and index_now
  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      std::string line;
      std::getline(*load_file,line);
      std::stringstream ss(line);
      double temp;
      ss >> temp;
      if (temp != limited_history.getSize()) 
	throw std::runtime_error("reloaded Limited History is of incorrect size");
      ss >> temp;
      limited_history.setNbHistoryPoints((int)temp);
      ss >> temp;
      limited_history.setIndexNow((int)temp);
      break;
    }
    case Format::Binary: {
      float temp;
      (*load_file).read((char *)&temp, sizeof(float));
      if (temp != limited_history.getSize()) 
	throw std::runtime_error("reloaded Limited History is of incorrect size");
      (*load_file).read((char *)&temp, sizeof(float));
      limited_history.setNbHistoryPoints((int)temp);
      (*load_file).read((char *)&temp, sizeof(float));
      limited_history.setIndexNow((int)temp);
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }
  
  int size = limited_history.getSize();
  double * lh_data = limited_history.getValues();
  this->read(load_file, lh_data, size);
}

/* -------------------------------------------------------------------------- */
void BaseIO::read(std::fstream * load_file,
		  double * data,
		  int size) {
  
  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      std::string line;
      std::getline(*load_file,line);
      std::stringstream ss(line);
      // count number of entries
      int count = 0;
      double temp = 0.;
      while (ss >> temp) {
	count++;
      }
      if (count != size) 
	throw std::runtime_error("BaseIO read: wrong number of entries "+std::to_string(count)+"!="+std::to_string(size));
      // get values
      ss = std::stringstream(line); // rewind line
      for (int n = 0; n < size; ++n) {
	ss >> data[n];
      }
      break;
    }
    case Format::Binary: {
      float temp = 0.;
      std::fstream * to_count = load_file;
      int count = 0;
      while ((*to_count).read((char *)&temp, sizeof(float))) {
	count++;
      }
      if (count != size) 
	throw std::runtime_error("BaseIO read: wrong number of entries "+std::to_string(count)+"!="+std::to_string(size));
      // get values
      for (int n = 0; n < size; ++n) {
	(*load_file).read((char *)&temp, sizeof(float));
	data[n] = (float)(temp);
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

  FieldMap::iterator f_it = this->registered_fields.begin();
  FieldMap::iterator f_end = this->registered_fields.end();

  for (; f_it!=f_end; ++f_it) {
    this->dumpField(this->open_files[f_it->first], *(f_it->second));
  }

  HistoryMap::iterator h_it = this->registered_histories.begin();
  HistoryMap::iterator h_end = this->registered_histories.end();

  for (; h_it!=h_end; ++h_it) {
    this->dumpHistory(this->open_files[h_it->first], *(h_it->second));
  }
}

/* -------------------------------------------------------------------------- */
void BaseIO::load(unsigned int) {

  if (!this->initiated) return;

  FieldMap::iterator f_it = this->registered_fields.begin();
  FieldMap::iterator f_end = this->registered_fields.end();

  for (; f_it!=f_end; ++f_it) {
    this->loadField(this->open_files[f_it->first], *(f_it->second));
  }

  HistoryMap::iterator h_it = this->registered_histories.begin();
  HistoryMap::iterator h_end = this->registered_histories.end();

  for (; h_it!=h_end; ++h_it) {
    this->loadHistory(this->open_files[h_it->first], *(h_it->second));
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
