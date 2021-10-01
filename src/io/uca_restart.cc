/**
 * @file   uca_restart.cc
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
#include "uca_restart.hh"
#include "static_communicator_mpi.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
Restart::Restart(const std::string &bname,
		 const std::string &path,
		 const Format format) :
  BaseIO() {
  this->initIO(bname, path, format);
  this->precision = 12; // restart should be precise
}

/* -------------------------------------------------------------------------- */
Restart::~Restart() {

}

/* -------------------------------------------------------------------------- */
void Restart::setBaseName(const std::string & bname) {
  this->base_name = bname;
  this->folder_name = this->base_name + "-restart";
}

/* -------------------------------------------------------------------------- */
std::string Restart::getFilePath(const std::string & name, unsigned int number) {

  std::string path_to_file = this->path + BaseIO::directorySeparator()
    + this->folder_name + BaseIO::directorySeparator()
    + name
    + this->rank_str + std::to_string(StaticCommunicatorMPI::getInstance()->whoAmI())
    + ".s" + std::to_string(number)
    + this->file_extension;
  
  return path_to_file;
}

/* -------------------------------------------------------------------------- */
void Restart::dump(unsigned int step) {

  if (!this->initiated) return;

  // open files
  FieldMap::iterator f_it = this->registered_fields.begin();
  FieldMap::iterator f_end = this->registered_fields.end();
  for (; f_it!=f_end; ++f_it) {
    this->open_files[f_it->first] = this->openFile(this->getFilePath(f_it->first, step),
						   std::ios::out);
  }
  HistoryMap::iterator h_it = this->registered_histories.begin();
  HistoryMap::iterator h_end = this->registered_histories.end();
  for (; h_it!=h_end; ++h_it) {
    this->open_files[h_it->first] = this->openFile(this->getFilePath(h_it->first, step),
						   std::ios::out);
  }
  
  BaseIO::dump(step);

  this->closeFiles(true);
}

/* -------------------------------------------------------------------------- */
void Restart::load(unsigned int step) {

  if (!this->initiated) return;

  // open files
  FieldMap::iterator f_it = this->registered_fields.begin();
  FieldMap::iterator f_end = this->registered_fields.end();
  for (; f_it!=f_end; ++f_it) {
    this->open_files[f_it->first] = this->openFile(this->getFilePath(f_it->first, step),
						   std::ios::in);
  }
  HistoryMap::iterator h_it = this->registered_histories.begin();
  HistoryMap::iterator h_end = this->registered_histories.end();
  for (; h_it!=h_end; ++h_it) {
    this->open_files[h_it->first] = this->openFile(this->getFilePath(h_it->first, step),
						   std::ios::in);
  }
    
  BaseIO::load(step);

  this->closeFiles(true);
}


__END_UGUCA__
