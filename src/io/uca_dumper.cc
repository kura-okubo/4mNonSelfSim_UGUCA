/**
 * @file   uca_dumper.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Fri Feb 5 2021
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
#include "uca_dumper.hh"
#include "nodal_field_component.hh"
#include "static_communicator_mpi.hh"
#include "uca_custom_mesh.hh"

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <typeinfo>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
Dumper::Dumper(BaseMesh & mesh) :
  BaseIO(),
  mesh(mesh) {
  
  // default name and path
  this->setBaseName("standard-bname");
  this->path = ".";
  
  this->time_file = NULL;
  this->field_file = NULL;
  
  this->initiated = false;
  this->parallel_dump = false;

  if (typeid(mesh) == typeid(CustomMesh))
    this->parallel_dump = true;
}

/* -------------------------------------------------------------------------- */
Dumper::~Dumper() {
  this->closeFiles(true);
}

/* -------------------------------------------------------------------------- */
void Dumper::closeFiles(bool release_memory) {

  int rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int root = this->mesh.getRoot();

  FileToFieldMap::iterator it = this->files_and_fields.begin();
  FileToFieldMap::iterator end = this->files_and_fields.end();
  
  for (; it != end; ++it) {
    it->first->close();
    if (release_memory) delete it->first;
  }
  
  if (this->initiated) {
    if (rank == root) {

      this->time_file->close();
      if (release_memory) delete this->time_file;
      
      this->coord_file->close();
      if (release_memory) delete this->coord_file;
      
      this->field_file->close();
      if (release_memory) delete this->field_file;
    }
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::setBaseName(const std::string & bname) {
  
  this->base_name = bname;
  
  this->info_file_name  = this->base_name + ".info";
  this->time_file_name  = this->base_name + ".time";
  this->coord_file_name = this->base_name + ".coords";
  this->field_file_name = this->base_name + ".fields";
  this->proc_file_name  = this->base_name + ".procs";
  this->folder_name     = this->base_name + "-DataFiles";
}

/* -------------------------------------------------------------------------- */
void Dumper::initDump(const std::string & bname,
		      const std::string & path,
                      const Format format) {

  this->initIO(bname, path, format);
  
  // only root dumps global data
  int rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int root = this->mesh.getRoot();
  
  this->initiated = true;

  // < --------------------------------------------------- only need of_string
  std::string of_string = "output_format undefined";
  switch (this->dump_format) {
  case Format::ASCII: {
    of_string = "output_format ascii";
    break;
  }
  case Format::CSV: {
    of_string = "output_format csv";
    break;
  }
  case Format::Binary: {
    of_string = "output_format binary";
    break;
  }
  default:
    throw std::runtime_error("Unsupported output format.");
  }
  
  if (rank == root) {
    
    // info file
    std::string path_to_info_file =
      this->path + BaseIO::directorySeparator() + this->info_file_name;
    std::ofstream info_file(path_to_info_file, std::ios::out);
    info_file << "field_description " << this->field_file_name << std::endl;
    info_file << "time_description " << this->time_file_name << std::endl;
    info_file << "coord_description " << this->coord_file_name << std::endl;
    if (this->parallel_dump)
      info_file << "rank_description " << this->proc_file_name << std::endl;
    info_file << "folder_name " << this->folder_name << std::endl;
    info_file << of_string << std::endl;
    info_file.close();

    // time file
    std::string path_to_time_file =
      this->path + BaseIO::directorySeparator() + this->time_file_name;
    
    this->time_file = new std::ofstream(path_to_time_file, std::ios::out);
    
    (*this->time_file) << std::scientific << std::setprecision(this->precision);

    // coord file
    std::string path_to_coord_file =
      this->path + BaseIO::directorySeparator() + this->coord_file_name;
    
    this->coord_file = new std::ofstream(path_to_coord_file, std::ios::out);
    
    (*this->coord_file) << std::scientific << std::setprecision(10);

    // field file
    std::string path_to_field_file =
      this->path + BaseIO::directorySeparator() + this->field_file_name;
    
    this->field_file = new std::ofstream(path_to_field_file, std::ios::out);

    // write coords file
    this->setCoords(this->coord_file);
  } // end of rank == root

  // if parallel dump: dump nb procs and additional coords
  if (this->parallel_dump) {

    std::string parallel_coord_file_name = "coords";
    
    // write rank string and number of procs to file
    if (rank == root) {
      std::string path_to_proc_file =
	this->path + BaseIO::directorySeparator() + this->proc_file_name;
      std::ofstream proc_file(path_to_proc_file, std::ios::out);
      proc_file << "coords_name " << parallel_coord_file_name << std::endl;
      proc_file << "rank_string " << this->rank_str << std::endl;
      proc_file << "psize " << StaticCommunicatorMPI::getInstance()->getNbProc() << std::endl;
      proc_file.close();
    } // end root
  
    // only if there are nodes
    if (this->mesh.getNbLocalNodes() > 0) {

      // file name: coords.prank
      std::string path_to_coord_file =
	this->path + BaseIO::directorySeparator() + this->folder_name
	+ BaseIO::directorySeparator() + parallel_coord_file_name
	+ this->rank_str + std::to_string(rank);
      std::ofstream p_coord_file(path_to_coord_file, std::ios::out);
      p_coord_file << std::scientific << std::setprecision(this->precision);
      
      // write coords file
      this->setCoords(&p_coord_file);
    } // end: nb_nodes > 0
  } // end: parallel dump
}

/* -------------------------------------------------------------------------- */
void Dumper::registerForDump(const std::string & field_name,
			     const NodalFieldComponent & nodal_field) {

  int rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int root = this->mesh.getRoot();

  std::string rank_name = "";
  if (this->parallel_dump)
    rank_name = this->rank_str + std::to_string(rank);

  std::string file_name = field_name + this->file_extension + rank_name;
  std::string path_to_file = this->path + BaseIO::directorySeparator() +
                             this->folder_name + BaseIO::directorySeparator() +
                             file_name;

  // open file
  std::ofstream * new_file = new std::ofstream();

  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      new_file->open(path_to_file, std::ios::out);
      (*new_file) << std::scientific << std::setprecision(this->precision);
      break;
    }
    case Format::Binary: {
      new_file->open(path_to_file, std::ios::out | std::ios::binary);  // open as binary file
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }

  // keep reference to file (only if there are nodes)
  if (this->mesh.getNbLocalNodes() > 0)
    this->files_and_fields[new_file] = (&nodal_field);

  // put info into field file
  if (rank == root) {
    std::string fname = field_name + this->file_extension;
    (*this->field_file) << field_name << " " << fname << std::endl;
  }
}

/* -------------------------------------------------------------------------- */
// coords has all coordinates of all points
// prints coords to file
void Dumper::setCoords(std::ofstream * cfile) {

  if (!this->initiated) return;
  double ** coords = this->mesh.getLocalCoords();

  for (int n=0; n<this->mesh.getNbLocalNodes(); ++n) {
    for (int d=0; d<this->mesh.getDim(); ++d) {
      if (d != 0) {
	(*cfile) << this->separator;
      }
      (*cfile) << coords[d][n];
    }
    (*cfile) << std::endl;
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::registerDumpFields(const std::string & field_names,
				char delim) {

  // separate string by delim
  std::vector<std::string> field_name_list;
  std::stringstream ss(field_names);
  std::string substring;
  while (std::getline(ss, substring, delim)) {
    field_name_list.push_back(substring);
  }

  // register dump fields
  for (std::vector<std::string>::iterator it = field_name_list.begin();
       it != field_name_list.end();
       ++it) {
    this->registerDumpField(*it);
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::dump(unsigned int step, double time) {

  if (!this->initiated) return;

  FileToFieldMap::iterator it = this->files_and_fields.begin();
  FileToFieldMap::iterator end = this->files_and_fields.end();

  for (; it!=end; ++it) {
    this->dumpField(it->first, it->second);
  }

  int rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int root = this->mesh.getRoot();
  if (rank == root) {
    (*this->time_file) << step << this->separator << time << std::endl;
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::dumpField(std::ofstream * dump_file,
		       const NodalFieldComponent * nodal_field) {
  if (!this->initiated) return;

  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      for (int n = 0; n < this->mesh.getNbLocalNodes(); ++n) {
        if (n != 0) (*dump_file) << this->separator;
        (*dump_file) << nodal_field->at(n);
      }
      (*dump_file) << std::endl;
      break;
    }
    case Format::Binary: {
      float temp = 0.0;
      for (int n = 0; n < this->mesh.getNbLocalNodes(); ++n) {
        temp = (float)(nodal_field->at(n));
        (*dump_file).write((char *)&temp, sizeof(float));
      }
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }
}

__END_UGUCA__
