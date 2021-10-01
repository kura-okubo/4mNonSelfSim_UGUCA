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
#include "uca_custom_mesh.hh"

#include <iomanip>
#include <typeinfo>
#include <sstream>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
Dumper::Dumper(BaseMesh & mesh) :
  BaseIO(),
  mesh(mesh) {
  
  this->time_file = NULL;
  this->field_file = NULL;
  
  this->parallel_dump = false;

  if (typeid(mesh) == typeid(CustomMesh))
    this->parallel_dump = true;
}

/* -------------------------------------------------------------------------- */
Dumper::~Dumper() {}

/* -------------------------------------------------------------------------- */
void Dumper::closeFiles(bool release_memory) {

  BaseIO::closeFiles(release_memory);
  
  int rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int root = this->mesh.getRoot();
  
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
// registers field and opens file
void Dumper::registerIO(const std::string & field_name,
			NodalFieldComponent & nodal_field) {

  // define path to file
  int rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::string rank_name = "";
  if (this->parallel_dump)
    rank_name = this->rank_str + std::to_string(rank);

  std::string file_name = field_name + this->file_extension + rank_name;
  std::string path_to_file = this->path + BaseIO::directorySeparator() +
                             this->folder_name + BaseIO::directorySeparator() +
                             file_name;

  // open file and keep reference to it (only if there are nodes)
  if (this->mesh.getNbLocalNodes() > 0) {
    BaseIO::registerIO(field_name, nodal_field);
    this->open_files[field_name] = this->openFile(path_to_file, std::ios::out);
  }

  // put info into field file
  int root = this->mesh.getRoot();
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

  BaseIO::dump(step,time);
  
  if (!this->initiated) return;

  int rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int root = this->mesh.getRoot();
  if (rank == root) {
    (*this->time_file) << step << this->separator << time << std::endl;
  }
}

__END_UGUCA__
