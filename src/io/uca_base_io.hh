/**
 * @file   uca_base_io.hh
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
#ifndef __BASE_IO_H__
#define __BASE_IO_H__

/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "nodal_field_component.hh"

#include <fstream>
#include <map>
//#include <sstream>
//#include <string>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class BaseIO {
public:
  enum class Format { ASCII, CSV, Binary };
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  typedef std::map<const std::string, NodalFieldComponent *> FieldMap;
  typedef std::map<const std::string, std::fstream *> FileMap;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  BaseIO();
  virtual ~BaseIO();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // sets variables and creates folder
  virtual void initIO(const std::string &bname,
		      const std::string &path,
		      const Format format = Format::ASCII);

  virtual void registerIO(const std::string & name,
			  NodalFieldComponent & nodal_field);

  virtual void dump(unsigned int step, double time = 0.);
  virtual void load(unsigned int step);

protected:

  virtual void setBaseName(const std::string & bname);
  std::fstream * openFile(const std::string & path_to_file,
			  std::fstream::openmode mode); 
  virtual void closeFiles(bool release_memory);

  void dumpField(std::fstream * dump_file,
		 const NodalFieldComponent & nodal_field);
  void loadField(std::fstream * load_file,
		 NodalFieldComponent & nodal_field);

  
  /* ------------------------------------------------------------------------ */
  /* File system related methods                                              */
  /* ------------------------------------------------------------------------ */
public:
  static std::string directorySeparator();
  static void createDirectory(std::string path);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
 /* ------------------------------------------------------------------------ */
 /* Class Members                                                            */
 /* ------------------------------------------------------------------------ */
protected:
  // has dump been initiated?
  bool initiated;

  // base name
  std::string base_name;
  std::string folder_name;

  // path to dumped files
  std::string path;

  // dump format
  Format dump_format;

  // characteristics of dumper
  std::string separator;

  // files extention
  std::string file_extension;

  // rank string
  std::string rank_str = ".proc";
  
  // precision of dump to text
  int precision = 6;

  // registered nodal fields
  FieldMap registered_fields;

  // open files
  FileMap open_files;
};

__END_UGUCA__

//#include "uca_base_io_impl.cc"

#endif /* __BASE_IO_H__ */
