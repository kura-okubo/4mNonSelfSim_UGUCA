/**
 * @file   uca_dumper.hh
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
#ifndef __DUMPER_H__
#define __DUMPER_H__

/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_base_io.hh"
#include "uca_base_mesh.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class Dumper : public BaseIO {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Dumper(BaseMesh & mesh);
  virtual ~Dumper();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initDump(const std::string &bname,
			const std::string &path,
                        const Format format = Format::ASCII);

  void registerIO(const std::string & field_name,
		  NodalField & nodal_field);

  void dump(unsigned int step, double time);

  // specific for dumper
  virtual void registerDumpField(const std::string & field_name) = 0;
  void registerDumpFields(const std::string & field_names,
			  char delim = ',');
  
 protected:
  virtual void setBaseName(const std::string & bname);

  void setCoords(std::ofstream * cfile);


protected:
  void closeFiles(bool release_memory);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
 /* ------------------------------------------------------------------------ */
 /* Class Members                                                            */
 /* ------------------------------------------------------------------------ */
protected:
  BaseMesh & mesh;

  bool parallel_dump = false;
  
private:
  // information based on base_name
  std::string info_file_name;
  std::string time_file_name;
  std::string coord_file_name;
  std::string field_file_name;
  std::string proc_file_name;

  // files corresponding to field
  // FileToFieldMap files_and_fields;

  // file with time stamps
  std::ofstream * time_file;

  // file with coord
  std::ofstream * coord_file;

  // file with field infos
  std::ofstream * field_file;
};

__END_UGUCA__

//#include "uca_dumper_impl.cc"

#endif /* __DUMPER_H__ */
