###
# @file   CMakeLists.txt
#
# @author David S. Kammer <dkammer@ethz.ch>
# @author Gabriele Albertini <ga288@cornell.edu>
# @author Chun-Yu Ke <ck659@cornell.edu>
#
# @date creation: Fri Feb 5 2021
# @date last modification: Fri Feb 5 2021
#
# @brief  TODO
#
#
# Copyright (C) 2021 ETH Zurich (David S. Kammer)
#
# This file is part of uguca.
#
# uguca is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# uguca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with uguca.  If not, see <https://www.gnu.org/licenses/>.
####
set(UGUCA_HEADER
  uca_common.hh

  # MPI static communicator
  static_communicator_mpi.hh

  # io
  io/uca_input_section.hh
  io/uca_parameter_reader.hh
  io/uca_base_io.hh
  io/uca_dumper.hh
  io/uca_restart.hh
  
  # meshes
  mesh/uca_base_mesh.hh
  mesh/uca_fftable_mesh.hh
  mesh/uca_distributed_fftable_mesh.hh
  mesh/uca_simple_mesh.hh
  mesh/uca_custom_mesh.hh

  # nodal fields
  mesh/nodal_field_component.hh
  mesh/nodal_field.hh
  mesh/fftable_nodal_field_component.hh
  mesh/fftable_nodal_field.hh
  
  # kernels and material
  half_space/kernel.hh
  half_space/precomputed_kernel.hh
  half_space/kernel_collection.hh
  half_space/material.hh

  # half spaces
  half_space/limited_history.hh
  half_space/preint_kernel.hh
  half_space/half_space.hh
  half_space/half_space_quasi_dynamic.hh
  half_space/half_space_dynamic.hh

  # interfaces
  interface/interface.hh
  interface/bimat_interface.hh
  interface/defrig_interface.hh
  interface/unimat_shear_interface.hh
  interface/unimat_normal_interface.hh

  # interface laws
  interface/laws/interface_law.hh
  interface/laws/linear_coulomb_friction_law.hh
  interface/laws/linear_shear_cohesive_law.hh
  interface/laws/linear_normal_cohesive_law.hh
  interface/laws/barras_law.hh
  interface/laws/rate_and_state_law.hh

  # use method as infinite boundary condition for FEM
  interface/infinite_boundary/infinite_boundary.hh
  )

set(UGUCA_SRC

  # MPI static communicator
  static_communicator_mpi.cc

  # io
  io/uca_input_section.cc
  io/uca_parameter_reader.cc
  io/uca_base_io.cc
  io/uca_dumper.cc  
  io/uca_restart.cc

  # meshes
  mesh/uca_base_mesh.cc
  mesh/uca_fftable_mesh.cc
  mesh/uca_distributed_fftable_mesh.cc
  mesh/uca_simple_mesh.cc
  mesh/uca_custom_mesh.cc

  # nodal fields
  mesh/nodal_field_component.cc
  mesh/nodal_field.cc
  mesh/fftable_nodal_field_component.cc
  mesh/fftable_nodal_field.cc

  # kernels and material
  half_space/kernel.cc
  half_space/precomputed_kernel.cc
  half_space/kernel_collection.cc
  half_space/material.cc

  # half spaces
  half_space/limited_history.cc
  half_space/preint_kernel.cc
  half_space/half_space.cc
  half_space/half_space_dynamic.cc
  half_space/half_space_quasi_dynamic.cc

  # interfaces
  interface/interface.cc
  interface/bimat_interface.cc
  interface/defrig_interface.cc
  interface/unimat_shear_interface.cc
  interface/unimat_normal_interface.cc

  # interface laws
  interface/laws/interface_law.cc
  interface/laws/linear_coulomb_friction_law.cc
  interface/laws/linear_shear_cohesive_law.cc
  interface/laws/linear_normal_cohesive_law.cc
  interface/laws/barras_law.cc
  interface/laws/rate_and_state_law.cc

  # use method as infinite boundary condition for FEM
  interface/infinite_boundary/infinite_boundary.cc
  )
set(UGUCA_DEPEND_FILES ${UGUCA_SRC} ${UGUCA_HEADER} PARENT_SCOPE)

add_library(uguca ${UGUCA_SRC})

target_link_libraries(uguca ${MPI_LIBRARIES} ${FFTW_LIBRARIES}
  ${BLAS_LIBRARIES})


list(APPEND UGUCA_EXPORT_LIST uguca)

if("${UGUCA_TARGETS_EXPORT}" STREQUAL "ugucaLibraryDepends")
  install(EXPORT ugucaLibraryDepends DESTINATION  lib/uguca
    COMPONENT dev)
endif()

#Export for build tree
export(TARGETS ${UGUCA_EXPORT_LIST}
  FILE "${CMAKE_BINARY_DIR}/ugucaLibraryDepends.cmake")
export(PACKAGE uguca)
