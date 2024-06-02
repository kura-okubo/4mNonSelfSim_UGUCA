#include "uca_common.hh"

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace uguca {

namespace uguca_wrappers {

void wrapMesh(py::module &mod);
void wrapHalfSpace(py::module &mod);
void wrapInterface(py::module &mod);
void wrapIo(py::module &mod);

} // namespace uguca_wrappers

PYBIND11_MODULE(puguca, mod) {
  mod.doc() = "python module of uguca";

  py::enum_<SolverMethod>(mod, "SolverMethod")
      .value("static", SolverMethod::_static)
      .value("quasi_dynamic", SolverMethod::_quasi_dynamic)
      .value("dynamic", SolverMethod::_dynamic)
      .value("adaptive", SolverMethod::_adaptive)
      .export_values();

  py::enum_<SpatialDirection>(mod, "SpatialDirection")
      .value("x", SpatialDirection::_x)
      .value("y", SpatialDirection::_y)
      .value("z", SpatialDirection::_z)
      .value("spatial_dir_count", SpatialDirection::_spatial_dir_count)
      .export_values();

  // Add submodules
  auto mesh_mod = mod.def_submodule("mesh", "mesh submodule");
  uguca_wrappers::wrapMesh(mesh_mod);

  auto module_half_space =
      mod.def_submodule("half_space", "half_space submodule");
  uguca_wrappers::wrapHalfSpace(module_half_space);

  auto module_interface = mod.def_submodule("interface", "interface submodule");
  uguca_wrappers::wrapInterface(module_interface);

  auto module_io = mod.def_submodule("io", "io submodule");
  uguca_wrappers::wrapIo(module_io);
}

} // namespace uguca
