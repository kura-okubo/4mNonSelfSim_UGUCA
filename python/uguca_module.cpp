#include "uca_common.hh"
#include "wrap.hh"

#include <pybind11/pybind11.h>

namespace uguca {

  namespace py = pybind11;

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


    wrap::wrapMesh(mod);
    wrap::wrapHalfSpace(mod);
    wrap::wrapInterface(mod);
    wrap::wrapIo(mod);
  }


}
