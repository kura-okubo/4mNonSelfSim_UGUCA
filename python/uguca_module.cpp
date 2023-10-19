#include "uca_common.hh"
#include "wrap.hh"

#include <pybind11/pybind11.h>

namespace uguca {

  namespace py = pybind11;

  PYBIND11_MODULE(puguca, mod) {
    mod.doc() = "python module of uguca";

    py::enum_<SolverMethod>(mod, "SolverMethod")
      .value("static", _static)
      .value("quasi_dynamic", _quasi_dynamic)
      .value("dynamic", _dynamic)
      .value("hybrid", _hybrid)
      .export_values();

    wrap::wrapMesh(mod);
  }


}
