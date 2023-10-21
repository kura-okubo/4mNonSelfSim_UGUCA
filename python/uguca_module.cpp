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
      .value("hybrid", SolverMethod::_hybrid)
      .export_values();

    wrap::wrapMesh(mod);
    wrap::wrapHalfSpace(mod);
    wrap::wrapInterface(mod);
    wrap::wrapIo(mod);
  }


}
