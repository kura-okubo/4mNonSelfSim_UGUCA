#ifndef __WRAP_HH__
#define __WRAP_HH__

#include "cast.hh"
#include "uca_common.hh"

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace uguca {

namespace wrap {

void wrapMesh(py::module &mod);
void wrapHalfSpace(py::module &mod);
void wrapInterface(py::module &mod);
void wrapIo(py::module &mod);

} // namespace wrap

} // namespace uguca

#endif
