#ifndef __WRAP_HH__
#define __WRAP_HH__

#include "uca_common.hh"

#include <pybind11/pybind11.h>

namespace uguca {

namespace wrap {

  namespace py = pybind11;

  void wrapMesh(py::module& mod);
  
}
  
}


#endif
