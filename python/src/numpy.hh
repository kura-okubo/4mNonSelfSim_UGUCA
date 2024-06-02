#ifndef __NUMPY_HH__
#define __NUMPY_HH__

#include "nodal_field_component.hh"
#include "uca_coommon.hh"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace uguca {

namespace wrap {
namespace py = pybind11;

using numpy = py::array_t<double, py::array::c_style | py::array::forecast>;

class NodalFieldBaseNumpy : public NodalField {
public:
  NodalFieldBaseNumpy(numpy &buffer) : NodalField {
    this->data.wrap(buffer.mutable_data(), buffer.size());
  }
};

} // namespace wrap

} // namespace uguca

#endif
