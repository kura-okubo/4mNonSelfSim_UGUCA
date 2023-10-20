#ifndef __NUMPY_HH__
#define __NUMPY_HH__

#include "uca_coommon.hh"
#include "nodal_field_component.hh"


#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace uguca {

namespace wrap {
namespace py = pybind11;

  
using numpy = py::array_t<double, py::array::c_style | py::array::forecast>;

  
class NodalFieldComponentNumpy : public NodalFieldComponent {
public:
  NodalFieldComponent(numpy & buffer) : NodalFieldComponent {
    this->nb_components = 1;
    this->data.wrapMemory(buffer.mutable_data(), buffer.size());
  }
};


  
}// namespace wrap


}

#endif
