#include "wrap.hh"
#include "uca_common.hh"
#include "uca_dumper.hh"

namespace uguca {

  namespace wrap {

    using namespace py::literals;

    void wrapIo(py::module& mod) {
      py::class_<Dumper, std::shared_ptr<Dumper>>(mod, "Dumper");
    }

  }
}
