#include <pybind11/pybind11.h>

#include "cast.hh"

#include "uca_base_io.hh"
#include "uca_common.hh"
#include "uca_dumper.hh"
#include "uca_restart.hh"

namespace py = pybind11;

namespace uguca {

namespace uguca_wrappers {

using namespace py::literals;

void wrapIo(py::module &mod) {

  py::enum_<BaseIO::Format>(mod, "Format")
      .value("ASCII", BaseIO::Format::ASCII)
      .value("CSV", BaseIO::Format::CSV)
      .value("Binary", BaseIO::Format::Binary);

  py::class_<BaseIO, std::shared_ptr<BaseIO>>(mod, "BaseIO")
      .def("initIO", &BaseIO::initIO)
      .def("load", &BaseIO::load)
      .def("dump", &BaseIO::dump);

  py::class_<Dumper, BaseIO, std::shared_ptr<Dumper>>(mod, "Dumper");

  py::class_<Restart, BaseIO, std::shared_ptr<Restart>>(mod, "Restart")
      .def(py::init<const std::string &, const std::string &,
                    const BaseIO::Format>(),
           py::arg("bname"), py::arg("path"),
           py::arg("format") = BaseIO::Format::ASCII);
}

} // namespace uguca_wrappers

} // namespace uguca
