#include <pybind11/pybind11.h>

#include "half_space.hh"
#include "material.hh"
#include "uca_fftable_mesh.hh"

namespace py = pybind11;

namespace uguca {

namespace uguca_wrappers {

using namespace py::literals;

class PyHalfSpace : public HalfSpace {
public:
  using HalfSpace::HalfSpace;

  void computeStressFourierCoeff(bool predicting = false,
                                 bool correcting = false,
                                 SolverMethod sm = _dynamic,
                                 unsigned int ts_factor = 1) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE_PURE(void, HalfSpace, computeStressFourierCoeff,
                           predicting, correcting, sm, ts_factor);
  }
};

void wrapHalfSpace(py::module &mod) {

  py::class_<Material, std::shared_ptr<Material>>(mod, "Material")
      .def(py::init<double, double, double, bool>(), "E"_a, "nu"_a, "rho"_a,
           "pstress"_a = false)
      .def("readPrecomputedKernels", &Material::readPrecomputedKernels,
           "path"_a = global_kernel_path);

  py::class_<HalfSpace, std::shared_ptr<HalfSpace>, PyHalfSpace>(mod,
                                                                 "HalfSpace")
      .def(py::init<Material &, FFTableMesh &, int, SpatialDirectionSet,
                    const std::string &>(),
           py::arg("material"), py::arg("mesh"), py::arg("side_factor"),
           py::arg("components"), py::arg("name") = "half_space")
      .def("getDisp", py::overload_cast<bool>(&HalfSpace::getDisp, py::const_),
           py::arg("predicting") = false, py::return_value_policy::reference)
      .def("getVelo", py::overload_cast<bool>(&HalfSpace::getVelo, py::const_),
           py::arg("predicting") = false, py::return_value_policy::reference);
}

} // namespace uguca_wrappers

} // namespace uguca
