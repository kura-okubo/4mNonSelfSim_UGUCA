#include "wrap.hh"
#include "material.hh"



namespace uguca {

  namespace wrap {

    using namespace py::literals;

    void wrapHalfSpace(py::module& mod) {

      
      py::class_<Material, std::shared_ptr<Material>>(mod, "Material")
	.def(py::init<double, double, double, bool>(),
	     "E"_a, "nu"_a, "rho"_a, "pstress"_a=false)
	.def("readPrecomputedKernels",
	     &Material::readPrecomputedKernels, 
	     "path"_a=global_kernel_path);

    }


  }

}

