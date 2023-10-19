#include "wrap.hh"
#include "interface_law.hh"
#include "barras_law.hh"
#include "uca_base_mesh.hh"
#include "uca_fftable_mesh.hh"
#include "bimat_interface.hh"
#include "material.hh"

namespace uguca {

  namespace wrap {

    using namespace py::literals;

    void wrapInterface(py::module& mod) {
      py::class_<BimatInterface,
		 std::shared_ptr<BimatInterface>>(mod, "BimatInterface")
	.def(py::init<FFTableMesh&, Material&, Material&, InterfaceLaw&,
	     const SolverMethod&>(),
	     "mesh"_a, "top_material"_a, "bot_material"_a, "law"_a,
	     "method"_a=_dynamic);

      
      py::class_<InterfaceLaw,
		 std::shared_ptr<InterfaceLaw>>(mod, "InterfaceLaw");
      
      py::class_<BarrasLaw, InterfaceLaw,
		 std::shared_ptr<BarrasLaw>>(mod, "BarrasLaw")
	.def(py::init<BaseMesh&, double, double, std::string&>(),
	     "mesh"_a, "tau_max_default"_a, "delta_c_default"_a,
	     "name"_a="blaw");
    }


  }

}
