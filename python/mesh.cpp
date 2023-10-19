#include "wrap.hh"
#include "uca_base_mesh.hh"
#include "uca_fftable_mesh.hh"
#include "uca_distributed_fftable_mesh.hh"
#include "uca_simple_mesh.hh"
#include "nodal_field_component.hh"
#include "nodal_field.hh"

namespace uguca {

  namespace wrap {

    using namespace py::literals;

    void wrapMesh(py::module& mod) {
      py::class_<BaseMesh, std::shared_ptr<BaseMesh>>(mod, "BaseMesh");

      py::class_<FFTableMesh, BaseMesh,
		 std::shared_ptr<FFTableMesh>>(mod, "FFTableMesh");
      
      py::class_<DistributedFFTableMesh, FFTableMesh,
	  std::shared_ptr<DistributedFFTableMesh>>(mod,
						   "DistributedFFTableMesh");

      py::class_<SimpleMesh, DistributedFFTableMesh,
		 std::shared_ptr<SimpleMesh>>(mod, "SimpleMesh")
	.def(py::init<double, int, bool>(),
	     py::arg("Lx"), py::arg("Nx"), py::arg("initialize")=true)
	.def(py::init<double, int, double, int, bool>(),
	     py::arg("Lx"), py::arg("Nx"), py::arg("Lz"), py::arg("Nz"),
	     py::arg("intialize")=true);

      py::class_<NodalFieldComponent,
	  std::shared_ptr<NodalFieldComponent>>(mod, "NodalFieldComponent")
	.def(py::init<const std::string&>(),
	     py::arg("name")="unnamed")
	.def(py::init<BaseMesh&, int, const std::string&>(),
	     py::arg("mesh"), py::arg("direction")=0,
	     py::arg("name")="unnamed")
	.def("setAllValuesTo",
	     &NodalFieldComponent::setAllValuesTo)
	.def("zeros",
	     &NodalFieldComponent::zeros);

      py::class_<NodalField, std::shared_ptr<NodalField>>(mod, "NodalField")
	.def(py::init<const std::string&>(),
	     py::arg("name")="unnamed")
	.def("component",
	     &NodalField::component, py::return_value_policy::reference)
	.def("setAllValuesTo",
	     &NodalField::setAllValuesTo)
	.def("zeros",
	     &NodalField::zeros);

    }
    
  }


}
