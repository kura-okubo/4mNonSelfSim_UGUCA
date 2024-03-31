#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "uca_common.hh"
#include "uca_base_mesh.hh"
#include "uca_fftable_mesh.hh"
#include "uca_distributed_fftable_mesh.hh"
#include "uca_simple_mesh.hh"
#include "nodal_field.hh"

#include "wrap.hh"

namespace uguca
{

	namespace wrap
	{

		using namespace py::literals;

		void wrapMesh(py::module &mod)
		{
			py::class_<BaseMesh, std::shared_ptr<BaseMesh>>(mod, "BaseMesh")
				.def("getLocalCoords", 
				     [](BaseMesh &self)
					 { py::array ret =  py::cast(self.getLocalCoords().getStorage());
					   return ret;
					 });

			py::class_<FFTableMesh, BaseMesh,
					   std::shared_ptr<FFTableMesh>>(mod, "FFTableMesh");

			py::class_<DistributedFFTableMesh, FFTableMesh,
					   std::shared_ptr<DistributedFFTableMesh>>(mod,
																"DistributedFFTableMesh");

			py::class_<SimpleMesh, DistributedFFTableMesh,
					   std::shared_ptr<SimpleMesh>>(mod, "SimpleMesh")
				.def(py::init<double, int>(),
					 py::arg("Lx"), py::arg("Nx"))
				.def(py::init<double, int, double, int>(),
					 py::arg("Lx"), py::arg("Nx"), py::arg("Lz"), py::arg("Nz"))
				.def("getNbLocalNodesAlloc",
					 &SimpleMesh::getNbLocalNodesAlloc);

			/*py::class_<NodalField, std::shared_ptr<NodalField>>(mod, "NodalField")
				.def(py::init<const std::string &>(),
					 py::arg("name") = "unnamed")
				.def(py::init<BaseMesh &, std::set<int>, const std::string &>(),
					 py::arg("mesh"),
					 py::arg("components"),
					 py::arg("name") = "unnamed")
				.def("setAllValuesTo",
					 &NodalField::setAllValuesTo)
				.def("zeros",
					 &NodalField::zeros);*/
		}

	}

}
