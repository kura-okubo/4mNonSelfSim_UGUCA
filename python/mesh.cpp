#include "wrap.hh"
#include "uca_base_mesh.hh"
#include "uca_fftable_mesh.hh"
#include "uca_distributed_fftable_mesh.hh"
#include "uca_simple_mesh.hh"


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
	     "Lx"_a, "Nx"_a, "initialize"_a=true)
	.def(py::init<double, int, double, int, bool>(),
	     "Lx"_a, "Nx"_a, "Lz"_a, "Nz"_a, "intialize"_a=true);
    }
    
  }


}
