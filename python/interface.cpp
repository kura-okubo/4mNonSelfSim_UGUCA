#include "wrap.hh"
#include "interface_law.hh"
#include "barras_law.hh"
#include "linear_coulomb_friction_law.hh"
#include "uca_base_mesh.hh"
#include "uca_fftable_mesh.hh"
#include "interface.hh"
#include "bimat_interface.hh"
#include "unimat_shear_interface.hh"
#include "defrig_interface.hh"
#include "material.hh"
#include "nodal_field.hh"

#include "uca_dumper.hh"

namespace uguca
{

	namespace wrap
	{

		using namespace py::literals;

		class PyInterfaceLaw : public InterfaceLaw
		{
		public:
			using InterfaceLaw::InterfaceLaw;

			void computeCohesiveForces(bool predicting = false) override
			{
				// NOLINTNEXTLINE
				PYBIND11_OVERRIDE_PURE(void, InterfaceLaw, computeCohesiveForces, predicting);
			}
		};

		void wrapInterface(py::module &mod)
		{

			py::class_<Interface, std::shared_ptr<Interface>>(mod, "Interface")
				.def("getLoad",
					 &Interface::getLoad,
					 py::return_value_policy::reference)
				.def("getCohesion",
					 &Interface::getCohesion,
					 py::return_value_policy::reference);

			py::class_<BimatInterface,
					   std::shared_ptr<BimatInterface>, Interface>(mod, "BimatInterface")
				.def(py::init<FFTableMesh &, SpatialDirectionSet, Material &, Material &, InterfaceLaw &,
							  const SolverMethod &>(),
					 py::arg("mesh"), py::arg("components"), py::arg("top_material"), py::arg("bot_material"), py::arg("law"),
					 py::arg("method") = _dynamic)
				.def("init",
					 &BimatInterface::init)
				.def("setTimeStep",
					 &BimatInterface::setTimeStep)
				.def("getStableTimeStep",
					 &BimatInterface::getStableTimeStep,
					 py::return_value_policy::reference)
				.def("initDump", [](BimatInterface &self, const std::string &bname,
									const std::string &path)
					 { self.initDump(bname, path); })
				.def("computeGap", [](BimatInterface &self, NodalField &gap, bool predicting)
					 { self.computeGap(gap, predicting); })
				.def("closingNormalGapForce",
					 &BimatInterface::closingNormalGapForce)
				.def("maintainShearGapForce",
					 &BimatInterface::maintainShearGapForce)
				.def("registerDumpFields", [](BimatInterface &self, const std::string &field_names)
					 { self.registerDumpFields(field_names); })
				.def("dump", [](BimatInterface &self, unsigned int step, double time)
					 { self.dump(step, time); })
				.def("advanceTimeStep",
					 &BimatInterface::advanceTimeStep);

			py::class_<UnimatShearInterface,
					   std::shared_ptr<UnimatShearInterface>, Interface>(mod,
																		 "UnimatShearInterface")
				.def(py::init<FFTableMesh &, SpatialDirectionSet, Material &, InterfaceLaw &,
							  const SolverMethod &>(),
					 py::arg("mesh"), py::arg("components"), py::arg("top_material"), py::arg("law"),
					 py::arg("method") = _dynamic)
				.def("init",
					 &UnimatShearInterface::init)
				.def("setTimeStep",
					 &UnimatShearInterface::setTimeStep)
				.def("getStableTimeStep",
					 &UnimatShearInterface::getStableTimeStep,
					 py::return_value_policy::reference)
				.def("initDump", [](UnimatShearInterface &self,
									const std::string &bname,
									const std::string &path)
					 { self.initDump(bname, path); })
				.def("registerDumpFields", [](UnimatShearInterface &self,
											  const std::string &field_names)
					 { self.registerDumpFields(field_names); })
				.def("dump", [](UnimatShearInterface &self, unsigned int step,
								double time)
					 { self.dump(step, time); })
				.def("advanceTimeStep",
					 &UnimatShearInterface::advanceTimeStep)
				.def("getLoad",
					 &UnimatShearInterface::getLoad,
					 py::return_value_policy::reference);

			py::class_<DefRigInterface,
					   std::shared_ptr<DefRigInterface>>(mod,
														 "DefRigInterface")
				.def(py::init<FFTableMesh &, SpatialDirectionSet, Material &, InterfaceLaw &,
							  const SolverMethod &>(),
					 py::arg("mesh"), py::arg("components"), py::arg("top_material"), py::arg("law"),
					 py::arg("method") = _dynamic)
				.def("init",
					 &DefRigInterface::init)
				.def("setTimeStep",
					 &DefRigInterface::setTimeStep)
				.def("getStableTimeStep",
					 &DefRigInterface::getStableTimeStep,
					 py::return_value_policy::reference)
				.def("initDump", [](DefRigInterface &self,
									const std::string &bname,
									const std::string &path)
					 { self.initDump(bname, path); })
				.def("registerDumpFields", [](DefRigInterface &self,
											  const std::string &field_names)
					 { self.registerDumpFields(field_names); })
				.def("dump", [](DefRigInterface &self, unsigned int step,
								double time)
					 { self.dump(step, time); })
				.def("advanceTimeStep",
					 &DefRigInterface::advanceTimeStep)
				.def("getLoad",
					 &DefRigInterface::getLoad,
					 py::return_value_policy::reference);

			py::class_<InterfaceLaw,
					   std::shared_ptr<InterfaceLaw>, PyInterfaceLaw>(mod, "InterfaceLaw")
				.def(py::init<BaseMesh &>(),
					 py::arg("mesh"))
				.def("getCohesion",
					 &InterfaceLaw::getCohesion,
					 py::return_value_policy::reference)
				.def("getInterface",
					 &InterfaceLaw::getInterface)
				.def("getMesh",
					 &InterfaceLaw::getMesh,
					 py::return_value_policy::reference)
				.def("computeCohesiveForces",
					 [](InterfaceLaw &self, bool predicting = false)
					 {
						 self.computeCohesiveForces(predicting);
					 });

			py::class_<BarrasLaw, InterfaceLaw,
					   std::shared_ptr<BarrasLaw>>(mod, "BarrasLaw")
				.def(py::init<BaseMesh &, double, double, std::string &>(),
					 py::arg("mesh"), py::arg("tau_max_default"), py::arg("delta_c_default"),
					 py::arg("name") = "blaw")
				.def("getTauMax",
					 &BarrasLaw::getTauMax, py::return_value_policy::reference)
				.def("getDc",
					 &BarrasLaw::getDc, py::return_value_policy::reference);

			py::class_<LinearCoulombFrictionLaw, InterfaceLaw,
					   std::shared_ptr<LinearCoulombFrictionLaw>>(mod,
																  "LinearCoulombFrictionLaw")
				.def(py::init<BaseMesh &, double, double, double,
							  double, std::string &>(),
					 py::arg("mesh"), py::arg("mu_s_default"),
					 py::arg("mu_k_default"),
					 py::arg("d_c_default"),
					 py::arg("char_reg_time") = 0.,
					 py::arg("name") = "lcflaw")
				.def("getMuS",
					 &LinearCoulombFrictionLaw::getMuS,
					 py::return_value_policy::reference)
				.def("getMuK",
					 &LinearCoulombFrictionLaw::getMuK,
					 py::return_value_policy::reference)
				.def("getDc",
					 &LinearCoulombFrictionLaw::getDc,
					 py::return_value_policy::reference)
				.def("getCharacteristicTime",
					 &LinearCoulombFrictionLaw::getCharacteristicTime,
					 py::return_value_policy::reference);
		}

	}

}
