#include "example_law.hh"
#include "interface.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

ExampleLaw::ExampleLaw(BaseMesh & mesh,
		       double example_parameter_default) :
  InterfaceLaw(mesh),
  example_parameter(mesh)
{
  this->example_parameter.setAllValuesTo(example_parameter_default);
}

/* -------------------------------------------------------------------------- */
void ExampleLaw::computeCohesiveForces(NodalField & cohesion,
				       bool predicting) {
  // YOUR CODE COMES HERE
  // YOU NEED TO FILL THE COHESION NODAL FIELDS
}

/* -------------------------------------------------------------------------- */
void ExampleLaw::registerDumpField(const std::string & field_name) {

  // example parameter
  if (field_name == "example_parameter") {
    this->interface->registerForDump(field_name,
				     this->example_parameter);
  }
  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }

}

__END_UGUCA__
