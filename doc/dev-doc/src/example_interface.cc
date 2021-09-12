#include "example_interface.hh"
#include "half_space.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
ExampleInterface::ExampleInterface(FFTableMesh & mesh,
				   Material & top_material,
				   InterfaceLaw & law,
				   const SolverMethod & method) :
  Interface(mesh, law)
{
  this->top = HalfSpace::newHalfSpace(mesh, 1, method);
  
  this->half_spaces.resize(1);
  this->half_spaces[0] = this->top;
  this->top->setMaterial(&top_material);
}

/* -------------------------------------------------------------------------- */
ExampleInterface::~ExampleInterface() {
  delete this->top;
}

/* -------------------------------------------------------------------------- */
void ExampleInterface::closingNormalGapForce(NodalFieldComponent & close_force,
					     bool predicting) {
  // YOUR CODE COMES HERE
  // YOU NEED TO FILL THE CLOSE_FORCE FIELD
}

/* -------------------------------------------------------------------------- */
void ExampleInterface::maintainShearGapForce(NodalField & maintain_force) {
  // YOUR CODE COMES HERE
  // YOU NEED TO FILL THE MAINTAIN_FORCE FIELD
}

/* -------------------------------------------------------------------------- */
void ExampleInterface::computeGap(NodalField & gap,
				  bool predicting) {
  // YOUR CODE COMES HERE
  // YOU NEED TO FILL THE GAP FIELD
}

/* -------------------------------------------------------------------------- */
void ExampleInterface::computeGapVelocity(NodalField & gap_velo,
					  bool predicting) {
  // YOUR CODE COMES HERE
  // YOU NEED TO FILL THE GAP_VELO FIELD
}

/* -------------------------------------------------------------------------- */
void ExampleInterface::registerDumpField(const std::string & field_name) {

  int d = std::atoi(&field_name[field_name.length() - 1]);

  if (d >= this->mesh.getDim())
    throw std::runtime_error("Field "+field_name
			     +" cannot be dumped, to high dimension");
  
  bool registered = false;
  // field_name starts with "top"
  if (field_name.rfind("top", 0) == 0) {
    // cut away "top_" from field_name and give interface as dumper
    registered = this->top->registerDumpFieldToDumper(field_name.substr(4),
						      field_name,
						      this);
  }

  // DO SAME WITH BOT IF NEEDED
  
  if (!registered) {
    // YOUR NEW FIELDS COME HERE
    if (field_name == "example_field" + std::to_string(d))
      this->registerForDump(field_name, this->example_field.component(d));
    else
      Interface::registerDumpField(field_name);
  }
}

__END_UGUCA__
