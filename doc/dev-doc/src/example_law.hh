#ifndef __EXAMPLE_LAW_H__
#define __EXAMPLE_LAW_H__

/* -------------------------------------------------------------------------- */
#include "interface_law.hh"

__BEGIN_UGUCA__

class ExampleLaw : public InterfaceLaw {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ExampleLaw(BaseMesh & mesh,
	     double example_parameter_default);

  virtual ~ExampleLaw() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void computeCohesiveForces(NodalField & cohesion,
			     bool predicting = false);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

 /* ------------------------------------------------------------------------ */
 /* Accessors                                                                */
 /* ------------------------------------------------------------------------ */
public:
  NodalFieldComponent & getExampleParameter() { return this->example_parameter; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  NodalFieldComponent example_parameter;
};

__END_UGUCA__

#endif /* __EXAMPLE_LAW_H__ */
