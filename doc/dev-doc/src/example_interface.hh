#ifndef __EXAMPLE_INTERFACE_H__
#define __EXAMPLE_INTERFACE_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "half_space.hh"
#include "interface_law.hh"
#include "interface.hh"

__BEGIN_UGUCA__

class ExampleInterface : public Interface {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ExampleInterface(FFTableMesh & mesh,
		   Material & top_material,
		   InterfaceLaw & law,
		   const SolverMethod & method = _dynamic);

  virtual ~ExampleInterface();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // compute force needed to close normal gap
  virtual void closingNormalGapForce(NodalFieldComponent & close_force,
				     bool predicting = false);

  // compute force needed to maintain current shear gap
  virtual void maintainShearGapForce(NodalField & maintain_force);

  // compute gap in displacement
  virtual void computeGap(NodalField & gap,
			  bool predicting = false);

  // compute gap relative velocity
  virtual void computeGapVelocity(NodalField & gap_velo,
				  bool predicting = false);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  virtual HalfSpace & getTop() { return *(this->top); }
  virtual HalfSpace & getBot() {
    throw "ExampleInterface::getBot() is not implemented.";
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // half spaces
  HalfSpace * top;
  
  // YOU CAN ADD NEW FIELDS HERE
  NodalField & example_field;
};

__END_UGUCA__

//#include "example_interface_impl.cc"

#endif /* __EXAMPLE_INTERFACE_H__ */
