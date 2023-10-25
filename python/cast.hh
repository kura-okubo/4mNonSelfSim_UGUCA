#ifndef __CAST_HH__
#define __CAST_HH__

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "nodal_field_component.hh"

namespace py = pybind11;

namespace pybind11{
  namespace detail{

    static inline handle policy_switch(return_value_policy policy, handle parent) {
      switch (policy) {
      case return_value_policy::copy:
      case return_value_policy::move:
	return handle();
      case return_value_policy::automatic_reference:  // happens in python-derived
                                                  // classes
      case return_value_policy::reference:
	return none();
      case return_value_policy::reference_internal:
	return parent;
      default:
	throw std::invalid_argument("return-value-policy cannot be handled");
      }
    }


    /**
     * Type caster for grid classes
     * inspired by https://tinyurl.com/y8m47qh3 from T. De Geus
     * and pybind11/eigen.h
     */

  
    /**
     * Type Caster : NodalFieldComponent <-> NumPy-array
     */ 

    template<>
    struct type_caster<uguca::NodalFieldComponent> {
      using type = uguca::NodalFieldComponent;
      template<typename U>
      using array_type = py::array_t<U, array::c_style | array::forcecast>;
      
    public:
      /**
       * This macro establishes the name 'NodalFieldComponentWrap' in
       * function signatures and declares a local variable
       * 'value' of type NodalFieldComponentWrap
       */
      PYBIND11_TYPE_CASTER(type, _("NodalFieldComponentWrap"));

      /// Conversion part 1 ---- (Python -> C++)
      bool load(py::handle src, bool /*convert*/) {

	auto buf = array_type<double>::ensure(src);
	if (!buf)
	  return false;

	auto dims = buf.ndim();
	if (dims != 1)
	  return false;

        value.setField(buf.mutable_data());
	return true;
      }

      
      /// Conversion part 2 --- (C++ -> python)
      static py::handle cast(const type & src,
			     py::return_value_policy policy, py::handle parent) {
	parent = policy_switch(policy, parent);
    
	py::array a(std::move(src.getNbNodes()),
		    src.storage(),
		    parent);
      
	return a.release();
      }
    };
    
 }
}

#endif
