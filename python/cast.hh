#ifndef __CAST_HH__
#define __CAST_HH__

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>

#include "nodal_field.hh"

namespace py = pybind11;

namespace pybind11
{
  namespace detail
  {

    static inline handle policy_switch(return_value_policy policy, handle parent)
    {
      switch (policy)
      {
      case return_value_policy::copy:
      case return_value_policy::move:
        return handle();
      case return_value_policy::automatic_reference: // happens in python-derived classes
      case return_value_policy::reference:
        return none();
      case return_value_policy::reference_internal:
        return parent;
      default:
        throw std::invalid_argument("return-value-policy cannot be handled");
      }
    }

    /**
     * Type caster for nodal field
     * inspired by https://tinyurl.com/y8m47qh3 from T. De Geus
     * and pybind11/eigen.h
     */

    /**
     * Type Caster : NodalField <-> NumPy-array
     */

   /* template <>
    struct type_caster<uguca::NodalField>
    {
      using type = uguca::NodalField;
      template <typename U>
      using array_type = py::array_t<U, array::c_style | array::forcecast>;

    public:
      PYBIND11_TYPE_CASTER(type, _("NodalFieldWrap"));

      /// Conversion part 1 ---- (Python -> C++)
      bool load(py::handle src, bool convert)
      {

        auto buf = array_type<double>::ensure(src);
        if (!buf)
          return false;

        auto dims = buf.ndim();
        std::vector<size_t> shape(dims);

        for (unsigned int i = 0; i < dims; ++i)
        {
          shape[i] = buf.shape()[i];
        }

        // value = type(shape[0], buf.data());

        // if (dims != 1)
        //   return false;

        // value.setField(buf.mutable_data());
        return true;
      }

      /// Conversion part 2 --- (C++ -> python)
      static py::handle cast(const type &field,
                             py::return_value_policy policy, py::handle parent)
      {
        parent = policy_switch(policy, parent);

        auto dims = field.getDataSize().size();

        std::vector<size_t> shapes(dims);
        std::vector<size_t> strides(dims);

        switch (dims)
        {
        case 1:
          shapes[0] = field.getDataSize()[0];
          strides[0] = 1. * sizeof(double);
          break;
        case 2:
          shapes[1] = field.getDataSize()[0];
          shapes[0] = 2;

          strides[0] = shapes[1] * sizeof(double);
          strides[1] = 1. * sizeof(double);
          break;
        default:
          std::cout << "Not implemented yet" << std::endl;
          break;
        }

        py::array a(std::move(shapes),
                    std::move(strides),
                    field.getInternalData(),
                    parent);

        return a.release();
      }
    };
    */

  }
}

#endif
