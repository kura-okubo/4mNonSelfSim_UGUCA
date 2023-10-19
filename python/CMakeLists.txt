pybind11_add_module(puguca MODULE
  uguca_module.cpp
)

target_link_libraries(puguca PRIVATE uguca)
