pybind11_add_module(puguca MODULE
  uguca_module.cpp
  mesh.cpp
  half_space.cpp
  interface.cpp
)

target_link_libraries(puguca PRIVATE uguca)
