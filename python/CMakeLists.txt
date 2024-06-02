set(PUGUCA_SRCS
  src/uguca.cpp
  src/half_space.cpp
  src/interface.cpp
  src/io.cpp
  src/mesh.cpp
)
  
pybind11_add_module(puguca SHARED ${PUGUCA_SRCS})

target_link_libraries(puguca PUBLIC uguca)

set_target_properties(puguca PROPERTIES
  DEBUG_POSTFIX ""
  LIBRARY_OUTPUT_DIRECTORY uguca)

file(COPY uguca DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY pyproject.toml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
