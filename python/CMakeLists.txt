set(PUGUCA_SRCS
  uguca_module.cpp
  mesh.cpp
  half_space.cpp
  interface.cpp
  io.cpp
)
  
pybind11_add_module(puguca ${PUGUCA_SRCS})

target_link_libraries(puguca PUBLIC uguca)

set_target_properties(puguca PROPERTIES
  DEBUG_POSTFIX ""
  LIBRARY_OUTPUT_DIRECTORY uguca)

file(COPY uguca DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY pyproject.toml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(TARGET puguca POST_BUILD
    COMMAND cd ${CMAKE_CURRENT_BINARY_DIR}/uguca && pybind11-stubgen --output-dir . puguca && cd ${CMAKE_CURRENT_BINARY_DIR})
