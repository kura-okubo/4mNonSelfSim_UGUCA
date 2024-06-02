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

# Generate static python stubs using pybind11-stubgen (python package)
add_custom_command(TARGET puguca POST_BUILD COMMAND 
    cd ${CMAKE_CURRENT_BINARY_DIR}/uguca && pybind11-stubgen --output-dir . puguca && cd ${CMAKE_CURRENT_BINARY_DIR}
)
