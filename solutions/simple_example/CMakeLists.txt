

add_executable(simple_example simple_example.cpp)
target_link_libraries(simple_example ${ROOT_LIBRARIES} ccm)
ccm_configure_target(simple_example)
set_target_properties(simple_example PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
ccm_install_executable(simple_example)
