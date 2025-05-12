

add_executable(simple_example simple_example.cpp)
target_link_libraries(simple_example ${ROOT_LIBRARIES} ccm)
set_target_properties(simple_example PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

