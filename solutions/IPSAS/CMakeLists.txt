

add_executable(simple_ipsas simple_ipsas.cpp)
target_link_libraries(simple_ipsas ${ROOT_LIBRARIES} ccm)
ccm_configure_target(simple_ipsas)
set_target_properties(simple_ipsas PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
ccm_install_executable(simple_ipsas)
