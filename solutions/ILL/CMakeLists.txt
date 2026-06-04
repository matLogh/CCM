add_executable(solveTimeEvo_ILL solveTimeEvo_ILL.cpp)
target_link_libraries(solveTimeEvo_ILL theuerkauf_fitter ${ROOT_LIBRARIES} ccm)
ccm_configure_target(solveTimeEvo_ILL)
set_target_properties(solveTimeEvo_ILL PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(detectTimeEvo_ILL detectTimeEvo_ILL.cpp)
target_link_libraries(detectTimeEvo_ILL ${ROOT_LIBRARIES} ccm)
ccm_configure_target(detectTimeEvo_ILL)
set_target_properties(detectTimeEvo_ILL PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(simple_ILL simple_ILL.cpp)
target_link_libraries(simple_ILL ${ROOT_LIBRARIES} ccm)
ccm_configure_target(simple_ILL)
set_target_properties(simple_ILL PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(simple_ILL_fw simple_ILL_fw.cpp)
target_link_libraries(simple_ILL_fw ${ROOT_LIBRARIES} ccm)
ccm_configure_target(simple_ILL_fw)
set_target_properties(simple_ILL_fw PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
