add_executable(solveTimeEvo_ILL solveTimeEvo_ILL.cpp)
target_link_libraries(solveTimeEvo_ILL theuerkauf_fitter ${ROOT_LIBRARIES} ccm)
set_target_properties(solveTimeEvo_ILL PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(detectTimeEvo_ILL detectTimeEvo_ILL.cpp)
target_link_libraries(detectTimeEvo_ILL ${ROOT_LIBRARIES} ccm)
set_target_properties(detectTimeEvo_ILL PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
