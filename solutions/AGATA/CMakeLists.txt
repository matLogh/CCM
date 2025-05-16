

add_executable(solveTimeEvo solveTimeEvo_AGATA.cpp)
target_link_libraries(solveTimeEvo theuerkauf_fitter ${ROOT_LIBRARIES} ccm)
set_target_properties(solveTimeEvo PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(matTimeEvo_cores matTimeEvo_AGATA.cpp)
target_link_libraries(matTimeEvo_cores  ${ROOT_LIBRARIES})
set_target_properties(matTimeEvo_cores PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(detectTimeEvo detectTimeEvo_AGATA.cpp)
target_link_libraries(detectTimeEvo  ${ROOT_LIBRARIES} ccm)
set_target_properties(detectTimeEvo PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(simple_AGATA simple_agata.cpp)
target_link_libraries(simple_AGATA  ${ROOT_LIBRARIES} ccm)
set_target_properties(simple_AGATA PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})