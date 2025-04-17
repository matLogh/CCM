

add_executable(calcTimeEvo calcTimeEvo_AGATA.cpp)
target_link_libraries(calcTimeEvo theuerkauf_fitter ${ROOT_LIBRARIES} ccm)
set_target_properties(calcTimeEvo PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(matTimeEvo_cores matTimeEvo_AGATA.cpp)
target_link_libraries(matTimeEvo_cores  ${ROOT_LIBRARIES})
set_target_properties(matTimeEvo_cores PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

