


add_library(theuerkauf_fitter TheuerkaufPeak.hpp Util.hpp)
set_target_properties(theuerkauf_fitter PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(theuerkauf_fitter ${ROOT_LIBRARIES} )

add_executable(fitter_test test.cpp)
# target_include_directories(fitter_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fitter_test theuerkauf_fitter ${ROOT_LIBRARIES})
set_target_properties(fitter_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)