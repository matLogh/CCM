
add_library(ccm CCM.cpp Cross_correlation.cpp CCMInterpolator.cpp RegionOfInterest.cpp)
set_target_properties(ccm PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ccm ${ROOT_LIBRARIES} )
