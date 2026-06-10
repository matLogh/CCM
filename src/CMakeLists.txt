
add_library(ccm CCM.cpp Cross_correlation.cpp CCMInterpolator.cpp RegionOfInterest.cpp)
set_target_properties(ccm PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ccm ${ROOT_LIBRARIES} )
ccm_configure_target(ccm)
ccm_install_library(ccm)

install(FILES
    CCM.h
    CCMInterpolator.h
    Cross_correlation.h
    RegionOfInterest.h
    variables.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ccm
)
