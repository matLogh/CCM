cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
   project(ccm)


add_compile_options(-pthread -std=c++14 -g -fmax-errors=5 )

find_package(ROOT REQUIRED COMPONENTS RIO Tree Graf)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})



# add_library(ccmlib SHARED  CheckCCM.h CCM.h Cross_correlation.h)
# target_link_libraries(ccmlib ${ROOT_LIBRARIES} )
message("   Project dir is : ${PROJECT_SOURCE_DIR}")


add_definitions("-DDATA_PATH=\"${PROJECT_SOURCE_DIR}\"")

add_executable(ccm main_example.cpp CheckCCM.h CCM.h Cross_correlation.h CheckCCM.cpp CCM.cpp Cross_correlation.cpp)
target_link_libraries(ccm ${ROOT_LIBRARIES} )

