cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
   project(ccm)


add_compile_options(-pthread -std=c++17 -g -fmax-errors=5 )

find_package(ROOT REQUIRED COMPONENTS RIO Tree Graf)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})



add_definitions("-DDATA_PATH=\"${PROJECT_SOURCE_DIR}\"")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fitter/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fitter)


add_executable(optimizer optimizer.cpp CCM.cpp Cross_correlation.cpp)
target_link_libraries(optimizer theuerkauf_fitter ${ROOT_LIBRARIES} )

add_executable(simple_example simple_example.cpp CCM.cpp Cross_correlation.cpp)
target_link_libraries(simple_example ${ROOT_LIBRARIES} )

