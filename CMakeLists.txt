cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
   project(ccm)


add_compile_options(-pthread -std=c++17 -g -fmax-errors=5 -O3)

find_package(ROOT REQUIRED COMPONENTS RIO Tree Graf MathMore)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})



add_definitions("-DDATA_PATH=\"${PROJECT_SOURCE_DIR}\"")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/modules/fitter/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/modules/fitter)

add_library(ccm CCM.cpp Cross_correlation.cpp CCMInterpolator.cpp)
set_target_properties(ccm PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ccm ${ROOT_LIBRARIES} )

# add_executable(optimizer optimizer.cpp CCM.cpp Cross_correlation.cpp)
# target_link_libraries(optimizer theuerkauf_fitter ${ROOT_LIBRARIES} )

# add_executable(simple_example simple_example.cpp CCM.cpp Cross_correlation.cpp)
# target_link_libraries(simple_example ${ROOT_LIBRARIES} )

add_executable(agata test_agata.cpp)
target_link_libraries(agata theuerkauf_fitter ${ROOT_LIBRARIES} ccm)


# add_executable(prisma_abberations prisma_abberations.cpp)
# target_link_libraries(prisma_abberations ${ROOT_LIBRARIES} ccm)