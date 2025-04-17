cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
   project(ccm)


# add_compile_options(-pthread -std=c++17 -g -Werror -fmax-errors=3 -O3 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 -Wformat-security -fstack-protector-strong)
add_compile_options(-pthread -std=c++17 -g -Werror -fmax-errors=3 -O3)
# add_compile_options(-fsanitize=address -fno-omit-frame-pointer -g)
# add_link_options(-fsanitize=address)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS RIO Tree Graf MathMore)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})



add_definitions("-DDATA_PATH=\"${PROJECT_SOURCE_DIR}\"")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/modules/fitter/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/modules/fitter)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/)

# add_executable(optimizer optimizer.cpp CCM.cpp Cross_correlation.cpp)
# target_link_libraries(optimizer theuerkauf_fitter ${ROOT_LIBRARIES} )

# add_executable(simple_example simple_example.cpp CCM.cpp Cross_correlation.cpp)
# target_link_libraries(simple_example ${ROOT_LIBRARIES} )

# add_executable(agata test_agata.cpp)
# target_link_libraries(agata theuerkauf_fitter ${ROOT_LIBRARIES} ccm)

# add_executable(test test_interpolator.cpp)
# target_link_libraries(test ${ROOT_LIBRARIES} ccm)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/solutions)


# add_executable(oscar_fix oscar_fix.cpp)
# target_link_libraries(oscar_fix theuerkauf_fitter ${ROOT_LIBRARIES} ccm)

# add_executable(prisma_abberations prisma_abberations.cpp)
# target_link_libraries(prisma_abberations ${ROOT_LIBRARIES} ccm)