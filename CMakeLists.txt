cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(ccm VERSION 1.0.0 LANGUAGES CXX)

# ============================================================================
# Project Configuration
# ============================================================================

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

# Build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")

# ============================================================================
# Compilation Flags
# ============================================================================

# Max number of diagnostics shown by the compiler
set(MAX_DIAGNOSTICS 3 CACHE STRING "Maximum number of compiler diagnostics to display")

# Common flags for all build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++17")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-fmax-errors=${MAX_DIAGNOSTICS})
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-ferror-limit=${MAX_DIAGNOSTICS})
endif()

# Debug build flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG" CACHE STRING "Debug flags" FORCE)

# Release build flags  
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Release flags" FORCE)

# RelWithDebInfo flags (release with debug info)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG" CACHE STRING "RelWithDebInfo flags" FORCE)

# MinSizeRel flags (minimal size release)
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG" CACHE STRING "MinSizeRel flags" FORCE)

# Enable additional warnings in Debug mode
# Note: -Werror is NOT added here globally to avoid breaking on external library warnings
# Apply -Werror selectively only to code targets
if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options(-Wall -Wextra -Werror -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 -Wformat-security -fstack-protector-strong)
endif()

# ============================================================================
# User Options
# ============================================================================

option(BUILD_EXAMPLES "Build example applications" ON)
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)

# Helper function to apply strict warnings to CCM code
# Usage: add_strict_warnings(my_target)
function(add_strict_warnings TARGET)
    target_compile_options(${TARGET} PRIVATE 
        -Werror 
        -Wall 
        -Wextra 
        -Wpedantic 
        -Wshadow 
        -Wconversion 
        -Wsign-conversion 
        -Wnull-dereference 
        -Wdouble-promotion 
        -Wformat=2 
        -Wformat-security 
        -fstack-protector-strong
    )
endfunction()

# AddressSanitizer support
if(ENABLE_ASAN)
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer -g)
    add_link_options(-fsanitize=address)
    message(STATUS "AddressSanitizer enabled")
endif()

# ============================================================================
# Report Configuration
# ============================================================================

message(STATUS "CCM Project Configuration:")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  CXX Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  Max Diagnostics: ${MAX_DIAGNOSTICS}")
message(STATUS "  Build Examples: ${BUILD_EXAMPLES}")
message(STATUS "  AddressSanitizer: ${ENABLE_ASAN}")

find_package(ROOT REQUIRED COMPONENTS RIO Tree Graf MathMore)

include(${ROOT_USE_FILE})
# Treat ROOT headers as system includes to suppress warnings from external library
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

link_libraries(stdc++fs)



add_definitions("-DDATA_PATH=\"${PROJECT_SOURCE_DIR}\"")

# Update git submodules
execute_process(
    COMMAND git submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "Failed to update git submodules. Please check your repository.")
endif()



add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/modules/fitter/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/modules/fitter)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/)

# ============================================================================
# Example Applications
# ============================================================================

if(BUILD_EXAMPLES)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/solutions)
    message(STATUS "Building example applications")
else()
    message(STATUS "Example applications disabled (use -DBUILD_EXAMPLES=ON to enable)")
endif()
