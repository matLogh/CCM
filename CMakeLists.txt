cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(ccm VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

# ============================================================================
# Project Configuration
# ============================================================================

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif()

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
option(BUILD_SOLUTION_AGATA "Build AGATA solution applications" ON)
option(BUILD_SOLUTION_ILL "Build ILL solution applications" ON)
option(BUILD_SOLUTION_PRISMA "Build PRISMA solution applications" ON)
option(BUILD_SOLUTION_OSCAR "Build OSCAR solution applications" ON)
option(BUILD_SOLUTION_SIMPLE_EXAMPLE "Build simple example solution application" ON)
option(INSTALL_EXAMPLES "Install example and solution applications" ON)
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
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")
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
if(BUILD_EXAMPLES)
    message(STATUS "  Install Examples: ${INSTALL_EXAMPLES}")
    message(STATUS "  Build Solution AGATA: ${BUILD_SOLUTION_AGATA}")
    message(STATUS "  Build Solution ILL: ${BUILD_SOLUTION_ILL}")
    message(STATUS "  Build Solution PRISMA: ${BUILD_SOLUTION_PRISMA}")
    message(STATUS "  Build Solution OSCAR: ${BUILD_SOLUTION_OSCAR}")
    message(STATUS "  Build Solution Simple Example: ${BUILD_SOLUTION_SIMPLE_EXAMPLE}")
endif()
message(STATUS "  AddressSanitizer: ${ENABLE_ASAN}")

find_package(Threads REQUIRED)
find_package(ROOT REQUIRED COMPONENTS RIO Tree Graf MathMore)

include(${ROOT_USE_FILE})
# Treat ROOT headers as system includes to suppress warnings from external library
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CCM_PLATFORM_LIBS Threads::Threads)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    list(APPEND CCM_PLATFORM_LIBS stdc++fs)
endif()

function(ccm_configure_target TARGET)
    target_compile_features(${TARGET} PUBLIC cxx_std_17)
    target_link_libraries(${TARGET} ${CCM_PLATFORM_LIBS})
    if(APPLE)
        set_target_properties(${TARGET} PROPERTIES
            MACOSX_RPATH ON
            BUILD_RPATH "${ROOT_LIBRARY_DIR};${ROOT_LIBRARY_DIRS}"
        )
    endif()
endfunction()

function(ccm_install_library TARGET)
    install(TARGETS ${TARGET}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endfunction()

function(ccm_install_executable TARGET)
    if(INSTALL_EXAMPLES)
        install(TARGETS ${TARGET}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
    endif()
endfunction()



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
