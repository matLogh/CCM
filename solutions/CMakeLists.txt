if(BUILD_SOLUTION_AGATA)
    add_subdirectory(AGATA)
endif()

if(BUILD_SOLUTION_ILL)
    add_subdirectory(ILL)
endif()

if(BUILD_SOLUTION_PRISMA)
    add_subdirectory(PRISMA)
endif()

if(BUILD_SOLUTION_OSCAR)
    add_subdirectory(OSCAR)
endif()

if(BUILD_SOLUTION_SIMPLE_EXAMPLE)
    add_subdirectory(simple_example)
endif()
