

add_executable(prisma_abberations prisma_abberations.cpp)
target_link_libraries(prisma_abberations  ${ROOT_LIBRARIES} ccm)
set_target_properties(prisma_abberations PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


