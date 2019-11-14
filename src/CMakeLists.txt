target_sources(${MPAS_EXE_NAME} PUBLIC main.F)

add_subdirectory(external)
add_dependencies(${MPAS_EXE_NAME} ezxml)

add_subdirectory(framework)
