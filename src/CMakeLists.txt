target_sources(${MPAS_EXE_NAME} PUBLIC main.F)

add_subdirectory(external)
add_subdirectory(tools)
add_subdirectory(framework)
add_subdirectory(core_atmosphere)
