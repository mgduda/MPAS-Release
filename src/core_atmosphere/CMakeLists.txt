add_custom_command(OUTPUT
                      block_dimension_routines.inc
                      core_variables.inc
                      define_packages.inc
                      domain_variables.inc
                      namelist_call.inc
                      namelist_defines.inc
                      setup_immutable_streams.inc
                      structs_and_variables.inc
                   COMMAND
                      ${CMAKE_BINARY_DIR}/src/tools/registry/parse ${CMAKE_SOURCE_DIR}/src/core_atmosphere/Registry_temp.xml
                   DEPENDS
                      parse Registry_temp.xml)

add_custom_target(atm_core
                  DEPENDS
                     block_dimension_routines.inc
                     core_variables.inc
                     define_packages.inc
                     domain_variables.inc
                     namelist_call.inc
                     namelist_defines.inc
                     setup_immutable_streams.inc
                     structs_and_variables.inc)

add_dependencies(${MPAS_EXE_NAME} atm_core)

set_property(GLOBAL APPEND PROPERTY MPAS_INTERNAL_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR})

target_sources(${MPAS_EXE_NAME} PUBLIC atm_core.F)
