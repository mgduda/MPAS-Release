add_executable(parse)

target_sources(parse
               PUBLIC
                  dictionary.c
                  fortprintf.c
                  gen_inc.c
                  parse.c
                  utility.c)

#
# Add include paths and libraries from packages that are internally built by MPAS
# In particular, the ezxml library is needed by the parse tool
#
get_property(MPAS_INTERNAL_INCLUDE_DIRS GLOBAL PROPERTY MPAS_INTERNAL_INCLUDE_DIRS)
get_property(MPAS_INTERNAL_LIBRARIES GLOBAL PROPERTY MPAS_INTERNAL_LIBRARIES)
target_include_directories(parse PRIVATE ${MPAS_INTERNAL_INCLUDE_DIRS})
target_link_libraries(parse PRIVATE ${MPAS_INTERNAL_LIBRARIES})
