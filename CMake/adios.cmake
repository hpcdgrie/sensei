if(ENABLE_ADIOS)
  find_package(ADIOS REQUIRED)
  add_library(adios INTERFACE)
  target_link_libraries(adios INTERFACE ${ADIOS_LIBRARIES})
  target_include_directories(adios SYSTEM INTERFACE ${ADIOS_INCLUDE_DIRS})
  if(ADIOS_DEFINITIONS)
    target_compile_definitions(adios INTERFACE ${ADIOS_DEFINITIONS})
  endif()
  install(TARGETS adios EXPORT adios)
  install(EXPORT adios DESTINATION lib/cmake EXPORT_LINK_INTERFACE_LIBRARIES)
endif()