if(ENABLE_VISTLE)
  find_package(VISTLE REQUIRED)
  add_library(sVistle INTERFACE)
  target_link_libraries(sLibsim INTERFACE ${VISTLE_LIBRARY})
  target_include_directories(sLibsim SYSTEM INTERFACE ${VISTLE_INCLUDE_DIR})
  install(TARGETS sVistle EXPORT sVistle)
  install(EXPORT sVistle DESTINATION lib/cmake EXPORT_LINK_INTERFACE_LIBRARIES)
endif()


