if(ENABLE_VISTLE)
find_package(Vistle REQUIRED COMPONENTS sensei sensei_vtk boost_mpi
)
  #find_package(vistle_sensei REQUIRED)
  #find_package(vistle_sensei_vtk REQUIRED)
  #find_package(vistle_boost_mpi REQUIRED)
  add_library(sVistle INTERFACE)
  target_link_libraries(sVistle
    INTERFACE
      Vistle::vistle_sensei
      Vistle::vistle_sensei_vtk
  )
  get_target_property(vistle_boost_mpi_includes Vistle::vistle_boost_mpi INTERFACE_INCLUDE_DIRECTORIES) 
  target_include_directories(sVistle INTERFACE ${vistle_boost_mpi_includes})      
  #get_target_property(vistle_sensei_includes Vistle::vistle_sensei INTERFACE_INCLUDE_DIRECTORIES)
  #target_include_directories(sVistle SYSTEM INTERFACE ${vistle_sensei_includes})
  install(TARGETS sVistle EXPORT sVistle)
  install(EXPORT sVistle DESTINATION lib/cmake EXPORT_LINK_INTERFACE_LIBRARIES)
endif()
