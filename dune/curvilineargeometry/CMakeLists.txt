add_subdirectory(integration)
add_subdirectory(interpolation)
add_subdirectory(test)
add_subdirectory(utility)

#install headers
set(HEADERS
  curvilineargeometry.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargeometry)
