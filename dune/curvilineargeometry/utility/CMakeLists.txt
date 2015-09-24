#install headers
set(HEADERS
  curvilinearvaliditycheck.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargeometry/utility)
