#install headers
set(HEADERS
  numericalrecursiveinterpolationintegrator.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargeometry/integration)
