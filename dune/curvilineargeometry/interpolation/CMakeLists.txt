#install headers
set(HEADERS
  polynomial.hh
  curvilinearelementinterpolator.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargeometry/interpolation)
