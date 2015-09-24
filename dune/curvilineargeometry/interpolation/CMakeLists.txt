#install headers
set(HEADERS
  polynomial.hh
  curvilineargeometryhelper.hh
  lagrangeinterpolator.hh
  differentialhelper.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargeometry/interpolation)
