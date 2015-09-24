#install headers
set(HEADERS
  adaptiveintegrator.hh
  integrationhelper.hh
  quadratureintegrator.hh
  quadraturerelativeerror.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargeometry/integration)
