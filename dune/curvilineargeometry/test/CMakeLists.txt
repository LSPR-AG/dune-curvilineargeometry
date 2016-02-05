set(TESTS
  test_polynomial
  test_lagrangeinterpolation
  test_adaptiveintegration
  test_quadratureintegration
  test_quadratureintegration-matrix
  test_quadratureintegration-speed
  test_curvilineargeometry)

add_executable(test_polynomial test-polynomial.cc)
add_executable(test_polynomialvector test-polynomialvector.cc)
add_executable(test_lagrangeinterpolation test-lagrangeinterpolation.cc)
add_executable(test_adaptiveintegration test-adaptiveintegration.cc)
add_executable(test_quadratureintegration test-quadratureintegration.cc)
add_executable(test_quadratureintegration-matrix test-quadratureintegration-matrix.cc)
add_executable(test_quadratureintegration-speed test-quadratureintegration-speed.cc)
add_executable(test_curvilineargeometry test-curvilineargeometry.cc)


target_link_libraries(test_lagrangeinterpolation dunecommon dunegeometry)
target_link_libraries(test_adaptiveintegration dunecommon dunegeometry)
target_link_libraries(test_quadratureintegration dunecommon dunegeometry)
target_link_libraries(test_quadratureintegration-matrix dunecommon dunegeometry)
target_link_libraries(test_quadratureintegration-speed dunecommon dunegeometry)
target_link_libraries(test_polynomial dunecommon dunegeometry)
target_link_libraries(test_polynomialvector dunecommon dunegeometry)
target_link_libraries(test_curvilineargeometry dunecommon dunegeometry)

# We do not want want to build the tests during make all,
# but just build them on demand
## add_directory_test_target(_test_target)
## add_dependencies(${_test_target} ${TESTS})

##foreach(_test ${TESTS})
##  add_test(${_test} ${_test})
##endforeach(_test ${TESTS})
