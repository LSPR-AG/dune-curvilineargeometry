# which tests to run
TESTS = \
	test-affinegeometry \
	test-axisalignedcubegeometry \
	test-geometrytype \
	test-multilineargeometry \
    test-polynomial \
	test-nonetype \
	test-quadrature \
	test-referenceelements \
	test-refinement \
	test-curvilinearelementinterpolation \
	test-lagrangegeometry

#\
#
#
#	test-interpolationintegratornumerical \
#


# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

geometrytestdir =  $(includedir)/dune/geometry/test
geometrytest_HEADERS = checkgeometry.hh

# output code coverage
#AM_CXXFLAGS = @AM_CXXFLAGS@ -fprofile-arcs -ftest-coverage

# define the programs
test_affinegeometry_SOURCES = test-affinegeometry.cc

test_axisalignedcubegeometry_SOURCES = test-axisalignedcubegeometry.cc

test_geometrytype_SOURCES = test-geometrytype.cc

test_multilineargeometry_SOURCES = test-multilineargeometry.cc

test_lagrangegeometry_SOURCES = test-lagrangegeometry.cc

test_nonetype_SOURCES = test-nonetype.cc

test_referenceelements_SOURCES = test-referenceelements.cc

test_quadrature_SOURCES = test-quadrature.cc

test_refinement_SOURCES = test-refinement.cc

test_polynomial_SOURCES = test-polynomial.cc

#test_interpolationintegratornumerical_SOURCES = test-interpolationintegratornumerical.cc

test_curvilinearelementinterpolation_SOURCES = test-curvilinearelementinterpolation.cc


include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
