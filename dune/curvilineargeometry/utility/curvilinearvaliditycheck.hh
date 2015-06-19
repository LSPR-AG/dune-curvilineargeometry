#ifndef DUNE_CURVILINEAR_VALIDITY_CHECK_HH_
#define DUNE_CURVILINEAR_VALIDITY_CHECK_HH_

#include <iostream>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

namespace Dune
{

namespace CurvilinearValidityCheck
{

// [TODO] Move to standard constant file
const int     CURVGRID_DEFAULT_SELF_INTERSECTION_SAMPLE_ORDER = 7;
const double  CURVGRID_DEFAULT_DETJ_SMALLEST_VALUE            = 1.0e-10;


template <class CurvGeom>
bool SelfIntersection(
	const CurvGeom & curvgeom,
	int sampleOrder = CURVGRID_DEFAULT_SELF_INTERSECTION_SAMPLE_ORDER,
	double tolerance = CURVGRID_DEFAULT_DETJ_SMALLEST_VALUE)
{
	typedef typename CurvGeom::ctype  ctype;
	const int coorddimension = CurvGeom::coorddimension;
	const int mydimension    = CurvGeom::mydimension;

	typedef typename CurvGeom::LocalPolynomial     LocalPolynomial;
	typedef typename CurvGeom::LocalCoordinate     LocalCoordinate;
	typedef typename std::vector<LocalCoordinate>  LocalCoordinateSet;

	LocalCoordinateSet samplePoints = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<ctype, mydimension>(sampleOrder);

	for (unsigned int i = 0; i < samplePoints.size(); i++)
	{
		LocalPolynomial  detjacP = Dune::DifferentialHelper::JacobianDeterminantAnalytical<CurvGeom, coorddimension, mydimension>::eval(curvgeom);
		ctype            detjac = detjacP.evaluate(samplePoints[i]);


		if (fabs(detjac) < tolerance )
		{
			std::cout << "Warning: JacobianDeterminant has evaluated to a small value " << detjac << " at sample point " << samplePoints[i] << std::endl;
		} else if (detjac < 0)
		{
			std::cout << "Jacobian Determinant: " << detjacP.to_string() << std::endl;
			std::cout << "Has negative value " << detjac << " at sample point " << samplePoints[i] << std::endl;

			DUNE_THROW(Dune::IOError, "__ERROR: CheckSelfIntersection test suspects that the tested entity is self-intersecting");
		}
	}
}


}


}


#endif //DUNE_CURVILINEAR_VALIDITY_CHECK_HH_
