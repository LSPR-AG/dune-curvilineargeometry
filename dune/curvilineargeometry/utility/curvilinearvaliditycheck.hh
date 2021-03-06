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

	LocalPolynomial  detjacP = curvgeom.JacobianDeterminantAnalytical();
	LocalCoordinateSet samplePoints = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<ctype, mydimension>(sampleOrder);

	ctype detJacMin = 1.0e+20;
	ctype detJacMax = -1.0e+20;

	for (unsigned int i = 0; i < samplePoints.size(); i++)
	{
		ctype detjac   = detjacP.evaluate(samplePoints[i]);
		ctype detjacv2 = curvgeom.integrationElement (samplePoints[i]);
		ctype err = fabs(detjac - detjacv2);

		if (err > 1.0e-8)
		{
			std::cout << "DetJac evaluated analytically = " << detjac << " and the semi-numerical one = " << detjacv2 << " did not match! Error = " << err << std::endl;
			//DUNE_THROW(Dune::IOError, "__ERROR: CheckSelfIntersection test says that detJac evaluation unstable");
		}

		detJacMin = std::min(detJacMin, detjac);
		detJacMax = std::max(detJacMax, detjac);


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

	if (detJacMin / detJacMax < 0.5)  { std::cout << "SelfIntersectionTest: detjac min/max ratio = " << detJacMin / detJacMax << std::endl;
		if (curvgeom.order() == 1)  {
			std::cout << "Unexpected for linear jacobian determinant: " << detjacP.to_string() << std::endl;
			DUNE_THROW(Dune::IOError, "__ERROR: CheckSelfIntersection finds that detJ varies over a linear element");
		}
	}
}


}


}


#endif //DUNE_CURVILINEAR_VALIDITY_CHECK_HH_
