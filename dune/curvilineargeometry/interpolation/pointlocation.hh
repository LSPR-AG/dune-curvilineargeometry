#ifndef DUNE_CURVILINEAR_POINT_LOCATION_HH_
#define DUNE_CURVILINEAR_POINT_LOCATION_HH_


#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>


namespace Dune
{

namespace CurvilinearPointLocation
{

// [TODO] Move to common/curvconst.hh
const int DIM0D = 0;
const int DIM1D = 1;
const int DIM2D = 2;
const int DIM3D = 3;



/*********************************************************************************************************************/
/** \brief Checks if a global coordinate is inside the element or not (!imperfect!)                                 **/
/*  The current algorithm for Simplices computes the global barycentric coordinates and checks if they sum up to 1  **/
/*  For linear elements, normalized sum > 1 implies that element is outside and sum = 1 implies that it is inside   **/
/*  For nonlinear nonconvex elements, normalized sum > 1 does not imply anything, because the total barycentric     **/
/*  area may be larger than the area of the element even with the sample point inside.                              **/
/*********************************************************************************************************************/
template <class CurvGeom, int cdim, int mydim>
struct BarycentricTest
{
	typedef typename CurvGeom::ctype              ctype;
	typedef typename CurvGeom::GlobalCoordinate   GlobalCoordinate;

	static bool isInside(const CurvGeom & curvgeom, const GlobalCoordinate &globalC, ctype tolerance)
	{
		// [FIXME] Throw dune-error if this generic method is called
		return false;
	}

};


// Checking if a point is inside an edge by simply verifying it lies between its two corners. Only for edges in 1D
template <class CurvGeom>
struct BarycentricTest<CurvGeom, DIM1D, DIM1D>
{
	typedef typename CurvGeom::ctype              ctype;
	typedef typename CurvGeom::GlobalCoordinate   GlobalCoordinate;

	static bool isInside(const CurvGeom & curvgeom, const GlobalCoordinate &globalC, ctype tolerance)
	{
		return ((curvgeom.corner(0)[0] - globalC[0] <= tolerance) && (globalC[0]- curvgeom.corner(1)[0] <= tolerance));
	}
};


// Checking if a point is inside a triangle by calculating global simplex coordinates. Only for triangles in 2D
template <class CurvGeom>
struct BarycentricTest<CurvGeom, DIM2D, DIM2D>
{
	typedef typename CurvGeom::ctype              ctype;
	typedef typename CurvGeom::GlobalCoordinate   GlobalCoordinate;

    typedef Dune::CurvilinearElementInterpolator< ctype, DIM1D, DIM2D >  SubentityInterpolator;
    typedef std::vector< SubentityInterpolator >                         SubentityInterpolatorVector;
    typedef Dune::Polynomial<ctype, DIM1D>                               SubentityPolynomial;
    typedef std::vector<SubentityPolynomial>                             SubentityPolynomialVector;

	static bool isInside(const CurvGeom & curvgeom, const GlobalCoordinate &globalC, ctype tolerance)
	{
	    SubentityInterpolatorVector edgeInterpolatorSet;
	    edgeInterpolatorSet.push_back(curvgeom.interpolator().SubentityInterpolator<DIM1D>(0));
	    edgeInterpolatorSet.push_back(curvgeom.interpolator().SubentityInterpolator<DIM1D>(1));
	    edgeInterpolatorSet.push_back(curvgeom.interpolator().SubentityInterpolator<DIM1D>(2));

	    ctype triangleArea   = volume(curvgeom.tolerance);
	    ctype barycentricSum = 0.0;

	    //std::cout << "    # total area " << tri_area << std::endl;

	    for (int i = 0; i < edgeInterpolatorSet.size(); i++)
	    {
	        // Obtains analytical vector function describing the global coordinate of the edge
	        SubentityPolynomialVector pv = edgeInterpolatorSet[i].interpolatoryVectorAnalytical();

	        // Calculates the area of barycentric triangle created by the corners of the curved edge and the point globalC
	        // In this calculation the only the curved edge is curved, the other two edges of this barycentric triangle are straight-sided
	        // We integrate over the generalized surface area, which is the cross product between the sweeping point and its derivative
	        ctype barycentricArea = 0.5 * ((pv[0] - globalC[0]) * pv[1].derivative(0) - (pv[1] - globalC[1]) * pv[0].derivative(0)).integrateRefSimplex();

	        //std::cout << "    # barycentric area " << barycentric_area << std::endl;

	        barycentricSum += barycentricArea;
	    }
	    return (barycentricSum / triangleArea - 1 < tolerance);
	}
};


// Checking if a point is inside a triangle by calculating global simplex coordinates. Only for triangles in 3D
template <class CurvGeom>
struct BarycentricTest<CurvGeom, DIM3D, DIM3D>
{
	typedef typename CurvGeom::ctype              ctype;
	typedef typename CurvGeom::GlobalCoordinate   GlobalCoordinate;

    typedef Dune::CurvilinearElementInterpolator< ctype, DIM2D, DIM3D >  SubentityInterpolator;
    typedef std::vector< SubentityInterpolator >                         SubentityInterpolatorVector;
    typedef Dune::Polynomial<ctype, DIM2D>                               SubentityPolynomial;
    typedef std::vector<SubentityPolynomial>                             SubentityPolynomialVector;

	static bool isInside(const CurvGeom & curvgeom, const GlobalCoordinate &globalC, ctype tolerance)
	{
	    SubentityInterpolatorVector faceInterpolatorSet;
	    faceInterpolatorSet.push_back(curvgeom.interpolator().SubentityInterpolator<DIM2D>(0));
	    faceInterpolatorSet.push_back(curvgeom.interpolator().SubentityInterpolator<DIM2D>(1));
	    faceInterpolatorSet.push_back(curvgeom.interpolator().SubentityInterpolator<DIM2D>(2));
	    faceInterpolatorSet.push_back(curvgeom.interpolator().SubentityInterpolator<DIM2D>(3));

	    ctype tetrahedronVolume = volume(tolerance);
	    ctype barycentricSum    = 0.0;

	    //std::cout << "    # total volume " << tet_vol << std::endl;

	    for (int i = 0; i < faceInterpolatorSet.size(); i++)
	    {
	        // Obtains analytical vector function describing the global coordinate of the face
	        SubentityPolynomialVector pv = faceInterpolatorSet[i].interpolatoryVectorAnalytical();

	        //for (int piter = 0; piter < 3; piter++) { pv[piter].print(); }

	        // Calculates the volume of barycentric tetrahedron created by the corners of the curved face and the point globalC
	        // In this calculation the only the curved face is curved, the other three faces of this barycentric tetrahedron are straight-sided
	        // We integrate over the generalized volume, which is the dot-cross product between the sweeping point and its two derivatives
	        ctype barycentricVolume = 0;
	        barycentricVolume += ((pv[0] - globalC[0]) * (pv[1].derivative(0) * pv[2].derivative(1) - pv[1].derivative(1) * pv[2].derivative(0))).integrateRefSimplex();
	        barycentricVolume += ((pv[1] - globalC[1]) * (pv[2].derivative(0) * pv[0].derivative(1) - pv[2].derivative(1) * pv[0].derivative(0))).integrateRefSimplex();
	        barycentricVolume += ((pv[2] - globalC[2]) * (pv[0].derivative(0) * pv[1].derivative(1) - pv[0].derivative(1) * pv[1].derivative(0))).integrateRefSimplex();


	        //std::cout << "    # barycentric volume " << barycentric_volume << std::endl;

	        barycentricSum += -barycentricVolume / 3.0;
	    }

	    //std::cout << "    # barycentric sum : " << barycentric_sum / tet_vol - 1 << std::endl;

	    return (barycentricSum / tetrahedronVolume - 1 < tolerance);
	}
};



/*********************************************************************************************************************/
/** \brief If the global point is too far away from the global centre of the element return false.                  **/
/** Then it can not be inside, because elements with surface curvature more than the internal radius are unexpected **/
/*********************************************************************************************************************/
// [TODO] Can replace CoM-radius by smallest containing radius. However, algorithm for that is sophisticated, so perhaps not much gain
// [TODO] Can consider surface curvature to have more accurate zealous prefactor
template <class CurvGeom, int cdim, int mydim>
struct FarPointTest
{
	typedef typename CurvGeom::ctype              ctype;
	typedef typename CurvGeom::GlobalCoordinate   GlobalCoordinate;

	static bool isInside(const CurvGeom & curvgeom, const GlobalCoordinate &globalC, ctype tolerance)
	{
		// [FIXME] Throw dune-error if this generic method is called
		return false;
	}

};


// Case of mydim==cdim, aka for 2D triangles and 3D tetrahedra
// Tolerance is technically unnecessary here for Zealousness being much more conservative.
// But if we get the exact estimate from surface curvature, then tolerance will become relevant again
template <class CurvGeom, int cdim>
struct FarPointTest<CurvGeom, cdim, cdim>
{
	typedef typename CurvGeom::ctype              ctype;
	typedef typename CurvGeom::GlobalCoordinate   GlobalCoordinate;

	static bool isInside(const CurvGeom & curvgeom, const GlobalCoordinate &globalC, ctype tolerance)
	{
	    GlobalCoordinate com = curvgeom.center();

	    // Find distance from center of element to the point
	    ctype d = (com - globalC).two_norm();

	    // Find corner to which the distance from center is largest, find that distance, call it CoM-radius
	    ctype comRadius = 0;
	    for (int i = 0; i < curvgeom.nCorner(); i++) { comRadius = std::max(comRadius, (com - curvgeom.corner(i)).two_norm()); }

	    // Assume that curvilinear element is entirely contained within 2 * its CoM-radius
	    const ctype ZEALOUS_PREFACTOR = 2.0;
	    return (d < ZEALOUS_PREFACTOR * comRadius);
	}

};


// For 1D this test is much simpler - just check if globalC is inbetween the two corners of the edge
template <class CurvGeom>
struct FarPointTest<CurvGeom, DIM1D, DIM1D>
{
	typedef typename CurvGeom::ctype              ctype;
	typedef typename CurvGeom::GlobalCoordinate   GlobalCoordinate;

	static bool isInside(const CurvGeom & curvgeom, const GlobalCoordinate &globalC, ctype tolerance)
	{
	    return (curvgeom.corner(0)[0] - globalC[0] < tolerance)&&(globalC[0] - curvgeom.corner(1)[0] < tolerance);
	}

};




} // namespace CurvilinearPointLocation

} // namespace Dune




#endif // DUNE_CURVILINEAR_POINT_LOCATION_HH_
