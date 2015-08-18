#ifndef DUNE_DIFFERENTIAL_HELPER_HH_
#define DUNE_DIFFERENTIAL_HELPER_HH_

#include <dune/curvilineargeometry/interpolation/polynomial.hh>


namespace Dune
{


// Forwards-declaration of the curvilinear Geometry Traits
//template< class ct >
//struct CurvilinearGeometryTraits;

// Forwards-declaration of the curvilinear Geometry
//template< class ct, int mydim, int cdim, class Traits = CurvilinearGeometryTraits<ct> >
//class CurvilinearGeometry;


// [TODO] Change template parameters to reuse CurvGeom-intrinsic cdim and mydim
namespace DifferentialHelper
{

const int DIM0D = 0;
const int DIM1D = 1;
const int DIM2D = 2;
const int DIM3D = 3;





/*****************************************************************************************************/
/** \brief Calculates Jacobian Determinant analytically as a polynomial, defined for mydim == cdim   */
/** Performs sign-correction by sampling the resulting map in the middle of the element              */
/*****************************************************************************************************/
template <class CurvGeom, int cdim, int mydim>
struct JacobianDeterminantAnalytical
{
	typedef typename CurvGeom::LocalPolynomial  LocalPolynomial;
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
		DUNE_THROW(NotImplemented, "Called generic method of DifferentialHelper::JacobianDeterminantAnalytical");
	    return LocalPolynomial();
	}
};


template <class CurvGeom>
struct JacobianDeterminantAnalytical<CurvGeom, DIM1D, DIM1D>
{
	typedef typename CurvGeom::LocalPolynomial   LocalPolynomial;
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
	    LocalPolynomial rez = analyticalMap[0].derivative(0);

	    // Clean-up in case some summands summed up to 0
	    rez.cleanUp();

	    // Change sign if determinant is negative
	    if (rez.evaluate(curvgeom.refElement().position( 0, 0 )) < 0) { rez *= -1; }

	    return rez;
	}
};


template <class CurvGeom>
struct JacobianDeterminantAnalytical<CurvGeom, DIM2D, DIM2D>
{
	typedef typename CurvGeom::LocalPolynomial   LocalPolynomial;
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
	    LocalPolynomial rez = analyticalMap[0].derivative(0) * analyticalMap[1].derivative(1) - analyticalMap[0].derivative(1) * analyticalMap[1].derivative(0);

	    // Clean-up in case some summands summed up to 0
	    rez.cleanUp();

	    // Change sign if determinant is negative
	    if (rez.evaluate(curvgeom.refElement().position( 0, 0 )) < 0) { rez *= -1; }

	    return rez;
	}
};


template <class CurvGeom>
struct JacobianDeterminantAnalytical<CurvGeom, DIM3D, DIM3D>
{
	typedef typename CurvGeom::LocalPolynomial   LocalPolynomial;
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
	    LocalPolynomial rez;

	    rez += analyticalMap[0].derivative(0) * ( analyticalMap[1].derivative(1) * analyticalMap[2].derivative(2) - analyticalMap[1].derivative(2) * analyticalMap[2].derivative(1) );
	    rez += analyticalMap[0].derivative(1) * ( analyticalMap[1].derivative(2) * analyticalMap[2].derivative(0) - analyticalMap[1].derivative(0) * analyticalMap[2].derivative(2) );
	    rez += analyticalMap[0].derivative(2) * ( analyticalMap[1].derivative(0) * analyticalMap[2].derivative(1) - analyticalMap[1].derivative(1) * analyticalMap[2].derivative(0) );

	    // Clean-up in case some summands summed up to 0
	    rez.cleanUp();

	    // Change sign if determinant is negative
	    if (rez.evaluate(curvgeom.refElement().position( 0, 0 )) < 0) { rez *= -1; }

	    return rez;
	}
};



/*****************************************************************************************************/
/** \brief Calculates the surface normal integration element (n * dJ), defined for mydim < cdim.     */
/** This quantity is always polynomial                                                               */
/*****************************************************************************************************/
template <class CurvGeom, int cdim, int mydim>
struct NormalIntegrationElementAnalytical
{
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static PolynomialVector eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
		DUNE_THROW(NotImplemented, "Called generic method of DifferentialHelper::NormalIntegrationElementAnalytical");
	    PolynomialVector rez;
	    return rez;
	}
};


// Case of edge in 2D
template <class CurvGeom>
struct NormalIntegrationElementAnalytical<CurvGeom, DIM2D, DIM1D>
{
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static PolynomialVector eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
	    PolynomialVector rez;

	    rez.push_back(analyticalMap[1].derivative(0));
	    rez.push_back(analyticalMap[0].derivative(0) * (-1));

	    return rez;
	}
};


// Case of face in 3D
template <class CurvGeom>
struct NormalIntegrationElementAnalytical<CurvGeom, DIM3D, DIM2D>
{
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static PolynomialVector eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
	    PolynomialVector rez;

	    rez.push_back(analyticalMap[2].derivative(0) * analyticalMap[1].derivative(1) - analyticalMap[1].derivative(0) * analyticalMap[2].derivative(1));
	    rez.push_back(analyticalMap[0].derivative(0) * analyticalMap[2].derivative(1) - analyticalMap[2].derivative(0) * analyticalMap[0].derivative(1));
	    rez.push_back(analyticalMap[1].derivative(0) * analyticalMap[0].derivative(1) - analyticalMap[0].derivative(0) * analyticalMap[1].derivative(1));

	    return rez;
	}
};



/*****************************************************************************************************/
/** \brief Calculates the integration element squared as a polynomial, defined for mydim < cdim      */
/** This quantity is always polynomial                                                               */
/*****************************************************************************************************/
// [TODO] In principle this quantity can be extended to all cdim-mydim pairs for completeness, but it is unnecessary for computation
template <class CurvGeom, int cdim, int mydim>
struct IntegrationElementSquaredAnalytical
{
	typedef typename CurvGeom::LocalPolynomial  LocalPolynomial;
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
		DUNE_THROW(NotImplemented, "Called generic method of DifferentialHelper::IntegrationElementSquaredAnalytical");
	    LocalPolynomial rez;
	    return rez;
	}
};


// Specialization for edges in any dimension
template <class CurvGeom, int cdim>
struct IntegrationElementSquaredAnalytical<CurvGeom, cdim, DIM1D>
{
	typedef typename CurvGeom::LocalPolynomial  LocalPolynomial;
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
		LocalPolynomial rez;

	    for (int i = 0; i < cdim; i++) {
	    	LocalPolynomial tmp = analyticalMap[i].derivative(0);
	        rez += tmp * tmp;
	    }

	    return rez;
	}
};


template <class CurvGeom>
struct IntegrationElementSquaredAnalytical<CurvGeom, DIM3D, DIM2D>
{
	typedef typename CurvGeom::LocalPolynomial   LocalPolynomial;
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialVector & analyticalMap)
	{
		LocalPolynomial rez;

	    PolynomialVector normalIntegrationElement = NormalIntegrationElementAnalytical<CurvGeom, DIM3D, DIM2D>::eval(curvgeom, analyticalMap);
	    for (int i = 0; i < DIM3D; i++) {
	        rez += normalIntegrationElement[i] * normalIntegrationElement[i];
	    }

	    return rez;
	}
};





} // namespace DifferentialHelper

} // namespace Dune




#endif
