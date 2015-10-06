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


/***********************************************************************************************************/
/** \brief Calculates Determinant of a matrix                                                              */
/** FIXME This functionality sort of available in MatrixHelper, but it throws asserts I can not deal with  */
/***********************************************************************************************************/
template <class ctype, int dim>
struct MatrixDeterminant
{
	template <class Matrix>
	static ctype eval(const Matrix & M)  { DUNE_THROW(NotImplemented, "Called generic method of DifferentialHelper::MatrixDeterminant"); }
};

template <class ctype>
struct MatrixDeterminant<ctype, DIM1D>
{
	template <class Matrix>
	static ctype eval(const Matrix & M)  { return M[0][0]; }
};

template <class ctype>
struct MatrixDeterminant<ctype, DIM2D>
{
	template <class Matrix>
	static ctype eval(const Matrix & M)  { return M[0][0] * M[1][1] - M[0][1] * M[1][0]; }
};

template <class ctype>
struct MatrixDeterminant<ctype, DIM3D>
{
	template <class Matrix>
	static ctype eval(const Matrix & M)  {
        return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
             + M[0][1] * (M[1][2] * M[2][0] - M[1][0] * M[2][2])
             + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
	}
};





/*****************************************************************************************************/
/** \brief Calculates Jacobian of element mapping Analytically for all combinations of mydim, codim  */
/*****************************************************************************************************/

template <class CurvGeom>
struct JacobianTransposeAnalytical
{
	static const int mydimension    = CurvGeom::mydimension;
	static const int coorddimension = CurvGeom::coorddimension;

	typedef typename CurvGeom::PolynomialGlobalCoordinate     PolynomialGlobalCoordinate;
	typedef typename CurvGeom::PolynomialJacobianTransposed   PolynomialJacobianTransposed;

	static PolynomialJacobianTransposed eval( const PolynomialGlobalCoordinate & analyticalMap )
    {
		PolynomialJacobianTransposed polyjt;

        for (int i = 0; i < coorddimension; i++)
        {
            for (int j = 0; j < mydimension; j++)
            {
            	polyjt[j][i] = (analyticalMap[i].derivative(j));
            }
        }

      return polyjt;
    }
};


/*****************************************************************************************************/
/** \brief Calculates Hessian of a given polynomial analytically for any local dimension mydim       */
/*****************************************************************************************************/
template <class CurvGeom>
struct HessianTransposeAnalytical
{
	static const int mydimension    = CurvGeom::mydimension;

	typedef typename CurvGeom::LocalPolynomial                LocalPolynomial;
	typedef typename CurvGeom::PolynomialGlobalCoordinate     PolynomialLocalCoordinate;
	typedef typename CurvGeom::PolynomialHessianTransposed    PolynomialHessianTransposed;

	static PolynomialHessianTransposed eval( const LocalPolynomial & f )
    {
		PolynomialLocalCoordinate grad;
		PolynomialHessianTransposed polyht;

		// Pre-compute gradient to reduce number of derivatives
		for (int j = 0; j < mydimension; j++)  { grad[j] = f.derivative(j); }

		// Compute Hessian Transposed
        for (int i = 0; i < mydimension; i++)
        {
            for (int j = 0; j < mydimension; j++)
            {
            	polyht[j][i] = (grad[i].derivative(j));
            }
        }

      return polyht;
    }
};


/*****************************************************************************************************/
/** \brief Calculates Jacobian Determinant analytically as a polynomial, defined for mydim == cdim   */
/** Performs sign-correction by sampling the resulting map in the middle of the element              */
/*****************************************************************************************************/
template <class CurvGeom, int cdim, int mydim>
struct JacobianDeterminantAnalytical
{
	typedef typename CurvGeom::LocalPolynomial             LocalPolynomial;
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
	{
		DUNE_THROW(NotImplemented, "Called generic method of DifferentialHelper::JacobianDeterminantAnalytical");
	    return LocalPolynomial();
	}
};


template <class CurvGeom>
struct JacobianDeterminantAnalytical<CurvGeom, DIM1D, DIM1D>
{
	typedef typename CurvGeom::LocalPolynomial   LocalPolynomial;
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
	{
	    LocalPolynomial rez = analyticalMap[0].derivative(0);

	    // Change sign if determinant is negative
	    if (rez.evaluate(curvgeom.refElement().position( 0, 0 )) < 0) { rez *= -1; }

	    return rez;
	}
};


template <class CurvGeom>
struct JacobianDeterminantAnalytical<CurvGeom, DIM2D, DIM2D>
{
	typedef typename CurvGeom::LocalPolynomial   LocalPolynomial;
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
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
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
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
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static PolynomialGlobalCoordinate eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
	{
		DUNE_THROW(NotImplemented, "Called generic method of DifferentialHelper::NormalIntegrationElementAnalytical");
		PolynomialGlobalCoordinate rez;
	    return rez;
	}
};


// Case of edge in 2D
template <class CurvGeom>
struct NormalIntegrationElementAnalytical<CurvGeom, DIM2D, DIM1D>
{
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static PolynomialGlobalCoordinate eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
	{
	    PolynomialGlobalCoordinate rez;

	    rez[0] = analyticalMap[1].derivative(0);
	    rez[1] = analyticalMap[0].derivative(0) * (-1);

	    return rez;
	}
};


// Case of face in 3D
template <class CurvGeom>
struct NormalIntegrationElementAnalytical<CurvGeom, DIM3D, DIM2D>
{
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static PolynomialGlobalCoordinate eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
	{
	    PolynomialGlobalCoordinate rez;

	    rez[0] = analyticalMap[2].derivative(0) * analyticalMap[1].derivative(1) - analyticalMap[1].derivative(0) * analyticalMap[2].derivative(1);
	    rez[1] = analyticalMap[0].derivative(0) * analyticalMap[2].derivative(1) - analyticalMap[2].derivative(0) * analyticalMap[0].derivative(1);
	    rez[2] = analyticalMap[1].derivative(0) * analyticalMap[0].derivative(1) - analyticalMap[0].derivative(0) * analyticalMap[1].derivative(1);

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
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
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
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
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
	typedef typename CurvGeom::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;

	static LocalPolynomial eval(const CurvGeom & curvgeom, const PolynomialGlobalCoordinate & analyticalMap)
	{
		LocalPolynomial rez;

	    PolynomialGlobalCoordinate normalIntegrationElement = NormalIntegrationElementAnalytical<CurvGeom, DIM3D, DIM2D>::eval(curvgeom, analyticalMap);
	    for (int i = 0; i < DIM3D; i++) {
	        rez += normalIntegrationElement[i] * normalIntegrationElement[i];
	    }

	    return rez;
	}
};





} // namespace DifferentialHelper

} // namespace Dune




#endif
