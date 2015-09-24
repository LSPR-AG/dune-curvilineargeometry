#ifndef DUNE_INTEGRATION_HELPER_HH_
#define DUNE_INTEGRATION_HELPER_HH_

#include <assert.h>

#include <dune/curvilineargeometry/integration/quadratureintegrator.hh>
#include <dune/curvilineargeometry/integration/adaptiveintegrator.hh>

namespace Dune
{


template<class CurvGeom, int NORM_TYPE>
class IntegralHelper
{

	typedef typename CurvGeom::ctype   ctype;
    static const int mydimension    = CurvGeom::mydimension;  //! geometry dimension
    static const int coorddimension = CurvGeom::coorddimension;       //! coordinate dimension

    typedef typename CurvGeom::LocalCoordinate   LocalCoordinate;
    typedef typename CurvGeom::GlobalCoordinate  GlobalCoordinate;

	typedef typename CurvGeom::LocalPolynomial   LocalPolynomial;
	typedef typename CurvGeom::PolynomialVector  PolynomialVector;




	  // Evaluates the generalized Jacobian determinant using pre-computed polynomial integration elements
	  // If the integral is over codim-0 element, the integration element is just the jacobian determinant
	  // Otherwise sqrt(det(J^T * J))
	  // [TODO] Template operator to avoid the unnecessary if-statement inside
	  struct JacobianFunctor
	  {
	      static const unsigned int RETURN_SIZE = 1;
	      typedef ctype                        ResultValue;
	      typedef typename std::vector<ctype>  ResultType;

	      LocalPolynomial integration_element_generalised_;

	      JacobianFunctor(const LocalPolynomial & integration_element_g) :
	          integration_element_generalised_ (integration_element_g)
	      {
	    	  // Optimize polynomial for faster evaluation
	    	  integration_element_generalised_.cache();
	      }

	      ResultType     operator()(const LocalCoordinate & x) const {
	          ctype rez;

	          if (mydimension == coorddimension)  { rez = integration_element_generalised_.evaluate(x); }
	          else                                { rez = sqrt(integration_element_generalised_.evaluate(x)); }

	          return ResultType(1, rez);
	      }

	      bool isPolynomial() const { return (mydimension == coorddimension)||(integration_element_generalised_.order() <= 1); }

	      unsigned int expectedOrder() const {
	          return (mydimension == coorddimension) ? integration_element_generalised_.order() : round(sqrt(integration_element_generalised_.order()));
	      }

	      ResultValue zeroValue(unsigned int rezIndex) const { return 0.0; }
	  };


	  // Functor that evaluates a polynomial that it carries inside
	  template<class ctype, int mydim>
	  struct PolynomialFunctor
	  {
	      typedef Polynomial<ctype, mydim> LocalPolynomial;
	      typedef FieldVector< ctype, mydim > LocalCoordinate;

	      static const unsigned int RETURN_SIZE = 1;
	      typedef ctype                        ResultValue;
	      typedef typename std::vector<ctype>  ResultType;


	      PolynomialFunctor(const LocalPolynomial & p) : p_(p) {
	    	  p_.cache();
	      }

	      ResultType operator()(const LocalCoordinate & x) const { return ResultType(1, p_.evaluate(x)); }

	      bool isPolynomial() const { return true; }

	      unsigned int expectedOrder() const  { return p_.order(); }

	      ResultValue zeroValue(unsigned int rezIndex) const { return 0.0; }

	      LocalPolynomial p_;
	  };


public:


	/** \brief Integrate the given polynomial over the element
	 *
	 *  If the dimension of the element and world match, the analytical integration will be performed.
	 *  Otherwise, analytical integration is not possible, and numerical integration will be performed.
	 *
	 *  \param[in]  P  Polynomial to integrate (for example, a basis function)
	 *  \param[in]  tolerance  the acceptable relative error for numerical integration
	 *
	 *  \returns the result of the integral
	 */
	static ctype integrateScalar(const CurvGeom & curvgeom, const LocalPolynomial & P, ctype RELATIVE_TOLERANCE, ctype ACCURACY_GOAL)
	{
	    if (mydimension == coorddimension)  { return integrateAnalyticalScalar(curvgeom, P); }
	    else                                { return integrateNumerical(curvgeom, PolynomialFunctor<ctype, mydimension>(P), RELATIVE_TOLERANCE, ACCURACY_GOAL); }
	}


	/** \brief Integrates the given functor numerically over the element
	 *
	 *
	 *  \param[in]  f  Functor which maps a local coordinate to a real value
	 *  \param[in]  tolerance  the acceptable relative error for numerical integration
	 *
	 *  \returns the result of the integral
	 */
	template <typename Functor>
	static typename Functor::ResultType::value_type integrateNumerical(
		const CurvGeom & curvgeom,
		const Functor & f,
		ctype RELATIVE_TOLERANCE,
	    ctype ACCURACY_GOAL)
	{
	    if (mydimension == coorddimension)
	    {
	    	LocalPolynomial detJ = curvgeom.JacobianDeterminantAnalytical();
	        JacobianFunctor g(detJ);
	        return integrateNumericalRecursive(curvgeom, f, g, RELATIVE_TOLERANCE, ACCURACY_GOAL)[0];
	    }
	    else
	    {
	    	LocalPolynomial normInt2 = curvgeom.IntegrationElementSquaredAnalytical();
	        JacobianFunctor g(normInt2);
	        return integrateNumericalRecursive(curvgeom, f, g, RELATIVE_TOLERANCE, ACCURACY_GOAL)[0];
	    }
	}

	/** \brief Integrate given polynomial analytically over the element
	 *
	 *  \param[in]  P    Polynomial to integrate over the element
	 *  \returns an integral of P over the element
	 *  \note (!) Analytical integration is only possible when mydimension == coorddimension
	 *  \note (!) At the moment only capable of integrating over a Simplex
	 */
	static ctype integrateAnalyticalScalar(const CurvGeom & curvgeom, const LocalPolynomial & P)
	{
		LocalPolynomial detJ = curvgeom.JacobianDeterminantAnalytical();
		return integrateAnalyticalScalar(P, detJ);
	}

	/** \brief Integrate given polynomial vector normal projection analytically over the element
	 *
	 *  \param[in]  PVec    Polynomial vector to project-integrate over the element
	 *  \returns an integral of PVec over the element boundary, projected
	 *  \note (!) This operation only has meaning to elements which have normals - edges in 2D and faces in 3D
	 *  \note (!) At the moment only capable of integrating over a Simplex
	 */
	static ctype integrateAnalyticalDot(const CurvGeom & curvgeom, const PolynomialVector & PVec)
	{
		PolynomialVector nintelem = curvgeom.NormalIntegrationElementAnalytical();
	    return integrateAnalyticalDot(PVec, nintelem);
	}


	/** \brief Integrate given polynomial vector normal cross product analytically over the element
	 *
	 *  \param[in]  PVec    Polynomial vector to cross-product-integrate over the element
	 *  \returns a coordinate - integral of PVec cross-product with normal over the element boundary
	 *  \note (!) This operation only has meaning to elements which have normals - edges in 2D and faces in 3D
	 *  \note (!) At the moment only capable of integrating over a Simplex
	 *
	 *  TODO: NOT IMPLEMENTED YET
	 */

	static GlobalCoordinate integrateAnalyticalTimes(const CurvGeom & curvgeom, const PolynomialVector & PVec)
	{
	    GlobalCoordinate rez;

	    // For 2D must return a scalar
	    // For 3D must return a 3D vector
	    // Maybe first enquire if this functionality is really necessary, then implement?

	    return rez;
	}



protected:


    /** \brief performs numerical integration of polynomial and non-polynomial functions. This functionality is necessary
     * for the CurvilinearGeometry to be able to calculate volumes of its own elements, which may have non-polynomial DetJac */
    template <typename Functor, typename JacobiFunctor>
    static typename Functor::ResultType integrateNumericalRecursive(
    		const CurvGeom & curvgeom,
            const Functor & f,
            const JacobiFunctor & jacobiDet,
            ctype RELATIVE_TOLERANCE,
            ctype ACCURACY_GOAL)
    {
        assert(mydimension > 0);                // We can not currently integrate over vertices, in principle this could be replaced by Dirac evaluation
        const int suggestedOrder = f.expectedOrder() + jacobiDet.expectedOrder();



        typedef Dune::QuadratureIntegrator<ctype, mydimension>  QuadratureIntegrator;

        if (f.isPolynomial() && jacobiDet.isPolynomial()) {  // If the integrand is polynomial, it can be integrated using fixed order
        	return QuadratureIntegrator::integrate(curvgeom.type(), f, suggestedOrder, jacobiDet);
        } else {                                             // Otherwise, the order has to be determined recursively
        	return QuadratureIntegrator::template integrateRecursive<JacobiFunctor, Functor, NORM_TYPE>(curvgeom.type(), f, jacobiDet, RELATIVE_TOLERANCE, ACCURACY_GOAL, suggestedOrder).second;
        }

        //AdaptiveIntegrator<ct, mydimension> NInt( type() );
        //return NInt.integrate( g, tolerance);
    }


    /** \brief performs analytical volume integral of a scalar function given that coorddimension==mydimension */
    static ctype integrateAnalyticalScalar(const LocalPolynomial & P, const LocalPolynomial & jacobianDeterminant)
    {
        assert(mydimension > 0);
        return (P * jacobianDeterminant).integrateRefSimplex();
    }


    /** \brief performs analytical surface integral of dot(f, n) */
    static ctype integrateAnalyticalDot(const PolynomialVector & PVec, const PolynomialVector & normalIntegrationElement )
    {
        assert(mydimension > 0);
        // Construct boundary integration element normal polynomial vector
        LocalPolynomial integrand;
        for (int i = 0; i < coorddimension; i++) { integrand += PVec[i] * normalIntegrationElement[i]; }
        return integrand.integrateRefSimplex();
    }




};


} // namespace Dune


#endif  //DUNE_INTEGRATION_HELPER_HH_
