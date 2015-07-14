/*******************************************************************
 * Numerical Recursive Interpolation Integrator
 * 
 * author: Aleksejs Fomins
 * date: 01.08.2014 - created
 * 
 * description:
 * Integrates arbitrary function over reference simplex (edge or triangle) using piecewise interpolation.
 * - Function is given in terms of a functor object. It need not be continuous
 * - The algorithm uses lagrange interpolation to sample the function at the interpolation points for 2nd and 4th order
 * - The integral of the real function over a simplex is approximated by an analytical integral of the interpolated function over this simplex
 * - If the total expected error is too high, the 4th order simplex with highest expected error is split into 2nd order simplices, which are then refined to 4th order and dealt with recursively. The algorithm re-uses previous sampling points when splitting
 * - The running error of a simplex is given by the difference between the approximated integrals for orders 2 and 4. All sub-simplices are always sorted by expected integration error regardless of their size
 * 
 * 
 * TODO: The algorithm is rather slow. It is possible that a higher interpolation order may be beneficial.
 * 
 *******************************************************************/


#ifndef DUNE_QUADRATUREINTEGRATOR_HH
#define DUNE_QUADRATUREINTEGRATOR_HH


#include <iostream>
#include <fstream>
#include <vector>
#include <queue>          // std::queue
#include <algorithm>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/function.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>





namespace Dune {




template<class ctype, int mydim>
class QuadratureIntegrator {
    typedef Dune::QuadratureRule<ctype, mydim>    QRule;
    typedef Dune::QuadratureRules<ctype, mydim>   QRules;

    static const int MAX_INT_ORDER = 60;


    template <class Geometry>
    struct IntegrationElementFunctor
    {
    	static const int mydimension = Geometry::mydimension;

    	typedef FieldVector< ctype, mydimension > LocalCoordinate;

        static const unsigned int RETURN_SIZE = 1;
        typedef ctype                        ResultValue;
        typedef typename std::vector<ctype>  ResultType;


        IntegrationElementFunctor(const Geometry & geom) : geom_(geom) {}

        ResultType operator()(const LocalCoordinate & x) const { return ResultType(1, geom_.integrationElement(x)); }

        ResultValue zeroValue(unsigned int rezIndex) const { return 0.0; }


        const Geometry & geom_;
    };









public:

    QuadratureIntegrator()  {}


    template<class IntegrandFunctor>
    struct Traits
    {
    	typedef typename IntegrandFunctor::ResultType         ResultType;
    	typedef typename IntegrandFunctor::ResultValue        ResultValue;
    	typedef typename std::pair<unsigned int, ResultType>  StatInfo;
    	typedef typename std::vector<StatInfo>                StatInfoVec;

    	static const unsigned int RETURN_SIZE = IntegrandFunctor::RETURN_SIZE;
    };


    // Computes integral of a functor over an element using the provided Geometry to calculate integration elements
    template<class Geometry, class IntegrandFunctor>
    static typename Traits<IntegrandFunctor>::ResultType integrate(
        const Geometry & geometry,
        const IntegrandFunctor & f,
        unsigned int integrOrder
    )
    {
        IntegrationElementFunctor<Geometry> detJ(geometry);
        return integrateImpl(geometry.type(), f, integrOrder, detJ);
    }


    // Computes integral of a functor over an element using the provided JacobianFunctor to calculate integration elements
    template<class JacobiFunctor, class IntegrandFunctor>
    static typename Traits<IntegrandFunctor>::ResultType integrate(
        Dune::GeometryType gt,
        const IntegrandFunctor & f,
        unsigned int integrOrder,
        const JacobiFunctor detJ)
    {
        return integrateImpl(gt, f, integrOrder, detJ);
    }


    // Computes integral of a functor over an element recursively using the provided Geometry to calculate integration elements
    template<class Geometry, class IntegrandFunctor>
    static typename Traits<IntegrandFunctor>::StatInfo integrateRecursive(
        const Geometry & geometry,
        const IntegrandFunctor & f,
        ctype rel_tol,
        unsigned int suggestedOrder = 1)
    {
        IntegrationElementFunctor<Geometry> detJ(geometry);
        return integrateRecursiveImpl(geometry.type(), f, rel_tol, suggestedOrder, detJ);
    }


    // Computes integral of a functor over an element recursively using the provided JacobianFunctor to calculate integration elements
    template<class JacobiFunctor, class IntegrandFunctor>
    static typename Traits<IntegrandFunctor>::StatInfo integrateRecursive(
        Dune::GeometryType gt,
        const IntegrandFunctor & f,
        ctype rel_tol,
        const JacobiFunctor & detJ,
        unsigned int suggestedOrder = 1)
    {
        return integrateRecursiveImpl(gt, f, rel_tol, suggestedOrder, detJ);
    }


    // Integrates functor using all specified quadrature orders, returns vector of values
    template<class Geometry, class IntegrandFunctor>
    static typename Traits<IntegrandFunctor>::StatInfoVec integrateStat(
    		const Geometry & geometry,
    		const IntegrandFunctor & f,
    		unsigned int integrOrderMax)
    {
    	typedef typename Traits<IntegrandFunctor>::ResultType   ResultType;
    	typedef typename Traits<IntegrandFunctor>::StatInfo     StatInfo;

    	std::vector<StatInfo> rez;

        for (unsigned int i = 1; i <= integrOrderMax; i++)
        {
            rez.push_back(StatInfo(
                QRules::rule(geometry.type(), i).size(),
                integrate(geometry, f, i)
            ));
        }

        return rez;
    }


protected:


    // Integrates functor over reference element using quadrature of a given order
    template<class JacobiFunctor, class IntegrandFunctor>
    static typename Traits<IntegrandFunctor>::ResultType integrateImpl(
        Dune::GeometryType gt,
        const IntegrandFunctor & f,
        unsigned int integrOrder,
        const JacobiFunctor & detJ
    )
    {
    	typedef typename Traits<IntegrandFunctor>::ResultType   ResultType;
    	const unsigned int nResult = Traits<IntegrandFunctor>::RETURN_SIZE;

        const QRule & rule = QRules::rule(gt, integrOrder);
        if (rule.order() < integrOrder)  { DUNE_THROW(Dune::Exception,"order not available"); }

        // Initialise a zero result vector
        // In dynamic case, different parts of the result could be of different size
        ResultType result;
        for (unsigned int i = 0; i < nResult; i++)  { result.push_back( f.zeroValue(i) ); }

        for (typename QRule::const_iterator i=rule.begin(); i!=rule.end(); ++i)
        {
            ResultType fval = f(i->position());
            ctype weight = i->weight();
            ctype detjac = detJ(i->position())[0];

            //if (integrOrder > 6)  { std::cout << "eval f(x)=" << fval[0] << " detjac=" << detjac << std::endl; }

            for (unsigned int iResult = 0; iResult < nResult; iResult++)
            {
                fval[iResult] *= weight * detjac;
                result[iResult] += fval[iResult];
            }
        }

        return result;
    }


    // Integrates using gradually increasing quadrature order until estimated smooth relative tolerance achieved
    template<class JacobiFunctor, class IntegrandFunctor>
    static typename Traits<IntegrandFunctor>::StatInfo integrateRecursiveImpl(
            Dune::GeometryType gt,
            const IntegrandFunctor & f,
            ctype rel_tol,
            unsigned int suggestedOrder,
            const JacobiFunctor & jacobiFunctor)
    {
    	typedef typename Traits<IntegrandFunctor>::ResultType   ResultType;
    	typedef typename Traits<IntegrandFunctor>::ResultValue  ResultValue;
    	typedef typename Traits<IntegrandFunctor>::StatInfo     StatInfo;
    	const unsigned int nResult = Traits<IntegrandFunctor>::RETURN_SIZE;


    	//std::cout << "Started recursive integral over geometry " << gt << "with relative tolerance "<< rel_tol << " suggester order " << suggestedOrder << std::endl;


        // Initialise a zero result vector
        // In dynamic case, different parts of the result could be of different size
        ResultType zeroresult(nResult);
        for (unsigned int i = 0; i < nResult; i++)  { zeroresult[i] = f.zeroValue(i); }

        ctype SMOOTH_FACTOR = 0.15;

        // Setting the initial estimated error of the integration to a reasonable value basically forces the routine to at least make two iterations
        // This avoids the case where the first guess for the integral is wrongly estimated as 0 due to unlucky symmetry of the integrand
        std::vector<ctype> deltaThis(nResult, 1.0);

        ctype relError = 1.0;                  // Set the initial relative error to non-zero such that there is at least one iteration
        ResultType         resultThis(zeroresult);  // Initialize the integral result to 0. This is purely conventional, it could be anything

        // Start at an order slightly earlier than the suggested one, to use it as an error reference
        unsigned int order = (suggestedOrder > 1) ? suggestedOrder - 2 : 0;

        // Initialize the the current size of quadrature set to zero
        int prevQuadSize = 0;
        int thisQuadSize = 0;

        while (relError > rel_tol)
        {
            // Increase quadrature order. In Dune for some magical reason consecutive orders sometimes have the same quadrature
            // Keep increasing the order if the number of quadrature points did not change
            do
            {
                order++;
                if (order >= MAX_INT_ORDER)  { DUNE_THROW(Dune::IOError, "QUAD_INTEGRATOR failed to converge to required accuracy"); }
                prevQuadSize = thisQuadSize;
                thisQuadSize = QRules::rule(gt, order).size();
            } while (prevQuadSize == thisQuadSize);

            // Integrate
            ResultType resultNew = integrate(gt, f, order, jacobiFunctor);

            // Refresh the relative error to calculate it anew
            relError = 0;

            for (unsigned int iResult = 0; iResult < nResult; iResult++)
            {
            	// The effective absolute error is defined as the difference between integrals at this and previous quadrature
            	ResultValue deltaNewVec  = resultNew[iResult];
                            deltaNewVec -= resultThis[iResult];
                ctype deltaNew      = absoluteValue(deltaNewVec);

                // Define smoothened absolute error as a linear combination between this and previous errors
                // This avoids the case when two consecutive samplings give same wrong answer due to unlucky symmetry of the integrand wrt quadrature points
                ctype deltaSmooth   = SMOOTH_FACTOR * deltaThis[iResult] +  (1 - SMOOTH_FACTOR) * deltaNew;


            	// [TODO] If the previous error is small enough, re-computing the magnitude is unnecessary computational effort
            	ctype magnitudeThis = absoluteValue(resultNew[iResult]);       // Find the magnitude of this result
                ctype errorThis = relativeError(deltaSmooth, magnitudeThis);   // Calculate the relative error

                // If several quantities are integrated at the same time, only consider the worst error among them
                // This is because we want all of the quantities to integrate to at least the requested precision
                relError = std::max(relError, errorThis);

                // Update values of result and error with new ones
                resultThis[iResult] = resultNew[iResult];
                deltaThis[iResult] = deltaNew;
            }
            //std::cout << "--- processed order=" << order << " quadrature size=" << thisQuadSize << " estimated relative error=" << relError << " desired error=" << rel_tol << std::endl;

            // Write a matrix to a file for debugging purposes
            //if (nResult == 6)  {
            //if (order > 8)
            //{
            //	writeMatrix(resultThis[0], "/home/fomins/Documents/test/duneMatrix" + std::to_string(order) + ".txt");
            //}
            //            }
        }

        std::cout << "Finished recursive integral over geometry " << gt << "suggested order: " << suggestedOrder << " final order " << order << std::endl;

        return StatInfo(order, resultThis);
    }


    template <class VT>           static void writeMatrix(VT & mat, std::string filename) {}
    template <class VT, int dim>  static void writeMatrix(Dune::FieldVector<VT, dim> & mat, std::string filename) {}
    template <class VT>           static void writeMatrix(Dune::DynamicVector<VT> & mat,    std::string filename) {}

    template <class VT>
    static void writeMatrix (Dune::DynamicMatrix<VT> & mat, std::string filename)
    {
    	std::ofstream myfile (filename);
    	if (myfile.is_open())
    	{
    		for (int i = 0; i < mat.N(); i++)
    		{
        		for (int j = 0; j < mat.M(); j++)
        		{
        			myfile << mat[i][j] << std::endl;
        		}
    		}
    		myfile.close();
    	}
    }


    static ctype relativeError(ctype abs_err, ctype abs_mag)
    {
		// Only calculate relative error if the absolute value of the error is large enough
		// Otherwise just return absolute error
    	return (abs_err > 1.0e-15) ? abs_err / abs_mag : abs_err;
    }


    template<class ValueType>
    static ctype absoluteValue(ValueType & v)
    {
    	return std::abs(v);
    }


    template<class ValueType, int RESULT_DIM>
    static ctype absoluteValue(Dune::FieldVector<ValueType, RESULT_DIM> & v)
    {
    	return v.two_norm();
    }


    template<class ValueType>
    static ctype absoluteValue(Dune::DynamicVector<ValueType> & v)
    {
    	return v.two_norm();
    }


    template<class ValueType>
    static ctype absoluteValue(Dune::DynamicMatrix<ValueType> & v)
    {
    	return v.frobenius_norm();
    }







};

} // Namespace Dune

#endif /** DUNE_QUADRATUREINTEGRATOR_HH **/
