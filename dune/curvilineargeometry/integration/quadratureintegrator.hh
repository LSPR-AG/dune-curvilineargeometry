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
#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/curvilineargeometry/integration/quadraturerelativeerror.hh>
#include <dune/curvilineargeometry/integration/quadratureabsolutevalue.hh>



namespace Dune {




template<class ctype, int mydim>
class QuadratureIntegrator {
    typedef Dune::QuadratureRule<ctype, mydim>    QRule;
    typedef Dune::QuadratureRules<ctype, mydim>   QRules;

    static const int MAX_INT_ORDER = 60;   // Maximal Quadrature order currently available in Dune


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

    struct IdentityFunctor
    {
    	static const int mydimension = mydim;
    	typedef FieldVector< ctype, mydimension > LocalCoordinate;

        typedef ctype                        ResultValue;
        typedef typename std::vector<ctype>  ResultType;

        static const unsigned int RETURN_SIZE = 1;

        IdentityFunctor() {}

        ResultType operator()(const LocalCoordinate & x) const { return ResultType(1, ResultValue(1.0)); }

        ResultValue zeroValue(unsigned int rezIndex) const { return ResultValue(0.0); }
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
    template<class Geometry, class IntegrandFunctor, int NormType>
    static typename Traits<IntegrandFunctor>::StatInfo integrateRecursive(
        const Geometry & geometry,
        const IntegrandFunctor & f,
        ctype RELATIVE_TOLERANCE,
        ctype ACCURACY_GOAL,
        unsigned int suggestedOrder = 1)
    {
    	typedef IntegrationElementFunctor<Geometry>  JacobiFunctor;

    	JacobiFunctor detJ(geometry);
        return integrateRecursiveImpl<JacobiFunctor, IntegrandFunctor, NormType>(geometry.type(), f, detJ, RELATIVE_TOLERANCE, ACCURACY_GOAL, suggestedOrder);
    }


    // Computes integral of a functor over an element recursively using the provided JacobianFunctor to calculate integration elements
    template<class JacobiFunctor, class IntegrandFunctor, int NormType>
    static typename Traits<IntegrandFunctor>::StatInfo integrateRecursive(
        Dune::GeometryType gt,
        const IntegrandFunctor & f,
        const JacobiFunctor & detJ,
        ctype RELATIVE_TOLERANCE,
        ctype ACCURACY_GOAL,
        unsigned int suggestedOrder = 1)
    {
        return integrateRecursiveImpl<JacobiFunctor, IntegrandFunctor, NormType>(gt, f, detJ, RELATIVE_TOLERANCE, ACCURACY_GOAL, suggestedOrder);
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
            rez.push_back(StatInfo( QRules::rule(geometry.type(), i).size(), integrate(geometry, f, i) ));
        }

        return rez;
    }


    // Integrates functor using all specified quadrature orders, returns vector of values
    template<class JacobiFunctor, class IntegrandFunctor>
    static typename Traits<IntegrandFunctor>::StatInfoVec integrateStat(
            Dune::GeometryType gt,
    		const IntegrandFunctor & f,
            const JacobiFunctor & detJ,
    		unsigned int integrOrderMax)
    {
    	typedef typename Traits<IntegrandFunctor>::ResultType   ResultType;
    	typedef typename Traits<IntegrandFunctor>::StatInfo     StatInfo;

    	std::vector<StatInfo> rez;

        for (unsigned int i = 1; i <= integrOrderMax; i++)
        {
            rez.push_back(StatInfo( QRules::rule(gt, i).size(), integrate(gt, f, i, detJ) ));
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
    template<class JacobiFunctor, class IntegrandFunctor, int NormType>
    static typename Traits<IntegrandFunctor>::StatInfo integrateRecursiveImpl(
            Dune::GeometryType gt,
            const IntegrandFunctor & f,
            const JacobiFunctor & jacobiFunctor,
            ctype RELATIVE_TOLERANCE,
            ctype ACCURACY_GOAL,
            unsigned int suggestedOrder)
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
        std::vector<ctype> relErrorThis(nResult, 1.0);

        ctype relErrorRez = 1.0;                  // Set the initial relative error to non-zero such that there is at least one iteration
        ResultType         resultThis(zeroresult);  // Initialize the integral result to 0. This is purely conventional, it could be anything

        // Start at an order slightly earlier than the suggested one, to use it as an error reference
        unsigned int order = (suggestedOrder > 1) ? suggestedOrder - 2 : 0;

        // Initialize the the current size of quadrature set to zero
        int prevQuadSize = 0;
        int thisQuadSize = 0;

        while (relErrorRez > RELATIVE_TOLERANCE)
        {
            // Increase quadrature order. In Dune for some magical reason consecutive orders sometimes have the same quadrature
            // Keep increasing the order if the number of quadrature points did not change
            do
            {
                order++;
                if (order >= MAX_INT_ORDER)  {
                	// [FIXME] Debug. One should really make the user aware if there is poor convergence
                	std::cout << "!!!Warning!!!: " << "QUAD_INTEGRATOR failed to converge to required accuracy" << std::endl;
                	std::cout << "Integral=" << writeVector(resultThis) << std::endl;
                	std::cout << "Error=" << writeVector(relErrorThis) << std::endl;

                	return StatInfo(order, resultThis);

                	//DUNE_THROW(Dune::IOError, "QUAD_INTEGRATOR failed to converge to required accuracy");
                }
                prevQuadSize = thisQuadSize;
                thisQuadSize = QRules::rule(gt, order).size();
            } while (prevQuadSize == thisQuadSize);

            // Integrate
            ResultType resultNew = integrate(gt, f, order, jacobiFunctor);

            // Refresh the relative error to calculate it anew
            relErrorRez = 0;

            for (unsigned int iResult = 0; iResult < nResult; iResult++)
            {
            	// The effective absolute error is defined as the difference between integrals at this and previous quadrature
            	ResultValue delta  = resultNew[iResult];
                            delta -= resultThis[iResult];

                ctype relErrorNew = Dune::QuadratureRelativeError<ctype, NormType>::eval(delta, resultNew[iResult], ACCURACY_GOAL);

                // Define smoothened absolute error as a linear combination between this and previous errors
                // This avoids the case when two consecutive samplings give same wrong answer due to unlucky symmetry of the integrand wrt quadrature points
                ctype relErrorSmooth   = SMOOTH_FACTOR * relErrorThis[iResult] +  (1 - SMOOTH_FACTOR) * relErrorNew;

                // If several quantities are integrated at the same time, only consider the worst error among them
                // This is because we want all of the quantities to integrate to at least the requested precision
                relErrorRez = std::max(relErrorRez, relErrorSmooth);

                // Update values of result and error with new ones
                resultThis[iResult]   = resultNew[iResult];
                relErrorThis[iResult] = relErrorNew;
            }



            // [FIXME] DEBUG
            /**
            std::vector<ctype> normRez(nResult);
            for (unsigned int iResult = 0; iResult < nResult; iResult++)  {
            	normRez[iResult] = Dune::QuadratureAbsoluteValue<ctype, QUADRATURE_NORM_L2>::eval(resultThis[iResult]);
            }

            std::cout << "--- processed order=" << order
            		  << " quadrature size="  << thisQuadSize
					  << " result absolute value=" << writeVector<ctype>(normRez)
            		  << " estimated relative error=" << writeVector<ctype>(relErrorThis)
            		  << " desired error=" << RELATIVE_TOLERANCE
            		  << " accuracy goal=" << ACCURACY_GOAL << std::endl;
            **/



            // Write a matrix to a file for debugging purposes
            // FIXME DEBUG
            /*
            if (nResult == 2)  {
            if (order > 15)
            {
            	writeMatrix(resultThis[1], "/home/fomins/Documents/test/duneMatrix" + std::to_string(order) + ".txt");
            }
                        }
            */
        }

        /** If suggested order significantly underestimates the actual order, the user should know about it  */
        /*
        if (order - suggestedOrder > 2) {
        	std::cout << "Warning: Integral over " << gt << " converged at order " << order << ", when suggested order is " << suggestedOrder << std::endl;
        }
        */



        return StatInfo(order, resultThis);
    }

    template<typename Val>
    static std::string writeVector(std::vector<Val>&  a)
    {
    	std::stringstream aaa;

    	aaa << "[";
    	for (int i = 0; i < a.size(); i++)  { aaa << a[i] << ","; }
    	aaa << "]";
    	return aaa.str();
    }



    template <class VT>           static void writeMatrix(VT & mat, std::string filename) {}
    template <class VT>           static void writeMatrix(Dune::DynamicVector<VT> & mat,    std::string filename)
    {
    	unsigned int L = mat.N();
    	unsigned int nDof = floor(0.5 * (sqrt(1 + 8 * L) - 1));
    	unsigned int pos = 0;

    	std::cout << "Writing matrix to matlab " << L << " " << nDof << " " << std::endl;

    	std::ofstream myfile (filename);
    	if (myfile.is_open())
    	{
        	for (int i = 0; i < nDof; i++)
        	{
            	for (int j = 0; j <= i; j++)
            	{
            		myfile << mat[pos++] << " ";
            	}
            	myfile << std::endl;
        	}
    		myfile.close();
    	}
    }



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









};

} // Namespace Dune

#endif /** DUNE_QUADRATUREINTEGRATOR_HH **/
