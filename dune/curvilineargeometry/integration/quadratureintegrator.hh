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
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>



namespace Dune {




template<class ctype, int mydim, class IntegrandFunctor>
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
        typedef typename std::vector<ctype>  ResultType;


        IntegrationElementFunctor(const Geometry & geom) : geom_(geom) {}

        ResultType operator()(const LocalCoordinate & x) const { return ResultType(1, geom_.integrationElement(x)); }

        const Geometry & geom_;
    };








public:

    static const unsigned int nResult = IntegrandFunctor::RETURN_SIZE;
	typedef typename IntegrandFunctor::ResultType   ResultType;
	typedef typename ResultType::value_type         ResultValue;

    typedef std::pair<int, ResultType>              StatInfo;       // To store pair (integrOrder, result)
    typedef typename std::vector<StatInfo>          StatInfoVec;






    QuadratureIntegrator()  {}


    // Computes integral of a functor over an element using the provided Geometry to calculate integration elements
    template<class Geometry>
    static ResultType integrate(
        const Geometry & geometry,
        const IntegrandFunctor & f,
        int integrOrder
    )
    {
        IntegrationElementFunctor<Geometry> detJ(geometry);
        return integrateImpl(geometry.type(), f, integrOrder, detJ);
    }


    // Computes integral of a functor over an element using the provided JacobianFunctor to calculate integration elements
    template<class JacobiFunctor>
    static ResultType integrate(
        Dune::GeometryType gt,
        const IntegrandFunctor & f,
        int integrOrder,
        const JacobiFunctor detJ
    )
    {
        return integrateImpl(gt, f, integrOrder, detJ);
    }


    // Computes integral of a functor over an element recursively using the provided Geometry to calculate integration elements
    template<class Geometry>
    static StatInfo integrateRecursive(
        const Geometry & geometry,
        const IntegrandFunctor & f,
        ctype rel_tol,
        unsigned int suggestedOrder = 1
    )
    {
        IntegrationElementFunctor<Geometry> detJ(geometry);
        return integrateRecursiveImpl(geometry.type(), f, rel_tol, suggestedOrder, detJ);
    }


    // Computes integral of a functor over an element recursively using the provided JacobianFunctor to calculate integration elements
    template<class JacobiFunctor>
    static StatInfo integrateRecursive(
        Dune::GeometryType gt,
        const IntegrandFunctor & f,
        ctype rel_tol,
        const JacobiFunctor & detJ,
        unsigned int suggestedOrder = 1
    )
    {
        return integrateRecursiveImpl(gt, f, rel_tol, suggestedOrder, detJ);
    }


    // Integrates functor using all specified quadrature orders, returns vector of values
    template<class Geometry>
    static StatInfoVec integrateStat(const Geometry & geometry, const IntegrandFunctor & f, int integrOrderMax)
    {
        StatInfoVec rez;

        for (int i = 1; i <= integrOrderMax; i++)
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
    template<class JacobiFunctor>
    static ResultType integrateImpl(
        Dune::GeometryType gt,
        const IntegrandFunctor & f,
        int integrOrder,
        const JacobiFunctor & detJ
    )
    {
        const QRule & rule = QRules::rule(gt, integrOrder);
        if (rule.order() < integrOrder)  { DUNE_THROW(Dune::Exception,"order not available"); }

        ResultType result(nResult, ResultValue(0.0));  // Assume automatic init by zero
        for (typename QRule::const_iterator i=rule.begin(); i!=rule.end(); ++i)
        {
            ResultType fval = f(i->position());

            ctype weight = i->weight();
            ctype detjac = detJ(i->position())[0];

            for (int iResult = 0; iResult < nResult; iResult++)
            {
                fval[iResult] *= weight * detjac;
                result[iResult] += fval[iResult];
            }
        }

        return result;
    }


    // Integrates using gradually increasing quadrature order until estimated smooth relative tolerance achieved
    template<class JacobiFunctor>
    static StatInfo integrateRecursiveImpl(
            Dune::GeometryType gt,
            const IntegrandFunctor & f,
            ctype rel_tol,
            unsigned int suggestedOrder,
            const JacobiFunctor & jacobiFunctor)
    {
    	std::cout << "Started recursive integral over geometry " << gt << "with relative tolerance "<< rel_tol << " suggester order " << suggestedOrder << std::endl;


        ctype SMOOTH_FACTOR = 0.15;

        ctype error = 1.0;
        ResultType rez_prev(nResult, ResultValue(0.0));
        ResultType rez_this(nResult, ResultValue(0.0));
        unsigned int order = suggestedOrder - 1;

        while (error > rel_tol)
        {
            // Increase order.
            // Keep increasing the order if the number of quadrature points did not increase
            int size_this = QRules::rule(gt, order).size();
            int size_next = size_this;

            while (size_this == size_next)
            {
                order++;
                if (order >= MAX_INT_ORDER)  { DUNE_THROW(Dune::IOError, "QUAD_INTEGRATOR failed to converge to required accuracy"); }
                size_this = size_next;
                size_next = QRules::rule(gt, order).size();
            }

            // Integrate
            ResultType rez_new = integrate(gt, f, order, jacobiFunctor);

            ResultType rez_smooth(nResult, ResultValue(0.0));
            ResultType rez_delta(nResult, ResultValue(0.0));

            for (int iResult = 0; iResult < nResult; iResult++)
            {
                ResultValue rez_smooth_part1 = rez_prev[iResult] + rez_new[iResult];
                ResultValue rez_smooth_part2 = rez_this[iResult];
                rez_smooth_part1 *= SMOOTH_FACTOR;
                rez_smooth_part2 *= (1 - 2 * SMOOTH_FACTOR);
                rez_smooth[iResult] = rez_smooth_part1 + rez_smooth_part2;
                rez_delta[iResult] = rez_new[iResult] - rez_smooth[iResult];
            }
            rez_prev = rez_this;
            rez_this = rez_new;

            // Compute error - compute smoothened error to avoid cases when two consecutive samplings give same wrong answer
            // If a vector of several quantities is integrated simultaneously, compare the relative error of the largest component
            //  with the requested relative tolerance
            error = absoluteValueResult(rez_delta, rez_this);
        }

        std::cout << "Finished recursive integral over geometry " << gt << std::endl;

        return StatInfo(order, rez_this);
    }


    template<typename ResultType>
    static ctype maxAbsoluteValueResult(ResultType & v, ResultType vmag)  {
    	ctype rez = 0;
    	for (int i = 0; i < v.size(); i++)  {
    		ctype abs_val = absoluteValue(vmag[i]);
    		ctype abs_err = absoluteValue(v[i]);

    		// Only calculate relative error if the absolute value of the error is large enough
    		// Otherwise just return absolute error
    		ctype rel_err = (abs_val > 1.0e-15) ? abs_err / abs_val : abs_err;

    		rez = std::max(rez, rel_err);
    	}
    	return rez;
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
    static ctype absoluteValue(Dune::DynamicMatrix<ValueType> & v)
    {
    	return v.frobenius_norm();
    }







};

} // Namespace Dune

#endif /** DUNE_QUADRATUREINTEGRATOR_HH **/
