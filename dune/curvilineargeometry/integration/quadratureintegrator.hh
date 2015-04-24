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




template<class ctype, int dim, int integranddim>
class QuadratureIntegrator {
    typedef Dune::QuadratureRule<ctype, dim>    QRule;
    typedef Dune::QuadratureRules<ctype, dim>   QRules;

    static const int MAX_INT_ORDER = 60;


    template <class Geometry>
    struct IntegrationElementFunctor
    {
    	static const int mydimension = Geometry::mydimension;

    	typedef FieldVector< ctype, mydimension > LocalCoordinate;

        IntegrationElementFunctor(const Geometry & geom) : geom_(geom) {}

        double operator()(const LocalCoordinate & x) const { return geom_.integrationElement(x); }

        const Geometry & geom_;
    };








public:

    typedef Dune::FieldVector<ctype, integranddim>  IntegrandType;
    typedef std::pair<int, IntegrandType>           StatInfo;       // To store pair (integrOrder, result)
    typedef typename std::vector<StatInfo>          StatInfoVec;

    QuadratureIntegrator()  {}


    // Computes integral of a functor over an element using the provided Geometry to calculate integration elements
    template<class Geometry, class Functor>
    static IntegrandType integrate(
        const Geometry & geometry,
        const Functor & f,
        int integrOrder
    )
    {
        IntegrationElementFunctor<Geometry> detJ(geometry);
        return integrateImpl(geometry.type(), f, integrOrder, detJ);
    }


    // Computes integral of a functor over an element using the provided JacobianFunctor to calculate integration elements
    template<class Functor, class JacobiFunctor>
    static IntegrandType integrate(
        Dune::GeometryType gt,
        const Functor & f,
        int integrOrder,
        const JacobiFunctor detJ
    )
    {
        return integrateImpl(gt, f, integrOrder, detJ);
    }


    // Computes integral of a functor over an element recursively using the provided Geometry to calculate integration elements
    template<class Geometry, class Functor>
    static StatInfo integrateRecursive(
        const Geometry & geometry,
        const Functor & f,
        ctype rel_tol,
        unsigned int suggestedOrder = 1
    )
    {
        IntegrationElementFunctor<Geometry> detJ(geometry);
        return integrateRecursiveImpl(geometry.type(), f, rel_tol, suggestedOrder, detJ);
    }


    // Computes integral of a functor over an element recursively using the provided JacobianFunctor to calculate integration elements
    template<class Functor, class JacobiFunctor>
    static StatInfo integrateRecursive(
        Dune::GeometryType gt,
        const Functor & f,
        ctype rel_tol,
        const JacobiFunctor & detJ,
        unsigned int suggestedOrder = 1
    )
    {
        return integrateRecursiveImpl(gt, f, rel_tol, suggestedOrder, detJ);
    }


    // Integrates functor using all specified quadrature orders, returns vector of values
    template<class Geometry, class Functor>
    static StatInfoVec integrateStat(const Geometry & geometry, const Functor & f, int integrOrderMax)
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
    template<class Functor, class JacobiFunctor>
    static IntegrandType integrateImpl(
        Dune::GeometryType gt,
        const Functor & f,
        int integrOrder,
        const JacobiFunctor & detJ
    )
    {
          const QRule & rule = QRules::rule(gt, integrOrder);
          if (rule.order() < integrOrder)  { DUNE_THROW(Dune::Exception,"order not available"); }

          IntegrandType result(0.0);  // Assume automatic init by zero
          for (typename QRule::const_iterator i=rule.begin(); i!=rule.end(); ++i)
          {
            IntegrandType fval = f(i->position());
            ctype weight = i->weight();
            ctype detjac = detJ(i->position());
            fval *= weight * detjac;
            result += fval;
          }

          return result;
    }


    // Integrates using gradually increasing quadrature order until estimated smooth relative tolerance achieved
    template<class Functor, class JacobiFunctor>
    static StatInfo integrateRecursiveImpl(
            Dune::GeometryType gt,
            const Functor & f,
            ctype rel_tol,
            unsigned int suggestedOrder,
            const JacobiFunctor & jacobiFunctor)
    {
        ctype SMOOTH_FACTOR = 0.15;

        ctype error = 1.0;
        IntegrandType rez_prev(0.0);
        IntegrandType rez_this(0.0);
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
            IntegrandType rez_new = integrate(gt, f, order, jacobiFunctor);

            IntegrandType rez_smooth;
            {
                IntegrandType rez_smooth_part1 = rez_prev + rez_new;
                IntegrandType rez_smooth_part2 = rez_this;
                rez_smooth_part1 *= SMOOTH_FACTOR;
                rez_smooth_part2 *= (1 - 2 * SMOOTH_FACTOR);
                rez_smooth = rez_smooth_part1 + rez_smooth_part2;
            }
            rez_prev = rez_this;
            rez_this = rez_new;

            // Compute error - compute smoothened error to avoid cases when two consecutive samplings give same wrong answer
            // Compute relative error unless result is close to zero, otherwise compute absolute error
            ctype abs_val = rez_this.one_norm();
            error = (rez_this - rez_smooth).one_norm();
            if (abs_val > 1.0e-15)  { error /= abs_val; }
        }
        return StatInfo(order, rez_this);
    }









};

} // Namespace Dune

#endif /** DUNE_QUADRATUREINTEGRATOR_HH **/
