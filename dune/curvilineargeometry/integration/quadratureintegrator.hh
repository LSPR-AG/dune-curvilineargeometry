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




template<class ctype, int dim>
class QuadratureIntegrator {
	typedef Dune::QuadratureRule<ctype, dim>    QRule;
	typedef Dune::QuadratureRules<ctype, dim>   QRules;

	static const int MAX_INT_ORDER = 60;

public:

	typedef std::pair<int, ctype> StatInfo;  // To store pair (integrOrder, result)
	typedef typename std::vector<StatInfo>  StatInfoVec;

	QuadratureIntegrator()  {}


	// Integrates functor over reference element using quadrature of a given order
	template<class Functor>
	static ctype integrate(Dune::GeometryType gt, Functor f, int integrOrder)
	{
		  const QRule & rule = QRules::rule(gt, integrOrder);
		  if (rule.order() < integrOrder)  { DUNE_THROW(Dune::Exception,"order not available"); }

		  ctype result = 0;
		  for (typename QRule::const_iterator i=rule.begin(); i!=rule.end(); ++i)
		  {
		    double fval = f(i->position());
		    double weight = i->weight();
		    result += fval * weight;
		  }

		  return result;
	}

	// Integrates using gradually increasing quadrature order until estimated smooth relative tolerance achieved
	template<class Functor>
	static StatInfo integrateRecursive(Dune::GeometryType gt, Functor f, ctype rel_tol)
	{
		ctype SMOOTH_FACTOR = 0.15;

		ctype error = 1.0;
		ctype rez_prev = 0.0;
		ctype rez_this = 0.0;
		int order = 0;

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
			ctype rez_new = integrate(gt, f, order);
			ctype rez_smooth = SMOOTH_FACTOR * (rez_prev + rez_new) + (1 - 2 * SMOOTH_FACTOR) * rez_this;
			rez_prev = rez_this;
			rez_this = rez_new;

			// Compute error - compute smoothened error to avoid cases when two consecutive samplings give same wrong answer
			// Compute relative error unless result is close to zero, otherwise compute absolute error
			error = fabs(rez_this - rez_smooth);
			if (fabs(rez_this) > 1.0e-15)  { error /= fabs(rez_this); }
		}
		return StatInfo(order, rez_this);
	}


	// Integrates functor using all specified quadrature orders, returns vector of values
	template<class Functor>
	static StatInfoVec integrateStat(Dune::GeometryType gt, Functor f, int integrOrderMax)
	{
		StatInfoVec rez;

		for (int i = 1; i <= integrOrderMax; i++)
		{
			rez.push_back(StatInfo(
				QRules::rule(gt, i).size(),
				integrate(gt, f, i)
			));
		}

		return rez;
	}


};

} // Namespace Dune

#endif /** DUNE_QUADRATUREINTEGRATOR_HH **/