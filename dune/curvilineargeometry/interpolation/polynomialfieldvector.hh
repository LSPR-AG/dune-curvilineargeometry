/*******************************************************************
 * Polynomial Field Vector Class
 * 
 * author: Aleksejs Fomins
 * date: 18.12.2014 - created.
 * 
 * description:
 * Provides the implementation of a polynomial vector,
 * 
 * 
 * [FIXME] Throw errors if the dimensions or lengths of polynomials/monomials mismatch
 * 
 *******************************************************************/



#ifndef DUNE_POLYNOMIAL_HH
#define DUNE_POLYNOMIAL_HH


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <config.h>
#include <dune/common/fvector.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>



namespace Dune {


template<class ctype, int mydim, int cdim>
class PolynomialFieldVector {

	typedef PolynomialFieldVector<ctype, mydim, cdim>  This;

  typedef Polynomial<ctype, mydim> LocalPolynomial;
  typedef FieldVector< ctype, mydim > LocalCoordinate;
  typedef FieldVector< ctype, cdim >  GlobalCoordinate;

public:
  LocalPolynomial poly_[cdim];


  // Operator Brackets

  // Operator Brackets-Assign

  // Operator Add

  // Operator Subtract

  // Operator Multiply by scalar

  GlobalCoordinate evaluate(LocalCoordinate x)
  {
	  GlobalCoordinate rez;
	  for (int i = 0; i < cdim; i++)  { rez[i] = poly_[i].evaluate(x); }
	  return rez;
  }

  LocalPolynomial mvDot(This p)
  {
	  LocalPolynomial rez;
	  for (int i = 0; i < cdim; i++)  { rez[i].mergeTo(poly_[i] * p[i]); }
	  return rez;
  }

  This mvCross(This p)  {}

  LocalPolynomial div()
  {
	  LocalPolynomial rez;
	  for (int i = 0; i < cdim; i++)  { rez[i].mergeTo(poly_[i].derivative(i)); }
	  return rez;
  }

  This rot()  {}







};

} // namespace Dune

#endif // #ifndef DUNE_POLYNOMIAL_HH
