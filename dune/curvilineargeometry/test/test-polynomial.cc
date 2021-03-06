/*******************************************************************
 * Test of Curvilinear Element Interpolator
 * 
 * author: Aleksejs Fomins
 * date: 01.07.2014 - created
 * 
 * description:
 * Tests the functionality of the Polynomial class and compares it with the expected analytic results
 * 
 * TODO: Implement automatic testing by recording analytic results, thus providing tests of the form PASS/FAIL
 * 
 *******************************************************************/


#include <config.h>
#include <dune/common/fvector.hh>
#include <dune/curvilineargeometry/interpolation/polynomial.hh>


using namespace Dune;


// [FIXME] Rewrite all tests in quantitative form, report pass/fail
int main ()
{
  typedef PolynomialTraits<double>::Monomial  Monomial;

  std::cout << "*********************************************************************" << std::endl;
  std::cout << "Started testing 1D polynomials " << std::endl;
  std::cout << "*********************************************************************" << std::endl;

  Polynomial<double, 1> poly1D_test;
  poly1D_test += Monomial(2, 0);
  poly1D_test += Monomial(-3, 1);
  poly1D_test += Monomial(4, 2);
  poly1D_test.compactify();

  Polynomial<double, 1> poly1D_test2 = poly1D_test;
  Polynomial<double, 1> poly1D_test3 = poly1D_test.derivative(0);
  Polynomial<double, 1> poly1D_test4 = poly1D_test;
  Polynomial<double, 1> poly1D_test5 = poly1D_test;
  poly1D_test2 *= poly1D_test;
  poly1D_test4 += poly1D_test3;
  poly1D_test5 *= -5;
  Polynomial<double, 1> poly1D_test6;
  poly1D_test6 += Monomial(3,1);
  poly1D_test6 = (poly1D_test6 * poly1D_test6 * poly1D_test6 + poly1D_test6 * poly1D_test6 * 2 - poly1D_test6 * 5 + 16.5) - 3;
  poly1D_test6.compactify();


  FieldVector< double, 1 > testEval1D; testEval1D[0] = 5.0;


  std::cout << "expecting Polynomial 2 x^0 - 3 x^1 + 4 x^2, result is: " << poly1D_test.to_string() << std::endl;
  std::cout << "expecting value p(x = 5) = 87, result is: " << poly1D_test.evaluate(testEval1D) << std::endl;
  std::cout << "expecting integral over x=[0,1] to be 1.83333, result is: " << poly1D_test.integrateRefSimplex() << std::endl;
  std::cout << "expecting product with itself +4 x^0  -12 x^1  +25 x^2  -24 x^3  +16 x^4, result is: " << poly1D_test2.to_string() << std::endl;
  std::cout << "expecting derivative Polynomial -3 x^0 + 8 x^1, result is: " << poly1D_test3.to_string() << std::endl;
  std::cout << "adding poly and its derivative, expected -1 x^0 + 5 x^1 + 4 x^2, result is: " << poly1D_test4.to_string() << std::endl;
  std::cout << "multiplying poly by -5, expected -10 x^0  +15 x^1  -20 x^2, result is: " << poly1D_test5.to_string() << std::endl;
  std::cout << "algebraic abilities test, expected 13.5 -15x + 18x^2 + 27 x^3, result is: " << poly1D_test6.to_string() << std::endl;


    // Testing 2D Polynomial
  std::cout << "*********************************************************************" << std::endl;
  std::cout << "Started testing 2D polynomials " << std::endl;
  std::cout << "*********************************************************************" << std::endl;

  Polynomial<double, 2> poly2D_test;

  poly2D_test += Monomial(2.0, 0, 0);
  poly2D_test += Monomial(-3.0, 1, 0);
  poly2D_test += Monomial(5.0, 0, 1);
  poly2D_test += Monomial(2.0, 1, 1);
  poly2D_test += Monomial(-1.0, 2, 0);
  poly2D_test += Monomial(3.0, 0, 2);
  poly2D_test.compactify();

  FieldVector< double, 2 > testEval2D; testEval2D[0] = 3; testEval2D[1] = 2;
  Polynomial<double, 2> poly2D_test2 = poly2D_test;
  Polynomial<double, 2> poly2D_test3_1 = poly2D_test.derivative(0);
  Polynomial<double, 2> poly2D_test3_2 = poly2D_test.derivative(1);
  Polynomial<double, 2> poly2D_test4 = poly2D_test3_1;
  Polynomial<double, 2> poly2D_test5 = poly2D_test;
  poly2D_test2 *= poly2D_test;
  poly2D_test4 += poly2D_test3_2;
  poly2D_test5 *= -5;
  Polynomial<double, 2> poly2D_test6;
  poly2D_test6 += Monomial(3,1,2);
  poly2D_test6 = (poly2D_test6 * poly2D_test6 * poly2D_test6 + poly2D_test6 * poly2D_test6 * 2 - poly2D_test6 * 5 + 16.5) - 3;
  poly2D_test6.compactify();


  std::cout << "expecting Polynomial +2 x^0 y^0  +5 x^0 y^1  +3 x^0 y^2  -3 x^1 y^0  +2 x^1 y^1  -1 x^2 y^0, result is: " << poly2D_test.to_string() << std::endl;
  std::cout << "expecting value p(x=3, y=2) = 18, result is: " << poly2D_test.evaluate(testEval2D) << std::endl;
  std::cout << "expecting integral over triangle to be 1.58333, result is: " << poly2D_test.integrateRefSimplex() << std::endl;
  std::cout << "expecting product with itself  , result is: " << poly2D_test2.to_string() << std::endl;
  std::cout << "x-derivative of the Polynomial, expected    -3 x^0 y^0  +2 x^0 y^1  -2 x^1 y^0, result is: " << poly2D_test3_1.to_string() << std::endl;
  std::cout << "y-derivative of the Polynomial, expected    +5 x^0 y^0  +6 x^0 y^1  +2 x^1 y^0, result is: " << poly2D_test3_2.to_string() << std::endl;
  std::cout << "adding both derivatives together, expected  +2 x^0 y^0  +8 x^0 y^1  +0 x^1 y^0, result is: " << poly2D_test4.to_string() << std::endl;
  std::cout << "multiplying poly by -5, expected            -10 x^0 y^0  -25 x^0 y^1  -15 x^0 y^2  +15 x^1 y^0  -10 x^1 y^1  +5 x^2 y^0, result is: " << poly2D_test5.to_string() << std::endl;
  std::cout << "algebraic abilities test, expected 13.5 -15x^0 y^2 + 18x^2 y^4 + 27 x^3 y^6, result is: " << poly2D_test6.to_string() << std::endl;

  // Testing 3D Polynomial
  std::cout << "*********************************************************************" << std::endl;
  std::cout << "Started testing 3D polynomials " << std::endl;
  std::cout << "*********************************************************************" << std::endl;

  Polynomial<double, 3> poly3D_test;

  poly3D_test += Monomial(2.0, 0, 0, 0);
  poly3D_test += Monomial(-3.0, 1, 0, 0);
  poly3D_test += Monomial(5.0, 0, 1, 0);
  poly3D_test += Monomial(-4.0, 0, 0, 1);
  poly3D_test += Monomial(2.0, 1, 1, 1);
  poly3D_test += Monomial(-1.0, 2, 0, 1);
  poly3D_test += Monomial(2.0, 0, 1, 2);
  poly3D_test.compactify();


  FieldVector< double, 3 > testEval3D; testEval3D[0] = 3; testEval3D[1] = 2; testEval3D[2] = 5;
  Polynomial<double, 3> poly3D_test2 = poly3D_test;
  Polynomial<double, 3> poly3D_test3_1 = poly3D_test.derivative(0);
  Polynomial<double, 3> poly3D_test3_2 = poly3D_test.derivative(1);
  Polynomial<double, 3> poly3D_test3_3 = poly3D_test.derivative(2);
  Polynomial<double, 3> poly3D_test4 = poly3D_test3_1;
  Polynomial<double, 3> poly3D_test5 = poly3D_test;
  poly3D_test2 *= poly3D_test;
  poly3D_test4 += poly3D_test3_2;
  poly3D_test4 += poly3D_test3_3;
  poly3D_test5 *= -5;

  Polynomial<double, 3> poly3D_test6;
  poly3D_test6 += Monomial(3,1,2,3);
  poly3D_test6 = (poly3D_test6 * poly3D_test6 * poly3D_test6 + poly3D_test6 * poly3D_test6 * 2 - poly3D_test6 * 5 + 16.5) - 3;
  poly3D_test6.compactify();


  std::cout << "expecting Polynomial +2 x^0 y^0 z^0  -4 x^0 y^0 z^1  +5 x^0 y^1 z^0  +2 x^0 y^1 z^2  -3 x^1 y^0 z^0  +2 x^1 y^1 z^1  -1 x^2 y^0 z^1, result is: "<< poly3D_test.to_string() << std::endl;
  std::cout << "expecting value p(x=3, y=2, z=5) = 98, result is: " << poly3D_test.evaluate(testEval3D) << std::endl;
  std::cout << "expecting integral over tetrahedron to be 0.255556, result is: " << poly3D_test.integrateRefSimplex() << std::endl;
  std::cout << "expecting product with itself  , result is: " << poly3D_test2.to_string() << std::endl;
  std::cout << "x-derivative of the Polynomial, expected    -3 x^0 y^0 z^0  +2 x^0 y^1 z^1  -2 x^1 y^0 z^1, result is: " << poly3D_test3_1.to_string() << std::endl;
  std::cout << "y-derivative of the Polynomial, expected    +5 x^0 y^0 z^0  +2 x^0 y^0 z^2  +2 x^1 y^0 z^1, result is: " << poly3D_test3_2.to_string() << std::endl;
  std::cout << "z-derivative of the Polynomial, expected    -4 x^0 y^0 z^0  +4 x^0 y^1 z^1  +2 x^1 y^1 z^0  -1 x^2 y^0 z^0, result is: " << poly3D_test3_3.to_string() << std::endl;
  std::cout << "adding all derivatives together, expected   -2 x^0 y^0 z^0  +2 x^0 y^0 z^2  +6 x^0 y^1 z^1  +0 x^1 y^0 z^1  +2 x^1 y^1 z^0  -1 x^2 y^0 z^0, result is: " << poly3D_test4.to_string() << std::endl;
  std::cout << "multiplying poly by -5, expected            -10 x^0 y^0 z^0  +20 x^0 y^0 z^1  -25 x^0 y^1 z^0  -10 x^0 y^1 z^2  +15 x^1 y^0 z^0  -10 x^1 y^1 z^1  +5 x^2 y^0 z^1, result is: " << poly3D_test5.to_string() << std::endl;
  std::cout << "algebraic abilities test, expected 13.5 -15x y^2 z^3 + 18x^2 y^4 z^6 + 27 x^3 y^6 z^9, result is: " << poly3D_test6.to_string() << std::endl;

  std::cout << "*********************************************************************" << std::endl;
  std::cout << "*********************************************************************" << std::endl;

  return 0;
}
