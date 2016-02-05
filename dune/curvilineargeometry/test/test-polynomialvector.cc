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

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include <config.h>

#include <dune/common/fvector.hh>
#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/interpolation/polynomialfieldvector.hh>


using namespace Dune;


/*********************************************************************************************************/
/**    Auxiliary Functions
/*********************************************************************************************************/

/** \brief Generate random numbers */
double randomReal(double a, double b) { return a + (b - a)*(double(rand()) / RAND_MAX); }
double randomInt(int a, int b)        { return a + rand() % (b - a + 1); }

/** \brief Generate random vector over a cube */
template <class ctype, int dim>
FieldVector<ctype, dim> randomVec(ctype a, ctype b) {
	FieldVector<ctype, dim> rez;
	for (int i = 0; i < dim; i++)  { rez[i] = randomReal(a, b); }
	return rez;
}

/** \brief Generate random unsigned integer vector with bound sum
 *   [TODO] There is likely a much faster way to do this, but for small powers its not too bad
 * */
template <int dim>
std::vector<int> randomIntVec(int maxSum) {
	int sum = randomInt(0, maxSum);
	std::vector<int> rez(dim, 0);
	for (int i = 0; i < sum; i++)  {  rez[randomInt(0, dim-1)]++; }
	return rez;
}


/** \brief Generate random polynomial
 * \param[in]  maxOrd      maximal power for each entire monomial (sum of powers)
 * \param[in]  len         Fixed number of terms in polynomial (unless some terms randomly have the same base)
 * */
template <class ctype, int mydim>
Polynomial<ctype, mydim>  randomPoly(int maxOrd, int len, ctype mag) {
	Polynomial<ctype, mydim> rez;

	for (int i = 0; i < len; i++) {
		typename Polynomial<ctype, mydim>::Monomial  mon(randomReal(-mag, mag), randomIntVec<mydim>(maxOrd));
		rez += mon;
	}

	// Note: In the polynomial class, the compactification does not automatically happen for +=monomial, because it is not worth it for 1 monomial,
	//       However, it is definitely worth it for an entire set of monomials
	rez.compactify();
	return rez;
}


/** \brief Generate random vector of polynomials
 * \param[in]  maxOrd      maximal power for each entire monomial (sum of powers)
 * \param[in]  len         Fixed number of terms in polynomial (unless some terms randomly have the same base)
 * */
template <class ctype, int mydim, int cdim>
PolynomialFieldVector<ctype, mydim, cdim>  randomPolyVec(int maxOrd, int len, ctype mag) {
	PolynomialFieldVector<ctype, mydim, cdim> rez;
	for (int i = 0; i < cdim; i++)  { rez[i] = randomPoly<ctype, mydim>(maxOrd, len, mag); }
	return rez;
}


// This structure is declared for convenience, such that the template arguments are used only once
template<class ctype, int mydim, int cdim>
struct MapWrapper {
	typedef Polynomial<ctype, mydim>                    LocalPoly;
	typedef typename LocalPoly::Monomial                Monomial;
	typedef FieldVector<ctype, mydim>                   LocalCoordinate;
	typedef FieldVector<ctype, cdim>                    GlobalCoordinate;
	typedef PolynomialFieldVector <ctype, mydim, cdim>  PolyVec;
	static  LocalPoly         randomPoly(int maxOrd, int len, ctype mag)     { return ::randomPoly<ctype, mydim>(maxOrd, len, mag); }
	static  PolyVec           randomPolyVec(int maxOrd, int len, ctype mag)  { return ::randomPolyVec<ctype, mydim, cdim>(maxOrd, len, mag); }
	static  LocalCoordinate   randomLocal(ctype a, ctype b)                  { return randomVec<ctype, mydim>(a, b); }
	static  GlobalCoordinate  randomGlobal(ctype a, ctype b)                 { return randomVec<ctype, cdim>(a, b); }
};






/*********************************************************************************************************/
/**    Tests
/*********************************************************************************************************/


/** \brief Test for nTest random polynomial vectors, that the Vector Calculus identity div(curl(x)) always holds
 *         The result should be zero-polynomial, thus is verified by considering the magnitude of the largest monomial prefactor */
template <class ctype, int mydim, int cdim>
bool testDivCurl(int nTest) {
	typedef MapWrapper<ctype, mydim, cdim>  Map;

	const int    MAX_POLY_ORDER = 10;
	const int    MAX_POLY_N_MONOMIAL = 10;
	const ctype  MAX_POLY_MAGNITUDE = 100.0;

	double maxErr = 0.0;
	std::cout << "[...] Running " << nTest << " tests of div(curl(x)) = 0";
	for (int i = 0; i < nTest; i++) {
		typename Map::PolyVec   myvec1  = Map::randomPolyVec(MAX_POLY_ORDER, MAX_POLY_N_MONOMIAL, MAX_POLY_MAGNITUDE);
		typename Map::PolyVec   rot1    = myvec1.rot();
		typename Map::LocalPoly divrot1 = rot1.div();
		//std::cout << std::endl << divrot1.to_string();

		maxErr = std::max(maxErr, divrot1.magnitude());
	}
	std::cout << "  ::: Max_error=" << maxErr << std::endl;

	return true;
}


/** \brief Test for nTest random polynomial vectors, that the Vector Calculus identity curl(grad(x)) always holds
 *         The result should be zero-polynomial, thus is verified by considering the magnitude of the largest monomial prefactor */
template <class ctype, int mydim, int cdim>
bool testCurlGrad(int nTest) {
	typedef MapWrapper<ctype, mydim, cdim>  Map;

	const int    MAX_POLY_ORDER = 10;
	const int    MAX_POLY_N_MONOMIAL = 10;
	const ctype  MAX_POLY_MAGNITUDE = 100.0;

	double maxErr = 0.0;
	std::cout << "[...] Running " << nTest << " tests of curl(grad(x)) = 0";
	for (int i = 0; i < nTest; i++) {
		typename Map::LocalPoly  p        = Map::randomPoly(MAX_POLY_ORDER, MAX_POLY_N_MONOMIAL, MAX_POLY_MAGNITUDE);
		typename Map::PolyVec    grad1    = Map::PolyVec::grad(p);
		typename Map::PolyVec    rotgrad1 = grad1.rot();

		for (int j = 0; j < cdim; j++)  { maxErr = std::max(maxErr, rotgrad1[j].magnitude()); }
		//std::cout << std::endl << rotgrad1.to_string();
	}
	std::cout << "  ::: Max_error=" << maxErr << std::endl;

	return true;
}


/** \brief Test for nTest random pairs of polynomial vectors the functionality of substituting one into another
 * The composite map is defined by
 *
 *    p'(x) = p(f(x)),
 *
 * The test computes p'(x) for composite map, also  r = f(x), and calculates the error of   |p'(x) - p(r)|
 * To do that, the error is evaluated for a set of random sample points
 *
 * [FIXME] Known bug:
 *   Assume have 2 random polynomial vectors of order 10. Plug one into another
 *   If all pre-factors are limited to fabs(a) <= 1, then the relative error is very small
 *   However, if all pre-factors are limited to fabs(a) <= 100, the error becomes very large
 *
 * */
template <class ctype, int mydim, int cdim, int mydim1>
bool testComposite(int nTest)
{
	typedef MapWrapper<ctype, mydim1, mydim>  Map1;
	typedef MapWrapper<ctype, mydim,  cdim>   Map2;
	typedef MapWrapper<ctype, mydim1, cdim>   MapComposite;

	const int    MAX_POLY_ORDER = 5;
	const int    MAX_POLY_N_MONOMIAL = 15;
	const ctype  MAX_POLY_MAGNITUDE = 1.0;

	std::cout << "[...] Running " << nTest << " tests of map composition for (mydim, cdim, mydim1)=(" << mydim << ", " << cdim << ", " << mydim1 << ")";

	double maxErr = 0.0;

	for (int iTest = 0; iTest < nTest; iTest++) {
		typename Map1::PolyVec         myvec1 = Map1::randomPolyVec(MAX_POLY_ORDER, MAX_POLY_N_MONOMIAL, MAX_POLY_MAGNITUDE);
		typename Map2::PolyVec         myvec2 = Map2::randomPolyVec(MAX_POLY_ORDER, MAX_POLY_N_MONOMIAL, MAX_POLY_MAGNITUDE);
		typename MapComposite::PolyVec myvec3 = myvec2.composite(myvec1);

		for (int i = 0; i < 100; i++) {
			typename Map1::LocalCoordinate           x       = Map1::randomLocal(-1.0, 1.0);   // NOTE: Choose a reasonable local coordinate, that does not grow with incr polynomial power
			typename Map1::GlobalCoordinate          r       = myvec1.evaluate(x);
			typename Map2::GlobalCoordinate          p       = myvec2.evaluate(r);
			typename MapComposite::GlobalCoordinate  p_prim  = myvec3.evaluate(x);

			ctype relErr = (p - p_prim).two_norm() / p.two_norm();

			//std::cout << "  ::: Rel_error=" << relErr << std::endl;

			maxErr = std::max(maxErr, relErr);
		}
	}
	std::cout << "  ::: Max_error=" << maxErr << std::endl;

	return true;
}








// [FIXME] Rewrite all tests in quantitative form, report pass/fail
int main ()
{
	srand (time(NULL));

	testDivCurl<double, 3, 3>(100);
	testCurlGrad<double, 3, 3>(100);

	testComposite<double, 3, 2, 1>(100);
	testComposite<double, 3, 2, 2>(100);
	testComposite<double, 3, 3, 1>(100);
	testComposite<double, 3, 3, 2>(100);

	return 0;
}
