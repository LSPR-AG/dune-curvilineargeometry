/*******************************************************************
 * Test of Numerical Recursive Interpolation Integrator
 * 
 * author: Aleksejs Fomins
 * date: 01.08.2014 - created
 * 
 * description:
 * Integrates a set of functions in 1D and 2D. Compares result to analytic or known numeric result
 * 
 * TODO: Record integration time
 * TODO: Implement automatic testing by evaluating analytic results, thus providing tests of the form PASS/FAIL
 * 
 *******************************************************************/

#include <config.h>
#include<iostream>
#include <math.h>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/integration/quadratureintegrator.hh>


using namespace Dune;

typedef FieldVector< double, 1 > GlobalVector1D;
typedef FieldVector< double, 2 > GlobalVector2D;
typedef FieldVector< double, 3 > GlobalVector3D;

typedef Function<GlobalVector1D, double> GlobalFunction1D;
typedef Function<GlobalVector2D, double> GlobalFunction2D;
typedef Function<GlobalVector3D, double> GlobalFunction3D;

typedef Dune::Polynomial<double, 1>  Polynomial1D;
typedef Dune::Polynomial<double, 2>  Polynomial2D;
typedef Dune::Polynomial<double, 3>  Polynomial3D;

typedef Dune::PolynomialTraits<double>::Monomial  Monomial;

struct function1d1 {   double operator()(const GlobalVector1D & in) const  { return 1; }  };
struct function1d2 {   double operator()(const GlobalVector1D & in) const  { return in[0]; }  };
struct function1d3 {   double operator()(const GlobalVector1D & in) const  { return in[0] * in[0] * in[0] - 3 * in[0] + 3; }  };
struct function1d4 {   double operator()(const GlobalVector1D & in) const  { return sqrt(in[0]); }  };

struct function2d1 {   double operator()(const GlobalVector2D & in) const  { return 1; }  };
struct function2d2 {   double operator()(const GlobalVector2D & in) const  { return 1 + in[0]; }  };
struct function2d3 {   double operator()(const GlobalVector2D & in) const  { return 1 + in[0] * in[0] + in[1] * in[1]; }  };
struct function2d4 {   double operator()(const GlobalVector2D & in) const  { return in[0] * in[1] * in[1]; }  };
struct function2d5 {   double operator()(const GlobalVector2D & in) const  { return sqrt(in[0] * in[1]); }  };
struct function2d6 {   double operator()(const GlobalVector2D & in) const  { return 2000 * in[0] * in[0] * in[0] * in[1] * in[1] * in[1]; }  };
struct function2d7 {   double operator()(const GlobalVector2D & in) const  { return 3628800 * pow(in[0], 7) * pow(in[1], 10); }  };
struct function2d8 {   double operator()(const GlobalVector2D & in) const  { return sqrt(pow(in[0], 7) * pow(in[1], 10) + 0.5); }  };

template<class ctype, int dim>
Dune::Polynomial<ctype, dim> NewtonPolynomial(int power)
{
	assert((dim > 0)&&(dim <= 3));
	assert((power > 0)&&(power <= 100));

	typedef Dune::Polynomial<ctype, dim>  Poly;

	// Construct base polynomial (1 + x + y + z) or its analogues in other dimensions
	std::vector<int> zeroVec(dim, 0);
	Poly base(Monomial(1.0, zeroVec));
	for (int i = 0; i < dim; i++)
	{
		std::vector<int> eigenVec(dim, 0);
		eigenVec[i] = 1.0;
		base += Monomial(1.0, eigenVec);
	}

	// Multiply polynomial to required power
	Poly rez = base;
	for (int i = 1; i < power; i++)  { rez *= base; }

	return rez;
}

template<class ct, int mydim>
struct PolynomialFunctor
{
    typedef Polynomial<ct, mydim> LocalPolynomial;
    typedef FieldVector< ct, mydim > LocalCoordinate;
    typedef FieldVector< ct, 1 >     ResultType;

    LocalPolynomial p_;

    PolynomialFunctor(const LocalPolynomial & p) : p_(p) {}

    ResultType operator()(const LocalCoordinate & x) const { return ResultType(p_.evaluate(x)); }
};


template <class StatInfo>
void recursiveWrite(StatInfo statInfo)
{
	std::cout << "converged to result " << statInfo.second << " at order " << statInfo.first << std::endl;
}

template <class StatInfoVec>
void statWrite(StatInfoVec statInfoVec)
{
	std::cout << "--------Writing statInfo---------" << std::endl;
	for (int i = 0; i < statInfoVec.size(); i++)
	{
		std::cout << "   ord=" << i + 1 << " np=" << statInfoVec[i].first << " rez=" << statInfoVec[i].second << std::endl;
	}
	std::cout << "--------Done writing statInfo---------" << std::endl;
}

template <int dim>
void depthTest(int maxdim, double rec_tol)
{
	typedef Dune::QuadratureIntegrator<double, dim, 1> Integrator;
	typedef typename Integrator::StatInfo              StatInfo;

	std::cout << "-----------Started " << dim << "D depth check-------------" << std::endl;
	Dune::GeometryType entityGeometry;   entityGeometry.makeSimplex(dim);
	Integrator funIntegrator;

	for (int i = 1; i <= maxdim; i++)
	{
		Dune::Polynomial<double, dim> thisPoly = NewtonPolynomial<double, dim>(i);

		double rez_poly = thisPoly.integrateRefSimplex();
		StatInfo rez_quad = funIntegrator.integrateRecursive(entityGeometry, PolynomialFunctor<double, dim>(thisPoly), rec_tol);

		double err = fabs(rez_poly - rez_quad.second[0]);
		if (fabs(rez_poly) > 1.0e-15)  { err /= rez_poly; }

		//std::cout << "Integrating polynomial " << thisPoly.to_string() << std::endl;
		std::cout << "effectve relative error = " << err;
		std::cout << " analytic integral=" << rez_poly;
		std::cout << " quadrature integral ";
		recursiveWrite(rez_quad);
	}
}


int main ()
{

  std::cout << "initialized" << std::endl;

  Dune::GeometryType edgeGeometry;   edgeGeometry.makeSimplex(1);
  Dune::GeometryType faceGeometry;   faceGeometry.makeSimplex(2);
  Dune::GeometryType elemGeometry;   elemGeometry.makeSimplex(3);

  Dune::QuadratureIntegrator<double, 1, 1> funIntegrator1DScalar;
  Dune::QuadratureIntegrator<double, 2, 1> funIntegrator2DScalar;
  Dune::QuadratureIntegrator<double, 3, 1> funIntegrator3DScalar;

  double rec_tol = 1.0e-5;

  recursiveWrite(funIntegrator1DScalar.integrateRecursive(edgeGeometry, function1d1(), rec_tol));
  recursiveWrite(funIntegrator1DScalar.integrateRecursive(edgeGeometry, function1d2(), rec_tol));
  recursiveWrite(funIntegrator1DScalar.integrateRecursive(edgeGeometry, function1d3(), rec_tol));
  recursiveWrite(funIntegrator1DScalar.integrateRecursive(edgeGeometry, function1d4(), rec_tol));

  recursiveWrite(funIntegrator2DScalar.integrateRecursive(faceGeometry, function2d1(), rec_tol));
  recursiveWrite(funIntegrator2DScalar.integrateRecursive(faceGeometry, function2d2(), rec_tol));
  recursiveWrite(funIntegrator2DScalar.integrateRecursive(faceGeometry, function2d3(), rec_tol));
  recursiveWrite(funIntegrator2DScalar.integrateRecursive(faceGeometry, function2d4(), rec_tol));
  recursiveWrite(funIntegrator2DScalar.integrateRecursive(faceGeometry, function2d5(), rec_tol));
  recursiveWrite(funIntegrator2DScalar.integrateRecursive(faceGeometry, function2d6(), rec_tol));
  recursiveWrite(funIntegrator2DScalar.integrateRecursive(faceGeometry, function2d7(), rec_tol));
  recursiveWrite(funIntegrator2DScalar.integrateRecursive(faceGeometry, function2d8(), rec_tol));

  depthTest<1>(25, rec_tol);
  depthTest<2>(25, rec_tol);
  depthTest<3>(25, rec_tol);


  /*
  statWrite(funIntegrator1D.integrateStat(edgeGeometry, function1d1(), 10));
  statWrite(funIntegrator1D.integrateStat(edgeGeometry, function1d3(), 10));
  statWrite(funIntegrator1D.integrateStat(edgeGeometry, function1d4(), 61));

  statWrite(funIntegrator2D.integrateStat(faceGeometry, function2d1(), 10));
  statWrite(funIntegrator2D.integrateStat(faceGeometry, function2d3(), 10));
  statWrite(funIntegrator2D.integrateStat(faceGeometry, function2d4(), 10));
  statWrite(funIntegrator2D.integrateStat(faceGeometry, function2d7(), 30));
  */



  return 0;
}
