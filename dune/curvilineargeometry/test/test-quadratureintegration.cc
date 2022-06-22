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



struct functionDouble
{
	typedef double                ResultValue;
	typedef std::vector<double>   ResultType;
	static const int RETURN_SIZE = 1;

	ResultValue zeroValue(unsigned int rezIndex) const { return 0.0; }
};

struct function1d1 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector1D & in) const  { return functionDouble::ResultType(1, 1.0); }  };
struct function1d2 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector1D & in) const  { return functionDouble::ResultType(1, in[0]); }  };
struct function1d3 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector1D & in) const  { return functionDouble::ResultType(1, in[0] * in[0] * in[0] - 3 * in[0] + 3); }  };
struct function1d4 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector1D & in) const  { return functionDouble::ResultType(1, sqrt(in[0])); }  };

struct function2d1 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, 1); }  };
struct function2d2 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, 1 + in[0]); }  };
struct function2d3 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, 1 + in[0] * in[0] + in[1] * in[1]); }  };
struct function2d4 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, in[0] * in[1] * in[1]); }  };
struct function2d5 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, sqrt(in[0] * in[1])); }  };
struct function2d6 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, 2000 * in[0] * in[0] * in[0] * in[1] * in[1] * in[1]); }  };
struct function2d7 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, 3628800 * pow(in[0], 7) * pow(in[1], 10)); }  };
struct function2d8 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, sqrt(pow(in[0], 7) * pow(in[1], 10) + 0.5)); }  };

struct function3d1 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector3D & in) const  { return functionDouble::ResultType(1, 1.0); }  };
struct function3d2 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector3D & in) const  { return functionDouble::ResultType(1, sqrt(in[0] * in[1] * in[2])); }  };
struct function3d3 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector3D & in) const  { return functionDouble::ResultType(1, sqrt(pow(in[0], 2) + pow(in[1], 2) + pow(in[2], 2))); }  };


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


template<class ctype, int mydim>
struct PolynomialFunctor
{
    typedef Polynomial<ctype, mydim> LocalPolynomial;
    typedef FieldVector< ctype, mydim > LocalCoordinate;

    static const unsigned int RETURN_SIZE = 1;
    typedef ctype                        ResultValue;
    typedef typename std::vector<ctype>  ResultType;

    LocalPolynomial p_;

    PolynomialFunctor(const LocalPolynomial & p) : p_(p) {}

    ResultType operator()(const LocalCoordinate & x) const { return ResultType(1, p_.evaluate(x)); }

    ResultValue zeroValue(unsigned int rezIndex) const { return 0.0; }
};


struct unityJacobianFunctor
{
	typedef double               ResultValue;
	typedef std::vector<double>  ResultType;

	template <typename Coordinate>
	ResultType operator()(const Coordinate & in) const  { return ResultType(1, 1.0); }

	ResultValue zeroValue(unsigned int rezIndex) const { return 0.0; }
};


template <class StatInfo>
void recursiveWrite(StatInfo statInfo)
{
	std::cout << "converged to result " << statInfo.second[0] << " at order " << statInfo.first << std::endl;
}

template <class StatInfoVec>
void statWrite(StatInfoVec statInfoVec)
{
	std::cout << "--------Writing statInfo---------" << std::endl;
	for (int i = 0; i < statInfoVec.size(); i++)
	{
		std::cout << "   ord=" << i + 1 << " np=" << statInfoVec[i].first << " rez=" << statInfoVec[i].second[0] << std::endl;
	}
	std::cout << "--------Done writing statInfo---------" << std::endl;
}

template <int dim, int NORM_TYPE>
void depthTest(int maxdim, double RELATIVE_TOLERANCE, double ACCURACY_GOAL)
{
	typedef PolynomialFunctor<double, dim>                               PolyFunctor;
	typedef Dune::QuadratureIntegrator<double, dim>                      Integrator;
	typedef typename Integrator::template Traits<PolyFunctor>::StatInfo  StatInfo;

	std::cout << "-----------Started " << dim << "D depth check-------------" << std::endl;
	Dune::GeometryType gt=Dune::GeometryTypes::simplex(dim);

	for (int i = 1; i <= maxdim; i++)
	{
		Dune::Polynomial<double, dim> thisPoly = NewtonPolynomial<double, dim>(i);

		double rez_poly = thisPoly.integrateRefSimplex();
		StatInfo rez_quad = Integrator::template integrateRecursive<unityJacobianFunctor, PolyFunctor, NORM_TYPE>(gt, PolyFunctor(thisPoly), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL);

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

  Dune::GeometryType edgeGeometry=Dune::GeometryTypes::simplex(1);
  Dune::GeometryType faceGeometry=Dune::GeometryTypes::simplex(2);
  Dune::GeometryType elemGeometry=Dune::GeometryTypes::simplex(3);

  typedef Dune::QuadratureIntegrator<double, 1> QuadIntegrator1D;
  typedef Dune::QuadratureIntegrator<double, 2> QuadIntegrator2D;
  typedef Dune::QuadratureIntegrator<double, 3> QuadIntegrator3D;
  //Dune::QuadratureIntegrator<double, 3, 1> funIntegrator3DScalar;

  double RELATIVE_TOLERANCE = 1.0e-5;
  double ACCURACY_GOAL = 1.0e-15;

  static const unsigned int NORM_TYPE = Dune::QUADRATURE_NORM_L2;


  for (int i = 0; i < 100; i++) {
	  // recursiveWrite(QuadIntegrator1D::template integrateRecursive<unityJacobianFunctor, function1d1, NORM_TYPE>(edgeGeometry, function1d1(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  // recursiveWrite(QuadIntegrator1D::template integrateRecursive<unityJacobianFunctor, function1d2, NORM_TYPE>(edgeGeometry, function1d2(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  // recursiveWrite(QuadIntegrator1D::template integrateRecursive<unityJacobianFunctor, function1d3, NORM_TYPE>(edgeGeometry, function1d3(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  // recursiveWrite(QuadIntegrator1D::template integrateRecursive<unityJacobianFunctor, function1d4, NORM_TYPE>(edgeGeometry, function1d4(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  //
	  //  recursiveWrite(QuadIntegrator2D::template integrateRecursive<unityJacobianFunctor, function2d1, NORM_TYPE>(faceGeometry, function2d1(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  //  recursiveWrite(QuadIntegrator2D::template integrateRecursive<unityJacobianFunctor, function2d2, NORM_TYPE>(faceGeometry, function2d2(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  //  recursiveWrite(QuadIntegrator2D::template integrateRecursive<unityJacobianFunctor, function2d3, NORM_TYPE>(faceGeometry, function2d3(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  //  recursiveWrite(QuadIntegrator2D::template integrateRecursive<unityJacobianFunctor, function2d4, NORM_TYPE>(faceGeometry, function2d4(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  //  recursiveWrite(QuadIntegrator2D::template integrateRecursive<unityJacobianFunctor, function2d5, NORM_TYPE>(faceGeometry, function2d5(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  //  recursiveWrite(QuadIntegrator2D::template integrateRecursive<unityJacobianFunctor, function2d6, NORM_TYPE>(faceGeometry, function2d6(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  //  recursiveWrite(QuadIntegrator2D::template integrateRecursive<unityJacobianFunctor, function2d7, NORM_TYPE>(faceGeometry, function2d7(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  //  recursiveWrite(QuadIntegrator2D::template integrateRecursive<unityJacobianFunctor, function2d8, NORM_TYPE>(faceGeometry, function2d8(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  //
	  //  recursiveWrite(QuadIntegrator3D::template integrateRecursive<unityJacobianFunctor, function3d1, NORM_TYPE>(elemGeometry, function3d1(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	  // recursiveWrite(QuadIntegrator3D::template integrateRecursive<unityJacobianFunctor, function3d2, NORM_TYPE>(elemGeometry, function3d2(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
	    recursiveWrite(QuadIntegrator3D::template integrateRecursive<unityJacobianFunctor, function3d3, NORM_TYPE>(elemGeometry, function3d3(), unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL));
  }

//  depthTest<1, NORM_TYPE>(25, RELATIVE_TOLERANCE, ACCURACY_GOAL);
//  depthTest<2, NORM_TYPE>(25, RELATIVE_TOLERANCE, ACCURACY_GOAL);
//  depthTest<3, NORM_TYPE>(25, RELATIVE_TOLERANCE, ACCURACY_GOAL);

  return 0;
}
