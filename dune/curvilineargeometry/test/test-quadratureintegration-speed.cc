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
#include <ctime>
#include <string>
#include <functional>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/integration/quadratureintegrator.hh>


using namespace Dune;

typedef FieldVector< double, 1 > GlobalVector1D;
typedef FieldVector< double, 2 > GlobalVector2D;
typedef FieldVector< double, 3 > GlobalVector3D;

typedef std::function<GlobalVector1D(double)> GlobalFunction1D;
typedef std::function<GlobalVector2D(double)> GlobalFunction2D;
typedef std::function<GlobalVector3D(double)> GlobalFunction3D;

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
struct function1d4 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector1D & in) const  { return functionDouble::ResultType(1, sqrt(1 + pow(in[0], 3))); }  };

struct function2d1 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, 1); }  };
struct function2d8 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector2D & in) const  { return functionDouble::ResultType(1, sqrt(pow(in[0], 3) + pow(in[1], 4))); }  };

struct function3d1 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector3D & in) const  { return functionDouble::ResultType(1, 1); }  };
struct function3d8 : public functionDouble {  functionDouble::ResultType operator()(const GlobalVector3D & in) const  { return functionDouble::ResultType(1, sqrt(pow(in[0], 3) + pow(in[1], 4) + pow(in[2], 5))); }  };


struct unityJacobianFunctor
{
	typedef double               ResultValue;
	typedef std::vector<double>  ResultType;

	template <typename Coordinate>
	ResultType operator()(const Coordinate & in) const  { return ResultType(1, 1.0); }

	ResultValue zeroValue(unsigned int rezIndex) const { return 0.0; }
};


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


std::string getTime()
{
	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	return std::string(asctime (timeinfo));
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

  unsigned int test_order = 20;

  std::cout << getTime() << " doing trivial integral" << std::endl;
  for (int i = 0; i < 1000; i++) { QuadIntegrator1D::integrateStat(edgeGeometry, function1d1(),  unityJacobianFunctor(), test_order); }
  for (int i = 0; i < 1000; i++) { QuadIntegrator2D::integrateStat(faceGeometry, function2d1(),  unityJacobianFunctor(), test_order); }
  for (int i = 0; i < 1000; i++) { QuadIntegrator3D::integrateStat(elemGeometry, function3d1(),  unityJacobianFunctor(), test_order); }

  std::cout << getTime() << " finished trivial integral" << std::endl;
  std::cout << getTime() << " started hard integral" << std::endl;

  for (int i = 0; i < 1000; i++) { QuadIntegrator1D::integrateStat(edgeGeometry, function1d4(),  unityJacobianFunctor(), test_order);  }
  for (int i = 0; i < 1000; i++) { QuadIntegrator2D::integrateStat(faceGeometry, function2d8(),  unityJacobianFunctor(), test_order);  }
  for (int i = 0; i < 1000; i++) { QuadIntegrator3D::integrateStat(elemGeometry, function3d8(),  unityJacobianFunctor(), test_order);  }

  std::cout << getTime() << " finished hard integral" << std::endl;

  return 0;
}
