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
#include <iostream>
#include <math.h>

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include <complex>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/integration/quadratureintegrator.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>


using namespace Dune;



template<class ctype>
ctype randomNumber()  { return ctype(double(rand()) / RAND_MAX * 2); }


template <class ctype, int dim>
Dune::Polynomial<ctype, dim>  randomPolynomial(int power)
{
	typedef typename Dune::PolynomialTraits<ctype>::Monomial  Monomial;

	Dune::Polynomial<ctype, dim> rez;

	std::vector<std::vector<int> > polyPower = Dune::CurvilinearGeometryHelper::template simplexGridEnumerate<dim>(power);

	for (int i = 0; i < polyPower.size(); i++)  { rez += Monomial(randomNumber<ctype>(), polyPower[i]); }

	return rez;
}


template <class ctype, int dim>
class matrixPolyRandomFunctor
{
protected:
	typedef Dune::Polynomial<ctype, dim>  Polynomial;
	typedef std::vector< Polynomial >     PolyVec;
	typedef std::vector< PolyVec >        PolyVecVec;

public:

	typedef Dune::FieldVector<ctype, dim> GlobalCoordinate;

	typedef Dune::DynamicMatrix<ctype>   ResultValue;
	typedef std::vector<ResultValue>     ResultType;
	static const int RETURN_SIZE = 1;


	matrixPolyRandomFunctor(int nRow, int nCol, int power) :
		nRow_(nRow),
		nCol_(nCol),
		power_(power)
	{
		polyMat_ = PolyVecVec(nRow, PolyVec(nCol, Polynomial()));

		for (int i = 0; i < nRow; i++)
		{
			for (int j = 0; j < nCol; j++)
			{
				polyMat_[i][j] = randomPolynomial<ctype, dim>(power);
				polyMat_[i][j].cache();
			}
		}

	}


	ResultValue zeroValue(unsigned int rezIndex) const { return ResultValue(nRow_, nCol_, 0.0); }


	// Evaluates all matrix elements at provided coordinate by evaluating each associated polynomial
	ResultType operator()(const GlobalCoordinate & in) const  {
		ResultType matrixVecRez = ResultType(1, ResultValue(nRow_, nCol_, ctype(0.0)));

		for (int i = 0; i < nRow_; i++)
		{
			for (int j = 0; j < nCol_; j++)
			{
				matrixVecRez[0][i][j] = polyMat_[i][j].evaluate(in);
			}
		}

		return matrixVecRez;
	}


	// Integrates all polynomial matrix elements over reference simplex using analytical formula
	ResultType integrateRef() const  {
		ResultType matrixVecRez = ResultType(1, ResultValue(nRow_, nCol_, ctype(0.0)));

		for (int i = 0; i < nRow_; i++)
		{
			for (int j = 0; j < nCol_; j++)
			{
				matrixVecRez[0][i][j] = polyMat_[i][j].integrateRefSimplex();
			}
		}

		return matrixVecRez;
	}



protected:

	int nRow_;
	int nCol_;
	int power_;

	std::vector<std::vector< Polynomial > >  polyMat_;
};


template <class ctype, int dim>
class matrixPolySqrtRandomFunctor
{
	typedef matrixPolyRandomFunctor<ctype, dim>    PolyFunctor;


	typedef typename PolyFunctor::GlobalCoordinate  GlobalCoordinate;

public:

	typedef typename PolyFunctor::ResultType        ResultType;
	typedef typename PolyFunctor::ResultValue       ResultValue;
	static const int RETURN_SIZE = 6;


	matrixPolySqrtRandomFunctor(int nRow, int nCol, int power) :
		nRow_(nRow),
		nCol_(nCol)
	{
		for (int i = 0; i < RETURN_SIZE; i++)  { A.push_back(PolyFunctor(nRow, nCol, power)); }
	}


	ResultValue zeroValue(unsigned int rezIndex) const { return A[0].zeroValue(rezIndex); }


	ResultType operator()(const GlobalCoordinate & in) const  {

		ResultType rez;

		for (int iSize = 0; iSize < RETURN_SIZE; iSize++)
		{
			rez.push_back(A[iSize](in)[0]);

			for (int i = 0; i < nRow_; i++)
			{
				for (int j = 0; j < nCol_; j++)
				{
					rez[iSize][i][j] = sqrt(rez[iSize][i][j]);
				}
			}
		}

		return rez;
	}

protected:
	std::vector<PolyFunctor> A;

	int nRow_;
	int nCol_;
};


struct unityJacobianFunctor
{
	typedef double               ResultValue;
	typedef std::vector<double>  ResultType;

	template <typename Coordinate>
	ResultType operator()(const Coordinate & in) const  { return ResultType(1, 1.0); }

	ResultValue zeroValue(unsigned int rezIndex) const { return 0.0; }
};


template <class Matrix>
double matrixError(Matrix & A, Matrix & B)
{
	Matrix C = A;
	C -= B;

	return C.frobenius_norm();
}


template<class ctype, int dim>
void powerTest(unsigned int N_TEST)
{
	const unsigned int MAX_MATRIX_DIM = 30;
//	Dune::GeometryType simplexGeometry;   simplexGeometry.makeSimplex(dim);
	Dune::GeometryType simplexGeometry=Dune::GeometryTypes::simplex(dim);

	typedef matrixPolyRandomFunctor<ctype, dim>      PolyFunctor;
	typedef Dune::QuadratureIntegrator<ctype, dim>   PolyQuadIntegrator;
	typedef typename PolyFunctor::ResultValue        Matrix;

	std::cout << "Performing test for dimension = " << dim << std::endl;

	for (int iPower = 1; iPower < 20; iPower++)
	{
		double err_max = 0;

		for (int iTest = 0; iTest < N_TEST; iTest++)
		{
			unsigned int nRow = 1 + (rand() % MAX_MATRIX_DIM);
			unsigned int nCol = 1 + (rand() % MAX_MATRIX_DIM);

			PolyFunctor         polyFunctor(nRow, nCol, iPower);
			PolyQuadIntegrator  integrator;

			Matrix rez  = integrator.integrate(simplexGeometry, polyFunctor, iPower, unityJacobianFunctor())[0];
			Matrix test = polyFunctor.integrateRef()[0];

			double err_this = matrixError(rez, test);
			err_max = std::max(err_max, err_this);
		}

		std::cout << "-- Computed test for power=" << iPower << " max error=" << err_max << std::endl;
	}
}


template<class ctype, int dim>
void powerTestRecursive(unsigned int N_TEST)
{
	const unsigned int MAX_MATRIX_DIM = 30;
	const double   RELATIVE_TOLERANCE = 1.0e-10;
	const double   ACCURACY_GOAL = 1.0e-20;
	const unsigned int NORM_TYPE = Dune::QUADRATURE_NORM_L2;

//	Dune::GeometryType simplexGeometry;   simplexGeometry.makeSimplex(dim);
	Dune::GeometryType simplexGeometry=Dune::GeometryTypes::simplex(dim);

	typedef matrixPolyRandomFunctor<ctype, dim>      PolyFunctor;
	typedef Dune::QuadratureIntegrator<ctype, dim>   PolyQuadIntegrator;
	typedef typename PolyFunctor::ResultValue        Matrix;

	std::cout << "Performing test for dimension = " << dim << std::endl;

	for (int iPower = 1; iPower < 20; iPower++)
	{
		double err_max = 0;

		for (int iTest = 0; iTest < N_TEST; iTest++)
		{
			unsigned int nRow = 1 + (rand() % MAX_MATRIX_DIM);
			unsigned int nCol = 1 + (rand() % MAX_MATRIX_DIM);

			PolyFunctor         polyFunctor(nRow, nCol, iPower);

			// Note that first we extract result vector from the return value, then the result from the vector
			Matrix rez  = PolyQuadIntegrator::template integrateRecursive<unityJacobianFunctor, polyFunctor, NORM_TYPE>(simplexGeometry, polyFunctor, unityJacobianFunctor(), RELATIVE_TOLERANCE, ACCURACY_GOAL, iPower).second[0];
			Matrix test = polyFunctor.integrateRef()[0];

			double err_this = matrixError(rez, test);
			err_max = std::max(err_max, err_this);
		}

		std::cout << "-- Computed test for power=" << iPower << " max error=" << err_max << std::endl;
	}
}


template<class ctype, int dim>
void powerTestRootRecursive(unsigned int N_TEST)
{
	const unsigned int MAX_MATRIX_DIM = 30;
	const double   REL_TOL = 1.0e-5;          // For non-polynomial integrals use higher tolerance

//	Dune::GeometryType simplexGeometry;   simplexGeometry.makeSimplex(dim);
	Dune::GeometryType simplexGeometry=Dune::GeometryTypes::simplex(dim);

	typedef matrixPolySqrtRandomFunctor<ctype, dim>  PolySqrtFunctor;
	typedef Dune::QuadratureIntegrator<ctype, dim>   PolyQuadIntegrator;
	typedef typename PolySqrtFunctor::ResultValue    Matrix;

	std::cout << "Performing test for dimension = " << dim << std::endl;

	for (int iPower = 1; iPower < 10; iPower++)
	{
		double err_max = 0;

		for (int iTest = 0; iTest < N_TEST; iTest++)
		{
			unsigned int nRow = 20; //1 + (rand() % MAX_MATRIX_DIM);
			unsigned int nCol = 20; //1 + (rand() % MAX_MATRIX_DIM);

			PolySqrtFunctor     polySqrtFunctor(nRow, nCol, iPower);
			PolyQuadIntegrator  integrator;

			integrator.integrateRecursive(simplexGeometry, polySqrtFunctor, REL_TOL, unityJacobianFunctor(), 1);  // Simply ignore the result for now
		}
	}
}





int main ()
{
	typedef std::complex<double> Complex;

	/* initialize random seed: */
	srand (time(NULL));

	std::cout << "initialized" << std::endl;
	const unsigned int N_TEST = 10;

	powerTest<double, 1>(N_TEST);
	powerTest<double, 2>(N_TEST);
	powerTest<double, 3>(N_TEST);

	//powerTestRecursive<double, 1>(N_TEST);
	//powerTestRecursive<double, 2>(N_TEST);
	//powerTestRecursive<double, 3>(N_TEST);

	//powerTestRootRecursive<double, 1>(1);
	//powerTestRootRecursive<double, 2>(1);
	//powerTestRootRecursive<double, 3>(1);

	// Still some trickery is needed to make the complex tests compile
	//powerTest<Complex, 1>(N_TEST);
	//powerTest<Complex, 2>(N_TEST);
	//powerTest<Complex, 3>(N_TEST);

	return 0;
}
