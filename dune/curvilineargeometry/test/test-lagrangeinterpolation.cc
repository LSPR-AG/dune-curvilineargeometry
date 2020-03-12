/*******************************************************************
 * Test of Curvilinear Element Interpolator
 * 
 * author: Aleksejs Fomins
 * date: 01.09.2014 - created
 * 
 * description:
 * Tests a set of analytic linear and non-linear geometries in 1D,2D and 3D. Compares real analytic, interpolated numeric and interpolated analytic values of the local-global map evaluated at interpolatory points and at some random points within the element. All three must match with very low tolerance, if the interpolatory order is sufficient for the requested geometry
 * 
 *******************************************************************/


#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <config.h>

#include <dune/common/fvector.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/interpolation/lagrangeinterpolator.hh>


using namespace Dune;

static const int DIM1D = 1;
static const int DIM2D = 2;
static const int DIM3D = 3;

typedef FieldVector< double, DIM1D > FieldVector1D;
typedef FieldVector< double, DIM2D > FieldVector2D;
typedef FieldVector< double, DIM3D > FieldVector3D;

typedef std::vector<FieldVector1D> FieldVectorVector1D;
typedef std::vector<FieldVector2D> FieldVectorVector2D;
typedef std::vector<FieldVector3D> FieldVectorVector3D;


// [TODO] Move all randomness to testhelperrandom.hh
double randomReal(double a, double b) { return a + (b - a)*(double(rand()) / RAND_MAX); }

// Constructs a random grid over a simplex with a given number of samples
template <int dim>
std::vector<FieldVector< double, dim > > sampleRandom(int sampleNo) {
    std::vector<FieldVector< double, dim > > tmprez;
    return tmprez;
}

// Generates set of coordinates randomly (uniformly) sampled over an edge
template <>
FieldVectorVector1D sampleRandom<DIM1D>(int sampleNo) {
    FieldVectorVector1D rez;

    for (int i = 0; i < sampleNo; i++)
    {
        FieldVector1D tmp;
        tmp[0] = randomReal(0,1);
        rez.push_back(tmp);
    }
    return rez;
}

// Generates set of coordinates randomly (uniformly) sampled over a triangle
template <>
FieldVectorVector2D sampleRandom<DIM2D>(int sampleNo) {
    FieldVectorVector2D rez;

    for (int i = 0; i < sampleNo; i++)
    {
        FieldVector2D tmp;
        tmp[0] = randomReal(0,1);
        tmp[1] = randomReal(0,1);

        // If outside unit triangle, map back in
        if (tmp[0] + tmp[1] > 1) { tmp[0] = 1 - tmp[0]; tmp[1] = 1 - tmp[1]; }

        rez.push_back(tmp);
    }
    return rez;
}

// Generates set of coordinates randomly (uniformly) sampled over a tetrahedron
template <>
FieldVectorVector3D sampleRandom<DIM3D>(int sampleNo) {
    FieldVectorVector3D rez;

    for (int i = 0; i < sampleNo; i++)
    {
        FieldVector3D tmp;
        tmp[0] = randomReal(0,1);
        tmp[1] = randomReal(0,1);
        tmp[1] = randomReal(0,1);

        // Redo while outside tetrahedron

        while (tmp[0] + tmp[1] + tmp[2] > 1) {
            tmp[0] = randomReal(0,1);
            tmp[2] = randomReal(0,1);
            tmp[3] = randomReal(0,1);
        }

        rez.push_back(tmp);
    }
    return rez;
}




// Evaluates a given vector function for a given vector of input coordinates, returns a vector of output coordinates
template <typename Interpolator>
struct FunctionEvaluator
{
	typedef std::vector<typename Interpolator::LocalCoordinate >   LocalCoordinateSet;
	typedef std::vector<typename Interpolator::GlobalCoordinate >  GlobalCoordinateSet;

	template<typename Functor>
	static GlobalCoordinateSet eval(const LocalCoordinateSet & localsample, const Functor & f)
	{
		GlobalCoordinateSet rez;
	    for (int i = 0; i < localsample.size(); i++) {rez.push_back(f(localsample[i])); }
	    return rez;
	}
};

// Evaluates an elementInterpolator for a given vector of input coordinates, returns a vector of output coordinates
template <typename Interpolator>
struct InterpolatorEvaluator
{
	typedef std::vector<typename Interpolator::LocalCoordinate >   LocalCoordinateSet;
	typedef std::vector<typename Interpolator::GlobalCoordinate >  GlobalCoordinateSet;

	static GlobalCoordinateSet eval(const LocalCoordinateSet & localsample, const Interpolator & interp)
	{
		GlobalCoordinateSet rez;
	    for (int i = 0; i < localsample.size(); i++) {rez.push_back(interp.global(localsample[i])); }
	    return rez;
	}
};

// Evaluates a vector of polynomials for a given vector of input coordinates, returns a vector of output coordinates
template <typename Interpolator>
struct PolynomialCoordEvaluator
{
	typedef typename  Interpolator::GlobalCoordinate               GlobalCoordinate;
	typedef std::vector<typename Interpolator::LocalCoordinate >   LocalCoordinateSet;
	typedef std::vector<typename Interpolator::GlobalCoordinate >  GlobalCoordinateSet;
	typedef typename Interpolator::PolynomialGlobalCoordinate      PolynomialGlobalCoordinate;

	static GlobalCoordinateSet eval(const LocalCoordinateSet & localsample, const PolynomialGlobalCoordinate & polynomialvector)
	{
		GlobalCoordinateSet rez;
	    for (int i = 0; i < localsample.size(); i++) {
	    	GlobalCoordinate tmp;
	        for (int j = 0; j < Interpolator::coorddimension; j++)  { tmp[j] = polynomialvector[j].evaluate(localsample[i]); }
	        rez.push_back(tmp);
	    }
	    return rez;
	}
};



// Initializes a field vector, such that it can be done in one line in main code
// [TODO] Move all fieldvector stuff to testhelpervector.hh
FieldVector1D initFieldVector(double x) {
    FieldVector1D rez;
    rez[0] = x;
    return rez;
}

FieldVector2D initFieldVector(double x, double y) {
    FieldVector2D rez;
    rez[0] = x;
    rez[1] = y;
    return rez;
}

FieldVector3D initFieldVector(double x, double y, double z) {
    FieldVector3D rez;
    rez[0] = x;
    rez[1] = y;
    rez[2] = z;
    return rez;
}

// Finds the 2-norm distance between two coordinate vectors
template<int dim>
double norm2sqVector(
        std::vector<FieldVector< double, dim >> & xVec1,
        std::vector<FieldVector< double, dim >> & xVec2
  ) {
    double rez = 0;
    for (int i = 0; i < xVec1.size(); i++) {
        rez += (xVec1[i] - xVec2[i]).two_norm2();
    }
    return rez;
}

// Prints out a given vector of coordinates
template<int dim >
void print(std::vector<FieldVector< double, dim >> xVec) {
    for (int i = 0; i < xVec.size(); i++)
    {
        std::cout << "(";
        for (int j = 0; j < xVec[i].size(); j++)
        {
            std::cout << xVec[i][j] << " ";
        }
        std::cout << ")" << std::endl;
    }
    std::cout << std::endl;
}





// Runs a test to compare how well the exact geometry given by f
// is approximated by Lagrange interpolation. Also checks if
// analytical and numerical lagrange interpolation give the same
// result which they should

// TODO: Implement Real-Numerical in LagrangeInterpolation, then add to this test
template <int mydim, int cdim, typename Functor>
void InterpolatorConsistencyTest(Functor f)
{
	static const int nTests = 5;
	static const int nPoints = 20;

	typedef LagrangeInterpolator<double, mydim, cdim>           Interpolator;
	typedef typename Interpolator::LocalCoordinate              LocalCoordinate;
	typedef typename Interpolator::GlobalCoordinate             GlobalCoordinate;
	typedef typename Interpolator::PolynomialGlobalCoordinate   PolynomialGlobalCoordinate;

	typedef typename std::vector<LocalCoordinate>               LocalCoordinateSet;
	typedef typename std::vector<GlobalCoordinate>              GlobalCoordinateSet;

	typedef FunctionEvaluator<Interpolator>         FunctionEvaluator;
	typedef InterpolatorEvaluator<Interpolator>     InterpolatorEvaluator;
	typedef PolynomialCoordEvaluator<Interpolator>  PolynomialCoordEvaluator;

	Interpolator                interps[nTests];
	PolynomialGlobalCoordinate  analytics[nTests];

//    Dune::GeometryType geomType;   geomType.makeSimplex(mydim);
    Dune::GeometryType geomType=Dune::GeometryTypes::simplex(mydim);


    std::cout << "*********************** Start " << geomType << " Test *************************" << std::endl;

    // Test mapping for actual interpolatory points
    for (int i = 0; i < nTests; i++)
    {
    	LocalCoordinateSet  regularRefGridVertexLocal  = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, mydim>(i+1);
    	GlobalCoordinateSet regularRefGridVertexGlobal = FunctionEvaluator::eval(regularRefGridVertexLocal, f);
        interps[i]   = Interpolator(geomType, regularRefGridVertexGlobal, i + 1);
        analytics[i] = interps[i].interpolatoryVectorAnalytical();

        GlobalCoordinateSet interpolatorPredictionGlobal = InterpolatorEvaluator::eval(regularRefGridVertexLocal, interps[i]);
        GlobalCoordinateSet polynomialPredictionGlobal   = PolynomialCoordEvaluator::eval(regularRefGridVertexLocal, analytics[i]);

        std::cout << "Gridsample> order: " << i + 1;
        std::cout << ", Real-Numerical error: " << norm2sqVector<cdim>(regularRefGridVertexGlobal, interpolatorPredictionGlobal);
        std::cout << ", Numerical-Analytical error: " << norm2sqVector<cdim>(interpolatorPredictionGlobal, polynomialPredictionGlobal) << std::endl;
    }

    // Test mapping for random points within reference element
    LocalCoordinateSet  randomPointSampleLocal  = sampleRandom<mydim>(nPoints);
    GlobalCoordinateSet randomPointSampleGlobal = FunctionEvaluator::eval(randomPointSampleLocal, f);
    for (int i = 0; i < nTests; i++)
    {
    	GlobalCoordinateSet interpolatorPredictionGlobal = InterpolatorEvaluator::eval(randomPointSampleLocal, interps[i]);
    	GlobalCoordinateSet polynomialPredictionGlobal   = PolynomialCoordEvaluator::eval(randomPointSampleLocal, analytics[i]);

        std::cout << "Random sample> pointnumber: " << nPoints;
        std::cout << ", Real-Numerical error: " << norm2sqVector<cdim>(randomPointSampleGlobal, interpolatorPredictionGlobal);
        std::cout << ", Numerical-Analytical error: " << norm2sqVector<cdim>(interpolatorPredictionGlobal, polynomialPredictionGlobal) << std::endl;
    }
    std::cout << "*********************** Finish " << geomType << " Test *************************" << std::endl << std::endl;
}



struct f3dLine1     { FieldVector3D operator()(const FieldVector1D & x) const { return initFieldVector(1, 0, 0); }  };
struct f3dLine2     { FieldVector3D operator()(const FieldVector1D & x) const { return initFieldVector(x[0], 0, 0); }  };
struct f3dLine3     { FieldVector3D operator()(const FieldVector1D & x) const { return initFieldVector(x[0] * x[0], x[0], 1); }  };
struct f3dLine4     { FieldVector3D operator()(const FieldVector1D & x) const { return initFieldVector(x[0] * x[0] * x[0], x[0], 1); }  };

struct f3dSurface1  { FieldVector3D operator()(const FieldVector2D & x) const { return initFieldVector(x[0], x[1], 0); }  };
struct f3dSurface2  { FieldVector3D operator()(const FieldVector2D & x) const { return initFieldVector(x[0], x[1], x[0]); }  };
struct f3dSurface3  { FieldVector3D operator()(const FieldVector2D & x) const { return initFieldVector(x[0], x[1], x[0] * x[1]); }  };

struct f3dVolume1   { FieldVector3D operator()(const FieldVector3D & x) const { return initFieldVector(x[0], x[1], x[0] * x[1] * x[2]); }  };
struct f3dVolume2   { FieldVector3D operator()(const FieldVector3D & x) const { return initFieldVector(4 * x[0] * x[0] * x[0] * x[0] , x[1] + 3*x[2], x[0] * x[1] * x[2]); }  };




int main() {
    srand (time(NULL));

    InterpolatorConsistencyTest<DIM1D, DIM3D, f3dLine1>(f3dLine1());
    InterpolatorConsistencyTest<DIM1D, DIM3D, f3dLine2>(f3dLine2());
    InterpolatorConsistencyTest<DIM1D, DIM3D, f3dLine3>(f3dLine3());
    InterpolatorConsistencyTest<DIM1D, DIM3D, f3dLine4>(f3dLine4());
    InterpolatorConsistencyTest<DIM2D, DIM3D, f3dSurface1>(f3dSurface1());
    InterpolatorConsistencyTest<DIM2D, DIM3D, f3dSurface2>(f3dSurface2());
    InterpolatorConsistencyTest<DIM2D, DIM3D, f3dSurface3>(f3dSurface3());
    InterpolatorConsistencyTest<DIM3D, DIM3D, f3dVolume1>(f3dVolume1());
    InterpolatorConsistencyTest<DIM3D, DIM3D, f3dVolume2>(f3dVolume2());
}
