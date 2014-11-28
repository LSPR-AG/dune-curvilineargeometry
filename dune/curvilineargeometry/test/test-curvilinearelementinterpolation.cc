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
#include <dune/curvilineargeometry/integration/numericalrecursiveinterpolationintegrator.hh>


using namespace Dune;

typedef FieldVector< double, 1 > FieldVector1D;
typedef FieldVector< double, 2 > FieldVector2D;
typedef FieldVector< double, 3 > FieldVector3D;

typedef std::vector<FieldVector1D> FieldVectorVector1D;
typedef std::vector<FieldVector2D> FieldVectorVector2D;
typedef std::vector<FieldVector3D> FieldVectorVector3D;

typedef polynomial<double, 1> Polynomial1D;
typedef polynomial<double, 2> Polynomial2D;
typedef polynomial<double, 3> Polynomial3D;
typedef std::vector<Polynomial1D> PolynomialVector1D;
typedef std::vector<Polynomial2D> PolynomialVector2D;
typedef std::vector<Polynomial3D> PolynomialVector3D;



double randomReal(double a, double b) { return a + (b - a)*(double(rand()) / RAND_MAX); }

// Constructs a random grid over a simplex with a given number of samples
template <int dim>
std::vector<FieldVector< double, dim > > sampleRandom(int sampleNo) {
    std::vector<FieldVector< double, dim > > tmprez;
    return tmprez;
}

// Generates set of coordinates randomly (uniformly) sampled over an edge
template <>
FieldVectorVector1D sampleRandom<1>(int sampleNo) {
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
FieldVectorVector2D sampleRandom<2>(int sampleNo) {
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
FieldVectorVector3D sampleRandom<3>(int sampleNo) {
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
template <int dim, int dimworld, typename Functor>
std::vector<FieldVector< double, dimworld > > evaluateFunction(
        std::vector<FieldVector< double, dim > > & localsample,
        Functor f )
{
    std::vector<FieldVector< double, dimworld > > rez;
    for (int i = 0; i < localsample.size(); i++) {rez.push_back(f(localsample[i])); }
    return rez;
}

// Evaluates an elementInterpolator for a given vector of input coordinates, returns a vector of output coordinates
template <int dim, int dimworld>
std::vector<FieldVector< double, dimworld > > evaluateInterpolatorNumerical(
        std::vector<FieldVector< double, dim > > & localsample,
        CurvilinearElementInterpolator<double, dim, dimworld> & interp )
{
    std::vector<FieldVector< double, dimworld > > rez;
    for (int i = 0; i < localsample.size(); i++) {rez.push_back(interp.realCoordinate(localsample[i])); }
    return rez;
}

// Evaluates a vector of polynomials for a given vector of input coordinates, returns a vector of output coordinates
template <int dim, int dimworld>
std::vector<FieldVector< double, dimworld > > evaluatePolynomialVector(
        std::vector< FieldVector< double, dim > > & localsample,
        std::vector< polynomial<double, dim> >& polynomialvector )
{
    std::vector<FieldVector< double, dimworld > > rez;
    for (int i = 0; i < localsample.size(); i++) {
        FieldVector< double, dimworld > tmp;
        for (int j = 0; j < dimworld; j++) {
            tmp[j] = polynomialvector[j].evaluate(localsample[i]);
        }
        rez.push_back(tmp);
    }
    return rez;
}

// Initializes a field vector, such that it can be done in one line in main code
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



// Runs a test to compare how well the exact edge given by f
// is approximated by Lagrange interpolation. Also checks if
// analytical and numerical lagrange interpolation give the same
// result which they should

// TODO: Implement Real-Numerical in LagrangeInterpolation, then add to this test
template <typename Functor>
void edgeTest3D(Functor f)
{
    CurvilinearElementInterpolator<double, 1, 3> interps[5];
    std::vector<polynomial<double, 1> > analytics[5];

    FieldVectorVector1D random1D_20 = sampleRandom<1>(20);
    FieldVectorVector3D sample_real_3D_line_p20_randomtest = evaluateFunction<1,3>(random1D_20, f);

    Dune::GeometryType edgeGeometry;   edgeGeometry.makeSimplex(1);

    std::cout << "*********************** Start Edge Test *************************" << std::endl;
    for (int i = 0; i < 5; i++)
    {
        FieldVectorVector1D grid1D = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, 1>(i+1);
        FieldVectorVector3D sample_real_3D_line_gridtest = evaluateFunction<1,3, Functor>(grid1D, f);
        interps[i] = CurvilinearElementInterpolator<double, 1, 3>(edgeGeometry, sample_real_3D_line_gridtest, i + 1);
        analytics[i] = interps[i].interpolatoryVectorAnalytical();

        FieldVectorVector3D sample_numerical_3D_line_gridtest = evaluateInterpolatorNumerical<1,3>(grid1D, interps[i]);
        FieldVectorVector3D sample_analytical_3D_line_gridtest = evaluatePolynomialVector<1,3>(grid1D, analytics[i]);

        std::cout << "Surface test - gridsample - order " << i + 1;
        std::cout << ", Real-Numerical error: " << norm2sqVector<3>(sample_real_3D_line_gridtest, sample_numerical_3D_line_gridtest);
        std::cout << ", Numerical-Analytical error: " << norm2sqVector<3>(sample_numerical_3D_line_gridtest, sample_analytical_3D_line_gridtest) << std::endl;
    }


    for (int i = 0; i < 5; i++)
    {
        FieldVectorVector3D sample_numerical_3D_line_randomtest = evaluateInterpolatorNumerical<1,3>(random1D_20, interps[i]);
        FieldVectorVector3D sample_analytical_3D_line_randomtest = evaluatePolynomialVector<1,3>(random1D_20, analytics[i]);

        std::cout << "Edge test - randomsample - pointnumber 20";
        std::cout << ", Real-Numerical error: " << norm2sqVector<3>(sample_real_3D_line_p20_randomtest, sample_numerical_3D_line_randomtest);
        std::cout << ", Numerical-Analytical error: " << norm2sqVector<3>(sample_numerical_3D_line_randomtest, sample_analytical_3D_line_randomtest) << std::endl;
    }
    std::cout << "*********************** Finish Edge Test *************************" << std::endl << std::endl;
}



// Runs a test to compare how well the exact surface given by f
// is approximated by Lagrange interpolation. Also checks if
// analytical and numerical lagrange interpolation give the same
// result which they should
template <typename Functor>
void surfaceTest3D(Functor f)
{
    CurvilinearElementInterpolator<double, 2, 3> interps[5];
    std::vector<polynomial<double, 2> > analytics[5];

    FieldVectorVector2D random2D_20 = sampleRandom<2>(20);
    FieldVectorVector3D sample_real_3D_surface_p20_randomtest = evaluateFunction<2,3, Functor>(random2D_20, f);

    Dune::GeometryType triangleGeometry;   triangleGeometry.makeSimplex(2);

    std::cout << "*********************** Start Surface Test *************************" << std::endl;
    for (int i = 0; i < 5; i++)
    {
        FieldVectorVector2D grid2D = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, 2>(i+1);
        FieldVectorVector3D sample_real_3D_surface_gridtest = evaluateFunction<2,3>(grid2D, f);
        interps[i] = CurvilinearElementInterpolator<double, 2, 3>(triangleGeometry, sample_real_3D_surface_gridtest, i + 1);
        analytics[i] = interps[i].interpolatoryVectorAnalytical();

        FieldVectorVector3D sample_numerical_3D_surface_gridtest = evaluateInterpolatorNumerical<2,3>(grid2D, interps[i]);
        FieldVectorVector3D sample_analytical_3D_surface_gridtest = evaluatePolynomialVector<2,3>(grid2D, analytics[i]);

        std::cout << "Surface test - gridsample - order " << i + 1;
        std::cout << ", Real-Numerical error: " << norm2sqVector<3>(sample_real_3D_surface_gridtest, sample_numerical_3D_surface_gridtest);
        std::cout << ", Numerical-Analytical error: " << norm2sqVector<3>(sample_numerical_3D_surface_gridtest, sample_analytical_3D_surface_gridtest) << std::endl;
    }


    for (int i = 0; i < 5; i++)
    {
        FieldVectorVector3D sample_numerical_3D_surface_randomtest = evaluateInterpolatorNumerical<2,3>(random2D_20, interps[i]);
        FieldVectorVector3D sample_analytical_3D_surface_randomtest = evaluatePolynomialVector<2,3>(random2D_20, analytics[i]);

        std::cout << "Surface test - randomsample - pointnumber 20";
        std::cout << ", Real-Numerical error: " << norm2sqVector<3>(sample_real_3D_surface_p20_randomtest, sample_numerical_3D_surface_randomtest);
        std::cout << ", Numerical-Analytical error: " << norm2sqVector<3>(sample_numerical_3D_surface_randomtest, sample_analytical_3D_surface_randomtest) << std::endl;
    }
    std::cout << "*********************** Finish Surface Test *************************" << std::endl << std::endl;
}

// Runs a test to compare how well the exact volume given by *f
// is approximated by Lagrange interpolation. Also checks if
// analytical and numerical lagrange interpolation give the same
// result which they should
template <typename Functor>
void volumeTest3D(Functor f)
{
    CurvilinearElementInterpolator<double, 3, 3> interps[5];
    std::vector<polynomial<double, 3> > analytics[5];

    FieldVectorVector3D random3D_20 = sampleRandom<3>(20);
    FieldVectorVector3D sample_real_3D_volume_p20_randomtest = evaluateFunction<3,3, Functor>(random3D_20, f);

    Dune::GeometryType tetrahedronGeometry;   tetrahedronGeometry.makeSimplex(3);

    std::cout << "*********************** Start Volume Test *************************" << std::endl;
    for (int i = 0; i < 5; i++)
    {
        FieldVectorVector3D grid3D = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, 3>(i+1);
        FieldVectorVector3D sample_real_3D_volume_gridtest = evaluateFunction<3,3>(grid3D, f);
        interps[i] = CurvilinearElementInterpolator<double, 3, 3>(tetrahedronGeometry, sample_real_3D_volume_gridtest, i + 1);
        analytics[i] = interps[i].interpolatoryVectorAnalytical();

        FieldVectorVector3D sample_numerical_3D_volume_gridtest = evaluateInterpolatorNumerical<3,3>(grid3D, interps[i]);
        FieldVectorVector3D sample_analytical_3D_volume_gridtest = evaluatePolynomialVector<3,3>(grid3D, analytics[i]);

        std::cout << "Volume test - gridsample - order " << i + 1;
        std::cout << ", Real-Numerical error: "       << norm2sqVector<3>(sample_real_3D_volume_gridtest, sample_numerical_3D_volume_gridtest);
        std::cout << ", Numerical-Analytical error: " << norm2sqVector<3>(sample_numerical_3D_volume_gridtest, sample_analytical_3D_volume_gridtest) << std::endl;
    }


    for (int i = 0; i < 5; i++)
    {
        FieldVectorVector3D sample_numerical_3D_volume_randomtest = evaluateInterpolatorNumerical<3,3>(random3D_20, interps[i]);
        FieldVectorVector3D sample_analytical_3D_volume_randomtest = evaluatePolynomialVector<3,3>(random3D_20, analytics[i]);

        std::cout << "Volume test - randomsample - pointnumber 20";
        std::cout << ", Real-Numerical error: "       << norm2sqVector<3>(sample_real_3D_volume_p20_randomtest, sample_numerical_3D_volume_randomtest);
        std::cout << ", Numerical-Analytical error: " << norm2sqVector<3>(sample_numerical_3D_volume_randomtest, sample_analytical_3D_volume_randomtest) << std::endl;
    }
    std::cout << "*********************** End Volume Test *************************" << std::endl << std::endl;
}



struct f3dLine1 {   FieldVector3D operator()(const FieldVector1D & x) { return initFieldVector(1, 0, 0); }  };
struct f3dLine2 {   FieldVector3D operator()(const FieldVector1D & x) { return initFieldVector(x[0], 0, 0); }  };
struct f3dLine3 {   FieldVector3D operator()(const FieldVector1D & x) { return initFieldVector(x[0] * x[0], x[0], 1); }  };
struct f3dLine4 {   FieldVector3D operator()(const FieldVector1D & x) { return initFieldVector(x[0] * x[0] * x[0], x[0], 1); }  };

struct f3dSurface1 {   FieldVector3D operator()(const FieldVector2D & x) { return initFieldVector(x[0], x[1], 0); }  };
struct f3dSurface2 {   FieldVector3D operator()(const FieldVector2D & x) { return initFieldVector(x[0], x[1], x[0]); }  };
struct f3dSurface3 {   FieldVector3D operator()(const FieldVector2D & x) { return initFieldVector(x[0], x[1], x[0] * x[1]); }  };

struct f3dVolume1 {   FieldVector3D operator()(const FieldVector3D & x) { return initFieldVector(x[0], x[1], x[0] * x[1] * x[2]); }  };
struct f3dVolume2 {   FieldVector3D operator()(const FieldVector3D & x) { return initFieldVector(4 * x[0] * x[0] * x[0] * x[0] , x[1] + 3*x[2], x[0] * x[1] * x[2]); }  };




int main() {
    srand (time(NULL));

    edgeTest3D(f3dLine1());
    edgeTest3D(f3dLine2());
    edgeTest3D(f3dLine3());
    edgeTest3D(f3dLine4());
    surfaceTest3D(f3dSurface1());
    surfaceTest3D(f3dSurface2());
    surfaceTest3D(f3dSurface3());
    volumeTest3D(f3dVolume1());
    volumeTest3D(f3dVolume2());
}
