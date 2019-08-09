/*******************************************************************
 * Test of Curvilinear Geometry
 * 
 * author: Aleksejs Fomins
 * date: 01.09.2014 - created
 * 
 * description:
 * Runs a set of tests checking local->global and global->local maps, isinside point location and integration
 * 
 * 
 *******************************************************************/


#include <config.h>

#include <vector>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/interpolation/differentialhelper.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/interpolation/lagrangeinterpolator.hh>
#include <dune/curvilineargeometry/integration/adaptiveintegrator.hh>

#include <dune/curvilineargeometry/curvilineargeometry.hh>





using namespace Dune;

typedef Dune::FieldVector< double, 1 > FieldVector1D;
typedef Dune::FieldVector< double, 2 > FieldVector2D;
typedef Dune::FieldVector< double, 3 > FieldVector3D;

typedef std::vector<FieldVector1D> FieldVectorVector1D;
typedef std::vector<FieldVector2D> FieldVectorVector2D;
typedef std::vector<FieldVector3D> FieldVectorVector3D;

typedef PolynomialTraits<double>::Monomial  Monomial;

struct TestStruct
{
	int nFail_;
	int nTot_;
	double worst_;

	TestStruct(int nFail, int nTot, double worst) : nFail_(nFail), nTot_(nTot), worst_(worst) {}
};

typedef std::pair<TestStruct, TestStruct>  TestPair;


// ************************************************************************************************
// IMPLEMENTING AUXILIARY METHODS
// ************************************************************************************************

// Generates uniform random numbers in interval [a,b]
// [TODO] Move all randomness to testhelperrandom.hh
double randomReal(double a, double b) { return a + (b - a)*(double(rand()) / RAND_MAX); }

// Constructs a random grid over a box, where each dimension is uniformly in the interval [a,b]
template <int dim>
std::vector<FieldVector< double, dim > > sampleBox(int sampleNo, double a, double b) {
    std::vector<FieldVector< double, dim > > rez;

    for (int i = 0; i < sampleNo; i++) {
        FieldVector< double, dim > tmp;
        for (int j = 0; j < dim; j++) { tmp[j] = randomReal(a,b); }
        rez.push_back(tmp);
    }
    return rez;
}

// Constructs a random grid over a simplex with a given number of samples
template <int dim>
std::vector<FieldVector< double, dim > > sampleRandom(int sampleNo) {
	DUNE_THROW(NotImplemented, "CurvilinearGeometryTest - using generic method");
    return std::vector<FieldVector< double, dim > >();
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


// ************************************************************************************************
// IMPLEMENTING FUNCTORS - EXAMPLE LOCAL-TO-GLOBAL MAPS
// ************************************************************************************************


// Generic Functor
template<class ctype, int mydim, int cdim> struct myFieldVectorFunctor
{
	typedef Dune::FieldVector<ctype, mydim>  InputValue;
	typedef Dune::FieldVector<ctype, cdim>   ResultValue;
	typedef std::vector<ResultValue>         ResultType;
	static const int RETURN_SIZE = 1;

	ResultValue zeroValue(unsigned int rezIndex) const { return ResultValue(0.0); }
};


// Generic Functors of different dimensions
typedef myFieldVectorFunctor<double, 1, 1>  Base11;
typedef myFieldVectorFunctor<double, 1, 2>  Base12;
typedef myFieldVectorFunctor<double, 1, 3>  Base13;
typedef myFieldVectorFunctor<double, 2, 2>  Base22;
typedef myFieldVectorFunctor<double, 2, 3>  Base23;
typedef myFieldVectorFunctor<double, 3, 3>  Base33;

template<class ctype, int mydim, int cdim> struct myFunctorIdentity    : public myFieldVectorFunctor<ctype, mydim, cdim> { };
template<class ctype, int mydim, int cdim> struct myFunctorLinear      : public myFieldVectorFunctor<ctype, mydim, cdim> { };
template<class ctype, int mydim, int cdim> struct myFunctorNonlinear1  : public myFieldVectorFunctor<ctype, mydim, cdim> { };
template<class ctype, int mydim, int cdim> struct myFunctorNonlinear2  : public myFieldVectorFunctor<ctype, mydim, cdim> { };

// Actual functors that evaluate some Local-To-Global transform
template<> struct myFunctorIdentity<double, 1, 1> : public Base11  { Base11::ResultType operator()(const Base11::InputValue & in) const { return Base11::ResultType(1, initFieldVector(in[0])); }  };
template<> struct myFunctorIdentity<double, 1, 2> : public Base12  { Base12::ResultType operator()(const Base12::InputValue & in) const { return Base12::ResultType(1, initFieldVector(in[0], 0)); }  };
template<> struct myFunctorIdentity<double, 1, 3> : public Base13  { Base13::ResultType operator()(const Base13::InputValue & in) const { return Base13::ResultType(1, initFieldVector(in[0], 0, 0)); }  };
template<> struct myFunctorIdentity<double, 2, 2> : public Base22  { Base22::ResultType operator()(const Base22::InputValue & in) const { return Base22::ResultType(1, initFieldVector(in[0], in[1])); }  };
template<> struct myFunctorIdentity<double, 2, 3> : public Base23  { Base23::ResultType operator()(const Base23::InputValue & in) const { return Base23::ResultType(1, initFieldVector(in[0], in[1], 0)); }  };
template<> struct myFunctorIdentity<double, 3, 3> : public Base33  { Base33::ResultType operator()(const Base33::InputValue & in) const { return Base33::ResultType(1, initFieldVector(in[0], in[1], in[2])); }  };

template<> struct myFunctorLinear<double, 1, 1> : public Base11  { Base11::ResultType operator()(const Base11::InputValue & in) const { return Base11::ResultType(1, initFieldVector(1.0 + 2.0 * in[0])); }  };
template<> struct myFunctorLinear<double, 1, 2> : public Base12  { Base12::ResultType operator()(const Base12::InputValue & in) const { return Base12::ResultType(1, initFieldVector(2.0 * in[0], 3.0 * in[0])); }  };
template<> struct myFunctorLinear<double, 1, 3> : public Base13  { Base13::ResultType operator()(const Base13::InputValue & in) const { return Base13::ResultType(1, initFieldVector(2.0 * in[0], 0.5 + 3.0 * in[0], 5.0 * in[0])); }  };
template<> struct myFunctorLinear<double, 2, 2> : public Base22  { Base22::ResultType operator()(const Base22::InputValue & in) const { return Base22::ResultType(1, initFieldVector(1.0 + in[0], in[0] + in[1])); }  };
template<> struct myFunctorLinear<double, 2, 3> : public Base23  { Base23::ResultType operator()(const Base23::InputValue & in) const { return Base23::ResultType(1, initFieldVector(in[1], 3.0 * in[0], in[0] + in[1])); }  };
template<> struct myFunctorLinear<double, 3, 3> : public Base33  { Base33::ResultType operator()(const Base33::InputValue & in) const { return Base33::ResultType(1, initFieldVector(in[0] + in[1], in[1] + in[2], in[2] + in[0])); }  };

template<> struct myFunctorNonlinear1<double, 1, 1> : public Base11  { Base11::ResultType operator()(const Base11::InputValue & in) const { return Base11::ResultType(1, initFieldVector(in[0] * in[0])); }  };
template<> struct myFunctorNonlinear1<double, 1, 2> : public Base12  { Base12::ResultType operator()(const Base12::InputValue & in) const { return Base12::ResultType(1, initFieldVector(in[0], in[0] * in[0])); }  };
template<> struct myFunctorNonlinear1<double, 1, 3> : public Base13  { Base13::ResultType operator()(const Base13::InputValue & in) const { return Base13::ResultType(1, initFieldVector(in[0], in[0] * in[0], 2.0)); }  };
template<> struct myFunctorNonlinear1<double, 2, 2> : public Base22  { Base22::ResultType operator()(const Base22::InputValue & in) const { return Base22::ResultType(1, initFieldVector(in[0]*in[0], in[1] * in[1])); }  };
template<> struct myFunctorNonlinear1<double, 2, 3> : public Base23  { Base23::ResultType operator()(const Base23::InputValue & in) const { return Base23::ResultType(1, initFieldVector(in[1] * in[1], in[0] * in[0], in[0] * in[1])); }  };
template<> struct myFunctorNonlinear1<double, 3, 3> : public Base33  { Base33::ResultType operator()(const Base33::InputValue & in) const { return Base33::ResultType(1, initFieldVector(in[0] * in[0], in[1] * in[1], in[2] * in[2])); }  };



// ************************************************************************************************
// IMPLEMENTING TEST RESULTS
// ************************************************************************************************

// Returns one of 5 possible polynomials such that polynomial order = index
template <class CurvGeom>
typename CurvGeom::LocalPolynomial BasisPolynomial(int index)
{
	typedef typename CurvGeom::LocalPolynomial  LocalPolynomial;

    switch (CurvGeom::mydimension)
    {
    case 1:
    {
    	LocalPolynomial rez(Monomial(1, 0));
        if (index > 0) { rez += Monomial(2, 1); }
        if (index > 1) { rez += Monomial(3, 2); }
        if (index > 2) { rez += Monomial(4, 3); }
        if (index > 3) { rez += Monomial(5, 4); }
        if (index > 4) { rez += Monomial(6, 5); }
        return rez;
    }  break;
    case 2:
    {
    	LocalPolynomial rez(Monomial(1, 0, 0));

        if (index > 0) { rez += Monomial(2, 1, 0);  rez += Monomial(2, 0, 1); }
        if (index > 1) { rez += Monomial(3, 2, 0);  rez += Monomial(3, 0, 2); rez += Monomial(1, 1, 1); }
        if (index > 2) { rez += Monomial(4, 3, 0);  rez += Monomial(4, 0, 3); rez += Monomial(1, 1, 2); }
        if (index > 3) { rez += Monomial(5, 4, 0);  rez += Monomial(5, 0, 4); rez += Monomial(1, 1, 3); }
        if (index > 4) { rez += Monomial(6, 5, 0);  rez += Monomial(6, 0, 5); rez += Monomial(1, 1, 4); }

        return rez;
    }  break;
    case 3:
    {
    	LocalPolynomial rez(Monomial(1, 0, 0, 0));
        if (index > 0) { rez += Monomial(2, 1, 0, 0);  rez += Monomial(2, 0, 1, 0);  rez += Monomial(2, 0, 0, 1); }
        if (index > 1) { rez += Monomial(3, 2, 0, 0);  rez += Monomial(3, 0, 2, 0);  rez += Monomial(3, 0, 0, 2); rez += Monomial(1, 1, 1, 0); }
        if (index > 2) { rez += Monomial(4, 3, 0, 0);  rez += Monomial(4, 0, 3, 0);  rez += Monomial(4, 0, 0, 3); rez += Monomial(1, 1, 1, 1); }
        if (index > 3) { rez += Monomial(5, 4, 0, 0);  rez += Monomial(5, 0, 4, 0);  rez += Monomial(5, 0, 0, 4); rez += Monomial(1, 1, 1, 2); }
        if (index > 4) { rez += Monomial(6, 5, 0, 0);  rez += Monomial(6, 0, 5, 0);  rez += Monomial(6, 0, 0, 5); rez += Monomial(1, 1, 1, 3); }
        return rez;
    }  break;
    default:  DUNE_THROW(RangeError, "Unexpected mydimension");   break;
    }
}


// Construct a mydim+1 polynomial vector to test the Surface Dot Product Integral
// Only for edges in 2d and faces in 3d
template<class CurvGeom>
typename CurvGeom::PolynomialGlobalCoordinate BasisVectorDiagonal()
{
	typedef typename CurvGeom::LocalPolynomial            LocalPolynomial;
	typedef typename CurvGeom::PolynomialGlobalCoordinate PolyCoord;
	PolyCoord rez;

    switch (CurvGeom::mydimension)
    {
    case 1:
        rez[0] = LocalPolynomial (Monomial(1, 1));
        rez[1] = LocalPolynomial (Monomial(1, 1));
        break;
    case 2:
        rez[0] = LocalPolynomial (Monomial(1, 1, 0));
        rez[1] = LocalPolynomial (Monomial(1, 0, 1));
        rez[2] = LocalPolynomial (Monomial(1, 1, 1));
        break;
    }

    return rez;
}

// A big array of all results of the integral tests calculated by hand
// Exact integrals and results should be given in Doc.
template<int mydim, int cdim>
double integralResult(int bf_order, int el_type)
{
    double rez_1D_7[6] = { 1.0, 7.0/3, 23.0 / 6, 163.0/30, 71.0/10, 617.0/70 };
    double rez_1D_89[6] = { 1.47894286, 3.175666172, 4.994678155, 6.89140143, 8.84167808, 10.83102449 };
    double rez_2D_123[6] = {0.5, 7.0/6, 41.0/24, 17.0/8, 37.0/15, 2.75714};
    double rez_2D_5[6] = {1.0/6, 13.0/30, 59.0/90, 103.0/126, 0.94127, 1.03915};
    double rez_2D_6[6] = {0.360858, 0.938231, 1.47326, 1.93004, 2.33506, 2.70079};
    double rez_3D_1[6] = {1.0/6, 5.0/12, 23.0/40, 0.676389, 0.748214, 0.801935};
    double rez_3D_3[6] = {1.0/90, 0.0301587, 0.0416667, 0.0481922, 0.0522134, 0.05483};


    switch (mydim)
    {
        case 1:
        {
            switch (el_type)
            {
            case 1:   return bf_order + 1;  break;
            case 2:
                switch (cdim)
                {
                case 1:  return 2 * (bf_order + 1);  break;
                case 2:  return sqrt(13) * (bf_order + 1);  break;
                case 3:  return sqrt(38) * (bf_order + 1);  break;
                } break;
            case 3:
                if (cdim == 1) { return rez_1D_7[bf_order]; }
                else           { return rez_1D_89[bf_order]; }
            }
        } break;
        case 2:
        {
            switch(el_type)
            {
                case 1:  return rez_2D_123[bf_order];  break;
                case 2:
                {
                    if (cdim == 2) { return rez_2D_123[bf_order]; }
                    else           { return sqrt(19) * rez_2D_123[bf_order]; }
                }  break;
                case 3:
                {
                    if (cdim == 2) { return rez_2D_5[bf_order]; }
                    else           { return rez_2D_6[bf_order]; }
                } break ;
            }
        } break;
        case 3:
        {
            switch(el_type)
            {
            case 1:  return rez_3D_1[bf_order];      break;
            case 2:  return 2 * rez_3D_1[bf_order];  break;
            case 3:  return rez_3D_3[bf_order];      break;
            }
        } break;
    }
}

// A big array of all results of the Surface Dot integral tests calculated by hand
// Exact integrals and results should be given in Doc.
template<int mydim, int cdim>
double integralDotResult(int el_type)
{
    if (mydim == 1)
    {
        switch (el_type)
        {
        case 1: return -0.5;  break;
        case 2: return 0.5;  break;
        case 3: return 1.0 / 6;  break;
        }
    } else
    {
        switch (el_type)
        {
        case 1: return -1.0 / 24;  break;
        case 2: return -13.0 / 24;  break;
        case 3: return -8.0/45;  break;
        }
    }
}






// ************************************************************************************************
// IMPLEMENTING TESTS
// ************************************************************************************************


// ---------------------------------------------------------------------------------------------------------------------------
// 1. CORNER-TEST
// This test checks if the corner coordinates returned by CurvilinearGeometry correspond to the ones used in initializing it.
// ---------------------------------------------------------------------------------------------------------------------------
template<class Functor, class SimplexGeom, class ReferenceElement>
TestStruct test01_corner(
	const Functor & f,
	const SimplexGeom & simplexGeom,
	const ReferenceElement & refElement,
	double ABSOLUTE_TOLERANCE)
{
	static const int mydim = SimplexGeom::mydimension;
	typedef typename SimplexGeom::GlobalCoordinate   GlobalCoordinate;

	TestStruct rez(0, 0, 0.0);

    for (int i = 0; i < mydim+1; i++) {
        //GlobalCoordinate tmpCorner = f(localCorner<double, mydim>(i))[0];
    	GlobalCoordinate tmpCorner = f(refElement.position(i, mydim))[0];
    	double err = (simplexGeom.corner(i) - tmpCorner).two_norm();

    	rez.nTot_++;
    	if (err > ABSOLUTE_TOLERANCE)  { rez.nFail_++; }
    	rez.worst_ = std::max(rez.worst_, err);
    }

    return rez;
}


// ---------------------------------------------------------------------------------------------------------------------------
// 2. LOCAL->GLOBAL TEST
// This test samples the example global map over the simplex, and checks if these samples match with the CurvilinearGeometry.global()
// ----------------------------------------------------------------------------------------------------------------------------------
template<class Functor, class SimplexGeom, class LocalVectorVector>
TestStruct test02_local_to_global(
	bool verbose,
	const Functor & f,
	const SimplexGeom & simplexGeom,
	const LocalVectorVector & randomLocalSample,
	int interpOrder,
	int functionOrder,
	double ABSOLUTE_TOLERANCE)
{
	TestStruct rez(0, 0, 0.0);

	// Only count as error if interpolatory polynomial is of sufficient order
    if (interpOrder < functionOrder) { if (verbose) { std::cout << ": Local-To-Global-test: --Omitted because polynomial order too small" << std::endl; } }
    else {
        for (unsigned int i = 0; i < randomLocalSample.size(); i++)
        {
            double err = (simplexGeom.global(randomLocalSample[i]) - f(randomLocalSample[i])[0]).two_norm();

        	rez.nTot_++;
        	if (err > ABSOLUTE_TOLERANCE)  { rez.nFail_++; }
        	rez.worst_ = std::max(rez.worst_, err);
        }
    }

    return rez;
}


// ---------------------------------------------------------------------------------------------------------------------------
// 3. GLOBAL->LOCAL TEST
// This test checks if local coordinates of the global interpolation vertices are recovered
// ----------------------------------------------------------------------------------------------------------------------------------
template<class SimplexGeom, class GlobalCoordinate, class LocalCoordinate>
TestPair test03_global_to_local(
	bool verbose,
	const SimplexGeom & simplexGeom,
	const std::vector<GlobalCoordinate> & global_vertices,
	const std::vector<LocalCoordinate>  & local_vertices,
	double ABSOLUTE_TOLERANCE)

{
	static const int cdim  = SimplexGeom::coorddimension;
	static const int mydim = SimplexGeom::mydimension;

	TestStruct rezTmp(0, 0, 0.0);
	TestPair rez(rezTmp, rezTmp);
	LocalCoordinate L;

    if (mydim != cdim) { if (verbose) { std::cout << ": Global-to-Local functionality not available for mismatching mydim and cdim" << std::endl; } }
    else {
        for (unsigned int i = 0; i < global_vertices.size(); i++) {
        	rez.first.nTot_++;
        	rez.second.nTot_++;

        	// Order important, local() method initializes L
        	bool isInside = simplexGeom.local(global_vertices[i], L);
        	double err = (L - local_vertices[i]).two_norm();

        	if (!isInside)                { rez.first.nFail_++; }
        	if (err > ABSOLUTE_TOLERANCE) { rez.second.nFail_++; }
        	rez.second.worst_ = std::max(rez.second.worst_, err);
            //std::cout << "  -- recovered vertex " << L << " from real vertex " << local_vertices[i] << " via global " << global_vertices[i] << " positioning is " << is_inside << std::endl;
        }
    }
    return rez;
}


// ---------------------------------------------------------------------------------------------------------------------------
// 4. GLOBAL->LOCAL TEST 2
// Test if LOCAL->GLOBAL->LOCAL for random sample inside element is preserved
// ---------------------------------------------------------------------------------------------------------------------------
template<class SimplexGeom, class LocalCoordinate>
TestPair test04_global_to_local(
	bool verbose,
	const SimplexGeom & simplexGeom,
	const std::vector<LocalCoordinate> & randomLocalSample,
	double ABSOLUTE_TOLERANCE)

{
	static const int cdim  = SimplexGeom::coorddimension;
	static const int mydim = SimplexGeom::mydimension;

	TestStruct rezTmp(0, 0, 0.0);
	TestPair rez(rezTmp, rezTmp);
	LocalCoordinate L;

    if (mydim != cdim) { if (verbose) { std::cout << ": Global-to-Local functionality not available for mismatching mydim and cdim" << std::endl; } }
    else {
        for (unsigned int i = 0; i < randomLocalSample.size(); i++)
        {
        	rez.first.nTot_++;
        	rez.second.nTot_++;

        	// Order important, local() method initializes L
        	bool isInside = simplexGeom.local(simplexGeom.global(randomLocalSample[i]), L);
        	double err = (randomLocalSample[i] - L).two_norm();

        	if (!isInside)                { rez.first.nFail_++; }
        	if (err > ABSOLUTE_TOLERANCE) { rez.second.nFail_++; }
        	rez.second.worst_ = std::max(rez.second.worst_, err);
            //std::cout << "  -- recovered vertex " << L << " from real vertex " << randomLocalSample[i] << " via global " << SimplexGeom.global(randomLocalSample[i]) << " positioning is " << is_inside << std::endl;
        }
    }
    return rez;
}


// ---------------------------------------------------------------------------------------------------------------------------
// 5. IS_INSIDE TEST
// This test generates a grid of global points just outside the boundary of the element
// the is_inside method should evaluate to false for all of them
// ---------------------------------------------------------------------------------------------------------------------------
template<class SimplexGeom, class LocalCoordinate, class GlobalCoordinate>
TestStruct test05_isoutside(
	bool verbose,
	const SimplexGeom & simplexGeom,
	int nSubentities)
{
	typedef typename SimplexGeom::ctype  ctype;
	static const int cdim  = SimplexGeom::coorddimension;
	static const int mydim = SimplexGeom::mydimension;
	static const int subdim = mydim - 1;

    // Define subentity geometries
    typedef FieldVector< ctype, subdim> SubLocalCoordinate;

    TestStruct rez(0, 0, 0.0);
    LocalCoordinate L;

	if ((mydim != cdim) || (mydim < 1) || (mydim > 3))  {
	    if (verbose) { std::cout << ": Global-to-Local functionality only available for tests of simplex geometries with matching dimensions" << std::endl; }
	}
	else if (mydim == 1) {
        // Check of points close to the corners but outside are really outside
        GlobalCoordinate c1 = simplexGeom.corner(0);
        GlobalCoordinate c2 = simplexGeom.corner(1);

        GlobalCoordinate delta = c2 - c1;
        delta *= 0.01;

        rez.nTot_ = 2;

        if (simplexGeom.local(c1 - delta, L))  { rez.nFail_++; }   // Check vertex just outside left corner of the edge
        if (simplexGeom.local(c2 + delta, L))  { rez.nFail_++; }   // Check vertex just outside right corner of the edge
	} else {
		// Produce points surrounding the element from the outside
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    int surround_grid_order = 10;
	    double displacement = 0.1;

        std::vector<SubLocalCoordinate > sublocal_coordinates = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, subdim>(surround_grid_order);

        std::vector<GlobalCoordinate> element_surround_points;

        // For each point on of the local grid, evaluate each subentity and its normal for that point
        // Thus compute point outside the element for each local grid point
        for (int biter = 0; biter < nSubentities; biter++)
        {
            for (unsigned int i = 0; i < sublocal_coordinates.size(); i++)
            {
            	LocalCoordinate  local  = Dune::CurvilinearGeometryHelper::coordinateInParent<ctype, subdim, mydim>(simplexGeom.type(), biter, sublocal_coordinates[i]);
                GlobalCoordinate global = simplexGeom.global(local);
                GlobalCoordinate normal = simplexGeom.subentityUnitNormal(biter, local);

                //GlobalCoordinate gp = subElement.global(sublocal_coordinates[i]);
                //GlobalCoordinate np = subElement.normal(sublocal_coordinates[i]);

                normal *= displacement;

                element_surround_points.push_back(global + normal);
            }
        }

        // Check if each of the surrounding points is indeed outside
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for (unsigned int i = 0; i < element_surround_points.size(); i++)
        {
        	rez.nTot_++;
        	bool isInside = simplexGeom.local(element_surround_points[i], L);
        	if (isInside)  { rez.nFail_++; }
            //std::cout << "-- for global vertex " << element_surround_points[i] << " recovered local vertex" << L << "and its is_inside = " << is_inside << std::endl;
        }
	}

	return rez;
}


// ---------------------------------------------------------------------------------------------------------------------------
// 6. INTEGRAL_TEST
// This test integrates over the element the basis functions given as polynomials of local coordinates
// CurvilinearGeometry should choose the analytical integration for matching dimensions and numerical for mismatching
// [FIXME] Introduce tests for the case of integrals integrating to 0
// ---------------------------------------------------------------------------------------------------------------------------
template<class SimplexGeom, class SimplexIntegrationHelper>
TestStruct test06_integration(bool verbose, const SimplexGeom & simplexGeom, int interpOrder, int functionOrder, int functType, double RELATIVE_TOLERANCE, double ACCURACY_GOAL)
{
	static const int cdim  = SimplexGeom::coorddimension;
	static const int mydim = SimplexGeom::mydimension;

	typedef typename SimplexGeom::LocalPolynomial  LocalPolynomial;

	TestStruct rez(0, 0, 0.0);

    if (interpOrder < functionOrder) {
        if (verbose) { std::cout << ": Integral-test: --Omitted because polynomial order too small" << std::endl; }
    } else
    {
        for (int bf_ord = 0; bf_ord <= 5; bf_ord++)
        {
        	LocalPolynomial basisP = BasisPolynomial<SimplexGeom>(bf_ord);

            double int_rez  = SimplexIntegrationHelper::integrateScalar(simplexGeom, basisP, RELATIVE_TOLERANCE, ACCURACY_GOAL);
            double int_true = integralResult<mydim, cdim>(bf_ord, functType);

            rez.nTot_++;
            double err = fabs((int_rez - int_true) / int_true);  // NOTE: Make sure that none of the integrands amounts to zero
            if (err > RELATIVE_TOLERANCE)  {
            	//std::cout << "Integrating polynomial " << basisP.to_string() << std::endl;
            	//std::cout << "Function type "          << functType << std::endl;
            	//std::cout << "result=" << int_rez << ", expected=" << int_true << std::endl;
            	//DUNE_THROW(MathError, "Unexpected test fail");

            	rez.nFail_++;
            }     // NOTE: It is important relative error requested from the integrator is at least as good as the real relative error
            rez.worst_ = std::max(rez.worst_, err);
        }
    }

    return rez;
}


// ---------------------------------------------------------------------------------------------------------------------------
// 7. Dot product integral test
// This test integrates over the element the polynomial vector basis function cdot the surface normal
// Only available for integrals which have a surface normal
// ---------------------------------------------------------------------------------------------------------------------------
template<class SimplexGeom, class SimplexIntegrationHelper>
TestStruct test07_dot_integration(bool verbose, const SimplexGeom & simplexGeom, int interpOrder, int functionOrder, int functType, double ABSOLUTE_TOLERANCE)
{
	static const int cdim  = SimplexGeom::coorddimension;
	static const int mydim = SimplexGeom::mydimension;
	typedef typename SimplexGeom::PolynomialGlobalCoordinate   PolynomialGlobalCoordinate;

	TestStruct rez(0, 0, 0.0);

    if (interpOrder < functionOrder) {
        if (verbose)  { std::cout << ": Integral-test: --Omitted because polynomial order too small" << std::endl; }
    }
    else if (((mydim == 1)&&(cdim == 2))||((mydim == 2)&&(cdim == 3)))
    {
        PolynomialGlobalCoordinate basisVectorP = BasisVectorDiagonal<SimplexGeom>();

        // Need a basis function which is a vector in global coordinates
        double int_rez  = SimplexIntegrationHelper::integrateAnalyticalDot(simplexGeom, basisVectorP);
        double int_true = integralDotResult<mydim, cdim>(functType);

        rez.nTot_++;
        double err = fabs((int_rez - int_true) / int_true);  // NOTE: Make sure that none of the integrands amounts to zero
        if (err > ABSOLUTE_TOLERANCE)  { rez.nFail_++; }
        rez.worst_ = std::max(rez.worst_, err);
    } else
    {
        if (verbose) { std::cout << ": Integral - Surface Dot Product Not Available for these dimensions " << std::endl; }
    }

    return rez;
}


// ************************************************************************************************
// IMPLEMENTING TEST METHODS
// ************************************************************************************************

bool testReport(const TestStruct & rez, std::string testname, bool verbose) {
	bool pass = (rez.nFail_ == 0);
	//if ((!pass) || verbose)  {
		std::cout << "   *** Test " << testname << " Fails=" << rez.nFail_ << "/" << rez.nTot_ << ", worst result=" << rez.worst_ << std::endl;
	//}
	return pass;
}


// [FIXME] RandomSample size 20 should be a user controlled constant from args, with default value 20
template<typename CurvGeom, typename Functor>
bool SimplexTest(bool verbose, Functor f, int testIndex, int interpOrder, int functionOrder, int functType, std::string f_name)
{
	if (verbose) { std::cout << " * testing order " << interpOrder << " behaviour for " << f_name << std::endl; }

	// Extract geometry templates
	typedef typename CurvGeom::ctype   ctype;
	static const int mydim = CurvGeom::mydimension;


	// Define derived types
	typedef typename CurvGeom::LocalCoordinate    		LocalCoordinate;
	typedef typename CurvGeom::GlobalCoordinate   		GlobalCoordinate;
    typedef std::vector< LocalCoordinate > 				LocalVectorVector;
    typedef std::vector< GlobalCoordinate > 			GlobalVectorVector;
    typedef typename CurvGeom::IntegrationHelper        SimplexIntegrationHelper;

    typedef typename Dune::Geo::ReferenceElements<ctype,mydim>::ReferenceElement	ReferenceElement;

    // Construct reference element
    //Dune::GeometryType simplexGeometryType( Dune::GenericGeometry::SimplexTopology< mydim >::type::id, mydim );

//    Dune::GeometryType simplexGeometryType;
//    simplexGeometryType.makeSimplex(mydim);

    Dune::GeometryType simplexGeometryType=Dune::GeometryTypes::simplex(mydim);


//    const Dune::ReferenceElement< ctype, mydim > refElement = Dune::ReferenceElements< ctype, mydim >::general( simplexGeometryType );
    const ReferenceElement& refElement = Dune::ReferenceElements< ctype, mydim >::general( simplexGeometryType );

    int nSubentities = refElement.size(1);

    // Construct local simplex grid and map it to global using given Functor
    LocalVectorVector local_vertices = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, mydim>(interpOrder);
    typedef typename CurvGeom::Vertices Vertices;
    Vertices global_vertices;
//    GlobalVectorVector global_vertices;
    for (unsigned int i = 0; i < local_vertices.size(); i++) { global_vertices.push_back(f(local_vertices[i])[0]); }

    // Construct a Curvilinear Geometry
    CurvGeom SimplexGeom(refElement, global_vertices, interpOrder);
//    auto SimplexGeom=typename CurvGeom.template CurvGeom<GlobalVectorVector>(refElement, global_vertices, interpOrder));
//    auto SimplexGeom=typename CurvGeom.template CurvGeom<GlobalVectorVector>(refElement, global_vertices, interpOrder);

    // Produce a set of sample points randomly sampled over the entity
    LocalVectorVector randomLocalSample = sampleRandom<mydim>(20);



    // ------------------------ Started Tests ------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------

    // Define integration constants
    double RELATIVE_TOLERANCE = 1.0e-5;
    double ACCURACY_GOAL = 1.0e-15;

    // Define precision constants
    double ABSOLUTE_TOLERANCE_1 = 1.0e-15;
    double ABSOLUTE_TOLERANCE_2 = 1.0e-12;
    double ABSOLUTE_TOLERANCE_3 = 1.0e-7;
    double ABSOLUTE_TOLERANCE_4 = 1.0e-7;
    double ABSOLUTE_TOLERANCE_7 = 1.0e-10;

    // Perform the actual tests
    bool testResult = false;
    switch(testIndex)
    {
    case 1 :
    {
    	TestStruct test1_pass = test01_corner(f, SimplexGeom,       refElement, ABSOLUTE_TOLERANCE_1);
    	testResult = testReport(test1_pass, "1", verbose);
    } break;
    case 2 :
    {
    	TestStruct test2_pass = test02_local_to_global(verbose, f, SimplexGeom,       randomLocalSample, interpOrder, functionOrder, ABSOLUTE_TOLERANCE_2);
    	testResult = testReport(test2_pass, "2", verbose);
    } break;
    case 3 :
    {
    	TestPair   test3_pass = test03_global_to_local(verbose, SimplexGeom,       global_vertices, local_vertices, ABSOLUTE_TOLERANCE_3);
    	bool rez3  = testReport(test3_pass.first,  "3-isinside", verbose);
    	     rez3 &= testReport(test3_pass.second, "3-global_to_local", verbose);
    	testResult = rez3;
    } break;
    case 4 :
    {
    	TestPair   test4_pass = test04_global_to_local(verbose, SimplexGeom,       randomLocalSample, ABSOLUTE_TOLERANCE_4);
    	bool rez4  = testReport(test4_pass.first,  "4-isinside", verbose);
    	     rez4 &= testReport(test4_pass.second, "4-global_to_local", verbose);
    	testResult = rez4;
    } break;
    case 5 :
    {
    	TestStruct test5_pass = test05_isoutside<CurvGeom, LocalCoordinate, GlobalCoordinate>(verbose, SimplexGeom, nSubentities);
    	testResult = testReport(test5_pass, "5", verbose);
    } break;
    case 6 :
    {
    	TestStruct test6_pass = test06_integration<CurvGeom, SimplexIntegrationHelper>(verbose, SimplexGeom, interpOrder, functionOrder, functType, RELATIVE_TOLERANCE, ACCURACY_GOAL);
    	testResult = testReport(test6_pass, "6", verbose);
    } break;
    case 7 :
    {
    	TestStruct test7_pass = test07_dot_integration<CurvGeom, SimplexIntegrationHelper>(verbose, SimplexGeom, interpOrder, functionOrder, functType, ABSOLUTE_TOLERANCE_7);
    	testResult = testReport(test7_pass, "7", verbose);
    } break;

    }

    // [TODO] Test for normal direction
    //
    // ************************************************************
    // ALL TESTS FINISHED
    // ************************************************************

    return testResult;
}


template<class ctype, int mydim, int cdim>
bool SimplexTestWrapper(bool verbose)
{
    std::cout << "----- started testing simplex with mydim=" << mydim << " and world_dim=" << cdim << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------" << std::endl;

    typedef Dune::CurvilinearGeometry<ctype, mydim, cdim>        CurvGeom;
    typedef Dune::CachedCurvilinearGeometry<ctype, mydim, cdim>  CachedCurvGeom;

    typedef myFunctorIdentity  <double, mydim, cdim>  FiD;
    typedef myFunctorLinear    <double, mydim, cdim>  Flin;
    typedef myFunctorNonlinear1<double, mydim, cdim>  Fquad1;

    static const int MIN_TEST = 1;  // FIXME 1
    static const int MAX_TEST = 7;  // FIXME 7
    static const int N_INTERP_ORDER = 5;

    bool rez = true;

    /*
    std::cout << "[Running CurvilinearGeometry tests..." << std::endl;
    for (int iTest = MIN_TEST; iTest <= MAX_TEST; iTest++)
    {

        bool test_pass = true;
        for (int ord = 1; ord <= N_INTERP_ORDER; ord++)
        {
        	test_pass &= SimplexTest<CurvGeom, FiD>   (verbose, FiD(),    iTest, ord, 1, 1, "identity map");
        	test_pass &= SimplexTest<CurvGeom, Flin>  (verbose, Flin(),   iTest, ord, 1, 2, "linear map");
        	test_pass &= SimplexTest<CurvGeom, Fquad1>(verbose, Fquad1(), iTest, ord, 2, 3, "quadratic map");
        }

        if (!test_pass)  { std::cout << "--Warning: Have Failed tests" << std::endl; }
        rez &= test_pass;
    }
    std::cout << "... Finished CurvilinearGeometry tests]" << std::endl;

     */

    std::cout << "[Running Cached CurvilinearGeometry tests..." << std::endl;
    for (int iTest = MIN_TEST; iTest <= MAX_TEST; iTest++)
    {
        bool test_pass_cache = true;
        for (int ord = 1; ord <= N_INTERP_ORDER; ord++)
        {
        	test_pass_cache &= SimplexTest<CachedCurvGeom, FiD>   (verbose, FiD(),    iTest, ord, 1, 1, "identity map");
        	test_pass_cache &= SimplexTest<CachedCurvGeom, Flin>  (verbose, Flin(),   iTest, ord, 1, 2, "linear map");
        	test_pass_cache &= SimplexTest<CachedCurvGeom, Fquad1>(verbose, Fquad1(), iTest, ord, 2, 3, "quadratic map");
        }

        if (!test_pass_cache)  { std::cout << "--Warning: Have Failed tests" << std::endl; }
        rez &= test_pass_cache;
    }
    std::cout << "... Finished Cached CurvilinearGeometry tests]" << std::endl;

    return rez;
}





int main ( int argc, char **argv )
{
  srand (time(NULL));

  bool pass = true;
  bool verbose = false;

  std::cout << ">>> Checking ctype = double" << std::endl;
  pass &= SimplexTestWrapper<double, 1, 1>(verbose);
  pass &= SimplexTestWrapper<double, 1, 2>(verbose);
  pass &= SimplexTestWrapper<double, 1, 3>(verbose);
  pass &= SimplexTestWrapper<double, 2, 2>(verbose);
  pass &= SimplexTestWrapper<double, 2, 3>(verbose);
  pass &= SimplexTestWrapper<double, 3, 3>(verbose);


/*  std::cout << "-----Joke test------" << std::endl;
  Dune::FieldVector< double, 2 > JokeTestVec1;
  Dune::FieldVector< double, 2 > JokeTestVec2;

  JokeTestVec1[0] = 1.43573e-16;  JokeTestVec1[1] = 1;
  JokeTestVec2[0] = 1.43573e-14;  JokeTestVec2[1] = 1;


  Dune::GeometryType simplexG( Dune::GenericGeometry::SimplexTopology< 2 >::type::id, 2 );
  //Dune::ReferenceElement< double, 2 > simplexR = Dune::ReferenceElements<double, 2>::general( simplexG );
  std::cout << "test result: " << Dune::ReferenceElements<double, 2>::general( simplexG ).checkInside(JokeTestVec1) << std::endl;
  std::cout << "test result: " << Dune::ReferenceElements<double, 2>::general( simplexG ).checkInside(JokeTestVec2) << std::endl;*/


  std::cout << "Totall test pass = " << pass << std::endl;

  return 0;
}
