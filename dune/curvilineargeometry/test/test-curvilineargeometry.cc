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
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>
#include <dune/curvilineargeometry/integration/numericalrecursiveinterpolationintegrator.hh>

#include <dune/curvilineargeometry/curvilineargeometry.hh>





using namespace Dune;

typedef Dune::FieldVector< double, 1 > FieldVector1D;
typedef Dune::FieldVector< double, 2 > FieldVector2D;
typedef Dune::FieldVector< double, 3 > FieldVector3D;

typedef std::vector<FieldVector1D> FieldVectorVector1D;
typedef std::vector<FieldVector2D> FieldVectorVector2D;
typedef std::vector<FieldVector3D> FieldVectorVector3D;


// ************************************************************************************************
// IMPLEMENTING AUXILIARY METHODS
// ************************************************************************************************

// Generates uniform random numbers in interval [a,b]
double randomReal(double a, double b) { return a + (b - a)*(double(rand()) / RAND_MAX); }

// Constructs a random grid over a simplex with a given number of samples
template <int dim>
std::vector<FieldVector< double, dim > > sampleRandom(int sampleNo) {
    std::vector<FieldVector< double, dim > > tmprez;
    return tmprez;
}

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

template<class ctype, int mydim>
Dune::FieldVector< ctype, mydim > localCorner(int i)
{
    Dune::FieldVector< ctype, mydim > rez;

    switch (mydim)
    {
    case 1:
        switch (i)
        {
        case 0: rez[0] = 0.0;  break;
        case 1: rez[0] = 1.0;  break;
        } break;
    case 2:
        switch (i)
        {
        case 0: rez[0] = 0.0;  rez[1] = 0.0;  break;
        case 1: rez[0] = 1.0;  rez[1] = 0.0;  break;
        case 2: rez[0] = 0.0;  rez[1] = 1.0;  break;
        } break;
    case 3:
        switch (i)
        {
        case 0: rez[0] = 0.0;  rez[1] = 0.0;  rez[2] = 0.0;  break;
        case 1: rez[0] = 1.0;  rez[1] = 0.0;  rez[2] = 0.0;  break;
        case 2: rez[0] = 0.0;  rez[1] = 1.0;  rez[2] = 0.0;  break;
        case 3: rez[0] = 0.0;  rez[1] = 0.0;  rez[2] = 1.0;  break;
        } break;
    }
    return rez;
}


// ************************************************************************************************
// IMPLEMENTING FUNCTORS - EXAMPLE LOCAL-TO-GLOBAL MAPS
// ************************************************************************************************

template<class ctype, int mydim, int cdim> struct myFunctorIdentity { };
template<class ctype, int mydim, int cdim> struct myFunctorLinear { };
template<class ctype, int mydim, int cdim> struct myFunctorNonlinear1 { };
template<class ctype, int mydim, int cdim> struct myFunctorNonlinear2 { };

template<> struct myFunctorIdentity<double, 1, 1> { FieldVector1D operator()(const FieldVector1D & in) { return initFieldVector(in[0]); }  };
template<> struct myFunctorIdentity<double, 1, 2> { FieldVector2D operator()(const FieldVector1D & in) { return initFieldVector(in[0], 0); }  };
template<> struct myFunctorIdentity<double, 1, 3> { FieldVector3D operator()(const FieldVector1D & in) { return initFieldVector(in[0], 0, 0); }  };
template<> struct myFunctorIdentity<double, 2, 2> { FieldVector2D operator()(const FieldVector2D & in) { return initFieldVector(in[0], in[1]); }  };
template<> struct myFunctorIdentity<double, 2, 3> { FieldVector3D operator()(const FieldVector2D & in) { return initFieldVector(in[0], in[1], 0); }  };
template<> struct myFunctorIdentity<double, 3, 3> { FieldVector3D operator()(const FieldVector3D & in) { return initFieldVector(in[0], in[1], in[2]); }  };

template<> struct myFunctorLinear<double, 1, 1> { FieldVector1D operator()(const FieldVector1D & in) { return initFieldVector(1.0 + 2.0 * in[0]); }  };
template<> struct myFunctorLinear<double, 1, 2> { FieldVector2D operator()(const FieldVector1D & in) { return initFieldVector(2.0 * in[0], 3.0 * in[0]); }  };
template<> struct myFunctorLinear<double, 1, 3> { FieldVector3D operator()(const FieldVector1D & in) { return initFieldVector(2.0 * in[0], 0.5 + 3.0 * in[0], 5.0 * in[0]); }  };
template<> struct myFunctorLinear<double, 2, 2> { FieldVector2D operator()(const FieldVector2D & in) { return initFieldVector(1.0 + in[0], in[0] + in[1]); }  };
template<> struct myFunctorLinear<double, 2, 3> { FieldVector3D operator()(const FieldVector2D & in) { return initFieldVector(in[1], 3.0 * in[0], in[0] + in[1]); }  };
template<> struct myFunctorLinear<double, 3, 3> { FieldVector3D operator()(const FieldVector3D & in) { return initFieldVector(in[0] + in[1], in[1] + in[2], in[2] + in[0]); }  };

template<> struct myFunctorNonlinear1<double, 1, 1> { FieldVector1D operator()(const FieldVector1D & in) { return initFieldVector(in[0] * in[0]); }  };
template<> struct myFunctorNonlinear1<double, 1, 2> { FieldVector2D operator()(const FieldVector1D & in) { return initFieldVector(in[0], in[0] * in[0]); }  };
template<> struct myFunctorNonlinear1<double, 1, 3> { FieldVector3D operator()(const FieldVector1D & in) { return initFieldVector(in[0], in[0] * in[0], 2.0); }  };
template<> struct myFunctorNonlinear1<double, 2, 2> { FieldVector2D operator()(const FieldVector2D & in) { return initFieldVector(in[0]*in[0], in[1] * in[1]); }  };
template<> struct myFunctorNonlinear1<double, 2, 3> { FieldVector3D operator()(const FieldVector2D & in) { return initFieldVector(in[1] * in[1], in[0] * in[0], in[0] * in[1]); }  };
template<> struct myFunctorNonlinear1<double, 3, 3> { FieldVector3D operator()(const FieldVector3D & in) { return initFieldVector(in[0] * in[0], in[1] * in[1], in[2] * in[2]); }  };





// Returns one of 5 possible polynomials such that polynomial order = index
template<class ctype, int mydim>
polynomial<ctype, mydim> BasisPolynomial1D(int index)
{
    polynomial<ctype, mydim> rez(polySummand(1, 0));

    if (index > 0) { rez.append(polySummand(2, 1)); }
    if (index > 1) { rez.append(polySummand(3, 2)); }
    if (index > 2) { rez.append(polySummand(4, 3)); }
    if (index > 3) { rez.append(polySummand(5, 4)); }
    if (index > 4) { rez.append(polySummand(6, 5)); }

    return rez;
}

// Returns one of 5 possible polynomials such that polynomial order = index
template<class ctype, int mydim>
polynomial<ctype, mydim> BasisPolynomial2D(int index)
{
    polynomial<ctype, mydim> rez(polySummand(1, 0, 0));

    if (index > 0) { rez.append(polySummand(2, 1, 0));  rez.append(polySummand(2, 0, 1)); }
    if (index > 1) { rez.append(polySummand(3, 2, 0));  rez.append(polySummand(3, 0, 2)); rez.append(polySummand(1, 1, 1)); }
    if (index > 2) { rez.append(polySummand(4, 3, 0));  rez.append(polySummand(4, 0, 3)); rez.append(polySummand(1, 1, 2)); }
    if (index > 3) { rez.append(polySummand(5, 4, 0));  rez.append(polySummand(5, 0, 4)); rez.append(polySummand(1, 1, 3)); }
    if (index > 4) { rez.append(polySummand(6, 5, 0));  rez.append(polySummand(6, 0, 5)); rez.append(polySummand(1, 1, 4)); }

    return rez;
}

// Returns one of 5 possible polynomials such that polynomial order = index
template<class ctype, int mydim>
polynomial<ctype, mydim> BasisPolynomial3D(int index)
{
    polynomial<ctype, mydim> rez(polySummand(1, 0, 0, 0));

    if (index > 0) { rez.append(polySummand(2, 1, 0, 0));  rez.append(polySummand(2, 0, 1, 0));  rez.append(polySummand(2, 0, 0, 1)); }
    if (index > 1) { rez.append(polySummand(3, 2, 0, 0));  rez.append(polySummand(3, 0, 2, 0));  rez.append(polySummand(3, 0, 0, 2)); rez.append(polySummand(1, 1, 1, 0)); }
    if (index > 2) { rez.append(polySummand(4, 3, 0, 0));  rez.append(polySummand(4, 0, 3, 0));  rez.append(polySummand(4, 0, 0, 3)); rez.append(polySummand(1, 1, 1, 1)); }
    if (index > 3) { rez.append(polySummand(5, 4, 0, 0));  rez.append(polySummand(5, 0, 4, 0));  rez.append(polySummand(5, 0, 0, 4)); rez.append(polySummand(1, 1, 1, 2)); }
    if (index > 4) { rez.append(polySummand(6, 5, 0, 0));  rez.append(polySummand(6, 0, 5, 0));  rez.append(polySummand(6, 0, 0, 5)); rez.append(polySummand(1, 1, 1, 3)); }

    return rez;
}

template<class ctype, int mydim>
polynomial<ctype, mydim> BasisPolynomial(int index)
{
    switch (mydim)
    {
    case 1: return BasisPolynomial1D<ctype, mydim>(index);  break;
    case 2: return BasisPolynomial2D<ctype, mydim>(index);  break;
    case 3: return BasisPolynomial3D<ctype, mydim>(index);  break;
    }
}


// Construct a mydim+1 polynomial vector to test the Surface Dot Product Integral
template<class ctype, int mydim>
std::vector< polynomial<ctype, mydim> > BasisVectorDiagonal()
{
    std::vector< polynomial<ctype, mydim> > rez;

    switch (mydim)
    {
    case 1:
        rez.push_back(polynomial<ctype, mydim> (polySummand(1, 1)));
        rez.push_back(polynomial<ctype, mydim> (polySummand(1, 1)));
        break;
    case 2:
        rez.push_back(polynomial<ctype, mydim> (polySummand(1, 1, 0)));
        rez.push_back(polynomial<ctype, mydim> (polySummand(1, 0, 1)));
        rez.push_back(polynomial<ctype, mydim> (polySummand(1, 1, 1)));
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
// IMPLEMENTING TEST METHODS
// ************************************************************************************************

template<class ctype, int mydim, int cdim, typename Functor>
bool SimplexTest(Functor f, int ord, int f_type, int f_order, std::string f_name)
{
    typedef Dune::FieldVector< ctype, mydim > LocalVector;
    typedef Dune::FieldVector< ctype, cdim > GlobalVector;
    typedef std::vector< LocalVector > LocalVectorVector;
    typedef std::vector< GlobalVector > GlobalVectorVector;

    typedef CurvilinearGeometry< ctype, mydim - 1, cdim> SubentityGeometry;
    typedef CachedCurvilinearGeometry< ctype, mydim - 1, cdim> SubentityCachedGeometry;
    typedef std::vector<SubentityGeometry> SubentityGeometryVector;
    typedef std::vector<SubentityCachedGeometry> SubentityCachedGeometryVector;
    typedef FieldVector< ctype, mydim - 1 > SubLocalCoordinate;


    std::cout << " * testing order " << ord << " behaviour for " << f_name << std::endl;

    // ************************************************************
    // Produce Geometry Type
    // ************************************************************

    Dune::GeometryType simplexGeometryType( Dune::GenericGeometry::SimplexTopology< mydim >::type::id, mydim );
    const Dune::ReferenceElement< ctype, mydim > &refElement = Dune::ReferenceElements< ctype, mydim >::general( simplexGeometryType );

    int nSubentities = refElement.size(1);

    // ************************************************************
    // Produce local simplex grid and map it to global using given Functor
    // ************************************************************

    LocalVectorVector local_vertices = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, mydim>(ord);
    GlobalVectorVector global_vertices;
    for (int i = 0; i < local_vertices.size(); i++) {global_vertices.push_back(f(local_vertices[i])); }

    // ************************************************************
    // Produce a Curvilinear Geometry
    // ************************************************************

    Dune::CurvilinearGeometry< ctype, mydim, cdim >  SimplexGeom(refElement, global_vertices, ord);

    // ************************************************************
    // Produce a Cached Curvilinear Geometry
    // ************************************************************

    Dune::CachedCurvilinearGeometry< ctype, mydim, cdim >  SimplexGeomCached(refElement, global_vertices, ord);




    // ------------------------ Started Tests ------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------
    bool test1_pass = true;   bool test1_cached_pass = true;
    bool test2_pass = true;   bool test2_cached_pass = true;
    bool test3_pass = true;   bool test3_cached_pass = true;
    bool test4_pass = true;   bool test4_cached_pass = true;
    bool test5_pass = true;   bool test5_cached_pass = true;
    bool test6_pass = true;   bool test6_cached_pass = true;
    bool test7_pass = true;   bool test7_cached_pass = true;

    bool test3_isinside_pass = true;   bool test3_isinside_cached_pass = true;
    bool test4_isinside_pass = true;   bool test4_isinside_cached_pass = true;


    // ---------------------------------------------------------------------------------------------------------------------------
    // 1. CORNER-TEST
    // This test checks if the corner coordinates returned by CurvilinearGeometry correspond to the ones used in initializing it.
    // ---------------------------------------------------------------------------------------------------------------------------
    for (int i = 0; i < mydim+1; i++) {
        double tmpTolerance = 1.0e-15;
        GlobalVector tmpCorner = f(localCorner<double, mydim>(i));

        test1_pass &= (SimplexGeom.corner(i) - tmpCorner).two_norm() < tmpTolerance;
        test1_cached_pass &= (SimplexGeomCached.corner(i) - tmpCorner).two_norm() < tmpTolerance;
    }
    std::cout << ": Corners-test pass = " << test1_pass << " cached-pass = " << test1_cached_pass << std::endl;


    // ---------------------------------------------------------------------------------------------------------------------------
    // 2. LOCAL->GLOBAL TEST
    // This test samples the example global map over the simplex, and checks if these samples match with the CurvilinearGeometry.global()
    // ----------------------------------------------------------------------------------------------------------------------------------
    if (ord < f_order) {
        // Only count as error if interpolatory polynomial is of sufficient order
        std::cout << ": Local-To-Global-test: --Omitted because polynomial order too small" << std::endl;
    }
    else
    {
        double tmpTolerance = 1.0e-12;
        LocalVectorVector randomLocalSample = sampleRandom<mydim>(20);

        for (int i = 0; i < randomLocalSample.size(); i++)
        {
            double errTmp = (SimplexGeom.global(randomLocalSample[i]) - f(randomLocalSample[i])).two_norm();
            double errTmpCached = (SimplexGeomCached.global(randomLocalSample[i]) - f(randomLocalSample[i])).two_norm();

            test2_pass &= errTmp < tmpTolerance;
            test2_cached_pass &= errTmpCached < tmpTolerance;
        }
        std::cout << ": Local-To-Global-test pass = " << test2_pass << " cached-pass = " << test2_cached_pass << std::endl ;
    }


    // ---------------------------------------------------------------------------------------------------------------------------
    // 3. GLOBAL->LOCAL TEST
    // This test checks if local coordinates of the global interpolation vertices are recovered
    // ----------------------------------------------------------------------------------------------------------------------------------
    if (mydim != cdim) { std::cout << ": Global-to-Local functionality not available for mismatching mydim and cdim" << std::endl; }
    else
    {
        double tmpTolerance = 1.0e-7;

        for (int i = 0; i < global_vertices.size(); i++) {
            LocalVector L;
            LocalVector L_c;

            // Order important, local() method initializes L
            test3_isinside_pass &= SimplexGeom.local(global_vertices[i], L);
            test3_isinside_cached_pass &= SimplexGeomCached.local(global_vertices[i], L_c);

            test3_pass &= (L - local_vertices[i]).two_norm() < tmpTolerance;
            test3_cached_pass &= (L_c - local_vertices[i]).two_norm() < tmpTolerance;

            //std::cout << "  -- recovered vertex " << L << " from real vertex " << local_vertices[i] << " via global " << global_vertices[i] << " positioning is " << is_inside << std::endl;
        }
        std::cout << ": Global-To-Local-test 1 pass = " << test3_pass << " cached-pass = " << test3_cached_pass;
        std::cout << ", is_inside test pass = " << test3_isinside_pass << " cached-pass = " << test3_isinside_cached_pass << std::endl ;
    }


    // ---------------------------------------------------------------------------------------------------------------------------
    // 4. GLOBAL->LOCAL TEST 2
    // Test if LOCAL->GLOBAL->LOCAL for random sample inside element is preserved
    // ----------------------------------------------------------------------------------------------------------------------------------
    if (mydim != cdim) { std::cout << ": Global-to-Local functionality not available for mismatching mydim and cdim" << std::endl; }
    else
    {
        double tmpTolerance = 1.0e-7;
        LocalVectorVector randomLocalSample = sampleRandom< mydim >(20);

        for (int i = 0; i < randomLocalSample.size(); i++)
        {
            LocalVector L;
            LocalVector L_c;

            // Order important, local() method initializes L
            test4_isinside_pass &= SimplexGeom.local(SimplexGeom.global(randomLocalSample[i]), L);
            test4_isinside_cached_pass &= SimplexGeomCached.local(SimplexGeomCached.global(randomLocalSample[i]), L_c);

            test4_pass &= (randomLocalSample[i] - L).two_norm() < tmpTolerance;
            test4_cached_pass &= (randomLocalSample[i] - L_c).two_norm() < tmpTolerance;

            //std::cout << "  -- recovered vertex " << L << " from real vertex " << randomLocalSample[i] << " via global " << SimplexGeom.global(randomLocalSample[i]) << " positioning is " << is_inside << std::endl;
        }
        std::cout << ": Global-To-Local-test 2 pass = " << test4_pass << " cached-pass = " << test4_cached_pass;
        std::cout << ", is_inside test pass = " << test4_isinside_pass << " cached-pass = " << test4_isinside_cached_pass << std::endl ;
    }


    // ---------------------------------------------------------------------------------------------------------------------------
    // 5. IS_INSIDE TEST
    // This test generates a grid of global points just outside the boundary of the element
    // the is_inside method should evaluate to false for all of them
    // ----------------------------------------------------------------------------------------------------------------------------------

    if ((mydim != cdim) || (mydim < 1) || (mydim > 3))
    {
        std::cout << ": Global-to-Local functionality only available for tests of simplex geometries with matching dimensions" << std::endl;
    }
    else if (mydim == 1)
    {
        // Check of points close to the corners but outside are really outside
        GlobalVector c1 = SimplexGeom.corner(0);
        GlobalVector c2 = SimplexGeom.corner(1);

        GlobalVector delta = c2 - c1;
        delta *= 0.01;

        LocalVector L;
        test5_pass &= !SimplexGeom.local(c1 - delta, L);  test5_cached_pass &= !SimplexGeomCached.local(c1 - delta, L);
        test5_pass &= !SimplexGeom.local(c2 + delta, L);  test5_cached_pass &= !SimplexGeomCached.local(c2 + delta, L);

        std::cout << ": IsOutside test pass = " << test5_pass << " cached-pass = " << test5_cached_pass << std::endl;
    }
    else
    {
        // Produce points surrounding the element from the outside
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        int surround_grid_order = 10;
        double displacement = 0.1;

        std::vector<SubLocalCoordinate > sublocal_coordinates = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, mydim-1>(surround_grid_order);

        GlobalVectorVector element_surround_points;
        GlobalVectorVector element_surround_points_cached;

        // For each point on of the local grid, evaluate each subentity and its normal for that point
        // Thus compute point outside the element for each local grid point

        //std::cout << " creating normals " << std::endl;
        for (int biter = 0; biter < nSubentities; biter++)
        {
        	SubentityGeometry subElement = SimplexGeom.template subentityGeometry < mydim - 1>(biter);
        	SubentityCachedGeometry subCachedElement = SimplexGeomCached.template subentityCachedGeometry < mydim - 1>(biter);

            for (int i = 0; i < sublocal_coordinates.size(); i++)
            {
                GlobalVector gp = subElement.global(sublocal_coordinates[i]);
                GlobalVector gp_cached = subCachedElement.global(sublocal_coordinates[i]);

                GlobalVector np = subElement.normal(sublocal_coordinates[i]);
                GlobalVector np_cached = subCachedElement.normal(sublocal_coordinates[i]);

                np *= displacement;
                np_cached *= displacement;

                element_surround_points.push_back(gp + np);
                element_surround_points_cached.push_back(gp_cached + np_cached);
            }
        }

        // Check if each of the surrounding points is indeed outside
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for (int i = 0; i < element_surround_points.size(); i++)
        {
            LocalVector L;
            test5_pass &= !SimplexGeom.local(element_surround_points[i], L);
            test5_cached_pass &= !SimplexGeomCached.local(element_surround_points_cached[i], L);

            //std::cout << "-- for global vertex " << element_surround_points[i] << " recovered local vertex" << L << "and its is_inside = " << is_inside << std::endl;
        }
        std::cout << ": IsOutside test pass = " << test5_pass << " cached-pass = " << test5_cached_pass << std::endl;
    }

    // ---------------------------------------------------------------------------------------------------------------------------
    // INTEGRAL_TEST
    // This test integrates over the element the basis functions given as polynomials of local coordinates
    // CurvilinearGeometry should choose the analytical integration for matching dimensions and numerical for mismatching
    // ---------------------------------------------------------------------------------------------------------------------------
    if (ord < f_order) {
        std::cout << ": Integral-test: --Omitted because polynomial order too small" << std::endl;
    } else
    {
        double tmpTolerance = 1.0e-5;

        for (int bf_ord = 0; bf_ord <= 5; bf_ord++)
        {
            polynomial<ctype, mydim> basisP = BasisPolynomial<ctype, mydim>(bf_ord);

            double int_rez = SimplexGeom.integrateScalar(basisP, tmpTolerance);
            double int_rez_cached = SimplexGeomCached.integrateScalar(basisP, tmpTolerance);

            double int_true = integralResult<mydim, cdim>(bf_ord, f_type);

            test6_pass &= fabs((int_rez - int_true) / int_true) < tmpTolerance;
            test6_cached_pass &= fabs((int_rez_cached - int_true) / int_true) < tmpTolerance;
            //std::cout << ": Integral of order " << bf_ord << " over entity dim " << mydim << " and order " << ord << " evaluated to ";
            //std::cout << int_rez << " which has error: " << int_rez - integralResults<mydim, cdim>(bf_ord, f_type) << std::endl;
        }
    }
    std::cout << ": Integral test pass = " << test6_pass << " cached-pass = " << test6_cached_pass << std::endl;

    // ---------------------------------------------------------------------------------------------------------------------------
    // Dot product integral test
    // This test integrates over the element the polynomial vector basis function cdot the surface normal
    // Only available for integrals which have a surface normal
    // ---------------------------------------------------------------------------------------------------------------------------
    if (ord < f_order) {
        std::cout << ": Integral-test: --Omitted because polynomial order too small" << std::endl;
    }
    else if (((mydim == 1)&&(cdim == 2))||((mydim == 2)&&(cdim == 3)))
    {
        double tmpTolerance = 1.0e-10;

        std::vector<polynomial<ctype, mydim>> basisVectorP = BasisVectorDiagonal<ctype, mydim>();

        // Need a basis function which is a vector in global coordinates
        double int_rez = SimplexGeom.integrateAnalyticalDot(basisVectorP);
        double int_rez_cached = SimplexGeomCached.integrateAnalyticalDot(basisVectorP);

        double int_true = integralDotResult<mydim, cdim>(f_type);

        test7_pass &= fabs((int_rez - int_true) / int_true) < tmpTolerance;
        test7_cached_pass &= fabs((int_rez_cached - int_true) / int_true) < tmpTolerance;

        std::cout << ": Integral - Surface Dot Product - pass =  " << test7_pass << " cached-pass = " << test7_cached_pass << std::endl;
    } else
    {
        std::cout << ": Integral - Surface Dot Product Not Available for these dimensions " << std::endl;
    }



    // Test for normal direction
    //
    // ************************************************************
    // ALL TESTS FINISHED
    // ************************************************************


}


template<class ctype, int mydim, int cdim>
bool SimplexTestWrapper()
{
    std::cout << "----- started testing simplex with mydim=" << mydim << " and world_dim=" << cdim << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------" << std::endl;

    for (int ord = 1; ord <= 5; ord++)
    {
        SimplexTest<ctype, mydim, cdim>(myFunctorIdentity<double, mydim, cdim>(), ord, 1, 1, "identity map");
        SimplexTest<ctype, mydim, cdim>(myFunctorLinear<double, mydim, cdim>(), ord, 2, 1, "linear map");
        SimplexTest<ctype, mydim, cdim>(myFunctorNonlinear1<double, mydim, cdim>(), ord, 3, 2, "quadratic map");
    }
    return true;
}





int main ( int argc, char **argv )
{
  srand (time(NULL));

  bool pass = true;

  std::cout << ">>> Checking ctype = double" << std::endl;
  pass &= SimplexTestWrapper<double, 1, 1>();
  pass &= SimplexTestWrapper<double, 1, 2>();
  pass &= SimplexTestWrapper<double, 1, 3>();
  pass &= SimplexTestWrapper<double, 2, 2>();
  pass &= SimplexTestWrapper<double, 2, 3>();
  pass &= SimplexTestWrapper<double, 3, 3>();


/*  std::cout << "-----Joke test------" << std::endl;
  Dune::FieldVector< double, 2 > JokeTestVec1;
  Dune::FieldVector< double, 2 > JokeTestVec2;

  JokeTestVec1[0] = 1.43573e-16;  JokeTestVec1[1] = 1;
  JokeTestVec2[0] = 1.43573e-14;  JokeTestVec2[1] = 1;


  Dune::GeometryType simplexG( Dune::GenericGeometry::SimplexTopology< 2 >::type::id, 2 );
  //Dune::ReferenceElement< double, 2 > simplexR = Dune::ReferenceElements<double, 2>::general( simplexG );
  std::cout << "test result: " << Dune::ReferenceElements<double, 2>::general( simplexG ).checkInside(JokeTestVec1) << std::endl;
  std::cout << "test result: " << Dune::ReferenceElements<double, 2>::general( simplexG ).checkInside(JokeTestVec2) << std::endl;*/



  return (pass ? 0 : 1);
}
