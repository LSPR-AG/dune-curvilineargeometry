/*******************************************************************
 * Curvilinear Geometry
 * 
 * author: Aleksejs Fomins
 * date: 01.09.2014 - created
 * 
 * description:
 * Provides the geometric structure of mesh elements given by lagrange polynomial interpolation
 * - Each element is defined by GeometryType, Interpolation Order, and Interpolatory Vertex Vector.
 * - Currently capable of working with simplices of orders 1-5 in dim 1-3
 * 
 *******************************************************************/



#ifndef DUNE_GEOMETRY_CURVILINEARGEOMETRY_HH
#define DUNE_GEOMETRY_CURVILINEARGEOMETRY_HH

#include <cassert>
#include <limits>
#include <vector>
#include <math.h>
#include <config.h>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>
#include <dune/curvilineargeometry/interpolation/differentialhelper.hh>
#include <dune/curvilineargeometry/interpolation/pointlocation.hh>
#include <dune/curvilineargeometry/integration/integrationhelper.hh>



namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ctype, int dim >
  class ReferenceElement;

  template< class ctype, int dim >
  struct ReferenceElements;



  // CurvilinearGeometryTraits
  // -------------------------

  /** \brief default traits class for CurvilinearGeometry
   *
   *  The CurvilinearGeometry (and CachedCurvilinearGeometry) allow tweaking
   *  some implementation details through a traits class.
   *
   *  This structure provides the default values.
   *
   *  \tparam  ct  coordinate type
   */
  template< class ct >
  struct CurvilinearGeometryTraits
  {
    typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ct > > MatrixHelper;

    static const unsigned int QUADRATURE_NORM_TYPE = Dune::QUADRATURE_NORM_L2;

    /** \brief tolerance to numerical algorithms */
    static ct tolerance () { return ct( 16 ) * std::numeric_limits< ct >::epsilon(); }

    template< int mydim, int cdim >
    struct VertexStorage
    {
      typedef std::vector< FieldVector< ct, cdim > > Type;
    };

    template< int dim >
    struct hasSingleGeometryType
    {
      static const bool v = true;
      static const unsigned int topologyId = ~0u;
    };
  };


  // CurvilinearGeometry
  // -------------------

  /** \brief generic geometry implementation based on corner coordinates
   *
   *  Based on the recursive definition of the reference elements, the
   *  CurvilinearGeometry provides a generic implementation of a geometry given
   *  the corner coordinates.
   *
   *  The geometric mapping is multilinear in the classical sense only in the
   *  case of cubes; for simplices it is linear.
   *  The name is still justified, because the mapping satisfies the important
   *  property of begin linear along edges.
   *
   *  \tparam  ct      coordinate type
   *  \tparam  mydim   geometry dimension
   *  \tparam  cdim    coordinate dimension
   *  \tparam  Traits  traits allowing to tweak some implementation details
   *                   (optional)
   *
   *  The requirements on the traits are documented along with their default,
   *  CurvilinearGeometryTraits.
   */
  template< class ct, int mydim, int cdim, class Traits = CurvilinearGeometryTraits< ct > >
  class CurvilinearGeometry
  {
    typedef CurvilinearGeometry< ct, mydim, cdim, Traits > This;

  public:

    typedef ct ctype;                           //! coordinate type
    static const int mydimension= mydim;        //! geometry dimension
    static const int coorddimension = cdim;     //! coordinate dimension

    // Plain data types used
    // [TODO] Move to traits
    typedef typename Dune::CurvilinearGeometryHelper::InternalIndexType         InternalIndexType;
    typedef typename Dune::CurvilinearGeometryHelper::InterpolatoryOrderType    InterpolatoryOrderType;

    typedef FieldVector< ctype, mydimension > LocalCoordinate;        //! type of local coordinates
    typedef FieldVector< ctype, coorddimension > GlobalCoordinate;    //! type of global coordinates

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

    //! type of jacobian inverse transposed
    class JacobianInverseTransposed;

    //! type of reference element
    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;

    //! type of element interpolator
    typedef CurvilinearElementInterpolator <ct, mydim, cdim> ElementInterpolator;


    // Define Polynomial and associated classes
    typedef typename PolynomialTraits<ctype>::Monomial  Monomial;
    typedef Polynomial<ctype, mydimension>              LocalPolynomial;
    typedef std::vector<LocalPolynomial>                PolynomialVector;

    // Define analytic differential operators over this geometry
    typedef typename Dune::DifferentialHelper::JacobianDeterminantAnalytical<This, cdim, mydim>         JacDetAnalytical;
    typedef typename Dune::DifferentialHelper::NormalIntegrationElementAnalytical<This, cdim, mydim>    IntElemNormalAnalytical;
    typedef typename Dune::DifferentialHelper::IntegrationElementSquaredAnalytical <This, cdim, mydim>  IntElemSquaredAnalytical;

    // Define integration routines for this geometry
    typedef typename Dune::IntegralHelper<This, Traits::QUADRATURE_NORM_TYPE>  IntegrationHelper;


    // Define tests that help to determine if a point is inside an entity
    typedef typename Dune::CurvilinearPointLocation::FarPointTest<This, cdim, mydim>      FarPointTest;
    typedef typename Dune::CurvilinearPointLocation::BarycentricTest<This, cdim, mydim>   BarycentricTest;


  private:

    static const bool hasSingleGeometryType = Traits::template hasSingleGeometryType< mydimension >::v;
    typedef typename Traits::template VertexStorage< mydimension, coorddimension >::Type::const_iterator VertexIterator;

  protected:
    typedef typename Traits::MatrixHelper MatrixHelper;
    typedef Dune::ReferenceElements< ctype, mydimension > ReferenceElements;


  public:
    /** \brief constructor
     *
     *  \param[in]  refElement  reference element for the geometry
     *  \param[in]  vertices    interpolatory vertices to store internally
     *  \param[in]  order       interpolatory order of the geometry
     *
     *  \note The type of corners is actually a template argument.
     *        It is only required that the internal corner storage can be
     *        constructed from this object.
     */
    CurvilinearGeometry ( const ReferenceElement &refElement,
                          const std::vector<GlobalCoordinate> &vertices,
                          InterpolatoryOrderType order)
    {
        elementInterpolator_ = ElementInterpolator( refElement, vertices, order);
    }

    /** \brief constructor
     *
     *  \param[in]  gt          geometry type
     *  \param[in]  vertices    interpolatory vertices to store internally
     *  \param[in]  order       interpolatory order of the geometry
     *
     *  \note The type of corners is actually a template argument.
     *        It is only required that the internal corner storage can be
     *        constructed from this object.
     */
    CurvilinearGeometry ( Dune::GeometryType gt,
                          const std::vector<GlobalCoordinate> &vertices,
                          InterpolatoryOrderType order)
    {
        elementInterpolator_ = ElementInterpolator( gt, vertices, order);
    }

    /** \brief constructor
     *
     *  \param[in]  elemInterp  element interpolator for this element
     *
     *  \note Construct a geometry from an existing interpolator
     */
    CurvilinearGeometry ( const ElementInterpolator & elemInterp) : elementInterpolator_(elemInterp)  {  }



    /** \brief is this mapping affine? */
    bool affine () const
    {
        return false;
    }

    const ElementInterpolator & interpolator() const { return elementInterpolator_; }

    InterpolatoryOrderType order() const { return elementInterpolator_.order(); }

    /** \brief obtain the name of the reference element */
    Dune::GeometryType type () const { return elementInterpolator_.type(); }

    /** \brief obtain the polynomial vector mapping from local to global geometry of the element */
    PolynomialVector interpolatoryVectorAnalytical() const { return elementInterpolator_.interpolatoryVectorAnalytical(); }

    /** \brief obtain an interpolatory vertex of the geometry */
    GlobalCoordinate vertex (int i) const { return elementInterpolator_.vertex(i); }

    /** \brief obtain all interpolatory vertices of the geometry */
    std::vector<GlobalCoordinate> vertexSet() const { return elementInterpolator_.vertexSet(); }

    /** \brief obtain number of vertices associated with the interpolatory polynomial */
    int nVertex () const { return elementInterpolator_.dofPerOrder(); }

    /** \brief obtain number of corners of the corresponding reference element */
    int nCorner () const { return refElement().size( mydimension ); }

    /** \brief obtain coordinates of the i-th corner */
    GlobalCoordinate corner ( InternalIndexType cornerLinearIndex ) const
    {
      assert( (cornerLinearIndex >= 0) && (cornerLinearIndex < nCorner()) );
      return elementInterpolator_.corner(cornerLinearIndex);
    }

    /** \brief obtain a vector of coordinates of all corners */
    std::vector< GlobalCoordinate > cornerSet() const
    {
        std::vector< GlobalCoordinate > rez;
        for (InternalIndexType i = 0; i < nCorner(); i++) { rez.push_back(corner(i)); }
        return rez;
    }

    /** \brief obtain the centroid of the mapping's image */
    GlobalCoordinate center () const { return global( refElement().position( 0, 0 ) ); }

    /** \brief evaluate the mapping
     *
     *  \param[in]  local  local coordinate to map
     *
     *  \returns corresponding global coordinate
     */
    GlobalCoordinate global ( const LocalCoordinate &local ) const
    {
        //std::cout << " requested global from local " << local << " of dim " << local.size() << " for global size " << coorddimension << "at interpvert size " << vertexSet().size()  << std::endl;
        return elementInterpolator_.realCoordinate(local);
    }

    /** \brief Construct CurvilinearGeometry classes for all mydim-1 subentities of this element
     *
     *  \returns a vector of CurvilinearGeometry classes corresponding to mydim-1 subentity geometries
     */
    template<int subdim>
    CurvilinearGeometry< ctype, subdim, cdim>  subentityGeometry(InternalIndexType subentityIndex) const
    {
        CurvilinearElementInterpolator <ct, subdim, cdim> interp = elementInterpolator_.template SubentityInterpolator<subdim>(subentityIndex);
        return CurvilinearGeometry< ctype, subdim, cdim> (interp);
    }

    /** \brief Construct a global coordinate normal of the curvilinear element evaluated at a given local coordinate
     *
     *    \param[in] local    local coordinate to at which the normal is evaluated
     *  \returns            a global vector which is the normal to the surface
     *
     *  \note: This method will throw an error if normal is not defined for this element type/dimension
     *  \note: Sign of the normal is determined entirely by the order of the corners
     *
     */
    // [TODO] Specialize only for mydim = cdim - 1
    GlobalCoordinate normal(const LocalCoordinate &local ) const
    {
        PolynomialVector analyticalMap = elementInterpolator_.interpolatoryVectorAnalytical();

        if (!type().isSimplex()) { DUNE_THROW(Dune::IOError, "__ERROR: normal() method only implemented for simplex geometries at the moment"); }

             if ((mydim == 1) && (cdim == 2)) { return normalEdge(local, analyticalMap); }
        else if ((mydim == 2) && (cdim == 3)) { return normalTriangle(local, analyticalMap); }
        else
        {
            DUNE_THROW(Dune::IOError, "__ERROR: normal() method only defined for edges in 2D and triangles in 3D");
        }
    }

    /** \brief Constructs the normal to the desired subentity of the element
     *
     *  \param[in]  indexInInside   the order number of the subentity of interest
     *  \param[in]  local           local coordinate within the subentity, at which to evaluate the normal
     *
     *  \return vector normal to the subentity, length is arbitrary
     *
     *  \note This should in principle work for any 2D and 3D objects
     *
     *  TODO: Try understand if normal to 2D curved face in 3D is meaningful
     *
     */
    // [TODO] Specialize only for cdim==mydim
    GlobalCoordinate subentityNormal(InternalIndexType indexInInside, const LocalCoordinate &local ) const
    {
        PolynomialVector analyticalMap = elementInterpolator_.interpolatoryVectorAnalytical();
        return subentityNormal(indexInInside, local, false, false, analyticalMap);
    }

    /** \brief Same as subentityNormal, but normalizes the vector
     *
     *  \param[in]  indexInInside   the order number of the subentity of interest
     *  \param[in]  local           local coordinate within the subentity, at which to evaluate the normal
     *
     *  \return unit vector normal to the subentity
     *
     */
    // [TODO] Specialize only for cdim==mydim
    GlobalCoordinate subentityUnitNormal(InternalIndexType indexInInside, const LocalCoordinate &local ) const
    {
        PolynomialVector analyticalMap = elementInterpolator_.interpolatoryVectorAnalytical();
        return subentityNormal(indexInInside, local, true, false, analyticalMap);
    }

    /** \brief Same as subentityNormal, but normalizes the vector
     *
     *  \param[in]  indexInInside   the order number of the subentity of interest
     *  \param[in]  local           local coordinate within the subentity, at which to evaluate the normal
     *
     *  \return unit vector normal to the subentity
     *
     */
    // [TODO] Specialize only for cdim==mydim
    GlobalCoordinate subentityIntegrationNormal(InternalIndexType indexInInside, const LocalCoordinate &local ) const
    {
        PolynomialVector analyticalMap = elementInterpolator_.interpolatoryVectorAnalytical();
        return subentityNormal(indexInInside, local, false, true, analyticalMap);
    }


    /** \brief evaluate the inverse mapping
     *
     *  \param[in]  global  global coordinate to map
     *  \param[in]  local   local coordinate returned
     *
     *  \return whether the local coordinate is inside the element
     *
     *  \note IMPORTANT! For curvilinear geometries local coordinate outside the element is not defined.
     *  The provided algorithm is not designed to find a suitable local coordinate outside the element.
     *  Thus, if the return value of the function is false, the value of local has no physical meaning.
     *
     *  \note IMPORTANT! It is very difficult (and unreliable) to find if a global point is inside of the element, if the
     *  dimensions of the element and the world do not match. Therefore, the local method is also only available for matching
     *  element and world dimensions.
     *
     *  \note IMPORTANT! Is_Inside() can not be checked by using global-to-local mapping in curvilinear case,
     *  because the for many global coordinates outside the element a local coordinate is simply not defined.
     *  For example, a local-to-global map x^2 is valid for an edge [0,1]. However, there is no local coordinate corresponding to negative global coordinates
     *
     *  \note For given global coordinate y the returned local coordinate x that minimizes
     *  the following function over the local coordinate space spanned by the reference element.
     *  \code
     *  (global( x ) - y).two_norm()
     *  \endcode
     */
    bool local ( const GlobalCoordinate &globalC, LocalCoordinate & localC ) const
    {
        //std::cout << " requested local from global " << globalC << " of dim " << globalC.size() << std::endl;
        if (mydim == 0)  { return true; }

        if (!type().isSimplex()) { DUNE_THROW(Dune::IOError, "__ERROR: curvilinear local() method only available for Simplex geometries at the moment :("); }

        if (mydimension != coorddimension) {
            std::cout << "ERROR: CURVGEOMETRY: local() requested for mydim=" << mydim << " cdim=" << cdim << std::endl;
            DUNE_THROW(Dune::IOError, "__ERROR: local() method is only defined if dimension of the element and world match :(");
        }

        PolynomialVector analyticalMap = elementInterpolator_.interpolatoryVectorAnalytical();
        return local(globalC, localC, analyticalMap);
    }

    /** \brief obtain the integration element
     *
     *  If the Jacobian of the mapping is denoted by $J(x)$, the integration
     *  integration element \f$\mu(x)\f$ is given by
     *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
     *
     *  \param[in]  local  local coordinate to evaluate the integration element in
     *
     *  \returns the integration element \f$\mu(x)\f$.
     *
     *  \note For affine mappings, it is more efficient to call
     *        jacobianInverseTransposed before integrationElement, if both
     *        are required.
     */
    ctype integrationElement ( const LocalCoordinate &local ) const
    {
        PolynomialVector analyticalMap = elementInterpolator_.interpolatoryVectorAnalytical();
        return integrationElement ( local, analyticalMap );
    }


    /** \brief Constructs polynomial Jacobian determinant.
     * \note (!) Only to be used for mydim == codim
     * \note It is an important observation that det(J) is not allowed to change sign within the element,
     * because that is equivalent to having self-overlapping geometry, and must have been avoided
     * during the interpolation stage. Thus to calculate |det(J)|, it is enough to find the sign of
     * det(J) inside the element and multiply by it.
     */
    LocalPolynomial JacobianDeterminantAnalytical() const
    {
        PolynomialVector analyticalMap = interpolatoryVectorAnalytical();
        //return JacobianDeterminantAnalytical(analyticalMap);
    	return JacDetAnalytical::eval(*this, analyticalMap);
    }


    /** \brief Constructs polynomial vector = integration element * element normal
     * \note Since element interpolation is polynomial, this quantity is also polynomial and thus is given analytically
     * \note (!) Only works for elements which have normals - edges in 2D and faces in 3D. Should NOT be called for any other combination
     */
    PolynomialVector NormalIntegrationElementAnalytical() const
    {
        PolynomialVector analyticalMap = interpolatoryVectorAnalytical();
        //return NormalIntegrationElementAnalytical(analyticalMap);
        return IntElemNormalAnalytical::eval(*this, analyticalMap);
    }


    /** \brief Calculates analytically det|JJ^T|, without taking the square root. Re-uses NormalIntegrationElementAnalytical */
    LocalPolynomial IntegrationElementSquaredAnalytical() const
    {
        PolynomialVector analyticalMap = interpolatoryVectorAnalytical();
        //return IntegrationElementSquaredAnalytical(analyticalMap);
    	return IntElemSquaredAnalytical::eval(*this, analyticalMap);
    }


    /** \brief obtain the volume of the mapping's image
     *
     *  \note The current implementation just returns
     *  \code
     *  integrationElement( refElement().position( 0, 0 ) ) * refElement().volume()
     *  \endcode
     *  which is wrong for n-linear surface maps and other nonlinear maps.
     */
    ctype volume (ctype RELATIVE_TOLERANCE) const
    {
        if (mydim == 0)  { return 0; }

        // As long as ACCURACY_GOAL is small, its actual value is irrelevant, because it only matters if integral is close to 0, but element volumes should not be 0
        ctype ACCURACY_GOAL = Traits::tolerance();

        return IntegrationHelper::integrateScalar(*this, identityPolynomial<ctype, mydim>(), RELATIVE_TOLERANCE, ACCURACY_GOAL);
    }

    /** \brief obtain the transposed of the Jacobian
     *
     *  \param[in]  local  local coordinate to evaluate Jacobian in
     *
     *  \returns a reference to the transposed of the Jacobian
     *
     *  \note The returned reference is reused on the next call to
     *        JacobianTransposed, destroying the previous value.
     */
    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
    {
        PolynomialVector analyticalMap = elementInterpolator_.interpolatoryVectorAnalytical();
        return jacobianTransposed (local, analyticalMap) ;
    }

    /** \brief obtain the transposed of the Jacobian's inverse
     *
     *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
     *  the Jacobian by \f$J(x)\f$, the following condition holds:
     *  \f[J^{-1}(x) J(x) = I.\f]
     */
    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const;


    /** \brief Returns reference element  */
    const ReferenceElement &refElement () const { return elementInterpolator_.refElement(); }

  protected:

    /** \brief Computes JT by differentiating the analytical map  */
    // [TODO][Optimization] Compare performance with numerical derivatives. For that, need to hard-code numerical derivatives into interpolator. For now too much effort
    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local, const PolynomialVector & analyticalMap ) const
    {
        JacobianTransposed jt;

        for (int i = 0; i < coorddimension; i++)
        {
            for (int j = 0; j < mydimension; j++)
            {
                jt[j][i] = (analyticalMap[i].derivative(j)).evaluate(local);
            }
        }

      return jt;
    }


    /** \brief Forwards-declaration of Jit constructing function */
    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local, const PolynomialVector & analyticalMap ) const;


    /** \brief Calculate 1D unit normal by rotating the tangential vector by pi/2 clockwise */
    // [TODO] Rename to intrinsic direction.
    // [TODO] Implement normal sign-correction by providing this routine with the 3rd point.
    GlobalCoordinate normalEdge(const LocalCoordinate &local, const PolynomialVector & analyticalMap ) const
    {
        GlobalCoordinate rez;
        rez[0] =  (analyticalMap[1].derivative(0)).evaluate(local);
        rez[1] = -(analyticalMap[0].derivative(0)).evaluate(local);

        // Sometimes we are unlucky and both derivatives evaluate to 0 at the selected point
        // Need to evaluate derivatives numerically
        if (rez.two_norm() <= 1.0e-10)
        {
            LocalCoordinate tmp_local = local;

            // This effect can only happen on the corners, so we need to sample 2nd point for the derivative
            // inside the edge
            tmp_local[0] += 0.002 * (0.5 - local[0]);

            rez[0] =       analyticalMap[1].evaluate(tmp_local) - analyticalMap[1].evaluate(local);
            rez[1] = -1 * (analyticalMap[0].evaluate(tmp_local) - analyticalMap[0].evaluate(local));

            if (local[0] > 0.5) { rez *= -1; }
        }


        // Normalize the normal
        rez /= rez.two_norm();

        return rez;
    }


    /** \brief Calculate 2D unit normal by computing the negative left cross product of tangential surface vectors
     * Advantage:    Does not use Jit, so is not sensitive to det(J) = 0
     * Disadvantage: Currently does not implement correct normal direction */
    // [TODO] Rename to intrinsic direction.
    // [TODO] Implement normal sign-correction by providing this routine with the 4th point.
    // [TODO] Specialize only for the 2D case
    GlobalCoordinate normalTriangle(const LocalCoordinate &local, const PolynomialVector & analyticalMap ) const
    {
        GlobalCoordinate rez;
        rez[0] =  -(analyticalMap[1].derivative(0) * analyticalMap[2].derivative(1) - analyticalMap[1].derivative(1) * analyticalMap[2].derivative(0) ).evaluate(local);
        rez[1] =  -(analyticalMap[2].derivative(0) * analyticalMap[0].derivative(1) - analyticalMap[2].derivative(1) * analyticalMap[0].derivative(0) ).evaluate(local);
        rez[2] =  -(analyticalMap[0].derivative(0) * analyticalMap[1].derivative(1) - analyticalMap[0].derivative(1) * analyticalMap[1].derivative(0) ).evaluate(local);

        // Normalize the normal
        rez /= rez.two_norm();

        return rez;
    }


    /** \brief Finds the normal of a (codim=1) subentity of a volume. Defined for cdim == mydim
     * Constructs normal in local coordinates using referenceElement, then maps it to global using inverse Jacobi transform  */
    // [TODO] Specialize only for the volume case
    // [FIXME] Defend against the case of det(jit) == 0
    GlobalCoordinate subentityNormal(InternalIndexType indexInInside, const LocalCoordinate &local, bool is_normalized, bool is_integrationelement, const PolynomialVector & analyticalMap ) const
    {
        JacobianInverseTransposed jit = jacobianInverseTransposed(local, analyticalMap);
        LocalCoordinate refNormal = refElement().integrationOuterNormal(indexInInside);

        GlobalCoordinate normal;
        jit.mv( refNormal, normal );

        if (is_normalized)               { normal *= 1.0 / normal.two_norm(); }
        else if (is_integrationelement)  { normal *= jit.detInv(); }


        /*
        CurvilinearGeometry< ctype, mydim-1, cdim>  triGeom = subentityGeometry<mydim-1>(indexInInside);
        GlobalCoordinate normal = triGeom.normalTriangle(local, triGeom.interpolatoryVectorAnalytical());
        if (is_normalized)               { normal *= 1.0 / normal.two_norm(); }
        else if (is_integrationelement)  {
            normal *= 1.0 / normal.two_norm();
            normal *= triGeom.integrationElement(local); }

        */

        return normal;
    }


    /** \brief Implements global->local mapping using amortized Newton's algorithm.
     * Returns false if the global coordinate is not inside the element. In this case the resulting local coordinate is not defined,
     * because Lagrange Polynomials are only bijective inside of the element. */
    // TODO: refElement.checkInside returns false if point very close to boundary, but outside, which is undesirable behaviour
    bool local ( const GlobalCoordinate &globalC, LocalCoordinate & localC, const PolynomialVector & analyticalMap ) const
    {
        const ctype tolerance = Traits::tolerance();

        // Check if the point is too far away from the element to be inside it
        if (!FarPointTest::isInside(*this, globalC, tolerance)) { return false; }
        else
        {
            // May calculate barycentric coordinate to have a chance of guessing if the point is inside
            // However it is not necessary, because regardless of its result we have to run the Newton Method
            // bool barycentric_test_rez = isInsideTestBarycentric(globalC, tolerance);


            LocalCoordinate c = refElement().position( 0, 0 );
            LocalCoordinate x = c;
            LocalCoordinate dx;

            // If the algorithm has low convergence, it must be stuck in some weird geometry,
            // which can only happen outside the element, because inside geometry is nice
            bool low_convergence = false;
            ctype running_error = (global( x ) - globalC).two_norm();
            int iterNo = 0;

            // If the local point is very far outside of the element, it is likely that the point is not inside
            bool far_point_local = false;

            do
            {
                iterNo++;

                // Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
                const GlobalCoordinate dglobal = global( x ) - globalC;

                // Calculate convergence rate and local distance every 10 steps
                if ( iterNo % 10 == 0 )
                {
                    ctype new_error = dglobal.two_norm();
                    if ( 2.0 * new_error > running_error ) { low_convergence = true; }
                    running_error = new_error;

                    if ((x - c).two_norm() > 4) { far_point_local = true; }
                }

                MatrixHelper::template xTRightInvA< mydimension, coorddimension >( jacobianTransposed( x, analyticalMap), dglobal, dx );
                x -= dx;
            } while ((dx.two_norm2() > tolerance)&&(!low_convergence)&&(!far_point_local));

            // Return the value of the local coordinate found
            localC = x;

            // Point is inside if convergence was reached and if the local point is inside refElement
            return ((dx.two_norm2() <= tolerance) && refElement().checkInside(x));
        }
    }


    /** \brief Implements generalized integration element I = sqrt(det(J^T J))  */
    ctype integrationElement ( const LocalCoordinate &local, const PolynomialVector & analyticalMap ) const
    {
      //assert(mydim > 0);
      return MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( jacobianTransposed( local, analyticalMap ) );
    }


  protected:
    ElementInterpolator elementInterpolator_;
  };



  // CurvilinearGeometry::JacobianInverseTransposed
  // ----------------------------------------------


  /*********************************************************/
  /* Implementation of CurvilinearGeometry                **/
  /* *******************************************************/

  /** \brief Definition of Jacobian Inverse class */
  // [FIXME] Defend against the case of det(jit) == 0
  template< class ct, int mydim, int cdim, class Traits >
  class CurvilinearGeometry< ct, mydim, cdim, Traits >::JacobianInverseTransposed
    : public FieldMatrix< ctype, coorddimension, mydimension >
  {
    typedef FieldMatrix< ctype, coorddimension, mydimension > Base;

  public:
    void setup ( const JacobianTransposed &jt )
    {
      detInv_ = MatrixHelper::template rightInvA< mydimension, coorddimension >( jt, static_cast< Base & >( *this ) );
    }

    void setupDeterminant ( const JacobianTransposed &jt )
    {
      detInv_ = MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( jt );
    }

    ctype det () const { return ctype( 1 ) / detInv_; }
    ctype detInv () const { return detInv_; }

  private:
    ctype detInv_;
  };


  /** \brief Specialized constructor of Jit for non-cached geometry */
  template< class ct, int mydim, int cdim, class Traits >
  inline typename CurvilinearGeometry< ct, mydim, cdim, Traits >::JacobianInverseTransposed
  CurvilinearGeometry< ct, mydim, cdim, Traits >::jacobianInverseTransposed ( const LocalCoordinate &local ) const
  {
      PolynomialVector analyticalMap = elementInterpolator_.interpolatoryVectorAnalytical();
      return jacobianInverseTransposed(local, analyticalMap);
  }


  /** \brief Generic constructor of Jit for curvilinear geometry, that allows specifying pre-computed analytical map */
  // [TODO] Unnecessary when we merge cached and non-cached geometries
  template< class ct, int mydim, int cdim, class Traits >
  inline typename CurvilinearGeometry< ct, mydim, cdim, Traits >::JacobianInverseTransposed
  CurvilinearGeometry< ct, mydim, cdim, Traits >::jacobianInverseTransposed ( const LocalCoordinate &local, const PolynomialVector & analyticalMap ) const
  {
    JacobianInverseTransposed jit;
    jit.setup( jacobianTransposed( local, analyticalMap ) );
    return jit;
  }





  /** \brief Implement a CurvilinearGeometry with additional caching
   *
   * This class implements the same interface and functionality as CurvilinearGeometry.
   * However, it additionally implements caching for various results.
   *
   *  \tparam  ct      coordinate type
   *  \tparam  mydim   geometry dimension
   *  \tparam  cdim    coordinate dimension
   *  \tparam  Traits  traits allowing to tweak some implementation details
   *                   (optional)
   *
   */
  template< class ct, int mydim, int cdim, class Traits = CurvilinearGeometryTraits< ct > >
  class CachedCurvilinearGeometry
    : public CurvilinearGeometry< ct, mydim, cdim, Traits >
  {
    typedef CachedCurvilinearGeometry< ct, mydim, cdim, Traits > This;
    typedef CurvilinearGeometry< ct, mydim, cdim, Traits > Base;

  public:

    typedef CurvilinearElementInterpolator <ct, mydim, cdim> ElementInterpolator;

    typedef typename Dune::CurvilinearGeometryHelper::InternalIndexType         InternalIndexType;
    typedef typename Dune::CurvilinearGeometryHelper::InterpolatoryOrderType    InterpolatoryOrderType;

    typedef ct ctype;                            //! coordinate type
    static const int mydimension= mydim;         //! geometry dimension
    static const int coorddimension = cdim;      //! coordinate dimension

    //! type of reference element
    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;

    typedef typename Base::LocalCoordinate LocalCoordinate;
    typedef typename Base::GlobalCoordinate GlobalCoordinate;

    typedef typename Base::JacobianTransposed JacobianTransposed;
    typedef typename Base::JacobianInverseTransposed JacobianInverseTransposed;

    typedef Polynomial<ctype, mydimension> LocalPolynomial;
    typedef std::vector<LocalPolynomial> PolynomialVector;

    typedef typename Dune::DifferentialHelper::JacobianDeterminantAnalytical<This, cdim, mydim>         JacDetAnalytical;
    typedef typename Dune::DifferentialHelper::NormalIntegrationElementAnalytical<This, cdim, mydim>    IntElemNormalAnalytical;
    typedef typename Dune::DifferentialHelper::IntegrationElementSquaredAnalytical <This, cdim, mydim>  IntElemSquaredAnalytical;

    typedef typename Dune::IntegralHelper<This, Traits::QUADRATURE_NORM_TYPE>  IntegrationHelper;


    template< class Vertices >
    CachedCurvilinearGeometry ( const ReferenceElement &refElement,
            const Vertices &vertices,
            int order)
       : Base(refElement, vertices, order)
    {
        init();
    }


    template< class Vertices >
    CachedCurvilinearGeometry (
            Dune::GeometryType gt,
            const Vertices &vertices,
            int order)
      : Base ( gt, vertices, order)
    {
        init();
    }


    CachedCurvilinearGeometry ( const ElementInterpolator & elemInterp)
      : Base (elemInterp)
    {
        init();
    }


    ~CachedCurvilinearGeometry()
    {
        if (JacobianDet_)  { delete JacobianDet_; };
        if (IntElem2_)     { delete IntElem2_; };
        if (NormIntElem_)  { delete NormIntElem_; };
    }

    /** \brief Construct CachedCurvilinearGeometry classes for all mydim-1 subentities of this element
     *
     *  \returns a vector of CachedCurvilinearGeometry classes corresponding to mydim-1 subentity geometries
     */
    template<int subdim>
    CachedCurvilinearGeometry< ctype, subdim, cdim>  subentityGeometry(InternalIndexType subentityIndex) const
    {
        return CachedCurvilinearGeometry< ctype, subdim, cdim> (elementInterpolator_. template SubentityInterpolator<subdim>(subentityIndex));
    }


    GlobalCoordinate normal(const LocalCoordinate &local ) const
    {
        if (!Base::type().isSimplex()) { DUNE_THROW(Dune::IOError, "__ERROR: normal() method only implemented for simplex geometries at the moment"); }

             if ((mydim == 1) && (cdim == 2)) { return Base::normalEdge(local, analyticalMap_); }
        else if ((mydim == 2) && (cdim == 3)) { return Base::normalTriangle(local, analyticalMap_); }
        else
        {
            DUNE_THROW(Dune::IOError, "__ERROR: normal() method only defined for edges in 2D and triangles in 3D");
        }
    }


    GlobalCoordinate subentityNormal(InternalIndexType indexInInside, const LocalCoordinate &local ) const
    {
        return Base::subentityNormal(indexInInside, local, false, false, analyticalMap_);
    }


    GlobalCoordinate subentityUnitNormal(InternalIndexType indexInInside, const LocalCoordinate &local ) const
    {
        return Base::subentityNormal(indexInInside, local, true, false, analyticalMap_);
    }


    GlobalCoordinate subentityIntegrationNormal(InternalIndexType indexInInside, const LocalCoordinate &local ) const
    {
        return Base::subentityNormal(indexInInside, local, false, true, analyticalMap_);
    }


    bool local ( const GlobalCoordinate &global, LocalCoordinate & local ) const
    {
        if (mydim == 0)  { return true; }
        return Base::local(global, local, analyticalMap_);
    }


    ctype integrationElement ( const LocalCoordinate &local ) const
    {
        return Base::integrationElement ( local, analyticalMap_ );
    }


    /** \brief Compute Jit using pre-computed analytical map */
    ctype volume (ctype RELATIVE_TOLERANCE) const
    {
        if (mydim == 0) { return 0; }

        // As long as ACCURACY_GOAL is small, its actual value is irrelevant, because it only matters if integral is close to 0, but element volumes should not be 0
        ctype ACCURACY_GOAL = Traits::tolerance();

        return IntegrationHelper::integrateScalar(*this, identityPolynomial<ctype, mydim>(), RELATIVE_TOLERANCE, ACCURACY_GOAL);
    }


    /** \brief Compute Jt using pre-computed analytical map */
    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
    {
        return Base::jacobianTransposed(local, analyticalMap_) ;
    }


    /** \brief Compute Jit using pre-computed analytical map */
    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
        return Base::jacobianInverseTransposed(local, analyticalMap_);
    }


    /** \brief Retrieve pre-computed analytical map from local to global coordinates */
    PolynomialVector interpolatoryVectorAnalytical() const       { return analyticalMap_; }


    /** \brief Retrieve pre-computed analytical Jacobian determinant. If it does not exist yet, pre-compute it */
    LocalPolynomial JacobianDeterminantAnalytical() const {
    	if (!JacobianDet_)  {
    		assert((mydim > 0)&&(mydim == cdim));  // Analytical Jacobian determinant only defined for volume-geometries
    		JacobianDet_ = new LocalPolynomial(JacDetAnalytical::eval(*this, analyticalMap_));
    	}

    	return *JacobianDet_;
    }


    /** \brief Retrieve pre-computed analytical normal integration element. If it does not exist yet, pre-compute it */
    PolynomialVector NormalIntegrationElementAnalytical() const  {
    	if (!NormIntElem_)  {
    		bool valid_dim = ((mydimension == 1) && (coorddimension == 2)) || ((mydimension == 2) && (coorddimension == 3));
    		assert(valid_dim);  // Analytical integration element only defined for surface-geometries and line-geometries
    		NormIntElem_ = new PolynomialVector(IntElemNormalAnalytical::eval(*this, analyticalMap_));
    	}

    	return *NormIntElem_;
    }


    /** \brief Retrieve pre-computed analytical generalized integration element squared. If it does not exist yet, pre-compute it */
    LocalPolynomial IntegrationElementSquaredAnalytical() const  {
    	if (!IntElem2_)  {
    		assert((mydim > 0)&&(mydim != cdim));  // Analytical integration element only defined for surface-geometries and line-geometries
    		IntElem2_ = new LocalPolynomial(IntElemSquaredAnalytical::eval(*this, analyticalMap_));
    	}

    	return *IntElem2_;
    }



  private:

    // Initialization runs at the constructor
    void init()  { analyticalMap_ = Base::interpolatoryVectorAnalytical();
    	JacobianDet_ = nullptr;
    	IntElem2_ = nullptr;
    	NormIntElem_ = nullptr;
    }


  protected:
    using Base::elementInterpolator_;


  private:
    mutable PolynomialVector analyticalMap_;
    mutable LocalPolynomial  * JacobianDet_;
    mutable LocalPolynomial  * IntElem2_;
    mutable PolynomialVector * NormIntElem_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_CURVILINEARGEOMETRY_HH
