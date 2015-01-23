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
#include <dune/curvilineargeometry/integration/numericalrecursiveinterpolationintegrator.hh>

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
    /** \brief helper structure containing some matrix routines
     *
     *  This helper allows exchanging the matrix inversion algorithms.
     *  It must provide the following static methods:
     *  \code
     *  template< int m, int n >
     *  static ctype sqrtDetAAT ( const FieldMatrix< ctype, m, n > &A );
     *
     *  template< int m, int n >
     *  static ctype rightInvA ( const FieldMatrix< ctype, m, n > &A,
     *                           FieldMatrix< ctype, n, m > &ret );
     *
     *  template< int m, int n >
     *  static void xTRightInvA ( const FieldMatrix< ctype, m, n > &A,
     *                            const FieldVector< ctype, n > &x,
     *                            FieldVector< ctype, m > &y );
     *  \endcode
     */
    typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ct > > MatrixHelper;

    /** \brief tolerance to numerical algorithms */
    static ct tolerance () { return ct( 16 ) * std::numeric_limits< ct >::epsilon(); }

    /** \brief template specifying the storage for the corners
     *
     *  Internally, the CurvilinearGeometry needs to store the corners of the
     *  geometry.
     *
     *  The corner storage may be chosen depending on geometry dimension and
     *  coordinate dimension. It is required to contain a type named Type, e.g.,
     *  \code
     *  template< int mydim, int cdim >
     *  struct CornerStorage
     *  {
     *    typedef std::vector< FieldVector< ctype, cdim > > Type;
     *  };
     *  \endcode
     *  By default, a std::vector of FieldVector is used.
     *
     *  Apart from being copy constructable and assignable, the corner storage
     *  must provide a constant input iterator, i.e., it must define a type
     *  const_iterator and a pair of constant begin / end methods.
     *
     *  \tparam  mydim  geometry dimension
     *  \tparam  cdim   coordinate dimension
     */
    template< int mydim, int cdim >
    struct VertexStorage
    {
      typedef std::vector< FieldVector< ct, cdim > > Type;
    };

    /** \brief will there be only one geometry type for a dimension?
     *
     *  If there is only a single geometry type for a certain dimension,
     *  <em>hasSingleGeometryType::v</em> can be set to true.
     *  Supporting only one geometry type might yield a gain in performance.
     *
     *  If <em>hasSingleGeometryType::v</em> is set to true, an additional
     *  parameter <em>topologyId</em> is required.
     *  Here's an example:
     *  \code
     *  static const unsigned int topologyId = SimplexTopology< dim >::type::id;
     *  \endcode
     */
    template< int dim >
    struct hasSingleGeometryType
    {
      static const bool v = true;
      static const unsigned int topologyId = ~0u;
    };
  };


  // Function that calculates the total integrand, by multiplying the passed integrand and the square root of the squared integration element.
  template<class ct, int mydim, typename Functor>
  struct BoundaryFunctor {
      typedef Polynomial<ct, mydim> LocalPolynomial;
      typedef FieldVector< ct, mydim > LocalCoordinate;

      Functor basis_function_;
      LocalPolynomial integration_element_squared_;

      BoundaryFunctor(const Functor & basis_function, const LocalPolynomial & integration_element_squared) :
          basis_function_ (basis_function),
          integration_element_squared_ (integration_element_squared)
      {}

      double operator()(const LocalCoordinate & x) const { return basis_function_(x) * sqrt(integration_element_squared_.evaluate(x)); }
  };

  // Functor that evaluates a polynomial that it carries inside
  template<class ct, int mydim>
  struct PolynomialFunctor
  {
      typedef Polynomial<ct, mydim> LocalPolynomial;
      typedef FieldVector< ct, mydim > LocalCoordinate;
      LocalPolynomial p_;

      PolynomialFunctor(const LocalPolynomial & p) : p_(p) {}

      double operator()(const LocalCoordinate & x) const { return p_.evaluate(x); }
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

    typedef Polynomial<ctype, mydimension> LocalPolynomial;
    typedef std::vector<LocalPolynomial> PolynomialVector;


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
    template< class Vertices>
    CurvilinearGeometry ( const ReferenceElement &refElement,
                          const Vertices &vertices,
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
    template< class Vertices>
    CurvilinearGeometry ( Dune::GeometryType gt,
                          const Vertices &vertices,
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
    //template< class Vertices>
    CurvilinearGeometry ( const ElementInterpolator & elemInterp) : elementInterpolator_(elemInterp)  {  }



    /** \brief is this mapping affine? */
    bool affine () const
    {
        return false;
    }

    ElementInterpolator interpolator() { return elementInterpolator_; }

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
        if (!type().isSimplex()) { DUNE_THROW(Dune::IOError, "__ERROR: curvilinear local() method only available for Simplex geometries at the moment :("); }

        if (mydimension != coorddimension) {
            std::cout << "ERROR_LOCAL" <<std::endl;
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
        return JacobianDeterminantAnalytical(analyticalMap);
    }

    /** \brief Constructs polynomial vector = integration element * element normal
     * \note Since element interpolation is polynomial, this quantity is also polynomial and thus is given analytically
     * \note (!) Only works for elements which have normals - edges in 2D and faces in 3D. Should NOT be called for any other combination
     */
    PolynomialVector NormalIntegrationElementAnalytical() const
    {
        PolynomialVector analyticalMap = interpolatoryVectorAnalytical();
        return NormalIntegrationElementAnalytical(analyticalMap);
    }

    // Calculates analytically det|JJ^T|, without taking the square root. Re-uses NormalIntegrationElementAnalytical
    LocalPolynomial IntegrationElementSquaredAnalytical() const
    {
        PolynomialVector analyticalMap = interpolatoryVectorAnalytical();
        return IntegrationElementSquaredAnalytical(analyticalMap);
    }


    /** \brief Integrate the given polynomial over the element
     *
     *  If the dimension of the element and world match, the analytical integration will be performed.
     *  Otherwise, analytical integration is not possible, and numerical integration will be performed.
     *
     *  \param[in]  P  Polynomial to integrate (for example, a basis function)
     *  \param[in]  tolerance  the acceptable relative error for numerical integration
     *
     *  \returns the result of the integral
     */
    ctype integrateScalar(const LocalPolynomial & P, double tolerance) const
    {
        if (mydimension == coorddimension) { return integrateAnalyticalScalar(P); }
        else                               { return integrateNumerical(PolynomialFunctor<ct, mydim>(P), tolerance); }
    }

    /** \brief Integrates the given functor numerically over the element
     *
     *
     *  \param[in]  f  Functor which maps a local coordinate to a real value
     *  \param[in]  tolerance  the acceptable relative error for numerical integration
     *
     *  \returns the result of the integral
     */
    template <typename Functor>
    ctype integrateNumerical(const Functor & f, double tolerance) const
    {
        LocalPolynomial integrationElementSquared = IntegrationElementSquaredAnalytical();
        return integrateNumerical(f, tolerance, integrationElementSquared);
    }

    /** \brief Integrate given polynomial analytically over the element
     *
     *  \param[in]  P    Polynomial to integrate over the element
     *  \returns an integral of P over the element
     *  \note (!) Analytical integration is only possible when mydim == cdim
     *  \note (!) At the moment only capable of integrating over a Simplex
     *
     *  TODO: Need to throw error if called with mydim != cdim
     */
    ctype integrateAnalyticalScalar(const LocalPolynomial & P) const
    {
        LocalPolynomial JDet = JacobianDeterminantAnalytical();
        return integrateAnalyticalScalar(P, JDet);
    }

    /** \brief Integrate given polynomial vector normal projection analytically over the element
     *
     *  \param[in]  PVec    Polynomial vector to project-integrate over the element
     *  \returns an integral of PVec over the element boundary, projected
     *  \note (!) This operation only has meaning to elements which have normals - edges in 2D and faces in 3D
     *  \note (!) At the moment only capable of integrating over a Simplex
     *
     *  TODO: Need to throw error if invalid mydim-cdim pair
     */
    ctype integrateAnalyticalDot(const PolynomialVector & PVec) const
    {
        // If the dimensionality makes sense for this integral
        bool valid_dim = ((mydimension == 1) && (coorddimension == 2)) || ((mydimension == 2) && (coorddimension == 3));

        PolynomialVector normalElement = NormalIntegrationElementAnalytical();

        return integrateAnalyticalDot(PVec, normalElement);
    }

    /** \brief Integrate given polynomial vector normal cross product analytically over the element
     *
     *  \param[in]  PVec    Polynomial vector to cross-product-integrate over the element
     *  \returns a coordinate - integral of PVec cross-product with normal over the element boundary
     *  \note (!) This operation only has meaning to elements which have normals - edges in 2D and faces in 3D
     *  \note (!) At the moment only capable of integrating over a Simplex
     *
     *  TODO: NOT IMPLEMENTED YET
     */
    GlobalCoordinate integrateAnalyticalTimes(const PolynomialVector & PVec) const
    {
        GlobalCoordinate rez;

        // For 2D must return a scalar
        // For 3D must return a 3D vector
        // Maybe first enquire if this functionality is really necessary, then implement?

        return rez;
    }


    /** \brief obtain the volume of the mapping's image
     *
     *  \note The current implementation just returns
     *  \code
     *  integrationElement( refElement().position( 0, 0 ) ) * refElement().volume()
     *  \endcode
     *  which is wrong for n-linear surface maps and other nonlinear maps.
     */
    ctype volume (double tolerance) const
    {
      return integrateScalar(identityPolynomial<ctype, mydim>(), tolerance);
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

  protected:

    const ReferenceElement &refElement () const { return elementInterpolator_.refElement(); }


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


    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local, const PolynomialVector & analyticalMap ) const;


    // If the global point is too far away from the global centre of the element
    // return false. Then it can not be inside, because elements with surface curvature
    // more than the internal radius are unexpected
    bool isInsideTestFarPoint( const GlobalCoordinate &globalC, ctype tolerance) const
    {
        // For 1D this test is much simpler - just check if globalC is inbetween the two corners of the edge
        if (cdim == 1) { return (corner(0)[0] - globalC[0] < tolerance)&&(globalC[0] - corner(1)[0] < tolerance); }


        GlobalCoordinate c = center();

        // Find distance from center of element to the point
        ctype d1 = (c - globalC).two_norm();

        // Find largest distance from center to corner
        ctype d2 = 0;
        for (int i = 0; i < nCorner(); i++) { d2 = std::max(d2, (c - corner(i)).two_norm()); }

        return (d1 < 2 * d2);
    }

    /** \brief Checks if a global coordinate is inside the element or not (!imperfect!)
     *
     *  The current algorithm for Simplices computes the global barycentric coordinates and checks if they sum up to 1
     *
     *  For linear elements, normalized sum > 1 implies that element is outside and sum = 1 implies that it is inside
     *
     *  For nonlinear nonconvex elements, normalized sum > 1 does not imply anything, because the total barycentric
     *  area may be larger than the area of the element even with the sample point inside.
     *
     */
    bool isInsideTestBarycentric(const GlobalCoordinate &globalC, ctype tolerance) const
    {
        switch (mydim)
        {
            case 1:  return(isInsideTestBarycentricEdge(globalC, tolerance));         break;
            case 2:  return(isInsideTestBarycentricTriangle(globalC, tolerance));     break;
            case 3:  return(isInsideTestBarycentricTetrahedron(globalC, tolerance));  break;
       }
    }


    // Checking if a point is inside an edge by simply verifying it lies between its two corners. Only for edges in 1D
    bool isInsideTestBarycentricEdge( const GlobalCoordinate &globalC, ctype tolerance) const
    {
        return ((corner(0)[0] - globalC[0] <= tolerance) && (globalC[0]- corner(1)[0] <= tolerance));
    }


    // Checking if a point is inside a triangle by calculating global simplex coordinates. Only for triangles in 2D
    bool isInsideTestBarycentricTriangle( const GlobalCoordinate &globalC, ctype tolerance) const
    {
        typedef CurvilinearElementInterpolator< ctype, mydim - 1, cdim > SubentityInterpolator;
        typedef std::vector< SubentityInterpolator > SubentityInterpolatorVector;
        typedef Polynomial<ctype, mydimension - 1> SubentityPolynomial;
        typedef std::vector<SubentityPolynomial> SubentityPolynomialVector;

        SubentityInterpolatorVector edgeInterpolatorSet;
        edgeInterpolatorSet.push_back(elementInterpolator_.Subentity<1>(0));
        edgeInterpolatorSet.push_back(elementInterpolator_.Subentity<1>(1));
        edgeInterpolatorSet.push_back(elementInterpolator_.Subentity<1>(2));


        ctype tri_area = volume(tolerance);
        ctype barycentric_sum = 0;

        //std::cout << "    # total area " << tri_area << std::endl;

        for (int i = 0; i < edgeInterpolatorSet.size(); i++)
        {
            // Obtains analytical vector function describing the global coordinate of the edge
            SubentityPolynomialVector pv = edgeInterpolatorSet[i].interpolatoryVectorAnalytical();

            // Calculates the area of barycentric triangle created by the corners of the curved edge and the point globalC
            // In this calculation the only the curved edge is curved, the other two edges of this barycentric triangle are straight-sided
            // We integrate over the generalized surface area, which is the cross product between the sweeping point and its derivative
            ctype barycentric_area = 0.5 * ((pv[0] - globalC[0]) * pv[1].derivative(0) - (pv[1] - globalC[1]) * pv[0].derivative(0)).integrateRefSimplex();

            //std::cout << "    # barycentric area " << barycentric_area << std::endl;

            barycentric_sum += barycentric_area;
        }
        return (barycentric_sum / tri_area - 1 < tolerance);
    }


    // Checking if a point is inside a triangle by calculating global simplex coordinates. Only for triangles in 3D
    bool isInsideTestBarycentricTetrahedron( const GlobalCoordinate &globalC, ctype tolerance) const
    {
    	typedef CurvilinearElementInterpolator< ctype, mydim - 1, cdim > SubentityInterpolator;
    	typedef std::vector< SubentityInterpolator > SubentityInterpolatorVector;
        typedef Polynomial<ctype, mydimension - 1> SubentityPolynomial;
        typedef std::vector<SubentityPolynomial> SubentityPolynomialVector;

        SubentityInterpolatorVector faceInterpolatorSet;
        faceInterpolatorSet.push_back(elementInterpolator_.Subentity<2>(0));
        faceInterpolatorSet.push_back(elementInterpolator_.Subentity<2>(1));
        faceInterpolatorSet.push_back(elementInterpolator_.Subentity<2>(2));
        faceInterpolatorSet.push_back(elementInterpolator_.Subentity<2>(3));

        ctype tet_vol = volume(tolerance);
        ctype barycentric_sum = 0;

        //std::cout << "    # total volume " << tet_vol << std::endl;

        for (int i = 0; i < faceInterpolatorSet.size(); i++)
        {
            // Obtains analytical vector function describing the global coordinate of the face
            SubentityPolynomialVector pv = faceInterpolatorSet[i].interpolatoryVectorAnalytical();

            //for (int piter = 0; piter < 3; piter++) { pv[piter].print(); }

            // Calculates the volume of barycentric tetrahedron created by the corners of the curved face and the point globalC
            // In this calculation the only the curved face is curved, the other three faces of this barycentric tetrahedron are straight-sided
            // We integrate over the generalized volume, which is the dot-cross product between the sweeping point and its two derivatives
            ctype barycentric_volume = 0;
            barycentric_volume += ((pv[0] - globalC[0]) * (pv[1].derivative(0) * pv[2].derivative(1) - pv[1].derivative(1) * pv[2].derivative(0))).integrateRefSimplex();
            barycentric_volume += ((pv[1] - globalC[1]) * (pv[2].derivative(0) * pv[0].derivative(1) - pv[2].derivative(1) * pv[0].derivative(0))).integrateRefSimplex();
            barycentric_volume += ((pv[2] - globalC[2]) * (pv[0].derivative(0) * pv[1].derivative(1) - pv[0].derivative(1) * pv[1].derivative(0))).integrateRefSimplex();


            //std::cout << "    # barycentric volume " << barycentric_volume << std::endl;

            barycentric_sum += -barycentric_volume / 3.0;
        }

        //std::cout << "    # barycentric sum : " << barycentric_sum / tet_vol - 1 << std::endl;

        return (barycentric_sum / tet_vol - 1 < tolerance);
    }


    // Calculate edge normal by rotating the tangential vector by pi/2 clockwise
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


    // Calculate triangle normal by computing the negative left cross product of the tangential vectors
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


    // Finds the normal in local coordinates using referenceElement, then maps it to global using inverse Jacobi transform
    GlobalCoordinate subentityNormal(InternalIndexType indexInInside, const LocalCoordinate &local, bool is_normalized, bool is_integrationelement, const PolynomialVector & analyticalMap ) const
    {
        JacobianInverseTransposed jit = jacobianInverseTransposed(local, analyticalMap);
        LocalCoordinate refNormal = refElement().integrationOuterNormal(indexInInside);

        GlobalCoordinate normal;
        jit.mv( refNormal, normal );

        if (is_normalized) { normal *= (ctype( 1 ) / normal.two_norm()); }
        // Constructs normal integration element. detInv(x) is exactly integrationElement(x), but faster
        else if (is_integrationelement) { normal *= jit.detInv(); }

        return normal;
    }


    // TODO: refElement.checkInside returns false if point very close to boundary, but outside, which is undesirable behaviour
    bool local ( const GlobalCoordinate &globalC, LocalCoordinate & localC, const PolynomialVector & analyticalMap ) const
    {
        const ctype tolerance = Traits::tolerance();

        // Check if the point is too far away to be inside
        if (!isInsideTestFarPoint(globalC, tolerance)) { return false; }
        else
        {
            // May calculate barycentric coordinate to have a chance of guessing if the point is inside
            // However it is not necessary, because regardless of its result we have to run the Newton Method
            // bool barycentric_test_rez = isInsideTestBarycentric(globalC, tolerance);


            LocalCoordinate c = refElement().position( 0, 0 );
            LocalCoordinate x = c;
            LocalCoordinate dx;

            //std::cout << "       @ Coordinate (" << x << std::endl;

            // If the algorithm has low convergence, it must be stuck in some weird geometry,
            // which can only happen outside the element, because inside geometry is nice
            bool low_convergence = false;
            ctype running_error = (global( x ) - globalC).two_norm();
            int iterNo = 0;

            // If the local point is very far outside of the element, it is likely that the point is not inside
            bool far_point_local = false;

            //std::cout << "searching for " << globalC << std::endl;

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

                //std::cout << "       @ Coordinate (" << x << "), using J^T = " << jacobianTransposed( x, analyticalMap) << std::endl;

                MatrixHelper::template xTRightInvA< mydimension, coorddimension >( jacobianTransposed( x, analyticalMap), dglobal, dx );
                x -= dx;

                //std::cout << "       @ Coordinate (" << x << std::endl;
            } while ((dx.two_norm2() > tolerance)&&(!low_convergence)&&(!far_point_local));


            //double sum_diff = 1;
            //for (int i = 0; i < mydim; i++) { sum_diff -= x[i]; }

            //std::cout << "   ===end tolerance " << dx.two_norm() << " of " << tolerance << " with total sum-1 " << sum_diff << " result " << x << " is " << refElement().checkInside(x) << std::endl;


            // Return the value of the local coordinate found
            localC = x;

            // Point is inside if convergence was reached and if the local point is inside refElement
            return ((dx.two_norm2() <= tolerance) && refElement().checkInside(x));
        }
    }

    ctype integrationElement ( const LocalCoordinate &local, const PolynomialVector & analyticalMap ) const
    {
      return MatrixHelper::template sqrtDetAAT< mydimension, coorddimension >( jacobianTransposed( local, analyticalMap ) );
    }

    LocalPolynomial JacobianDeterminantAnalytical(const PolynomialVector & analyticalMap) const
    {
        LocalPolynomial rez;

        switch(coorddimension)
        {
        case 1:  rez = analyticalMap[0].derivative(0);
            break;
        case 2:  rez = analyticalMap[0].derivative(0) * analyticalMap[1].derivative(1) - analyticalMap[0].derivative(1) * analyticalMap[1].derivative(0);
            break;
        case 3:
            rez.mergeTo( analyticalMap[0].derivative(0) * ( analyticalMap[1].derivative(1) * analyticalMap[2].derivative(2) - analyticalMap[1].derivative(2) * analyticalMap[2].derivative(1) ) );
            rez.mergeTo( analyticalMap[0].derivative(1) * ( analyticalMap[1].derivative(2) * analyticalMap[2].derivative(0) - analyticalMap[1].derivative(0) * analyticalMap[2].derivative(2) ) );
            rez.mergeTo( analyticalMap[0].derivative(2) * ( analyticalMap[1].derivative(0) * analyticalMap[2].derivative(1) - analyticalMap[1].derivative(1) * analyticalMap[2].derivative(0) ) );
            break;
        }
        // Change sign if determinant is negative
        if (rez.evaluate(refElement().position( 0, 0 )) < 0) { rez.multScalar(-1); }

        return rez;
    }

    PolynomialVector NormalIntegrationElementAnalytical(const PolynomialVector & analyticalMap) const
    {
        PolynomialVector rez;

        // Case of edge in 2D
        if ((mydimension == 1) && (coorddimension == 2))
        {
            rez.push_back(analyticalMap[1].derivative(0));
            rez.push_back(analyticalMap[0].derivative(0) * (-1));
        } else
        // Case of face in 3D
        if ((mydimension == 2) && (coorddimension == 3))
        {
            rez.push_back(analyticalMap[2].derivative(0) * analyticalMap[1].derivative(1) - analyticalMap[1].derivative(0) * analyticalMap[2].derivative(1));
            rez.push_back(analyticalMap[0].derivative(0) * analyticalMap[2].derivative(1) - analyticalMap[2].derivative(0) * analyticalMap[0].derivative(1));
            rez.push_back(analyticalMap[1].derivative(0) * analyticalMap[0].derivative(1) - analyticalMap[0].derivative(0) * analyticalMap[1].derivative(1));
        }

        return rez;
    }

    LocalPolynomial IntegrationElementSquaredAnalytical(const PolynomialVector & analyticalMap) const
    {
        LocalPolynomial rez;
        switch (mydimension)
        {
        case 1:
            for (int i = 0; i < coorddimension; i++) {
                LocalPolynomial tmp = analyticalMap[i].derivative(0);
                rez.mergeTo(tmp * tmp);
            }
            break;
        case 2:
            PolynomialVector normalIntegrationElement = NormalIntegrationElementAnalytical(analyticalMap);
            for (int i = 0; i < coorddimension; i++) {
                rez.mergeTo(normalIntegrationElement[i] * normalIntegrationElement[i]);
            }
            break;
        }
        return rez;
    }

    template <typename Functor>
    ctype integrateNumerical(Functor f, double tolerance, const LocalPolynomial & integrationElementSquared) const
    {
        BoundaryFunctor<ct, mydim, Functor> Integrand(f, integrationElementSquared);
        NumericalRecursiveInterpolationIntegrator<ct, mydim> NInt( type() );
        return NInt.integrate( Integrand, tolerance);
    }

    ctype integrateAnalyticalScalar(const LocalPolynomial & P, const LocalPolynomial & jacobianDeterminant) const
    {
        return (P * jacobianDeterminant).integrateRefSimplex();
    }

    ctype integrateAnalyticalDot(const PolynomialVector & PVec, const PolynomialVector & normalIntegrationElement ) const
    {
        // Construct boundary integration element normal polynomial vector
        LocalPolynomial integrand;
        for (int i = 0; i < coorddimension; i++) { integrand.mergeTo(PVec[i] * normalIntegrationElement[i]); }
        return integrand.integrateRefSimplex();
    }

  protected:
    ElementInterpolator elementInterpolator_;
  };



  // CurvilinearGeometry::JacobianInverseTransposed
  // ----------------------------------------------

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


  // Implementation of CurvilinearGeometry
  // -------------------------------------

  template< class ct, int mydim, int cdim, class Traits >
  inline typename CurvilinearGeometry< ct, mydim, cdim, Traits >::JacobianInverseTransposed
  CurvilinearGeometry< ct, mydim, cdim, Traits >::jacobianInverseTransposed ( const LocalCoordinate &local ) const
  {
      PolynomialVector analyticalMap = elementInterpolator_.interpolatoryVectorAnalytical();
      return jacobianInverseTransposed(local, analyticalMap);
  }

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

    typedef CurvilinearElementInterpolator <ct, mydim, cdim> ElementInterpolator;

    //! coordinate type
    typedef ct ctype;
    //! geometry dimension
    static const int mydimension= mydim;
    //! coordinate dimension
    static const int coorddimension = cdim;

    //! type of reference element
    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;

  public:

    typedef typename Base::LocalCoordinate LocalCoordinate;
    typedef typename Base::GlobalCoordinate GlobalCoordinate;

    typedef typename Base::JacobianTransposed JacobianTransposed;
    typedef typename Base::JacobianInverseTransposed JacobianInverseTransposed;

    typedef Polynomial<ctype, mydimension> LocalPolynomial;
    typedef std::vector<LocalPolynomial> PolynomialVector;


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

    //template< class Vertices >
    CachedCurvilinearGeometry ( const ElementInterpolator & elemInterp)
      : Base (elemInterp)
    {
        init();
    }


    /** \brief Construct CachedCurvilinearGeometry classes for all mydim-1 subentities of this element
     *
     *  \returns a vector of CachedCurvilinearGeometry classes corresponding to mydim-1 subentity geometries
     */
    template<int subdim>
    CachedCurvilinearGeometry< ctype, subdim, cdim>  subentityCachedGeometry(InternalIndexType subentityIndex) const
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
        return Base::local(global, local, analyticalMap_);
    }

    ctype integrationElement ( const LocalCoordinate &local ) const
    {
        return Base::integrationElement ( local, analyticalMap_ );
    }

    LocalPolynomial JacobianDeterminantAnalytical() const
    {
        return Base::JacobianDeterminantAnalytical(analyticalMap_);
    }

    PolynomialVector NormalIntegrationElementAnalytical() const
    {
        return Base::NormalIntegrationElementAnalytical(analyticalMap_);
    }

    LocalPolynomial IntegrationElementSquaredAnalytical() const
    {
        return Base::IntegrationElementSquaredAnalytical(analyticalMap_);
    }

    ctype integrateScalar(const LocalPolynomial & P, double tolerance) const
    {
        if (mydimension == coorddimension) { return integrateAnalyticalScalar(P); }
        else                               { return integrateNumerical(PolynomialFunctor<ct, mydim>(P), tolerance); }
    }

    template <typename Functor>
    ctype integrateNumerical(const Functor & f, double tolerance) const
    {
        return Base::integrateNumerical(f, tolerance, IntElem2_);
    }

    ctype integrateAnalyticalScalar(const LocalPolynomial & P) const
    {
        return Base::integrateAnalyticalScalar(P, JacobianDet_);
    }

    // TODO: Throw error if invalid dim-cdim pair
    ctype integrateAnalyticalDot(const PolynomialVector & PVec) const
    {
        // If the dimensionality makes sense for this integral
        bool valid_dim = ((mydim == 1) && (cdim == 2)) || ((mydim == 2) && (cdim == 3));

        return Base::integrateAnalyticalDot(PVec, NormIntElem_);
    }

    // TODO: Implement integrateAnalyticalTimes function for Cached

    ctype volume (double tolerance) const
    {
      return integrateScalar(identityPolynomial<ctype, mydim>(), tolerance);
    }

    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
    {
        return Base::jacobianTransposed(local, analyticalMap_) ;
    }

    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
        return Base::jacobianInverseTransposed(local, analyticalMap_);
    }

  private:

    // Initialization runs at the constructor
    void init()
    {
        analyticalMap_ = Base::interpolatoryVectorAnalytical();

        if (mydim == cdim) { JacobianDet_ = Base::JacobianDeterminantAnalytical(analyticalMap_); }
        else               { IntElem2_ = Base::IntegrationElementSquaredAnalytical(analyticalMap_); }

        bool valid_dim = ((mydimension == 1) && (coorddimension == 2)) || ((mydimension == 2) && (coorddimension == 3));
        if (valid_dim)     { NormIntElem_ = Base::NormalIntegrationElementAnalytical(analyticalMap_); }
    }



  protected:
    using Base::elementInterpolator_;


  private:
    mutable PolynomialVector analyticalMap_;
    mutable LocalPolynomial JacobianDet_;
    mutable LocalPolynomial IntElem2_;
    mutable PolynomialVector NormIntElem_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_CURVILINEARGEOMETRY_HH
