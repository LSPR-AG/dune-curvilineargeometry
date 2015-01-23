/*******************************************************************
 * Curvilinear Element Interpolator
 * 
 * author: Aleksejs Fomins
 * date: 01.08.2014 - created
 * 
 * description:
 * Provides the lagrange polynomial interpolation of mesh elements
 * - Stores interpolatory vertices of the element, given by its GeometryType and Interpolation Order
 * - Currently only manages Simplex geometries, error is thrown for any other geometry type
 * - Constructs analytic polynomial representing the map from reference to curvilinear element. Also provides direct numerical evaluation of this map for computational speedup
 * - Provides curvilinear interpolators for subentities of this element
 * 
 * 
 * 
 *******************************************************************/



#ifndef DUNE_CURVILINEARELEMENTINTERPOLATOR_HH
#define DUNE_CURVILINEARELEMENTINTERPOLATOR_HH

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <config.h>
#include <stdio.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>


namespace Dune {


template<class ctype, int mydim, int cdim>
class CurvilinearElementInterpolator {
  protected:
    static const int mydimension= mydim;
    static const int coorddimension = cdim;

    typedef typename Dune::CurvilinearGeometryHelper::InternalIndexType        InternalIndexType;
    typedef typename Dune::CurvilinearGeometryHelper::InterpolatoryOrderType   InterpolatoryOrderType;
    typedef typename Dune::CurvilinearGeometryHelper::IntegerCoordinateVector  IntegerCoordinateVector;

    typedef FieldVector<ctype, mydimension> LocalVector;
    typedef FieldVector< ctype, coorddimension > GlobalVector;
    typedef Polynomial<ctype, mydimension> LocalPolynomial;
    typedef std::vector<LocalPolynomial> PolynomialVector;

    typedef Dune::ReferenceElement< ctype, mydimension > ReferenceElement;
    typedef Dune::ReferenceElements< ctype, mydimension > ReferenceElements;

    const ReferenceElement *refElement_;
    InterpolatoryOrderType order_;

    std::vector<GlobalVector> point_;

  public:

    /** \brief Empty constructor. Used mostly as a patch, please use the other two constructors for actual interpolators */
    CurvilinearElementInterpolator() {}

    
    /** \brief Constructor
     *
     *  \param[in] refElement	provides reference element to obtain base functionality of an element and know its type
     *  \param[in] point	a vector of coordinates of interpolatory vertices using Dune convention
     *  \param[in] order	the interpolation order of this element
     *
     */
    CurvilinearElementInterpolator (
    		const ReferenceElement &refElement,
    		const std::vector<GlobalVector> & point,
    		InterpolatoryOrderType order) :
        refElement_( &refElement ),
        point_(point)
  {
    order_ = order;
  }

  
    /** \brief Constructor
     *
     *  \param[in] gt		GeometryType of this element, to know what it is and to construct a reference element.
     *  \param[in] point	a vector of coordinates of interpolatory vertices using Dune convention
     *  \param[in] order	the interpolation order of this element
     *
     */  
    CurvilinearElementInterpolator (
    		Dune::GeometryType gt,
    		const std::vector<GlobalVector> & point,
    		InterpolatoryOrderType order) :
        refElement_( &ReferenceElements::general( gt ) ),
        point_(point)
  {
    order_ = order;
  }


  


    /** \brief Interpolaton order of this element */
    InterpolatoryOrderType order() const { return order_; }

    /** \brief GeometryType of this element */
    Dune::GeometryType type () const { return refElement_->type(); }

    /** \brief Reference element associated to this element */
    const ReferenceElement &refElement () const { return *refElement_; }

    /** \brief Number of degrees of freedom of this element (number of interpolation points) */
    int dofPerOrder() const { return Dune::CurvilinearGeometryHelper::dofPerOrder(type(), order_); }

    /** \brief Coordinate of the interpolatory point_
     *  \param[in]  i	interpolatory vertex internal index
     */
    GlobalVector vertex(InternalIndexType vertexIndex) const { return point_[vertexIndex]; }

    /** \brief Coordinates of all interpolatory points */
    std::vector<GlobalVector> vertexSet() const { return point_; }

    /** \brief Internal vertex index of a corner in the interpolatory vertex vector
     *  \param[in]  cornerLinearIndex	index of a corner wrt set of corners of the entity
     */
    InternalIndexType cornerIndex(InternalIndexType cornerLinearIndex) const
    {
    	return Dune::CurvilinearGeometryHelper::cornerIndex(type(), order_, cornerLinearIndex);
    }

    /** \brief Coordinate of a corner, given corner index
     *  \param[in]  cornerLinearIndex	index of a corner wrt set of corners of the entity
     */
    GlobalVector corner(InternalIndexType cornerLinearIndex) const { return vertex(cornerIndex(cornerLinearIndex)); }

    /** \brief Number of corners of this element */
    int nCorner() const { return refElement_->size(mydim); }

    /** \brief Number of codim-subentities of this element */
    int nSubentity(int codim) const { return refElement_->size(codim); }


    /** \brief  Lagrange LocalPolynomial corresponding to a given vertex, evaluated at local reference coordinates (u,v,w)
     *  \param[in]  vertexIndex		The internal index of the associated interpolatory vertex
     *  \param[in]  local		The local coordinate to evaluate the polynomial
     * 
     *  note: there is one lagrange polynomial corresponding to each interpolatory vertex, therefore it makes sense to have the same index for polynomials and vertices
     */
    double lagrangePolynomial(InternalIndexType vertexIndex, const LocalVector &local) const
    {
        if (!type().isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: lagrangePolynomial() only implemented for Simplex geometries at the moment"); }

        double rez = 0;

        switch (mydim) {
        case 1:  rez = lagrangePolynomialEdge       (vertexIndex, local[0]                     );   break;
        case 2:  rez = lagrangePolynomialTriangle   (vertexIndex, local[0], local[1]           );   break;
        case 3:  rez = lagrangePolynomialTetrahedron(vertexIndex, local[0], local[1], local[2] );   break;
        }

        return rez;
    }

    /** \brief  Numerically evaluates the curvilinear map from local to global coordinates
     *  \param[in]  local		The local coordinate to evaluate the polynomial
     * 
     *  note: the map is given by linear superposition of lagrange polynomials with interpolatory vertex coordinates
     */
    GlobalVector realCoordinate(const LocalVector &local) const
    {
        GlobalVector rez;

        for (int j = 0; j < coorddimension; j++) { rez[j] = 0; }

        for (int i = 0; i < dofPerOrder(); i++)
        {
            rez.axpy(lagrangePolynomial(i, local), point_[i]);

            //double LP = lagrangePolynomial(i, local);
            //for (int j = 0; j < coorddimension; j++) { rez[j] += point_[i][j] * LP; }
        }

        return rez;
    }

    /** \brief  Analytic map from local to global coordinates, given explicitly by the polynomial class  */
    PolynomialVector interpolatoryVectorAnalytical() const {
        if (!type().isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: interpolatoryVectorAnalytical() only implemented for Simplex geometries at the moment"); }

        return interpolatoryVectorAnalyticalSimplex();
    }


    /** \brief  Generate a subentity interpolator, templated by the dimension of the subentity
     *  \param[in]  subentityNo		Subentity internal index inside the element
     */
    template<int subdim>
    CurvilinearElementInterpolator< ctype, subdim, cdim > SubentityInterpolator(InternalIndexType subentityNo) const
    {
        if (!type().isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: SubentityInterpolator() only implemented for Simplex geometries at the moment"); }
        if ((subentityNo < 0)||(subentityNo >= nSubentity(mydim - subdim)))
        {
        	DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: SubentityInterpolator() - Unexpected subentity index");
        }

        std::vector<InternalIndexType> subentityIndex = Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ctype, mydim>(type(), order_, mydim - subdim, subentityNo);
        std::vector<GlobalVector> subentityPoint;
        for (int i = 0; i < subentityIndex.size(); i++)  { subentityPoint.push_back(point_[subentityIndex[i]]); }

        Dune::GeometryType subentityType;
        subentityType.makeSimplex(subdim);
        return CurvilinearElementInterpolator< ctype, subdim, cdim > (subentityType, subentityPoint, order_);
    }

  private:

    // Implementation of the SimplexInterpolator Class
    // *****************************************************

    /** \brief  Lagrange polynomial for edge
     *  \param[in]  vertexIndex		The internal index of the associated interpolatory vertex
     *  \param[in]  u			first local coordinate
     */
    double lagrangePolynomialEdge(InternalIndexType vertexIndex, double u) const {
        double u2; double u3; double u4; double u5;

        switch (order_)
        {
        case 1:
            switch (vertexIndex)
            {
            case 0    :     return -u + 1;   break;
            case 1    :     return u;   break;
            } break;
        case 2:
            u2 = pow(u,2);

            switch (vertexIndex)
            {
            case 0    :     return 2*u2 - 3*u + 1;   break;
            case 1    :     return -4*u2 + 4*u;   break;
            case 2    :     return 2*u2 - u;   break;
            } break;
        case 3:
            u2 = pow(u,2);  u3 = pow(u,3);

              switch (vertexIndex)
               {
               case 0    :     return -9.0/2*u3 + 9*u2 - 11.0/2*u + 1;   break;
               case 1    :     return 27.0/2*u3 - 45.0/2*u2 + 9*u;   break;
               case 2    :     return -27.0/2*u3 + 18*u2 - 9.0/2*u;   break;
               case 3    :     return 9.0/2*u3 - 9.0/2*u2 + u;   break;
               } break;
        case 4:
            u2 = pow(u,2);  u3 = pow(u,3); u4 = pow(u,4);

              switch (vertexIndex)
            {
               case 0    :     return 32.0/3*u4 - 80.0/3*u3 + 70.0/3*u2 - 25.0/3*u + 1;   break;
               case 1    :     return -128.0/3*u4 + 96*u3 - 208.0/3*u2 + 16*u;    break;
               case 2    :     return 64*u4 - 128*u3 + 76*u2 - 12*u;    break;
               case 3    :     return -128.0/3*u4 + 224.0/3*u3 - 112.0/3*u2 + 16.0/3*u;   break;
               case 4    :     return 32.0/3*u4 - 16*u3 + 22.0/3*u2 - u;    break;
               } break;
        case 5:
            u2 = pow(u,2);  u3 = pow(u,3); u4 = pow(u,4); u5 = pow(u,5);

              switch (vertexIndex)
            {
               case 0    :     return -625.0/24*u5 + 625.0/8*u4 - 2125.0/24*u3 + 375.0/8*u2 - 137.0/12*u + 1;   break;
               case 1    :     return 3125.0/24*u5 - 4375.0/12*u4 + 8875.0/24*u3 - 1925.0/12*u2 + 25*u;    break;
               case 2    :     return -3125.0/12*u5 + 8125.0/12*u4 - 7375.0/12*u3 + 2675.0/12*u2 - 25*u;    break;
               case 3    :     return 3125.0/12*u5 - 625*u4 + 6125.0/12*u3 - 325.0/2*u2 + 50.0/3*u;    break;
               case 4    :     return -3125.0/24*u5 + 6875.0/24*u4 - 5125.0/24*u3 + 1525.0/24*u2 - 25.0/4*u;   break;
               case 5    :     return 625.0/24*u5 - 625.0/12*u4 + 875.0/24*u3 - 125.0/12*u2 + u;   break;
               } break;
        }
    }


    /** \brief  Lagrange polynomial for triangle
     *  \param[in]  vertexIndex		The internal index of the associated interpolatory vertex
     *  \param[in]  u			first local coordinate
     *  \param[in]  v			second local coordinate
     */
    double lagrangePolynomialTriangle(InternalIndexType vertexIndex, double u, double v) const {
        double u2; double u3; double u4; double u5;
        double v2; double v3; double v4; double v5;


      switch(order_)
      {
      case 1 :
          switch (vertexIndex)
          {
        case 0    :    return -u - v + 1;    break;
        case 1    :    return u;        break;
        case 2    :    return v;        break;
          }    break;
      case 2:
          u2 = pow(u,2);  v2 = pow(v,2);

          switch (vertexIndex)
          {
        case 0    :    return 2*u2 + 4*u*v + 2*v2 - 3*u - 3*v + 1;    break;
        case 1    :    return -4*u2 - 4*u*v + 4*u;            break;
        case 2    :    return 2*u2 - u;                break;
        case 3    :    return -4*u*v - 4*v2 + 4*v;            break;
        case 4    :    return 4*u*v;                    break;
        case 5    :    return 2*v2 - v;                break;
          }    break;
      case 3:
          u2 = pow(u,2);  v2 = pow(v,2);
          u3 = pow(u,3);  v3 = pow(v,3);

          switch (vertexIndex)
          {
        case 0    :    return -9.0/2*u3 - 27.0/2*u2*v - 27.0/2*u*v2 - 9.0/2*v3 + 9*u2 + 18*u*v + 9*v2 - 11.0/2*u - 11.0/2*v + 1;    break;
        case 1    :    return 27.0/2*u3 + 27*u2*v + 27.0/2*u*v2 - 45.0/2*u2 - 45.0/2*u*v + 9*u;                    break;
        case 2    :    return -27.0/2*u3 - 27.0/2*u2*v + 18*u2 + 9.0/2*u*v - 9.0/2*u;                            break;
        case 3    :    return 9.0/2*u3 - 9.0/2*u2 + u;                                            break;
        case 4    :    return 27.0/2*u2*v + 27*u*v2 + 27.0/2*v3 - 45.0/2*u*v - 45.0/2*v2 + 9*v;                    break;
        case 5    :    return -27*u2*v - 27*u*v2 + 27*u*v;                                        break;
        case 6    :    return 27.0/2*u2*v - 9.0/2*u*v;                                            break;
        case 7    :    return -27.0/2*u*v2 - 27.0/2*v3 + 9.0/2*u*v + 18*v2 - 9.0/2*v;                            break;
        case 8    :    return 27.0/2*u*v2 - 9.0/2*u*v;                                            break;
        case 9    :    return 9.0/2*v3 - 9.0/2*v2 + v;                                            break;
          }    break;
      case 4:
            u2 = pow(u,2);  v2 = pow(v,2);
            u3 = pow(u,3);  v3 = pow(v,3);
            u4 = pow(u,4);  v4 = pow(v,4);

          switch (vertexIndex)
          {
        case 0    :    return 32.0/3*u4 + 128.0/3*u3*v + 64*u2*v2 + 128.0/3*u*v3 + 32.0/3*v4 - 80.0/3*u3 - 80*u2*v - 80*u*v2 - 80.0/3*v3 + 70.0/3*u2 + 140.0/3*u*v + 70.0/3*v2 - 25.0/3*u - 25.0/3*v + 1;    break;
        case 1    :    return -128.0/3*u4 - 128*u3*v - 128*u2*v2 - 128.0/3*u*v3 + 96*u3 + 192*u2*v + 96*u*v2 - 208.0/3*u2 - 208.0/3*u*v + 16*u;                                break;
        case 2    :    return 64*u4 + 128*u3*v + 64*u2*v2 - 128*u3 - 144*u2*v - 16*u*v2 + 76*u2 + 28*u*v - 12*u;                                                break;
        case 3    :    return -128.0/3*u4 - 128.0/3*u3*v + 224.0/3*u3 + 32*u2*v - 112.0/3*u2 - 16.0/3*u*v + 16.0/3*u;                                                break;
        case 4    :    return 32.0/3*u4 - 16*u3 + 22.0/3*u2 - u;                                                                        break;
        case 5    :    return -128.0/3*u3*v - 128*u2*v2 - 128*u*v3 - 128.0/3*v4 + 96*u2*v + 192*u*v2 + 96*v3 - 208.0/3*u*v - 208.0/3*v2 + 16*v;                                break;
        case 6    :    return 128*u3*v + 256*u2*v2 + 128*u*v3 - 224*u2*v - 224*u*v2 + 96*u*v;                                                            break;
        case 7    :    return -128*u3*v - 128*u2*v2 + 160*u2*v + 32*u*v2 - 32*u*v;                                                                break;
        case 8    :    return 128.0/3*u3*v - 32*u2*v + 16.0/3*u*v;                                                                        break;
        case 9    :    return 64*u2*v2 + 128*u*v3 + 64*v4 - 16*u2*v - 144*u*v2 - 128*v3 + 28*u*v + 76*v2 - 12*v;                                                break;
        case 10   :    return -128*u2*v2 - 128*u*v3 + 32*u2*v + 160*u*v2 - 32*u*v;                                                                break;
        case 11   :    return 64*u2*v2 - 16*u2*v - 16*u*v2 + 4*u*v;                                                                        break;
        case 12   :    return -128.0/3*u*v3 - 128.0/3*v4 + 32*u*v2 + 224.0/3*v3 - 16.0/3*u*v - 112.0/3*v2 + 16.0/3*v;                                                break;
        case 13   :    return 128.0/3*u*v3 - 32*u*v2 + 16.0/3*u*v;                                                                        break;
        case 14   :    return 32.0/3*v4 - 16*v3 + 22.0/3*v2 - v;                                                                        break;
          }    break;
      case 5:
          u2 = pow(u,2);  v2 = pow(v,2);
          u3 = pow(u,3);  v3 = pow(v,3);
          u4 = pow(u,4);  v4 = pow(v,4);
          u5 = pow(u,5);  v5 = pow(v,5);

          switch (vertexIndex)
          {
        case 0    :    return -625.0/24*u5 - 3125.0/24*u4*v - 3125.0/12*u3*v2 - 3125.0/12*u2*v3 - 3125.0/24*u*v4 - 625.0/24*v5 + 625.0/8*u4 + 625.0/2*u3*v + 1875.0/4*u2*v2 + 625.0/2*u*v3 + 625.0/8*v4 - 2125.0/24*u3 - 2125.0/8*u2*v - 2125.0/8*u*v2 - 2125.0/24*v3 + 375.0/8*u2 + 375.0/4*u*v + 375.0/8*v2 - 137.0/12*u - 137.0/12*v + 1;    break;
        case 1    :    return 3125.0/24*u5 + 3125.0/6*u4*v + 3125.0/4*u3*v2 + 3125.0/6*u2*v3 + 3125.0/24*u*v4 - 4375.0/12*u4 - 4375.0/4*u3*v - 4375.0/4*u2*v2 - 4375.0/12*u*v3 + 8875.0/24*u3 + 8875.0/12*u2*v + 8875.0/24*u*v2 - 1925.0/12*u2 - 1925.0/12*u*v + 25*u;                                        break;
        case 2    :    return -3125.0/12*u5 - 3125.0/4*u4*v - 3125.0/4*u3*v2 - 3125.0/12*u2*v3 + 8125.0/12*u4 + 5625.0/4*u3*v + 3125.0/4*u2*v2 + 625.0/12*u*v3 - 7375.0/12*u3 - 8875.0/12*u2*v - 125*u*v2 + 2675.0/12*u2 + 1175.0/12*u*v - 25*u;                                                break;
        case 3    :    return 3125.0/12*u5 + 3125.0/6*u4*v + 3125.0/12*u3*v2 - 625*u4 - 3125.0/4*u3*v - 625.0/4*u2*v2 + 6125.0/12*u3 + 3875.0/12*u2*v + 125.0/6*u*v2 - 325.0/2*u2 - 75.0/2*u*v + 50.0/3*u;                                                                    break;
        case 4    :    return -3125.0/24*u5 - 3125.0/24*u4*v + 6875.0/24*u4 + 625.0/4*u3*v - 5125.0/24*u3 - 1375.0/24*u2*v + 1525.0/24*u2 + 25.0/4*u*v - 25.0/4*u;                                                                                        break;
        case 5    :    return 625.0/24*u5 - 625.0/12*u4 + 875.0/24*u3 - 125.0/12*u2 + u;                                                                                                                            break;
        case 6    :    return 3125.0/24*u4*v + 3125.0/6*u3*v2 + 3125.0/4*u2*v3 + 3125.0/6*u*v4 + 3125.0/24*v5 - 4375.0/12*u3*v - 4375.0/4*u2*v2 - 4375.0/4*u*v3 - 4375.0/12*v4 + 8875.0/24*u2*v + 8875.0/12*u*v2 + 8875.0/24*v3 - 1925.0/12*u*v - 1925.0/12*v2 + 25*v;                                        break;
        case 7    :    return -3125.0/6*u4*v - 3125.0/2*u3*v2 - 3125.0/2*u2*v3 - 3125.0/6*u*v4 + 1250*u3*v + 2500*u2*v2 + 1250*u*v3 - 5875.0/6*u2*v - 5875.0/6*u*v2 + 250*u*v;                                                                                    break;
        case 8    :    return 3125.0/4*u4*v + 3125.0/2*u3*v2 + 3125.0/4*u2*v3 - 3125.0/2*u3*v - 6875.0/4*u2*v2 - 625.0/4*u*v3 + 3625.0/4*u2*v + 1125.0/4*u*v2 - 125*u*v;                                                                                    break;
        case 9    :    return -3125.0/6*u4*v - 3125.0/6*u3*v2 + 2500.0/3*u3*v + 625.0/2*u2*v2 - 2125.0/6*u2*v - 125.0/3*u*v2 + 125.0/3*u*v;                                                                                                    break;
        case 10   :    return 3125.0/24*u4*v - 625.0/4*u3*v + 1375.0/24*u2*v - 25.0/4*u*v;                                                                                                                            break;
        case 11   :    return -3125.0/12*u3*v2 - 3125.0/4*u2*v3 - 3125.0/4*u*v4 - 3125.0/12*v5 + 625.0/12*u3*v + 3125.0/4*u2*v2 + 5625.0/4*u*v3 + 8125.0/12*v4 - 125*u2*v - 8875.0/12*u*v2 - 7375.0/12*v3 + 1175.0/12*u*v + 2675.0/12*v2 - 25*v;                                                break;
        case 12   :    return 3125.0/4*u3*v2 + 3125.0/2*u2*v3 + 3125.0/4*u*v4 - 625.0/4*u3*v - 6875.0/4*u2*v2 - 3125.0/2*u*v3 + 1125.0/4*u2*v + 3625.0/4*u*v2 - 125*u*v;                                                                                    break;
        case 13   :    return -3125.0/4*u3*v2 - 3125.0/4*u2*v3 + 625.0/4*u3*v + 4375.0/4*u2*v2 + 625.0/4*u*v3 - 375.0/2*u2*v - 375.0/2*u*v2 + 125.0/4*u*v;                                                                                            break;
        case 14   :    return 3125.0/12*u3*v2 - 625.0/12*u3*v - 625.0/4*u2*v2 + 125.0/4*u2*v + 125.0/6*u*v2 - 25.0/6*u*v;                                                                                                            break;
        case 15   :    return 3125.0/12*u2*v3 + 3125.0/6*u*v4 + 3125.0/12*v5 - 625.0/4*u2*v2 - 3125.0/4*u*v3 - 625*v4 + 125.0/6*u2*v + 3875.0/12*u*v2 + 6125.0/12*v3 - 75.0/2*u*v - 325.0/2*v2 + 50.0/3*v;                                                                    break;
        case 16   :    return -3125.0/6*u2*v3 - 3125.0/6*u*v4 + 625.0/2*u2*v2 + 2500.0/3*u*v3 - 125.0/3*u2*v - 2125.0/6*u*v2 + 125.0/3*u*v;                                                                                                    break;
        case 17   :    return 3125.0/12*u2*v3 - 625.0/4*u2*v2 - 625.0/12*u*v3 + 125.0/6*u2*v + 125.0/4*u*v2 - 25.0/6*u*v;                                                                                                            break;
        case 18   :    return -3125.0/24*u*v4 - 3125.0/24*v5 + 625.0/4*u*v3 + 6875.0/24*v4 - 1375.0/24*u*v2 - 5125.0/24*v3 + 25.0/4*u*v + 1525.0/24*v2 - 25.0/4*v;                                                                                        break;
        case 19   :    return 3125.0/24*u*v4 - 625.0/4*u*v3 + 1375.0/24*u*v2 - 25.0/4*u*v;                                                                                                                            break;
        case 20   :    return 625.0/24*v5 - 625.0/12*v4 + 875.0/24*v3 - 125.0/12*v2 + v;                                                                                                                            break;
          }    break;
      }
    }


    /** \brief  Lagrange polynomial for tetrahedron
     *  \param[in]  vertexIndex		The internal index of the associated interpolatory vertex
     *  \param[in]  u			first local coordinate
     *  \param[in]  v			second local coordinate
     *  \param[in]  w			third local coordinate
     */
    double lagrangePolynomialTetrahedron(InternalIndexType vertexIndex, double u, double v, double w) const {
      double u2; double u3; double u4; double u5;
      double v2; double v3; double v4; double v5;
      double w2; double w3; double w4; double w5;


    switch(order_)
    {
    case 1 :
        switch (vertexIndex)
        {
      case 0    :    return -u - v - w + 1;    break;
      case 1    :    return u;        break;
      case 2    :    return v;        break;
      case 3    :    return w;        break;
        }    break;
    case 2:
        u2 = pow(u,2);  v2 = pow(v,2);  w2 = pow(w,2);

        switch (vertexIndex)
        {
      case 0    :    return 2*u2 + 4*u*v + 2*v2 + 4*u*w + 4*v*w + 2*w2 - 3*u - 3*v - 3*w + 1;    break;
      case 1    :    return -4*u2 - 4*u*v - 4*u*w + 4*u;                        break;
      case 2    :    return 2*u2 - u;                                break;
      case 3    :    return -4*u*v - 4*v2 - 4*v*w + 4*v;                        break;
      case 4    :    return 4*u*v;                                    break;
      case 5    :    return 2*v2 - v;                                break;
      case 6    :    return -4*u*w - 4*v*w - 4*w2 + 4*w;                        break;
      case 7    :    return 4*u*w;                                    break;
      case 8    :    return 4*v*w;                                    break;
      case 9    :    return 2*w2 - w;                                break;
        }    break;
    case 3:
        u2 = pow(u,2);  v2 = pow(v,2);  w2 = pow(w,2);
        u3 = pow(u,3);  v3 = pow(v,3);  w3 = pow(w,3);

        switch (vertexIndex)
        {
      case 0    :    return -9.0/2*u3 - 27.0/2*u2*v - 27.0/2*u*v2 - 9.0/2*v3 - 27.0/2*u2*w - 27*u*v*w - 27.0/2*v2*w - 27.0/2*u*w2 - 27.0/2*v*w2 - 9.0/2*w3 + 9*u2 + 18*u*v + 9*v2 + 18*u*w + 18*v*w + 9*w2 - 11.0/2*u - 11.0/2*v - 11.0/2*w + 1;    break;
      case 1    :    return 27.0/2*u3 + 27*u2*v + 27.0/2*u*v2 + 27*u2*w + 27*u*v*w + 27.0/2*u*w2 - 45.0/2*u2 - 45.0/2*u*v - 45.0/2*u*w + 9*u;                                                    break;
      case 2    :    return -27.0/2*u3 - 27.0/2*u2*v - 27.0/2*u2*w + 18*u2 + 9.0/2*u*v + 9.0/2*u*w - 9.0/2*u;                                                                    break;
      case 3    :    return 9.0/2*u3 - 9.0/2*u2 + u;                                                                                                    break;
      case 4    :    return 27.0/2*u2*v + 27*u*v2 + 27.0/2*v3 + 27*u*v*w + 27*v2*w + 27.0/2*v*w2 - 45.0/2*u*v - 45.0/2*v2 - 45.0/2*v*w + 9*v;                                                    break;
      case 5    :    return -27*u2*v - 27*u*v2 - 27*u*v*w + 27*u*v;                                                                                            break;
      case 6    :    return 27.0/2*u2*v - 9.0/2*u*v;                                                                                                    break;
      case 7    :    return -27.0/2*u*v2 - 27.0/2*v3 - 27.0/2*v2*w + 9.0/2*u*v + 18*v2 + 9.0/2*v*w - 9.0/2*v;                                                                    break;
      case 8    :    return 27.0/2*u*v2 - 9.0/2*u*v;                                                                                                    break;
      case 9    :    return 9.0/2*v3 - 9.0/2*v2 + v;                                                                                                    break;
      case 10    :    return 27.0/2*u2*w + 27*u*v*w + 27.0/2*v2*w + 27*u*w2 + 27*v*w2 + 27.0/2*w3 - 45.0/2*u*w - 45.0/2*v*w - 45.0/2*w2 + 9*w;                                                    break;
      case 11    :    return -27*u2*w - 27*u*v*w - 27*u*w2 + 27*u*w;                                                                                            break;
      case 12    :    return 27.0/2*u2*w - 9.0/2*u*w;                                                                                                    break;
      case 13    :    return -27*u*v*w - 27*v2*w - 27*v*w2 + 27*v*w;                                                                                            break;
      case 14    :    return 27*u*v*w;                                                                                                        break;
      case 15    :    return 27.0/2*v2*w - 9.0/2*v*w;                                                                                                    break;
      case 16    :    return -27.0/2*u*w2 - 27.0/2*v*w2 - 27.0/2*w3 + 9.0/2*u*w + 9.0/2*v*w + 18*w2 - 9.0/2*w;                                                                    break;
      case 17    :    return 27.0/2*u*w2 - 9.0/2*u*w;                                                                                                    break;
      case 18    :    return 27.0/2*v*w2 - 9.0/2*v*w;                                                                                                    break;
      case 19    :    return 9.0/2*w3 - 9.0/2*w2 + w;                                                                                                    break;
        }    break;
    case 4:
        u2 = pow(u,2);  v2 = pow(v,2);  w2 = pow(w,2);
        u3 = pow(u,3);  v3 = pow(v,3);  w3 = pow(w,3);
        u4 = pow(u,4);  v4 = pow(v,4);  w4 = pow(w,4);

        switch (vertexIndex)
        {
      case 0    :    return 32.0/3*u4 + 128.0/3*u3*v + 64*u2*v2 + 128.0/3*u*v3 + 32.0/3*v4 + 128.0/3*u3*w + 128*u2*v*w + 128*u*v2*w + 128.0/3*v3*w + 64*u2*w2 + 128*u*v*w2 + 64*v2*w2 + 128.0/3*u*w3 + 128.0/3*v*w3 + 32.0/3*w4 - 80.0/3*u3 - 80*u2*v - 80*u*v2 - 80.0/3*v3 - 80*u2*w - 160*u*v*w - 80*v2*w - 80*u*w2 - 80*v*w2 - 80.0/3*w3 + 70.0/3*u2 + 140.0/3*u*v + 70.0/3*v2 + 140.0/3*u*w + 140.0/3*v*w + 70.0/3*w2 - 25.0/3*u - 25.0/3*v - 25.0/3*w + 1;    break;
      case 1    :    return -128.0/3*u4 - 128*u3*v - 128*u2*v2 - 128.0/3*u*v3 - 128*u3*w - 256*u2*v*w - 128*u*v2*w - 128*u2*w2 - 128*u*v*w2 - 128.0/3*u*w3 + 96*u3 + 192*u2*v + 96*u*v2 + 192*u2*w + 192*u*v*w + 96*u*w2 - 208.0/3*u2 - 208.0/3*u*v - 208.0/3*u*w + 16*u;                                                                                                break;
      case 2    :    return 64*u4 + 128*u3*v + 64*u2*v2 + 128*u3*w + 128*u2*v*w + 64*u2*w2 - 128*u3 - 144*u2*v - 16*u*v2 - 144*u2*w - 32*u*v*w - 16*u*w2 + 76*u2 + 28*u*v + 28*u*w - 12*u;                                                                                                                                        break;
      case 3    :    return -128.0/3*u4 - 128.0/3*u3*v - 128.0/3*u3*w + 224.0/3*u3 + 32*u2*v + 32*u2*w - 112.0/3*u2 - 16.0/3*u*v - 16.0/3*u*w + 16.0/3*u;                                                                                                                                                        break;
      case 4    :    return 32.0/3*u4 - 16*u3 + 22.0/3*u2 - u;                                                                                                                                                                                                    break;
      case 5    :    return -128.0/3*u3*v - 128*u2*v2 - 128*u*v3 - 128.0/3*v4 - 128*u2*v*w - 256*u*v2*w - 128*v3*w - 128*u*v*w2 - 128*v2*w2 - 128.0/3*v*w3 + 96*u2*v + 192*u*v2 + 96*v3 + 192*u*v*w + 192*v2*w + 96*v*w2 - 208.0/3*u*v - 208.0/3*v2 - 208.0/3*v*w + 16*v;                                                                                                break;
      case 6    :    return 128*u3*v + 256*u2*v2 + 128*u*v3 + 256*u2*v*w + 256*u*v2*w + 128*u*v*w2 - 224*u2*v - 224*u*v2 - 224*u*v*w + 96*u*v;                                                                                                                                                            break;
      case 7    :    return -128*u3*v - 128*u2*v2 - 128*u2*v*w + 160*u2*v + 32*u*v2 + 32*u*v*w - 32*u*v;                                                                                                                                                                                break;
      case 8    :    return 128.0/3*u3*v - 32*u2*v + 16.0/3*u*v;                                                                                                                                                                                                    break;
      case 9    :    return 64*u2*v2 + 128*u*v3 + 64*v4 + 128*u*v2*w + 128*v3*w + 64*v2*w2 - 16*u2*v - 144*u*v2 - 128*v3 - 32*u*v*w - 144*v2*w - 16*v*w2 + 28*u*v + 76*v2 + 28*v*w - 12*v;                                                                                                                                        break;
      case 10    :    return -128*u2*v2 - 128*u*v3 - 128*u*v2*w + 32*u2*v + 160*u*v2 + 32*u*v*w - 32*u*v;                                                                                                                                                                                break;
      case 11    :    return 64*u2*v2 - 16*u2*v - 16*u*v2 + 4*u*v;                                                                                                                                                                                                    break;
      case 12    :    return -128.0/3*u*v3 - 128.0/3*v4 - 128.0/3*v3*w + 32*u*v2 + 224.0/3*v3 + 32*v2*w - 16.0/3*u*v - 112.0/3*v2 - 16.0/3*v*w + 16.0/3*v;                                                                                                                                                        break;
      case 13    :    return 128.0/3*u*v3 - 32*u*v2 + 16.0/3*u*v;                                                                                                                                                                                                    break;
      case 14    :    return 32.0/3*v4 - 16*v3 + 22.0/3*v2 - v;                                                                                                                                                                                                    break;
      case 15    :    return -128.0/3*u3*w - 128*u2*v*w - 128*u*v2*w - 128.0/3*v3*w - 128*u2*w2 - 256*u*v*w2 - 128*v2*w2 - 128*u*w3 - 128*v*w3 - 128.0/3*w4 + 96*u2*w + 192*u*v*w + 96*v2*w + 192*u*w2 + 192*v*w2 + 96*w3 - 208.0/3*u*w - 208.0/3*v*w - 208.0/3*w2 + 16*w;                                                                                                break;
      case 16    :    return 128*u3*w + 256*u2*v*w + 128*u*v2*w + 256*u2*w2 + 256*u*v*w2 + 128*u*w3 - 224*u2*w - 224*u*v*w - 224*u*w2 + 96*u*w;                                                                                                                                                            break;
      case 17    :    return -128*u3*w - 128*u2*v*w - 128*u2*w2 + 160*u2*w + 32*u*v*w + 32*u*w2 - 32*u*w;                                                                                                                                                                                break;
      case 18    :    return 128.0/3*u3*w - 32*u2*w + 16.0/3*u*w;                                                                                                                                                                                                    break;
      case 19    :    return 128*u2*v*w + 256*u*v2*w + 128*v3*w + 256*u*v*w2 + 256*v2*w2 + 128*v*w3 - 224*u*v*w - 224*v2*w - 224*v*w2 + 96*v*w;                                                                                                                                                            break;
      case 20    :    return -256*u2*v*w - 256*u*v2*w - 256*u*v*w2 + 256*u*v*w;                                                                                                                                                                                            break;
      case 21    :    return 128*u2*v*w - 32*u*v*w;                                                                                                                                                                                                            break;
      case 22    :    return -128*u*v2*w - 128*v3*w - 128*v2*w2 + 32*u*v*w + 160*v2*w + 32*v*w2 - 32*v*w;                                                                                                                                                                                break;
      case 23 :    return 128*u*v2*w - 32*u*v*w;                                                                                                                                                                                                            break;
      case 24    :    return 128.0/3*v3*w - 32*v2*w + 16.0/3*v*w;                                                                                                                                                                                                    break;
      case 25    :    return 64*u2*w2 + 128*u*v*w2 + 64*v2*w2 + 128*u*w3 + 128*v*w3 + 64*w4 - 16*u2*w - 32*u*v*w - 16*v2*w - 144*u*w2 - 144*v*w2 - 128*w3 + 28*u*w + 28*v*w + 76*w2 - 12*w;                                                                                                                                        break;
      case 26    :    return -128*u2*w2 - 128*u*v*w2 - 128*u*w3 + 32*u2*w + 32*u*v*w + 160*u*w2 - 32*u*w;                                                                                                                                                                                break;
      case 27    :    return 64*u2*w2 - 16*u2*w - 16*u*w2 + 4*u*w;                                                                                                                                                                                                    break;
      case 28    :    return -128*u*v*w2 - 128*v2*w2 - 128*v*w3 + 32*u*v*w + 32*v2*w + 160*v*w2 - 32*v*w;                                                                                                                                                                                break;
      case 29    :    return 128*u*v*w2 - 32*u*v*w;    break;
      case 30    :    return 64*v2*w2 - 16*v2*w - 16*v*w2 + 4*v*w;    break;
      case 31    :    return -128.0/3*u*w3 - 128.0/3*v*w3 - 128.0/3*w4 + 32*u*w2 + 32*v*w2 + 224.0/3*w3 - 16.0/3*u*w - 16.0/3*v*w - 112.0/3*w2 + 16.0/3*w;    break;
      case 32    :    return 128.0/3*u*w3 - 32*u*w2 + 16.0/3*u*w;    break;
      case 33    :    return 128.0/3*v*w3 - 32*v*w2 + 16.0/3*v*w;    break;
      case 34    :    return 32.0/3*w4 - 16*w3 + 22.0/3*w2 - w;    break;
        }    break;
    case 5:
        u2 = pow(u,2);  v2 = pow(v,2);  w2 = pow(w,2);
        u3 = pow(u,3);  v3 = pow(v,3);  w3 = pow(w,3);
        u4 = pow(u,4);  v4 = pow(v,4);  w4 = pow(w,4);
        u5 = pow(u,5);  v5 = pow(v,5);  w5 = pow(w,5);

        switch (vertexIndex)
        {
      case 0    :    return -625.0/24*u5 - 3125.0/24*u4*v - 3125.0/12*u3*v2 - 3125.0/12*u2*v3 - 3125.0/24*u*v4 - 625.0/24*v5 - 3125.0/24*u4*w - 3125.0/6*u3*v*w - 3125.0/4*u2*v2*w - 3125.0/6*u*v3*w - 3125.0/24*v4*w - 3125.0/12*u3*w2 - 3125.0/4*u2*v*w2 - 3125.0/4*u*v2*w2 - 3125.0/12*v3*w2 - 3125.0/12*u2*w3 - 3125.0/6*u*v*w3 - 3125.0/12*v2*w3 - 3125.0/24*u*w4 - 3125.0/24*v*w4 - 625.0/24*w5 + 625.0/8*u4 + 625.0/2*u3*v + 1875.0/4*u2*v2 + 625.0/2*u*v3 + 625.0/8*v4 + 625.0/2*u3*w + 1875.0/2*u2*v*w + 1875.0/2*u*v2*w + 625.0/2*v3*w + 1875.0/4*u2*w2 + 1875.0/2*u*v*w2 + 1875.0/4*v2*w2 + 625.0/2*u*w3 + 625.0/2*v*w3 + 625.0/8*w4 - 2125.0/24*u3 - 2125.0/8*u2*v - 2125.0/8*u*v2 - 2125.0/24*v3 - 2125.0/8*u2*w - 2125.0/4*u*v*w - 2125.0/8*v2*w - 2125.0/8*u*w2 - 2125.0/8*v*w2 - 2125.0/24*w3 + 375.0/8*u2 + 375.0/4*u*v + 375.0/8*v2 + 375.0/4*u*w + 375.0/4*v*w + 375.0/8*w2 - 137.0/12*u - 137.0/12*v - 137.0/12*w + 1;    break;
      case 1    :    return 3125.0/24*u5 + 3125.0/6*u4*v + 3125.0/4*u3*v2 + 3125.0/6*u2*v3 + 3125.0/24*u*v4 + 3125.0/6*u4*w + 3125.0/2*u3*v*w + 3125.0/2*u2*v2*w + 3125.0/6*u*v3*w + 3125.0/4*u3*w2 + 3125.0/2*u2*v*w2 + 3125.0/4*u*v2*w2 + 3125.0/6*u2*w3 + 3125.0/6*u*v*w3 + 3125.0/24*u*w4 - 4375.0/12*u4 - 4375.0/4*u3*v - 4375.0/4*u2*v2 - 4375.0/12*u*v3 - 4375.0/4*u3*w - 4375.0/2*u2*v*w - 4375.0/4*u*v2*w - 4375.0/4*u2*w2 - 4375.0/4*u*v*w2 - 4375.0/12*u*w3 + 8875.0/24*u3 + 8875.0/12*u2*v + 8875.0/24*u*v2 + 8875.0/12*u2*w + 8875.0/12*u*v*w + 8875.0/24*u*w2 - 1925.0/12*u2 - 1925.0/12*u*v - 1925.0/12*u*w + 25*u;    break;
      case 2    :    return -3125.0/12*u5 - 3125.0/4*u4*v - 3125.0/4*u3*v2 - 3125.0/12*u2*v3 - 3125.0/4*u4*w - 3125.0/2*u3*v*w - 3125.0/4*u2*v2*w - 3125.0/4*u3*w2 - 3125.0/4*u2*v*w2 - 3125.0/12*u2*w3 + 8125.0/12*u4 + 5625.0/4*u3*v + 3125.0/4*u2*v2 + 625.0/12*u*v3 + 5625.0/4*u3*w + 3125.0/2*u2*v*w + 625.0/4*u*v2*w + 3125.0/4*u2*w2 + 625.0/4*u*v*w2 + 625.0/12*u*w3 - 7375.0/12*u3 - 8875.0/12*u2*v - 125*u*v2 - 8875.0/12*u2*w - 250*u*v*w - 125*u*w2 + 2675.0/12*u2 + 1175.0/12*u*v + 1175.0/12*u*w - 25*u;    break;
      case 3    :    return 3125.0/12*u5 + 3125.0/6*u4*v + 3125.0/12*u3*v2 + 3125.0/6*u4*w + 3125.0/6*u3*v*w + 3125.0/12*u3*w2 - 625*u4 - 3125.0/4*u3*v - 625.0/4*u2*v2 - 3125.0/4*u3*w - 625.0/2*u2*v*w - 625.0/4*u2*w2 + 6125.0/12*u3 + 3875.0/12*u2*v + 125.0/6*u*v2 + 3875.0/12*u2*w + 125.0/3*u*v*w + 125.0/6*u*w2 - 325.0/2*u2 - 75.0/2*u*v - 75.0/2*u*w + 50.0/3*u;    break;
      case 4    :    return -3125.0/24*u5 - 3125.0/24*u4*v - 3125.0/24*u4*w + 6875.0/24*u4 + 625.0/4*u3*v + 625.0/4*u3*w - 5125.0/24*u3 - 1375.0/24*u2*v - 1375.0/24*u2*w + 1525.0/24*u2 + 25.0/4*u*v + 25.0/4*u*w - 25.0/4*u;    break;
      case 5    :    return 625.0/24*u5 - 625.0/12*u4 + 875.0/24*u3 - 125.0/12*u2 + u;    break;
      case 6    :    return 3125.0/24*u4*v + 3125.0/6*u3*v2 + 3125.0/4*u2*v3 + 3125.0/6*u*v4 + 3125.0/24*v5 + 3125.0/6*u3*v*w + 3125.0/2*u2*v2*w + 3125.0/2*u*v3*w + 3125.0/6*v4*w + 3125.0/4*u2*v*w2 + 3125.0/2*u*v2*w2 + 3125.0/4*v3*w2 + 3125.0/6*u*v*w3 + 3125.0/6*v2*w3 + 3125.0/24*v*w4 - 4375.0/12*u3*v - 4375.0/4*u2*v2 - 4375.0/4*u*v3 - 4375.0/12*v4 - 4375.0/4*u2*v*w - 4375.0/2*u*v2*w - 4375.0/4*v3*w - 4375.0/4*u*v*w2 - 4375.0/4*v2*w2 - 4375.0/12*v*w3 + 8875.0/24*u2*v + 8875.0/12*u*v2 + 8875.0/24*v3 + 8875.0/12*u*v*w + 8875.0/12*v2*w + 8875.0/24*v*w2 - 1925.0/12*u*v - 1925.0/12*v2 - 1925.0/12*v*w + 25*v;    break;
      case 7    :    return -3125.0/6*u4*v - 3125.0/2*u3*v2 - 3125.0/2*u2*v3 - 3125.0/6*u*v4 - 3125.0/2*u3*v*w - 3125*u2*v2*w - 3125.0/2*u*v3*w - 3125.0/2*u2*v*w2 - 3125.0/2*u*v2*w2 - 3125.0/6*u*v*w3 + 1250*u3*v + 2500*u2*v2 + 1250*u*v3 + 2500*u2*v*w + 2500*u*v2*w + 1250*u*v*w2 - 5875.0/6*u2*v - 5875.0/6*u*v2 - 5875.0/6*u*v*w + 250*u*v;    break;
      case 8    :    return 3125.0/4*u4*v + 3125.0/2*u3*v2 + 3125.0/4*u2*v3 + 3125.0/2*u3*v*w + 3125.0/2*u2*v2*w + 3125.0/4*u2*v*w2 - 3125.0/2*u3*v - 6875.0/4*u2*v2 - 625.0/4*u*v3 - 6875.0/4*u2*v*w - 625.0/2*u*v2*w - 625.0/4*u*v*w2 + 3625.0/4*u2*v + 1125.0/4*u*v2 + 1125.0/4*u*v*w - 125*u*v;    break;
      case 9    :    return -3125.0/6*u4*v - 3125.0/6*u3*v2 - 3125.0/6*u3*v*w + 2500.0/3*u3*v + 625.0/2*u2*v2 + 625.0/2*u2*v*w - 2125.0/6*u2*v - 125.0/3*u*v2 - 125.0/3*u*v*w + 125.0/3*u*v;    break;
      case 10    :    return 3125.0/24*u4*v - 625.0/4*u3*v + 1375.0/24*u2*v - 25.0/4*u*v;    break;
      case 11    :    return -3125.0/12*u3*v2 - 3125.0/4*u2*v3 - 3125.0/4*u*v4 - 3125.0/12*v5 - 3125.0/4*u2*v2*w - 3125.0/2*u*v3*w - 3125.0/4*v4*w - 3125.0/4*u*v2*w2 - 3125.0/4*v3*w2 - 3125.0/12*v2*w3 + 625.0/12*u3*v + 3125.0/4*u2*v2 + 5625.0/4*u*v3 + 8125.0/12*v4 + 625.0/4*u2*v*w + 3125.0/2*u*v2*w + 5625.0/4*v3*w + 625.0/4*u*v*w2 + 3125.0/4*v2*w2 + 625.0/12*v*w3 - 125*u2*v - 8875.0/12*u*v2 - 7375.0/12*v3 - 250*u*v*w - 8875.0/12*v2*w - 125*v*w2 + 1175.0/12*u*v + 2675.0/12*v2 + 1175.0/12*v*w - 25*v;    break;
      case 12    :    return 3125.0/4*u3*v2 + 3125.0/2*u2*v3 + 3125.0/4*u*v4 + 3125.0/2*u2*v2*w + 3125.0/2*u*v3*w + 3125.0/4*u*v2*w2 - 625.0/4*u3*v - 6875.0/4*u2*v2 - 3125.0/2*u*v3 - 625.0/2*u2*v*w - 6875.0/4*u*v2*w - 625.0/4*u*v*w2 + 1125.0/4*u2*v + 3625.0/4*u*v2 + 1125.0/4*u*v*w - 125*u*v;    break;
      case 13    :    return -3125.0/4*u3*v2 - 3125.0/4*u2*v3 - 3125.0/4*u2*v2*w + 625.0/4*u3*v + 4375.0/4*u2*v2 + 625.0/4*u*v3 + 625.0/4*u2*v*w + 625.0/4*u*v2*w - 375.0/2*u2*v - 375.0/2*u*v2 - 125.0/4*u*v*w + 125.0/4*u*v;    break;
      case 14    :    return 3125.0/12*u3*v2 - 625.0/12*u3*v - 625.0/4*u2*v2 + 125.0/4*u2*v + 125.0/6*u*v2 - 25.0/6*u*v;    break;
      case 15    :    return 3125.0/12*u2*v3 + 3125.0/6*u*v4 + 3125.0/12*v5 + 3125.0/6*u*v3*w + 3125.0/6*v4*w + 3125.0/12*v3*w2 - 625.0/4*u2*v2 - 3125.0/4*u*v3 - 625*v4 - 625.0/2*u*v2*w - 3125.0/4*v3*w - 625.0/4*v2*w2 + 125.0/6*u2*v + 3875.0/12*u*v2 + 6125.0/12*v3 + 125.0/3*u*v*w + 3875.0/12*v2*w + 125.0/6*v*w2 - 75.0/2*u*v - 325.0/2*v2 - 75.0/2*v*w + 50.0/3*v;    break;
      case 16    :    return -3125.0/6*u2*v3 - 3125.0/6*u*v4 - 3125.0/6*u*v3*w + 625.0/2*u2*v2 + 2500.0/3*u*v3 + 625.0/2*u*v2*w - 125.0/3*u2*v - 2125.0/6*u*v2 - 125.0/3*u*v*w + 125.0/3*u*v;    break;
      case 17    :    return 3125.0/12*u2*v3 - 625.0/4*u2*v2 - 625.0/12*u*v3 + 125.0/6*u2*v + 125.0/4*u*v2 - 25.0/6*u*v;    break;
      case 18    :    return -3125.0/24*u*v4 - 3125.0/24*v5 - 3125.0/24*v4*w + 625.0/4*u*v3 + 6875.0/24*v4 + 625.0/4*v3*w - 1375.0/24*u*v2 - 5125.0/24*v3 - 1375.0/24*v2*w + 25.0/4*u*v + 1525.0/24*v2 + 25.0/4*v*w - 25.0/4*v;    break;
      case 19    :    return 3125.0/24*u*v4 - 625.0/4*u*v3 + 1375.0/24*u*v2 - 25.0/4*u*v;    break;
      case 20    :    return 625.0/24*v5 - 625.0/12*v4 + 875.0/24*v3 - 125.0/12*v2 + v;    break;
      case 21    :    return 3125.0/24*u4*w + 3125.0/6*u3*v*w + 3125.0/4*u2*v2*w + 3125.0/6*u*v3*w + 3125.0/24*v4*w + 3125.0/6*u3*w2 + 3125.0/2*u2*v*w2 + 3125.0/2*u*v2*w2 + 3125.0/6*v3*w2 + 3125.0/4*u2*w3 + 3125.0/2*u*v*w3 + 3125.0/4*v2*w3 + 3125.0/6*u*w4 + 3125.0/6*v*w4 + 3125.0/24*w5 - 4375.0/12*u3*w - 4375.0/4*u2*v*w - 4375.0/4*u*v2*w - 4375.0/12*v3*w - 4375.0/4*u2*w2 - 4375.0/2*u*v*w2 - 4375.0/4*v2*w2 - 4375.0/4*u*w3 - 4375.0/4*v*w3 - 4375.0/12*w4 + 8875.0/24*u2*w + 8875.0/12*u*v*w + 8875.0/24*v2*w + 8875.0/12*u*w2 + 8875.0/12*v*w2 + 8875.0/24*w3 - 1925.0/12*u*w - 1925.0/12*v*w - 1925.0/12*w2 + 25*w;    break;
      case 22    :    return -3125.0/6*u4*w - 3125.0/2*u3*v*w - 3125.0/2*u2*v2*w - 3125.0/6*u*v3*w - 3125.0/2*u3*w2 - 3125*u2*v*w2 - 3125.0/2*u*v2*w2 - 3125.0/2*u2*w3 - 3125.0/2*u*v*w3 - 3125.0/6*u*w4 + 1250*u3*w + 2500*u2*v*w + 1250*u*v2*w + 2500*u2*w2 + 2500*u*v*w2 + 1250*u*w3 - 5875.0/6*u2*w - 5875.0/6*u*v*w - 5875.0/6*u*w2 + 250*u*w;    break;
      case 23    :    return 3125.0/4*u4*w + 3125.0/2*u3*v*w + 3125.0/4*u2*v2*w + 3125.0/2*u3*w2 + 3125.0/2*u2*v*w2 + 3125.0/4*u2*w3 - 3125.0/2*u3*w - 6875.0/4*u2*v*w - 625.0/4*u*v2*w - 6875.0/4*u2*w2 - 625.0/2*u*v*w2 - 625.0/4*u*w3 + 3625.0/4*u2*w + 1125.0/4*u*v*w + 1125.0/4*u*w2 - 125*u*w;    break;
      case 24    :    return -3125.0/6*u4*w - 3125.0/6*u3*v*w - 3125.0/6*u3*w2 + 2500.0/3*u3*w + 625.0/2*u2*v*w + 625.0/2*u2*w2 - 2125.0/6*u2*w - 125.0/3*u*v*w - 125.0/3*u*w2 + 125.0/3*u*w;    break;
      case 25    :    return 3125.0/24*u4*w - 625.0/4*u3*w + 1375.0/24*u2*w - 25.0/4*u*w;    break;
      case 26    :    return -3125.0/6*u3*v*w - 3125.0/2*u2*v2*w - 3125.0/2*u*v3*w - 3125.0/6*v4*w - 3125.0/2*u2*v*w2 - 3125*u*v2*w2 - 3125.0/2*v3*w2 - 3125.0/2*u*v*w3 - 3125.0/2*v2*w3 - 3125.0/6*v*w4 + 1250*u2*v*w + 2500*u*v2*w + 1250*v3*w + 2500*u*v*w2 + 2500*v2*w2 + 1250*v*w3 - 5875.0/6*u*v*w - 5875.0/6*v2*w - 5875.0/6*v*w2 + 250*v*w;    break;
      case 27    :    return 3125.0/2*u3*v*w + 3125*u2*v2*w + 3125.0/2*u*v3*w + 3125*u2*v*w2 + 3125*u*v2*w2 + 3125.0/2*u*v*w3 - 5625.0/2*u2*v*w - 5625.0/2*u*v2*w - 5625.0/2*u*v*w2 + 1250*u*v*w;    break;
      case 28    :    return -3125.0/2*u3*v*w - 3125.0/2*u2*v2*w - 3125.0/2*u2*v*w2 + 1875*u2*v*w + 625.0/2*u*v2*w + 625.0/2*u*v*w2 - 625.0/2*u*v*w;    break;
      case 29    :    return 3125.0/6*u3*v*w - 625.0/2*u2*v*w + 125.0/3*u*v*w;    break;
      case 30    :    return 3125.0/4*u2*v2*w + 3125.0/2*u*v3*w + 3125.0/4*v4*w + 3125.0/2*u*v2*w2 + 3125.0/2*v3*w2 + 3125.0/4*v2*w3 - 625.0/4*u2*v*w - 6875.0/4*u*v2*w - 3125.0/2*v3*w - 625.0/2*u*v*w2 - 6875.0/4*v2*w2 - 625.0/4*v*w3 + 1125.0/4*u*v*w + 3625.0/4*v2*w + 1125.0/4*v*w2 - 125*v*w;    break;
      case 31    :    return -3125.0/2*u2*v2*w - 3125.0/2*u*v3*w - 3125.0/2*u*v2*w2 + 625.0/2*u2*v*w + 1875*u*v2*w + 625.0/2*u*v*w2 - 625.0/2*u*v*w;    break;
      case 32    :    return 3125.0/4*u2*v2*w - 625.0/4*u2*v*w - 625.0/4*u*v2*w + 125.0/4*u*v*w;    break;
      case 33    :    return -3125.0/6*u*v3*w - 3125.0/6*v4*w - 3125.0/6*v3*w2 + 625.0/2*u*v2*w + 2500.0/3*v3*w + 625.0/2*v2*w2 - 125.0/3*u*v*w - 2125.0/6*v2*w - 125.0/3*v*w2 + 125.0/3*v*w;    break;
      case 34    :    return 3125.0/6*u*v3*w - 625.0/2*u*v2*w + 125.0/3*u*v*w;    break;
      case 35    :    return 3125.0/24*v4*w - 625.0/4*v3*w + 1375.0/24*v2*w - 25.0/4*v*w;    break;
      case 36    :    return -3125.0/12*u3*w2 - 3125.0/4*u2*v*w2 - 3125.0/4*u*v2*w2 - 3125.0/12*v3*w2 - 3125.0/4*u2*w3 - 3125.0/2*u*v*w3 - 3125.0/4*v2*w3 - 3125.0/4*u*w4 - 3125.0/4*v*w4 - 3125.0/12*w5 + 625.0/12*u3*w + 625.0/4*u2*v*w + 625.0/4*u*v2*w + 625.0/12*v3*w + 3125.0/4*u2*w2 + 3125.0/2*u*v*w2 + 3125.0/4*v2*w2 + 5625.0/4*u*w3 + 5625.0/4*v*w3 + 8125.0/12*w4 - 125*u2*w - 250*u*v*w - 125*v2*w - 8875.0/12*u*w2 - 8875.0/12*v*w2 - 7375.0/12*w3 + 1175.0/12*u*w + 1175.0/12*v*w + 2675.0/12*w2 - 25*w;    break;
      case 37    :    return 3125.0/4*u3*w2 + 3125.0/2*u2*v*w2 + 3125.0/4*u*v2*w2 + 3125.0/2*u2*w3 + 3125.0/2*u*v*w3 + 3125.0/4*u*w4 - 625.0/4*u3*w - 625.0/2*u2*v*w - 625.0/4*u*v2*w - 6875.0/4*u2*w2 - 6875.0/4*u*v*w2 - 3125.0/2*u*w3 + 1125.0/4*u2*w + 1125.0/4*u*v*w + 3625.0/4*u*w2 - 125*u*w;    break;
      case 38    :    return -3125.0/4*u3*w2 - 3125.0/4*u2*v*w2 - 3125.0/4*u2*w3 + 625.0/4*u3*w + 625.0/4*u2*v*w + 4375.0/4*u2*w2 + 625.0/4*u*v*w2 + 625.0/4*u*w3 - 375.0/2*u2*w - 125.0/4*u*v*w - 375.0/2*u*w2 + 125.0/4*u*w;    break;
      case 39    :    return 3125.0/12*u3*w2 - 625.0/12*u3*w - 625.0/4*u2*w2 + 125.0/4*u2*w + 125.0/6*u*w2 - 25.0/6*u*w;    break;
      case 40    :    return 3125.0/4*u2*v*w2 + 3125.0/2*u*v2*w2 + 3125.0/4*v3*w2 + 3125.0/2*u*v*w3 + 3125.0/2*v2*w3 + 3125.0/4*v*w4 - 625.0/4*u2*v*w - 625.0/2*u*v2*w - 625.0/4*v3*w - 6875.0/4*u*v*w2 - 6875.0/4*v2*w2 - 3125.0/2*v*w3 + 1125.0/4*u*v*w + 1125.0/4*v2*w + 3625.0/4*v*w2 - 125*v*w;    break;
      case 41    :    return -3125.0/2*u2*v*w2 - 3125.0/2*u*v2*w2 - 3125.0/2*u*v*w3 + 625.0/2*u2*v*w + 625.0/2*u*v2*w + 1875*u*v*w2 - 625.0/2*u*v*w;    break;
      case 42    :    return 3125.0/4*u2*v*w2 - 625.0/4*u2*v*w - 625.0/4*u*v*w2 + 125.0/4*u*v*w;    break;
      case 43    :    return -3125.0/4*u*v2*w2 - 3125.0/4*v3*w2 - 3125.0/4*v2*w3 + 625.0/4*u*v2*w + 625.0/4*v3*w + 625.0/4*u*v*w2 + 4375.0/4*v2*w2 + 625.0/4*v*w3 - 125.0/4*u*v*w - 375.0/2*v2*w - 375.0/2*v*w2 + 125.0/4*v*w;    break;
      case 44    :    return 3125.0/4*u*v2*w2 - 625.0/4*u*v2*w - 625.0/4*u*v*w2 + 125.0/4*u*v*w;    break;
      case 45    :    return 3125.0/12*v3*w2 - 625.0/12*v3*w - 625.0/4*v2*w2 + 125.0/4*v2*w + 125.0/6*v*w2 - 25.0/6*v*w;    break;
      case 46    :    return 3125.0/12*u2*w3 + 3125.0/6*u*v*w3 + 3125.0/12*v2*w3 + 3125.0/6*u*w4 + 3125.0/6*v*w4 + 3125.0/12*w5 - 625.0/4*u2*w2 - 625.0/2*u*v*w2 - 625.0/4*v2*w2 - 3125.0/4*u*w3 - 3125.0/4*v*w3 - 625*w4 + 125.0/6*u2*w + 125.0/3*u*v*w + 125.0/6*v2*w + 3875.0/12*u*w2 + 3875.0/12*v*w2 + 6125.0/12*w3 - 75.0/2*u*w - 75.0/2*v*w - 325.0/2*w2 + 50.0/3*w;    break;
      case 47    :    return -3125.0/6*u2*w3 - 3125.0/6*u*v*w3 - 3125.0/6*u*w4 + 625.0/2*u2*w2 + 625.0/2*u*v*w2 + 2500.0/3*u*w3 - 125.0/3*u2*w - 125.0/3*u*v*w - 2125.0/6*u*w2 + 125.0/3*u*w;    break;
      case 48    :    return 3125.0/12*u2*w3 - 625.0/4*u2*w2 - 625.0/12*u*w3 + 125.0/6*u2*w + 125.0/4*u*w2 - 25.0/6*u*w;    break;
      case 49    :    return -3125.0/6*u*v*w3 - 3125.0/6*v2*w3 - 3125.0/6*v*w4 + 625.0/2*u*v*w2 + 625.0/2*v2*w2 + 2500.0/3*v*w3 - 125.0/3*u*v*w - 125.0/3*v2*w - 2125.0/6*v*w2 + 125.0/3*v*w;    break;
      case 50    :    return 3125.0/6*u*v*w3 - 625.0/2*u*v*w2 + 125.0/3*u*v*w;    break;
      case 51    :    return 3125.0/12*v2*w3 - 625.0/4*v2*w2 - 625.0/12*v*w3 + 125.0/6*v2*w + 125.0/4*v*w2 - 25.0/6*v*w;    break;
      case 52    :    return -3125.0/24*u*w4 - 3125.0/24*v*w4 - 3125.0/24*w5 + 625.0/4*u*w3 + 625.0/4*v*w3 + 6875.0/24*w4 - 1375.0/24*u*w2 - 1375.0/24*v*w2 - 5125.0/24*w3 + 25.0/4*u*w + 25.0/4*v*w + 1525.0/24*w2 - 25.0/4*w;    break;
      case 53    :    return 3125.0/24*u*w4 - 625.0/4*u*w3 + 1375.0/24*u*w2 - 25.0/4*u*w;    break;
      case 54    :    return 3125.0/24*v*w4 - 625.0/4*v*w3 + 1375.0/24*v*w2 - 25.0/4*v*w;    break;
      case 55    :    return 625.0/24*w5 - 625.0/12*w4 + 875.0/24*w3 - 125.0/12*w2 + w;    break;
        }    break;
    }
  }


    /** \brief  Analytic map from local to global coordinates, given explicitly by the polynomial class. Implementation for simplex  */
    PolynomialVector interpolatoryVectorAnalyticalSimplex() const {
        int dof = dofPerOrder();

        PolynomialVector rez;

        PolynomialVector monomial_basis;
        PolynomialVector lagrange_basis;

        // Step 0. Construct an interpolatory grid over the element
        IntegerCoordinateVector simplexPoints = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(order_);
        std::vector<LocalVector> localCoordinateSet = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<ctype, mydim>(simplexPoints, order_);

        // Step 1. Construct monomial basis and set of local coordinates
        for (int i = 0; i < simplexPoints.size(); i++)
        {
                monomial_basis.push_back(LocalPolynomial(PolynomialTraits::Monomial(1.0, simplexPoints[i])));
        }

        // Step 2. Evaluate monomial basis into matrix
        DynamicMatrix< ctype > lagrange_matrix (dof, dof, 0.0);

        for (int i = 0; i < dof; i++)
        {
            for (int j = 0; j < dof; j++)
            {
                lagrange_matrix[i][j] = monomial_basis[i].evaluate(localCoordinateSet[j]);
            }
        }

        // Step 3. Inverse matrix
        lagrange_matrix.invert();

        // Step 4. Construct Lagrange polynomials - linear combinations of monomials with matrix rows
        for (int i = 0; i < dof; i++)
        {
            LocalPolynomial tmpPoly;

            for (int j = 0; j < dof; j++)
            {
                tmpPoly.axpy(monomial_basis[j], lagrange_matrix[i][j]);
            }
            lagrange_basis.push_back(tmpPoly);
        }

        // Step 5. Construct interpolatoryVector - for each coordinate a linear combination of coordinate of points and lagrange polynomials
        for (int d = 0; d < coorddimension; d++)
        {
            LocalPolynomial tmpPoly;
            for (int i = 0; i < dof; i++)
            {
                tmpPoly.axpy(lagrange_basis[i], point_[i][d]);
            }

            rez.push_back(tmpPoly);
            // Remove all terms that are insignificantly small
            rez[d].cleanUp();
        }

        return rez;
    }


}; // SimplexInterpolator




} // Namespace Dune

#endif /** DUNE_CURVILINEARELEMENTINTERPOLATOR_HH **/
