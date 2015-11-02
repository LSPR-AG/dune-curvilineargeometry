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



#ifndef DUNE_LAGRANGE_INTERPOLATOR_HH
#define DUNE_LAGRANGE_INTERPOLATOR_HH

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
class LagrangeInterpolator {
public:

	/***********************************************************************************/
	/*  Type definitions                                                               */
	/***********************************************************************************/

    static const int mydimension= mydim;
    static const int coorddimension = cdim;

    typedef typename Dune::CurvilinearGeometryHelper::InternalIndexType        InternalIndexType;
    typedef typename Dune::CurvilinearGeometryHelper::InterpolatoryOrderType   InterpolatoryOrderType;
    typedef typename Dune::CurvilinearGeometryHelper::IntegerCoordinateVector  IntegerCoordinateVector;


    typedef std::vector<ctype>                                 NumVector1D;
    typedef std::vector<NumVector1D>                           NumVector2D;
    typedef std::vector<NumVector2D>                           NumVector3D;

    typedef std::vector<std::vector<ctype> >                   PowerVector;
    typedef FieldVector<ctype, mydimension>                    LocalCoordinate;
    typedef FieldVector< ctype, coorddimension >               GlobalCoordinate;
    typedef Polynomial<ctype, mydimension>                     LocalPolynomial;
    typedef std::vector<LocalPolynomial>                       PolynomialVector;
    typedef FieldVector<LocalPolynomial, mydimension>          PolynomialLocalCoordinate;
    typedef FieldVector<LocalPolynomial, coorddimension>       PolynomialGlobalCoordinate;
    typedef typename Dune::PolynomialTraits<ctype>::Monomial   Monomial;

    typedef Dune::ReferenceElement< ctype, mydimension >       ReferenceElement;
    typedef Dune::ReferenceElements< ctype, mydimension >      ReferenceElements;

    const ReferenceElement *refElement_;
    InterpolatoryOrderType order_;

    std::vector<GlobalCoordinate> point_;



	/***********************************************************************************/
	/*  Methods                                                                        */
	/***********************************************************************************/

    /** \brief Empty constructor. Used mostly as a patch, please use the other two constructors for actual interpolators */
    LagrangeInterpolator() {}

    
    /** \brief Constructor
     *
     *  \param[in] refElement	provides reference element to obtain base functionality of an element and know its type
     *  \param[in] point	a vector of coordinates of interpolatory vertices using Dune convention
     *  \param[in] order	the interpolation order of this element
     *
     */
    LagrangeInterpolator (
    		const ReferenceElement &refElement,
    		const std::vector<GlobalCoordinate> & point,
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
    LagrangeInterpolator (
    		Dune::GeometryType gt,
    		const std::vector<GlobalCoordinate> & point,
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
    GlobalCoordinate vertex(InternalIndexType vertexIndex) const {
    	assert(vertexIndex < dofPerOrder());
    	return point_[vertexIndex];
    }

    /** \brief Coordinates of all interpolatory points */
    std::vector<GlobalCoordinate> vertexSet() const { return point_; }

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
    GlobalCoordinate corner(InternalIndexType cornerLinearIndex) const {return vertex(cornerIndex(cornerLinearIndex)); }

    /** \brief Number of corners of this element */
    int nCorner() const { return refElement_->size(mydim); }

    /** \brief Number of codim-subentities of this element */
    int nSubentity(int codim) const { return refElement_->size(codim); }


    /** \brief  Lagrange LocalPolynomial corresponding to a given vertex, evaluated at local reference coordinates (U[0][1],U[1][1],U[2][1])
     *  \param[in]  vertexIndex		The internal index of the associated interpolatory vertex
     *  \param[in]  local		The local coordinate to evaluate the polynomial
     *
     *  IMPORTANT: This method is inefficient. When multiple Lagrange Polynomials are needed for a single coordinate, pre-compute power vector.
     * 
     *  note: there is one lagrange polynomial corresponding to each interpolatory vertex, therefore it makes sense to have the same index for polynomials and vertices
     */
    double lagrangePolynomial(InternalIndexType vertexIndex, const LocalCoordinate &local) const
    {
    	PowerVector pw = powerVector(local);
    	return lagrangePolynomial(vertexIndex, pw);
    }

    /** \brief  Numerically evaluates the curvilinear map from local to global coordinates
     *  \param[in]  local		The local coordinate to evaluate the polynomial
     * 
     *  note: the map is given by linear superposition of lagrange polynomials with interpolatory vertex coordinates
     */
    GlobalCoordinate global(const LocalCoordinate &local) const
    {
    	// In case of 0-dim geometry, a 0-dim vector is mapped to a single point defining the geometry
    	if (mydim == 0) { return point_[0]; }

    	// Pre-compute all monomials used in computation to save time
    	NumVector1D M = hierarchicMonomial(local);

        GlobalCoordinate rez(0.0);
        //for (int j = 0; j < coorddimension; j++) { rez[j] = 0; }

        for (int i = 0; i < dofPerOrder(); i++)
        {
            rez.axpy(lagrangePolynomial(i, M), point_[i]);
        }

        return rez;
    }

    /** \brief  Analytic map from local to global coordinates, given explicitly by the polynomial class  */
    PolynomialGlobalCoordinate interpolatoryVectorAnalytical() const {
        if (!type().isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: interpolatoryVectorAnalytical() only implemented for Simplex geometries at the moment"); }

        return interpolatoryVectorAnalyticalSimplex();
    }


    /** \brief  Generate a subentity interpolator, templated by the dimension of the subentity
     *  \param[in]  subentityNo		Subentity internal index inside the element
     */
    template<int subdim>
    LagrangeInterpolator< ctype, subdim, cdim > SubentityInterpolator(InternalIndexType subentityNo) const
    {
        if (!type().isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: SubentityInterpolator() only implemented for Simplex geometries at the moment"); }
        if ((subentityNo < 0)||(subentityNo >= nSubentity(mydim - subdim)))
        {
        	DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: SubentityInterpolator() - Unexpected subentity index");
        }

        // Vertex subentity can not be lower than a vertex itself
        if (mydim == 0)  {
        	assert(subdim == 0);
        	return LagrangeInterpolator< ctype, subdim, cdim > (type(), point_, order_);
        }

        std::vector<InternalIndexType> subentityIndex = Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ctype, mydim>(type(), order_, mydim - subdim, subentityNo);

        std::vector<GlobalCoordinate> subentityPoint;
        for (int i = 0; i < subentityIndex.size(); i++)  { subentityPoint.push_back(point_[subentityIndex[i]]); }

        Dune::GeometryType subentityType;
        subentityType.makeSimplex(subdim);
        return LagrangeInterpolator< ctype, subdim, cdim > (subentityType, subentityPoint, order_);
    }

  private:

    // Implementation of the SimplexInterpolator Class
    // *****************************************************

    // Generates monomials in hierarchically-increasing-order, since this way same monomial vector is applicable to all interpolation orders
    // Hierarchic means that first all monomials of order N appear, and then are followed by order N+1
    // 1D: 1, x, x^2, x^3, ...
    // 2D: 1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3, ...
    // 3D: 1, x, y, z, x^2, xy, xz, y^2, yz, z^2, x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z yz^2, z^3, ...
    NumVector1D hierarchicMonomial(const LocalCoordinate &local) const
	{
    	NumVector1D rez(1, 1.0);

    	switch(mydim)
    	{
    	case 1 :
    	{
    		for (int iOrd = 1; iOrd <= order_; iOrd++)  { rez.push_back(rez[iOrd - 1] * local[0]); }
    	} break;
    	case 2 :
    	{
    		// First generate all monomials in conventional order
    		NumVector2D tmp2d(order_ + 1, NumVector1D(order_ + 1, 1.0));
    		for (int ix = 1; ix <= order_; ix++)  { tmp2d[ix][0] = tmp2d[ix - 1][0] * local[0]; }
    		for (int ix = 0; ix <= order_; ix++)  {
        		for (int iy = 1; iy <= order_ - ix; iy++)  { tmp2d[ix][iy] = tmp2d[ix][iy - 1] * local[1]; }
    		}

    		// Then convert them to hierarchical order
    		for (int iOrd = 1; iOrd <= order_; iOrd++)  {
    			for (int ix = iOrd; ix >= 0; ix--)  {
    				int iy = iOrd - ix;
    				rez.push_back(tmp2d[ix][iy]);
    			}

    		}
    	} break;
    	case 3 :
    	{
    		// First generate all monomials in conventional order
    		NumVector3D tmp3d(order_ + 1, NumVector2D(order_ + 1, NumVector1D(order_ + 1, 1.0)));
    		for (int ix = 1; ix <= order_; ix++)  { tmp3d[ix][0][0] = tmp3d[ix - 1][0][0] * local[0]; }
    		for (int ix = 0; ix <= order_; ix++)  {
    			for (int iy = 1; iy <= order_ - ix; iy++)  { tmp3d[ix][iy][0] = tmp3d[ix][iy - 1][0] * local[1]; }
        		for (int iy = 0; iy <= order_ - ix; iy++)  {
            		for (int iz = 1; iz <= order_ - ix - iy; iz++)  { tmp3d[ix][iy][iz] = tmp3d[ix][iy][iz - 1] * local[2]; }
        		}
    		}

    		// Then convert them to hierarchical order
    		for (int iOrd = 1; iOrd <= order_; iOrd++)  {
    			for (int ix = iOrd; ix >= 0; ix--)  {
        			for (int iy = iOrd - ix; iy >= 0; iy--)  {
        				int iz = iOrd - ix - iy;
        				rez.push_back(tmp3d[ix][iy][iz]);
        			}
    			}
    		}
    	} break;

    	}

    	//std::cout << "order = " << order_ << "cdim = " << cdim << ";   expected = " << dofPerOrder() << " got " << rez.size() << std::endl;
    	assert(rez.size() == dofPerOrder());


    	return rez;
	}


    /** \brief  Lagrange LocalPolynomial corresponding to a given vertex, evaluated at local reference coordinates (U[0][1],U[1][1],U[2][1])
     *  \param[in]  vertexIndex		The internal index of the associated interpolatory vertex
     *  \param[in]  local		The local coordinate to evaluate the polynomial
     *
     *  IMPORTANT: This method is inefficient. When multiple Lagrange Polynomials are needed for a single coordinate, pre-compute power vector.
     *
     *  note: there is one lagrange polynomial corresponding to each interpolatory vertex, therefore it makes sense to have the same index for polynomials and vertices
     */
    double lagrangePolynomial(InternalIndexType vertexIndex, const NumVector1D & M) const
    {
    	assert(mydim > 0);
        if (!type().isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: lagrangePolynomial() only implemented for Simplex geometries at the moment"); }

        double rez = 0;

        switch (mydim) {
        case 1:  rez = lagrangePolynomialEdge       (vertexIndex, M);   break;
        case 2:  rez = lagrangePolynomialTriangle   (vertexIndex, M);   break;
        case 3:  rez = lagrangePolynomialTetrahedron(vertexIndex, M);   break;
        }

        return rez;
    }



    /** \brief  Lagrange polynomial for edge
     *  \param[in]  vertexIndex		The internal index of the associated interpolatory vertex
     *  \param[in]  U[0][1]			first local coordinate
     */
    double lagrangePolynomialEdge(InternalIndexType vertexIndex, const NumVector1D & M) const {
        switch (order_)
        {
        case 1:
            switch (vertexIndex)
            {
            case 0    :     return -M[1] + 1;   break;
            case 1    :     return M[1];   break;
            } break;
        case 2:
            switch (vertexIndex)
            {
            case 0    :     return 2*M[2] - 3*M[1] + 1;   break;
            case 1    :     return -4*M[2] + 4*M[1];   break;
            case 2    :     return 2*M[2] - M[1];   break;
            } break;
        case 3:
              switch (vertexIndex)
               {
               case 0    :     return -9.0/2*M[3] + 9*M[2] - 11.0/2*M[1] + 1;   break;
               case 1    :     return 27.0/2*M[3] - 45.0/2*M[2] + 9*M[1];   break;
               case 2    :     return -27.0/2*M[3] + 18*M[2] - 9.0/2*M[1];   break;
               case 3    :     return 9.0/2*M[3] - 9.0/2*M[2] + M[1];   break;
               } break;
        case 4:
              switch (vertexIndex)
            {
               case 0    :     return 32.0/3*M[4] - 80.0/3*M[3] + 70.0/3*M[2] - 25.0/3*M[1] + 1;   break;
               case 1    :     return -128.0/3*M[4] + 96*M[3] - 208.0/3*M[2] + 16*M[1];    break;
               case 2    :     return 64*M[4] - 128*M[3] + 76*M[2] - 12*M[1];    break;
               case 3    :     return -128.0/3*M[4] + 224.0/3*M[3] - 112.0/3*M[2] + 16.0/3*M[1];   break;
               case 4    :     return 32.0/3*M[4] - 16*M[3] + 22.0/3*M[2] - M[1];    break;
               } break;
        case 5:
              switch (vertexIndex)
            {
               case 0    :     return -625.0/24*M[5] + 625.0/8*M[4] - 2125.0/24*M[3] + 375.0/8*M[2] - 137.0/12*M[1] + 1;   break;
               case 1    :     return 3125.0/24*M[5] - 4375.0/12*M[4] + 8875.0/24*M[3] - 1925.0/12*M[2] + 25*M[1];    break;
               case 2    :     return -3125.0/12*M[5] + 8125.0/12*M[4] - 7375.0/12*M[3] + 2675.0/12*M[2] - 25*M[1];    break;
               case 3    :     return 3125.0/12*M[5] - 625*M[4] + 6125.0/12*M[3] - 325.0/2*M[2] + 50.0/3*M[1];    break;
               case 4    :     return -3125.0/24*M[5] + 6875.0/24*M[4] - 5125.0/24*M[3] + 1525.0/24*M[2] - 25.0/4*M[1];   break;
               case 5    :     return 625.0/24*M[5] - 625.0/12*M[4] + 875.0/24*M[3] - 125.0/12*M[2] + M[1];   break;
               } break;
        }
    }


    /** \brief  Lagrange polynomial for triangle
     *  \param[in]  vertexIndex		The internal index of the associated interpolatory vertex
     *  \param[in]  U[0][1]			first local coordinate
     *  \param[in]  U[1][1]			second local coordinate
     */
    double lagrangePolynomialTriangle(InternalIndexType vertexIndex, const NumVector1D & M) const {
    	switch(order_)
    	{
    	case 1 :
    		switch (vertexIndex)
    		{
    	   		case 0    :    return -M[1] - M[2] + 1;    break;
    	   		case 1    :    return M[1];        break;
    	   		case 2    :    return M[2];        break;
    		}    break;
    	case 2:
    		switch (vertexIndex)
    		{
    	   		case 0    :    return 2*M[3] + 4*M[4] + 2*M[5] - 3*M[1] - 3*M[2] + 1;    break;
    	   		case 1    :    return -4*M[3] - 4*M[4] + 4*M[1];            break;
    	   		case 2    :    return 2*M[3] - M[1];                break;
    	   		case 3    :    return -4*M[4] - 4*M[5] + 4*M[2];            break;
    	   		case 4    :    return 4*M[4];                    break;
    	   		case 5    :    return 2*M[5] - M[2];                break;
    	   	}    break;
    	case 3:
    		switch (vertexIndex)
    		{
    			case 0    :    return -9.0/2*M[6] - 27.0/2*M[7] - 27.0/2*M[8] - 9.0/2*M[9] + 9*M[3] + 18*M[4] + 9*M[5] - 11.0/2*M[1] - 11.0/2*M[2] + 1;    break;
    			case 1    :    return 27.0/2*M[6] + 27*M[7] + 27.0/2*M[8] - 45.0/2*M[3] - 45.0/2*M[4] + 9*M[1];                    break;
    			case 2    :    return -27.0/2*M[6] - 27.0/2*M[7] + 18*M[3] + 9.0/2*M[4] - 9.0/2*M[1];                            break;
    			case 3    :    return 9.0/2*M[6] - 9.0/2*M[3] + M[1];                                            break;
    			case 4    :    return 27.0/2*M[7] + 27*M[8] + 27.0/2*M[9] - 45.0/2*M[4] - 45.0/2*M[5] + 9*M[2];                    break;
    			case 5    :    return -27*M[7] - 27*M[8] + 27*M[4];                                        break;
    			case 6    :    return 27.0/2*M[7] - 9.0/2*M[4];                                            break;
    			case 7    :    return -27.0/2*M[8] - 27.0/2*M[9] + 9.0/2*M[4] + 18*M[5] - 9.0/2*M[2];                            break;
    			case 8    :    return 27.0/2*M[8] - 9.0/2*M[4];                                            break;
    			case 9    :    return 9.0/2*M[9] - 9.0/2*M[5] + M[2];                                            break;
    		}    break;
    	case 4:
    		switch (vertexIndex)
    		{
    		case 0    :    return 32.0/3*M[10] + 128.0/3*M[11] + 64*M[12] + 128.0/3*M[13] + 32.0/3*M[14] - 80.0/3*M[6] - 80*M[7] - 80*M[8] - 80.0/3*M[9] + 70.0/3*M[3] + 140.0/3*M[4] + 70.0/3*M[5] - 25.0/3*M[1] - 25.0/3*M[2] + 1;    break;
    		case 1    :    return -128.0/3*M[10] - 128*M[11] - 128*M[12] - 128.0/3*M[13] + 96*M[6] + 192*M[7] + 96*M[8] - 208.0/3*M[3] - 208.0/3*M[4] + 16*M[1];                                break;
    		case 2    :    return 64*M[10] + 128*M[11] + 64*M[12] - 128*M[6] - 144*M[7] - 16*M[8] + 76*M[3] + 28*M[4] - 12*M[1];                                                break;
    		case 3    :    return -128.0/3*M[10] - 128.0/3*M[11] + 224.0/3*M[6] + 32*M[7] - 112.0/3*M[3] - 16.0/3*M[4] + 16.0/3*M[1];                                                break;
    		case 4    :    return 32.0/3*M[10] - 16*M[6] + 22.0/3*M[3] - M[1];                                                                        break;
    		case 5    :    return -128.0/3*M[11] - 128*M[12] - 128*M[13] - 128.0/3*M[14] + 96*M[7] + 192*M[8] + 96*M[9] - 208.0/3*M[4] - 208.0/3*M[5] + 16*M[2];                                break;
    		case 6    :    return 128*M[11] + 256*M[12] + 128*M[13] - 224*M[7] - 224*M[8] + 96*M[4];                                                            break;
    		case 7    :    return -128*M[11] - 128*M[12] + 160*M[7] + 32*M[8] - 32*M[4];                                                                break;
    		case 8    :    return 128.0/3*M[11] - 32*M[7] + 16.0/3*M[4];                                                                        break;
    		case 9    :    return 64*M[12] + 128*M[13] + 64*M[14] - 16*M[7] - 144*M[8] - 128*M[9] + 28*M[4] + 76*M[5] - 12*M[2];                                                break;
    		case 10   :    return -128*M[12] - 128*M[13] + 32*M[7] + 160*M[8] - 32*M[4];                                                                break;
    		case 11   :    return 64*M[12] - 16*M[7] - 16*M[8] + 4*M[4];                                                                        break;
    		case 12   :    return -128.0/3*M[13] - 128.0/3*M[14] + 32*M[8] + 224.0/3*M[9] - 16.0/3*M[4] - 112.0/3*M[5] + 16.0/3*M[2];                                                break;
    		case 13   :    return 128.0/3*M[13] - 32*M[8] + 16.0/3*M[4];                                                                        break;
    		case 14   :    return 32.0/3*M[14] - 16*M[9] + 22.0/3*M[5] - M[2];                                                                        break;
    		}    break;
    	case 5:
    		switch (vertexIndex)
    		{
    		case 0    :    return -625.0/24*M[15] - 3125.0/24*M[16] - 3125.0/12*M[17] - 3125.0/12*M[18] - 3125.0/24*M[19] - 625.0/24*M[20] + 625.0/8*M[10] + 625.0/2*M[11] + 1875.0/4*M[12] + 625.0/2*M[13] + 625.0/8*M[14] - 2125.0/24*M[6] - 2125.0/8*M[7] - 2125.0/8*M[8] - 2125.0/24*M[9] + 375.0/8*M[3] + 375.0/4*M[4] + 375.0/8*M[5] - 137.0/12*M[1] - 137.0/12*M[2] + 1;    break;
    		case 1    :    return 3125.0/24*M[15] + 3125.0/6*M[16] + 3125.0/4*M[17] + 3125.0/6*M[18] + 3125.0/24*M[19] - 4375.0/12*M[10] - 4375.0/4*M[11] - 4375.0/4*M[12] - 4375.0/12*M[13] + 8875.0/24*M[6] + 8875.0/12*M[7] + 8875.0/24*M[8] - 1925.0/12*M[3] - 1925.0/12*M[4] + 25*M[1];                                        break;
    		case 2    :    return -3125.0/12*M[15] - 3125.0/4*M[16] - 3125.0/4*M[17] - 3125.0/12*M[18] + 8125.0/12*M[10] + 5625.0/4*M[11] + 3125.0/4*M[12] + 625.0/12*M[13] - 7375.0/12*M[6] - 8875.0/12*M[7] - 125*M[8] + 2675.0/12*M[3] + 1175.0/12*M[4] - 25*M[1];                                                break;
    		case 3    :    return 3125.0/12*M[15] + 3125.0/6*M[16] + 3125.0/12*M[17] - 625*M[10] - 3125.0/4*M[11] - 625.0/4*M[12] + 6125.0/12*M[6] + 3875.0/12*M[7] + 125.0/6*M[8] - 325.0/2*M[3] - 75.0/2*M[4] + 50.0/3*M[1];                                                                    break;
    		case 4    :    return -3125.0/24*M[15] - 3125.0/24*M[16] + 6875.0/24*M[10] + 625.0/4*M[11] - 5125.0/24*M[6] - 1375.0/24*M[7] + 1525.0/24*M[3] + 25.0/4*M[4] - 25.0/4*M[1];                                                                                        break;
    		case 5    :    return 625.0/24*M[15] - 625.0/12*M[10] + 875.0/24*M[6] - 125.0/12*M[3] + M[1];                                                                                                                            break;
    		case 6    :    return 3125.0/24*M[16] + 3125.0/6*M[17] + 3125.0/4*M[18] + 3125.0/6*M[19] + 3125.0/24*M[20] - 4375.0/12*M[11] - 4375.0/4*M[12] - 4375.0/4*M[13] - 4375.0/12*M[14] + 8875.0/24*M[7] + 8875.0/12*M[8] + 8875.0/24*M[9] - 1925.0/12*M[4] - 1925.0/12*M[5] + 25*M[2];                                        break;
    		case 7    :    return -3125.0/6*M[16] - 3125.0/2*M[17] - 3125.0/2*M[18] - 3125.0/6*M[19] + 1250*M[11] + 2500*M[12] + 1250*M[13] - 5875.0/6*M[7] - 5875.0/6*M[8] + 250*M[4];                                                                                    break;
    		case 8    :    return 3125.0/4*M[16] + 3125.0/2*M[17] + 3125.0/4*M[18] - 3125.0/2*M[11] - 6875.0/4*M[12] - 625.0/4*M[13] + 3625.0/4*M[7] + 1125.0/4*M[8] - 125*M[4];                                                                                    break;
    		case 9    :    return -3125.0/6*M[16] - 3125.0/6*M[17] + 2500.0/3*M[11] + 625.0/2*M[12] - 2125.0/6*M[7] - 125.0/3*M[8] + 125.0/3*M[4];                                                                                                    break;
    		case 10   :    return 3125.0/24*M[16] - 625.0/4*M[11] + 1375.0/24*M[7] - 25.0/4*M[4];                                                                                                                            break;
    		case 11   :    return -3125.0/12*M[17] - 3125.0/4*M[18] - 3125.0/4*M[19] - 3125.0/12*M[20] + 625.0/12*M[11] + 3125.0/4*M[12] + 5625.0/4*M[13] + 8125.0/12*M[14] - 125*M[7] - 8875.0/12*M[8] - 7375.0/12*M[9] + 1175.0/12*M[4] + 2675.0/12*M[5] - 25*M[2];                                                break;
    		case 12   :    return 3125.0/4*M[17] + 3125.0/2*M[18] + 3125.0/4*M[19] - 625.0/4*M[11] - 6875.0/4*M[12] - 3125.0/2*M[13] + 1125.0/4*M[7] + 3625.0/4*M[8] - 125*M[4];                                                                                    break;
    		case 13   :    return -3125.0/4*M[17] - 3125.0/4*M[18] + 625.0/4*M[11] + 4375.0/4*M[12] + 625.0/4*M[13] - 375.0/2*M[7] - 375.0/2*M[8] + 125.0/4*M[4];                                                                                            break;
    		case 14   :    return 3125.0/12*M[17] - 625.0/12*M[11] - 625.0/4*M[12] + 125.0/4*M[7] + 125.0/6*M[8] - 25.0/6*M[4];                                                                                                            break;
    		case 15   :    return 3125.0/12*M[18] + 3125.0/6*M[19] + 3125.0/12*M[20] - 625.0/4*M[12] - 3125.0/4*M[13] - 625*M[14] + 125.0/6*M[7] + 3875.0/12*M[8] + 6125.0/12*M[9] - 75.0/2*M[4] - 325.0/2*M[5] + 50.0/3*M[2];                                                                    break;
    		case 16   :    return -3125.0/6*M[18] - 3125.0/6*M[19] + 625.0/2*M[12] + 2500.0/3*M[13] - 125.0/3*M[7] - 2125.0/6*M[8] + 125.0/3*M[4];                                                                                                    break;
    		case 17   :    return 3125.0/12*M[18] - 625.0/4*M[12] - 625.0/12*M[13] + 125.0/6*M[7] + 125.0/4*M[8] - 25.0/6*M[4];                                                                                                            break;
    		case 18   :    return -3125.0/24*M[19] - 3125.0/24*M[20] + 625.0/4*M[13] + 6875.0/24*M[14] - 1375.0/24*M[8] - 5125.0/24*M[9] + 25.0/4*M[4] + 1525.0/24*M[5] - 25.0/4*M[2];                                                                                        break;
    		case 19   :    return 3125.0/24*M[19] - 625.0/4*M[13] + 1375.0/24*M[8] - 25.0/4*M[4];                                                                                                                            break;
    		case 20   :    return 625.0/24*M[20] - 625.0/12*M[14] + 875.0/24*M[9] - 125.0/12*M[5] + M[2];                                                                                                                            break;
    	   	}    break;
    	}
    }


    /** \brief  Lagrange polynomial for tetrahedron
     *  \param[in]  vertexIndex		The internal index of the associated interpolatory vertex
     *  \param[in]  u			first local coordinate
     *  \param[in]  U[1][1]			second local coordinate
     *  \param[in]  U[2][1]			third local coordinate
     */
    double lagrangePolynomialTetrahedron(InternalIndexType vertexIndex, const NumVector1D & M) const {
    	switch(order_)
    	{
    	    case 1 :
    	        switch (vertexIndex)
    	        {
    	      case 0    :    return -M[1] - M[2] - M[3] + 1;    break;
    	      case 1    :    return M[1];        break;
    	      case 2    :    return M[2];        break;
    	      case 3    :    return M[3];        break;
    	        }    break;
    	    case 2:
    	        switch (vertexIndex)
    	        {
    	      case 0    :    return 2*M[4] + 4*M[5] + 2*M[7] + 4*M[6] + 4*M[8] + 2*M[9] - 3*M[1] - 3*M[2] - 3*M[3] + 1;    break;
    	      case 1    :    return -4*M[4] - 4*M[5] - 4*M[6] + 4*M[1];                        break;
    	      case 2    :    return 2*M[4] - M[1];                                break;
    	      case 3    :    return -4*M[5] - 4*M[7] - 4*M[8] + 4*M[2];                        break;
    	      case 4    :    return 4*M[5];                                    break;
    	      case 5    :    return 2*M[7] - M[2];                                break;
    	      case 6    :    return -4*M[6] - 4*M[8] - 4*M[9] + 4*M[3];                        break;
    	      case 7    :    return 4*M[6];                                    break;
    	      case 8    :    return 4*M[8];                                    break;
    	      case 9    :    return 2*M[9] - M[3];                                break;
    	        }    break;
    	    case 3:
    	        switch (vertexIndex)
    	        {
    	      case 0    :    return -9.0/2*M[10] - 27.0/2*M[11] - 27.0/2*M[13] - 9.0/2*M[16] - 27.0/2*M[12] - 27*M[14] - 27.0/2*M[17] - 27.0/2*M[15] - 27.0/2*M[18] - 9.0/2*M[19] + 9*M[4] + 18*M[5] + 9*M[7] + 18*M[6] + 18*M[8] + 9*M[9] - 11.0/2*M[1] - 11.0/2*M[2] - 11.0/2*M[3] + 1;    break;
    	      case 1    :    return 27.0/2*M[10] + 27*M[11] + 27.0/2*M[13] + 27*M[12] + 27*M[14] + 27.0/2*M[15] - 45.0/2*M[4] - 45.0/2*M[5] - 45.0/2*M[6] + 9*M[1];                                                    break;
    	      case 2    :    return -27.0/2*M[10] - 27.0/2*M[11] - 27.0/2*M[12] + 18*M[4] + 9.0/2*M[5] + 9.0/2*M[6] - 9.0/2*M[1];                                                                    break;
    	      case 3    :    return 9.0/2*M[10] - 9.0/2*M[4] + M[1];                                                                                                    break;
    	      case 4    :    return 27.0/2*M[11] + 27*M[13] + 27.0/2*M[16] + 27*M[14] + 27*M[17] + 27.0/2*M[18] - 45.0/2*M[5] - 45.0/2*M[7] - 45.0/2*M[8] + 9*M[2];                                                    break;
    	      case 5    :    return -27*M[11] - 27*M[13] - 27*M[14] + 27*M[5];                                                                                            break;
    	      case 6    :    return 27.0/2*M[11] - 9.0/2*M[5];                                                                                                    break;
    	      case 7    :    return -27.0/2*M[13] - 27.0/2*M[16] - 27.0/2*M[17] + 9.0/2*M[5] + 18*M[7] + 9.0/2*M[8] - 9.0/2*M[2];                                                                    break;
    	      case 8    :    return 27.0/2*M[13] - 9.0/2*M[5];                                                                                                    break;
    	      case 9    :    return 9.0/2*M[16] - 9.0/2*M[7] + M[2];                                                                                                    break;
    	      case 10    :    return 27.0/2*M[12] + 27*M[14] + 27.0/2*M[17] + 27*M[15] + 27*M[18] + 27.0/2*M[19] - 45.0/2*M[6] - 45.0/2*M[8] - 45.0/2*M[9] + 9*M[3];                                                    break;
    	      case 11    :    return -27*M[12] - 27*M[14] - 27*M[15] + 27*M[6];                                                                                            break;
    	      case 12    :    return 27.0/2*M[12] - 9.0/2*M[6];                                                                                                    break;
    	      case 13    :    return -27*M[14] - 27*M[17] - 27*M[18] + 27*M[8];                                                                                            break;
    	      case 14    :    return 27*M[14];                                                                                                        break;
    	      case 15    :    return 27.0/2*M[17] - 9.0/2*M[8];                                                                                                    break;
    	      case 16    :    return -27.0/2*M[15] - 27.0/2*M[18] - 27.0/2*M[19] + 9.0/2*M[6] + 9.0/2*M[8] + 18*M[9] - 9.0/2*M[3];                                                                    break;
    	      case 17    :    return 27.0/2*M[15] - 9.0/2*M[6];                                                                                                    break;
    	      case 18    :    return 27.0/2*M[18] - 9.0/2*M[8];                                                                                                    break;
    	      case 19    :    return 9.0/2*M[19] - 9.0/2*M[9] + M[3];                                                                                                    break;
    	        }    break;
    	    case 4:
    	        switch (vertexIndex)
    	        {
    	      case 0    :    return 32.0/3*M[20] + 128.0/3*M[21] + 64*M[23] + 128.0/3*M[26] + 32.0/3*M[30] + 128.0/3*M[22] + 128*M[24] + 128*M[27] + 128.0/3*M[31] + 64*M[25] + 128*M[28] + 64*M[32] + 128.0/3*M[29] + 128.0/3*M[33] + 32.0/3*M[34] - 80.0/3*M[10] - 80*M[11] - 80*M[13] - 80.0/3*M[16] - 80*M[12] - 160*M[14] - 80*M[17] - 80*M[15] - 80*M[18] - 80.0/3*M[19] + 70.0/3*M[4] + 140.0/3*M[5] + 70.0/3*M[7] + 140.0/3*M[6] + 140.0/3*M[8] + 70.0/3*M[9] - 25.0/3*M[1] - 25.0/3*M[2] - 25.0/3*M[3] + 1;    break;
    	      case 1    :    return -128.0/3*M[20] - 128*M[21] - 128*M[23] - 128.0/3*M[26] - 128*M[22] - 256*M[24] - 128*M[27] - 128*M[25] - 128*M[28] - 128.0/3*M[29] + 96*M[10] + 192*M[11] + 96*M[13] + 192*M[12] + 192*M[14] + 96*M[15] - 208.0/3*M[4] - 208.0/3*M[5] - 208.0/3*M[6] + 16*M[1];                                                                                                break;
    	      case 2    :    return 64*M[20] + 128*M[21] + 64*M[23] + 128*M[22] + 128*M[24] + 64*M[25] - 128*M[10] - 144*M[11] - 16*M[13] - 144*M[12] - 32*M[14] - 16*M[15] + 76*M[4] + 28*M[5] + 28*M[6] - 12*M[1];                                                                                                                                        break;
    	      case 3    :    return -128.0/3*M[20] - 128.0/3*M[21] - 128.0/3*M[22] + 224.0/3*M[10] + 32*M[11] + 32*M[12] - 112.0/3*M[4] - 16.0/3*M[5] - 16.0/3*M[6] + 16.0/3*M[1];                                                                                                                                                        break;
    	      case 4    :    return 32.0/3*M[20] - 16*M[10] + 22.0/3*M[4] - M[1];                                                                                                                                                                                                    break;
    	      case 5    :    return -128.0/3*M[21] - 128*M[23] - 128*M[26] - 128.0/3*M[30] - 128*M[24] - 256*M[27] - 128*M[31] - 128*M[28] - 128*M[32] - 128.0/3*M[33] + 96*M[11] + 192*M[13] + 96*M[16] + 192*M[14] + 192*M[17] + 96*M[18] - 208.0/3*M[5] - 208.0/3*M[7] - 208.0/3*M[8] + 16*M[2];                                                                                                break;
    	      case 6    :    return 128*M[21] + 256*M[23] + 128*M[26] + 256*M[24] + 256*M[27] + 128*M[28] - 224*M[11] - 224*M[13] - 224*M[14] + 96*M[5];                                                                                                                                                            break;
    	      case 7    :    return -128*M[21] - 128*M[23] - 128*M[24] + 160*M[11] + 32*M[13] + 32*M[14] - 32*M[5];                                                                                                                                                                                break;
    	      case 8    :    return 128.0/3*M[21] - 32*M[11] + 16.0/3*M[5];                                                                                                                                                                                                    break;
    	      case 9    :    return 64*M[23] + 128*M[26] + 64*M[30] + 128*M[27] + 128*M[31] + 64*M[32] - 16*M[11] - 144*M[13] - 128*M[16] - 32*M[14] - 144*M[17] - 16*M[18] + 28*M[5] + 76*M[7] + 28*M[8] - 12*M[2];                                                                                                                                        break;
    	      case 10    :    return -128*M[23] - 128*M[26] - 128*M[27] + 32*M[11] + 160*M[13] + 32*M[14] - 32*M[5];                                                                                                                                                                                break;
    	      case 11    :    return 64*M[23] - 16*M[11] - 16*M[13] + 4*M[5];                                                                                                                                                                                                    break;
    	      case 12    :    return -128.0/3*M[26] - 128.0/3*M[30] - 128.0/3*M[31] + 32*M[13] + 224.0/3*M[16] + 32*M[17] - 16.0/3*M[5] - 112.0/3*M[7] - 16.0/3*M[8] + 16.0/3*M[2];                                                                                                                                                        break;
    	      case 13    :    return 128.0/3*M[26] - 32*M[13] + 16.0/3*M[5];                                                                                                                                                                                                    break;
    	      case 14    :    return 32.0/3*M[30] - 16*M[16] + 22.0/3*M[7] - M[2];                                                                                                                                                                                                    break;
    	      case 15    :    return -128.0/3*M[22] - 128*M[24] - 128*M[27] - 128.0/3*M[31] - 128*M[25] - 256*M[28] - 128*M[32] - 128*M[29] - 128*M[33] - 128.0/3*M[34] + 96*M[12] + 192*M[14] + 96*M[17] + 192*M[15] + 192*M[18] + 96*M[19] - 208.0/3*M[6] - 208.0/3*M[8] - 208.0/3*M[9] + 16*M[3];                                                                                                break;
    	      case 16    :    return 128*M[22] + 256*M[24] + 128*M[27] + 256*M[25] + 256*M[28] + 128*M[29] - 224*M[12] - 224*M[14] - 224*M[15] + 96*M[6];                                                                                                                                                            break;
    	      case 17    :    return -128*M[22] - 128*M[24] - 128*M[25] + 160*M[12] + 32*M[14] + 32*M[15] - 32*M[6];                                                                                                                                                                                break;
    	      case 18    :    return 128.0/3*M[22] - 32*M[12] + 16.0/3*M[6];                                                                                                                                                                                                    break;
    	      case 19    :    return 128*M[24] + 256*M[27] + 128*M[31] + 256*M[28] + 256*M[32] + 128*M[33] - 224*M[14] - 224*M[17] - 224*M[18] + 96*M[8];                                                                                                                                                            break;
    	      case 20    :    return -256*M[24] - 256*M[27] - 256*M[28] + 256*M[14];                                                                                                                                                                                            break;
    	      case 21    :    return 128*M[24] - 32*M[14];                                                                                                                                                                                                            break;
    	      case 22    :    return -128*M[27] - 128*M[31] - 128*M[32] + 32*M[14] + 160*M[17] + 32*M[18] - 32*M[8];                                                                                                                                                                                break;
    	      case 23 :    return 128*M[27] - 32*M[14];                                                                                                                                                                                                            break;
    	      case 24    :    return 128.0/3*M[31] - 32*M[17] + 16.0/3*M[8];                                                                                                                                                                                                    break;
    	      case 25    :    return 64*M[25] + 128*M[28] + 64*M[32] + 128*M[29] + 128*M[33] + 64*M[34] - 16*M[12] - 32*M[14] - 16*M[17] - 144*M[15] - 144*M[18] - 128*M[19] + 28*M[6] + 28*M[8] + 76*M[9] - 12*M[3];                                                                                                                                        break;
    	      case 26    :    return -128*M[25] - 128*M[28] - 128*M[29] + 32*M[12] + 32*M[14] + 160*M[15] - 32*M[6];                                                                                                                                                                                break;
    	      case 27    :    return 64*M[25] - 16*M[12] - 16*M[15] + 4*M[6];                                                                                                                                                                                                    break;
    	      case 28    :    return -128*M[28] - 128*M[32] - 128*M[33] + 32*M[14] + 32*M[17] + 160*M[18] - 32*M[8];                                                                                                                                                                                break;
    	      case 29    :    return 128*M[28] - 32*M[14];    break;
    	      case 30    :    return 64*M[32] - 16*M[17] - 16*M[18] + 4*M[8];    break;
    	      case 31    :    return -128.0/3*M[29] - 128.0/3*M[33] - 128.0/3*M[34] + 32*M[15] + 32*M[18] + 224.0/3*M[19] - 16.0/3*M[6] - 16.0/3*M[8] - 112.0/3*M[9] + 16.0/3*M[3];    break;
    	      case 32    :    return 128.0/3*M[29] - 32*M[15] + 16.0/3*M[6];    break;
    	      case 33    :    return 128.0/3*M[33] - 32*M[18] + 16.0/3*M[8];    break;
    	      case 34    :    return 32.0/3*M[34] - 16*M[19] + 22.0/3*M[9] - M[3];    break;
    	        }    break;
    	    case 5:
    	        switch (vertexIndex)
    	        {
    	      case 0    :    return -625.0/24*M[35] - 3125.0/24*M[36] - 3125.0/12*M[38] - 3125.0/12*M[41] - 3125.0/24*M[45] - 625.0/24*M[50] - 3125.0/24*M[37] - 3125.0/6*M[39] - 3125.0/4*M[42] - 3125.0/6*M[46] - 3125.0/24*M[51] - 3125.0/12*M[40] - 3125.0/4*M[43] - 3125.0/4*M[47] - 3125.0/12*M[52] - 3125.0/12*M[44] - 3125.0/6*M[48] - 3125.0/12*M[53] - 3125.0/24*M[49] - 3125.0/24*M[54] - 625.0/24*M[55] + 625.0/8*M[20] + 625.0/2*M[21] + 1875.0/4*M[23] + 625.0/2*M[26] + 625.0/8*M[30] + 625.0/2*M[22] + 1875.0/2*M[24] + 1875.0/2*M[27] + 625.0/2*M[31] + 1875.0/4*M[25] + 1875.0/2*M[28] + 1875.0/4*M[32] + 625.0/2*M[29] + 625.0/2*M[33] + 625.0/8*M[34] - 2125.0/24*M[10] - 2125.0/8*M[11] - 2125.0/8*M[13] - 2125.0/24*M[16] - 2125.0/8*M[12] - 2125.0/4*M[14] - 2125.0/8*M[17] - 2125.0/8*M[15] - 2125.0/8*M[18] - 2125.0/24*M[19] + 375.0/8*M[4] + 375.0/4*M[5] + 375.0/8*M[7] + 375.0/4*M[6] + 375.0/4*M[8] + 375.0/8*M[9] - 137.0/12*M[1] - 137.0/12*M[2] - 137.0/12*M[3] + 1;    break;
    	      case 1    :    return 3125.0/24*M[35] + 3125.0/6*M[36] + 3125.0/4*M[38] + 3125.0/6*M[41] + 3125.0/24*M[45] + 3125.0/6*M[37] + 3125.0/2*M[39] + 3125.0/2*M[42] + 3125.0/6*M[46] + 3125.0/4*M[40] + 3125.0/2*M[43] + 3125.0/4*M[47] + 3125.0/6*M[44] + 3125.0/6*M[48] + 3125.0/24*M[49] - 4375.0/12*M[20] - 4375.0/4*M[21] - 4375.0/4*M[23] - 4375.0/12*M[26] - 4375.0/4*M[22] - 4375.0/2*M[24] - 4375.0/4*M[27] - 4375.0/4*M[25] - 4375.0/4*M[28] - 4375.0/12*M[29] + 8875.0/24*M[10] + 8875.0/12*M[11] + 8875.0/24*M[13] + 8875.0/12*M[12] + 8875.0/12*M[14] + 8875.0/24*M[15] - 1925.0/12*M[4] - 1925.0/12*M[5] - 1925.0/12*M[6] + 25*M[1];    break;
    	      case 2    :    return -3125.0/12*M[35] - 3125.0/4*M[36] - 3125.0/4*M[38] - 3125.0/12*M[41] - 3125.0/4*M[37] - 3125.0/2*M[39] - 3125.0/4*M[42] - 3125.0/4*M[40] - 3125.0/4*M[43] - 3125.0/12*M[44] + 8125.0/12*M[20] + 5625.0/4*M[21] + 3125.0/4*M[23] + 625.0/12*M[26] + 5625.0/4*M[22] + 3125.0/2*M[24] + 625.0/4*M[27] + 3125.0/4*M[25] + 625.0/4*M[28] + 625.0/12*M[29] - 7375.0/12*M[10] - 8875.0/12*M[11] - 125*M[13] - 8875.0/12*M[12] - 250*M[14] - 125*M[15] + 2675.0/12*M[4] + 1175.0/12*M[5] + 1175.0/12*M[6] - 25*M[1];    break;
    	      case 3    :    return 3125.0/12*M[35] + 3125.0/6*M[36] + 3125.0/12*M[38] + 3125.0/6*M[37] + 3125.0/6*M[39] + 3125.0/12*M[40] - 625*M[20] - 3125.0/4*M[21] - 625.0/4*M[23] - 3125.0/4*M[22] - 625.0/2*M[24] - 625.0/4*M[25] + 6125.0/12*M[10] + 3875.0/12*M[11] + 125.0/6*M[13] + 3875.0/12*M[12] + 125.0/3*M[14] + 125.0/6*M[15] - 325.0/2*M[4] - 75.0/2*M[5] - 75.0/2*M[6] + 50.0/3*M[1];    break;
    	      case 4    :    return -3125.0/24*M[35] - 3125.0/24*M[36] - 3125.0/24*M[37] + 6875.0/24*M[20] + 625.0/4*M[21] + 625.0/4*M[22] - 5125.0/24*M[10] - 1375.0/24*M[11] - 1375.0/24*M[12] + 1525.0/24*M[4] + 25.0/4*M[5] + 25.0/4*M[6] - 25.0/4*M[1];    break;
    	      case 5    :    return 625.0/24*M[35] - 625.0/12*M[20] + 875.0/24*M[10] - 125.0/12*M[4] + M[1];    break;
    	      case 6    :    return 3125.0/24*M[36] + 3125.0/6*M[38] + 3125.0/4*M[41] + 3125.0/6*M[45] + 3125.0/24*M[50] + 3125.0/6*M[39] + 3125.0/2*M[42] + 3125.0/2*M[46] + 3125.0/6*M[51] + 3125.0/4*M[43] + 3125.0/2*M[47] + 3125.0/4*M[52] + 3125.0/6*M[48] + 3125.0/6*M[53] + 3125.0/24*M[54] - 4375.0/12*M[21] - 4375.0/4*M[23] - 4375.0/4*M[26] - 4375.0/12*M[30] - 4375.0/4*M[24] - 4375.0/2*M[27] - 4375.0/4*M[31] - 4375.0/4*M[28] - 4375.0/4*M[32] - 4375.0/12*M[33] + 8875.0/24*M[11] + 8875.0/12*M[13] + 8875.0/24*M[16] + 8875.0/12*M[14] + 8875.0/12*M[17] + 8875.0/24*M[18] - 1925.0/12*M[5] - 1925.0/12*M[7] - 1925.0/12*M[8] + 25*M[2];    break;
    	      case 7    :    return -3125.0/6*M[36] - 3125.0/2*M[38] - 3125.0/2*M[41] - 3125.0/6*M[45] - 3125.0/2*M[39] - 3125*M[42] - 3125.0/2*M[46] - 3125.0/2*M[43] - 3125.0/2*M[47] - 3125.0/6*M[48] + 1250*M[21] + 2500*M[23] + 1250*M[26] + 2500*M[24] + 2500*M[27] + 1250*M[28] - 5875.0/6*M[11] - 5875.0/6*M[13] - 5875.0/6*M[14] + 250*M[5];    break;
    	      case 8    :    return 3125.0/4*M[36] + 3125.0/2*M[38] + 3125.0/4*M[41] + 3125.0/2*M[39] + 3125.0/2*M[42] + 3125.0/4*M[43] - 3125.0/2*M[21] - 6875.0/4*M[23] - 625.0/4*M[26] - 6875.0/4*M[24] - 625.0/2*M[27] - 625.0/4*M[28] + 3625.0/4*M[11] + 1125.0/4*M[13] + 1125.0/4*M[14] - 125*M[5];    break;
    	      case 9    :    return -3125.0/6*M[36] - 3125.0/6*M[38] - 3125.0/6*M[39] + 2500.0/3*M[21] + 625.0/2*M[23] + 625.0/2*M[24] - 2125.0/6*M[11] - 125.0/3*M[13] - 125.0/3*M[14] + 125.0/3*M[5];    break;
    	      case 10    :    return 3125.0/24*M[36] - 625.0/4*M[21] + 1375.0/24*M[11] - 25.0/4*M[5];    break;
    	      case 11    :    return -3125.0/12*M[38] - 3125.0/4*M[41] - 3125.0/4*M[45] - 3125.0/12*M[50] - 3125.0/4*M[42] - 3125.0/2*M[46] - 3125.0/4*M[51] - 3125.0/4*M[47] - 3125.0/4*M[52] - 3125.0/12*M[53] + 625.0/12*M[21] + 3125.0/4*M[23] + 5625.0/4*M[26] + 8125.0/12*M[30] + 625.0/4*M[24] + 3125.0/2*M[27] + 5625.0/4*M[31] + 625.0/4*M[28] + 3125.0/4*M[32] + 625.0/12*M[33] - 125*M[11] - 8875.0/12*M[13] - 7375.0/12*M[16] - 250*M[14] - 8875.0/12*M[17] - 125*M[18] + 1175.0/12*M[5] + 2675.0/12*M[7] + 1175.0/12*M[8] - 25*M[2];    break;
    	      case 12    :    return 3125.0/4*M[38] + 3125.0/2*M[41] + 3125.0/4*M[45] + 3125.0/2*M[42] + 3125.0/2*M[46] + 3125.0/4*M[47] - 625.0/4*M[21] - 6875.0/4*M[23] - 3125.0/2*M[26] - 625.0/2*M[24] - 6875.0/4*M[27] - 625.0/4*M[28] + 1125.0/4*M[11] + 3625.0/4*M[13] + 1125.0/4*M[14] - 125*M[5];    break;
    	      case 13    :    return -3125.0/4*M[38] - 3125.0/4*M[41] - 3125.0/4*M[42] + 625.0/4*M[21] + 4375.0/4*M[23] + 625.0/4*M[26] + 625.0/4*M[24] + 625.0/4*M[27] - 375.0/2*M[11] - 375.0/2*M[13] - 125.0/4*M[14] + 125.0/4*M[5];    break;
    	      case 14    :    return 3125.0/12*M[38] - 625.0/12*M[21] - 625.0/4*M[23] + 125.0/4*M[11] + 125.0/6*M[13] - 25.0/6*M[5];    break;
    	      case 15    :    return 3125.0/12*M[41] + 3125.0/6*M[45] + 3125.0/12*M[50] + 3125.0/6*M[46] + 3125.0/6*M[51] + 3125.0/12*M[52] - 625.0/4*M[23] - 3125.0/4*M[26] - 625*M[30] - 625.0/2*M[27] - 3125.0/4*M[31] - 625.0/4*M[32] + 125.0/6*M[11] + 3875.0/12*M[13] + 6125.0/12*M[16] + 125.0/3*M[14] + 3875.0/12*M[17] + 125.0/6*M[18] - 75.0/2*M[5] - 325.0/2*M[7] - 75.0/2*M[8] + 50.0/3*M[2];    break;
    	      case 16    :    return -3125.0/6*M[41] - 3125.0/6*M[45] - 3125.0/6*M[46] + 625.0/2*M[23] + 2500.0/3*M[26] + 625.0/2*M[27] - 125.0/3*M[11] - 2125.0/6*M[13] - 125.0/3*M[14] + 125.0/3*M[5];    break;
    	      case 17    :    return 3125.0/12*M[41] - 625.0/4*M[23] - 625.0/12*M[26] + 125.0/6*M[11] + 125.0/4*M[13] - 25.0/6*M[5];    break;
    	      case 18    :    return -3125.0/24*M[45] - 3125.0/24*M[50] - 3125.0/24*M[51] + 625.0/4*M[26] + 6875.0/24*M[30] + 625.0/4*M[31] - 1375.0/24*M[13] - 5125.0/24*M[16] - 1375.0/24*M[17] + 25.0/4*M[5] + 1525.0/24*M[7] + 25.0/4*M[8] - 25.0/4*M[2];    break;
    	      case 19    :    return 3125.0/24*M[45] - 625.0/4*M[26] + 1375.0/24*M[13] - 25.0/4*M[5];    break;
    	      case 20    :    return 625.0/24*M[50] - 625.0/12*M[30] + 875.0/24*M[16] - 125.0/12*M[7] + M[2];    break;
    	      case 21    :    return 3125.0/24*M[37] + 3125.0/6*M[39] + 3125.0/4*M[42] + 3125.0/6*M[46] + 3125.0/24*M[51] + 3125.0/6*M[40] + 3125.0/2*M[43] + 3125.0/2*M[47] + 3125.0/6*M[52] + 3125.0/4*M[44] + 3125.0/2*M[48] + 3125.0/4*M[53] + 3125.0/6*M[49] + 3125.0/6*M[54] + 3125.0/24*M[55] - 4375.0/12*M[22] - 4375.0/4*M[24] - 4375.0/4*M[27] - 4375.0/12*M[31] - 4375.0/4*M[25] - 4375.0/2*M[28] - 4375.0/4*M[32] - 4375.0/4*M[29] - 4375.0/4*M[33] - 4375.0/12*M[34] + 8875.0/24*M[12] + 8875.0/12*M[14] + 8875.0/24*M[17] + 8875.0/12*M[15] + 8875.0/12*M[18] + 8875.0/24*M[19] - 1925.0/12*M[6] - 1925.0/12*M[8] - 1925.0/12*M[9] + 25*M[3];    break;
    	      case 22    :    return -3125.0/6*M[37] - 3125.0/2*M[39] - 3125.0/2*M[42] - 3125.0/6*M[46] - 3125.0/2*M[40] - 3125*M[43] - 3125.0/2*M[47] - 3125.0/2*M[44] - 3125.0/2*M[48] - 3125.0/6*M[49] + 1250*M[22] + 2500*M[24] + 1250*M[27] + 2500*M[25] + 2500*M[28] + 1250*M[29] - 5875.0/6*M[12] - 5875.0/6*M[14] - 5875.0/6*M[15] + 250*M[6];    break;
    	      case 23    :    return 3125.0/4*M[37] + 3125.0/2*M[39] + 3125.0/4*M[42] + 3125.0/2*M[40] + 3125.0/2*M[43] + 3125.0/4*M[44] - 3125.0/2*M[22] - 6875.0/4*M[24] - 625.0/4*M[27] - 6875.0/4*M[25] - 625.0/2*M[28] - 625.0/4*M[29] + 3625.0/4*M[12] + 1125.0/4*M[14] + 1125.0/4*M[15] - 125*M[6];    break;
    	      case 24    :    return -3125.0/6*M[37] - 3125.0/6*M[39] - 3125.0/6*M[40] + 2500.0/3*M[22] + 625.0/2*M[24] + 625.0/2*M[25] - 2125.0/6*M[12] - 125.0/3*M[14] - 125.0/3*M[15] + 125.0/3*M[6];    break;
    	      case 25    :    return 3125.0/24*M[37] - 625.0/4*M[22] + 1375.0/24*M[12] - 25.0/4*M[6];    break;
    	      case 26    :    return -3125.0/6*M[39] - 3125.0/2*M[42] - 3125.0/2*M[46] - 3125.0/6*M[51] - 3125.0/2*M[43] - 3125*M[47] - 3125.0/2*M[52] - 3125.0/2*M[48] - 3125.0/2*M[53] - 3125.0/6*M[54] + 1250*M[24] + 2500*M[27] + 1250*M[31] + 2500*M[28] + 2500*M[32] + 1250*M[33] - 5875.0/6*M[14] - 5875.0/6*M[17] - 5875.0/6*M[18] + 250*M[8];    break;
    	      case 27    :    return 3125.0/2*M[39] + 3125*M[42] + 3125.0/2*M[46] + 3125*M[43] + 3125*M[47] + 3125.0/2*M[48] - 5625.0/2*M[24] - 5625.0/2*M[27] - 5625.0/2*M[28] + 1250*M[14];    break;
    	      case 28    :    return -3125.0/2*M[39] - 3125.0/2*M[42] - 3125.0/2*M[43] + 1875*M[24] + 625.0/2*M[27] + 625.0/2*M[28] - 625.0/2*M[14];    break;
    	      case 29    :    return 3125.0/6*M[39] - 625.0/2*M[24] + 125.0/3*M[14];    break;
    	      case 30    :    return 3125.0/4*M[42] + 3125.0/2*M[46] + 3125.0/4*M[51] + 3125.0/2*M[47] + 3125.0/2*M[52] + 3125.0/4*M[53] - 625.0/4*M[24] - 6875.0/4*M[27] - 3125.0/2*M[31] - 625.0/2*M[28] - 6875.0/4*M[32] - 625.0/4*M[33] + 1125.0/4*M[14] + 3625.0/4*M[17] + 1125.0/4*M[18] - 125*M[8];    break;
    	      case 31    :    return -3125.0/2*M[42] - 3125.0/2*M[46] - 3125.0/2*M[47] + 625.0/2*M[24] + 1875*M[27] + 625.0/2*M[28] - 625.0/2*M[14];    break;
    	      case 32    :    return 3125.0/4*M[42] - 625.0/4*M[24] - 625.0/4*M[27] + 125.0/4*M[14];    break;
    	      case 33    :    return -3125.0/6*M[46] - 3125.0/6*M[51] - 3125.0/6*M[52] + 625.0/2*M[27] + 2500.0/3*M[31] + 625.0/2*M[32] - 125.0/3*M[14] - 2125.0/6*M[17] - 125.0/3*M[18] + 125.0/3*M[8];    break;
    	      case 34    :    return 3125.0/6*M[46] - 625.0/2*M[27] + 125.0/3*M[14];    break;
    	      case 35    :    return 3125.0/24*M[51] - 625.0/4*M[31] + 1375.0/24*M[17] - 25.0/4*M[8];    break;
    	      case 36    :    return -3125.0/12*M[40] - 3125.0/4*M[43] - 3125.0/4*M[47] - 3125.0/12*M[52] - 3125.0/4*M[44] - 3125.0/2*M[48] - 3125.0/4*M[53] - 3125.0/4*M[49] - 3125.0/4*M[54] - 3125.0/12*M[55] + 625.0/12*M[22] + 625.0/4*M[24] + 625.0/4*M[27] + 625.0/12*M[31] + 3125.0/4*M[25] + 3125.0/2*M[28] + 3125.0/4*M[32] + 5625.0/4*M[29] + 5625.0/4*M[33] + 8125.0/12*M[34] - 125*M[12] - 250*M[14] - 125*M[17] - 8875.0/12*M[15] - 8875.0/12*M[18] - 7375.0/12*M[19] + 1175.0/12*M[6] + 1175.0/12*M[8] + 2675.0/12*M[9] - 25*M[3];    break;
    	      case 37    :    return 3125.0/4*M[40] + 3125.0/2*M[43] + 3125.0/4*M[47] + 3125.0/2*M[44] + 3125.0/2*M[48] + 3125.0/4*M[49] - 625.0/4*M[22] - 625.0/2*M[24] - 625.0/4*M[27] - 6875.0/4*M[25] - 6875.0/4*M[28] - 3125.0/2*M[29] + 1125.0/4*M[12] + 1125.0/4*M[14] + 3625.0/4*M[15] - 125*M[6];    break;
    	      case 38    :    return -3125.0/4*M[40] - 3125.0/4*M[43] - 3125.0/4*M[44] + 625.0/4*M[22] + 625.0/4*M[24] + 4375.0/4*M[25] + 625.0/4*M[28] + 625.0/4*M[29] - 375.0/2*M[12] - 125.0/4*M[14] - 375.0/2*M[15] + 125.0/4*M[6];    break;
    	      case 39    :    return 3125.0/12*M[40] - 625.0/12*M[22] - 625.0/4*M[25] + 125.0/4*M[12] + 125.0/6*M[15] - 25.0/6*M[6];    break;
    	      case 40    :    return 3125.0/4*M[43] + 3125.0/2*M[47] + 3125.0/4*M[52] + 3125.0/2*M[48] + 3125.0/2*M[53] + 3125.0/4*M[54] - 625.0/4*M[24] - 625.0/2*M[27] - 625.0/4*M[31] - 6875.0/4*M[28] - 6875.0/4*M[32] - 3125.0/2*M[33] + 1125.0/4*M[14] + 1125.0/4*M[17] + 3625.0/4*M[18] - 125*M[8];    break;
    	      case 41    :    return -3125.0/2*M[43] - 3125.0/2*M[47] - 3125.0/2*M[48] + 625.0/2*M[24] + 625.0/2*M[27] + 1875*M[28] - 625.0/2*M[14];    break;
    	      case 42    :    return 3125.0/4*M[43] - 625.0/4*M[24] - 625.0/4*M[28] + 125.0/4*M[14];    break;
    	      case 43    :    return -3125.0/4*M[47] - 3125.0/4*M[52] - 3125.0/4*M[53] + 625.0/4*M[27] + 625.0/4*M[31] + 625.0/4*M[28] + 4375.0/4*M[32] + 625.0/4*M[33] - 125.0/4*M[14] - 375.0/2*M[17] - 375.0/2*M[18] + 125.0/4*M[8];    break;
    	      case 44    :    return 3125.0/4*M[47] - 625.0/4*M[27] - 625.0/4*M[28] + 125.0/4*M[14];    break;
    	      case 45    :    return 3125.0/12*M[52] - 625.0/12*M[31] - 625.0/4*M[32] + 125.0/4*M[17] + 125.0/6*M[18] - 25.0/6*M[8];    break;
    	      case 46    :    return 3125.0/12*M[44] + 3125.0/6*M[48] + 3125.0/12*M[53] + 3125.0/6*M[49] + 3125.0/6*M[54] + 3125.0/12*M[55] - 625.0/4*M[25] - 625.0/2*M[28] - 625.0/4*M[32] - 3125.0/4*M[29] - 3125.0/4*M[33] - 625*M[34] + 125.0/6*M[12] + 125.0/3*M[14] + 125.0/6*M[17] + 3875.0/12*M[15] + 3875.0/12*M[18] + 6125.0/12*M[19] - 75.0/2*M[6] - 75.0/2*M[8] - 325.0/2*M[9] + 50.0/3*M[3];    break;
    	      case 47    :    return -3125.0/6*M[44] - 3125.0/6*M[48] - 3125.0/6*M[49] + 625.0/2*M[25] + 625.0/2*M[28] + 2500.0/3*M[29] - 125.0/3*M[12] - 125.0/3*M[14] - 2125.0/6*M[15] + 125.0/3*M[6];    break;
    	      case 48    :    return 3125.0/12*M[44] - 625.0/4*M[25] - 625.0/12*M[29] + 125.0/6*M[12] + 125.0/4*M[15] - 25.0/6*M[6];    break;
    	      case 49    :    return -3125.0/6*M[48] - 3125.0/6*M[53] - 3125.0/6*M[54] + 625.0/2*M[28] + 625.0/2*M[32] + 2500.0/3*M[33] - 125.0/3*M[14] - 125.0/3*M[17] - 2125.0/6*M[18] + 125.0/3*M[8];    break;
    	      case 50    :    return 3125.0/6*M[48] - 625.0/2*M[28] + 125.0/3*M[14];    break;
    	      case 51    :    return 3125.0/12*M[53] - 625.0/4*M[32] - 625.0/12*M[33] + 125.0/6*M[17] + 125.0/4*M[18] - 25.0/6*M[8];    break;
    	      case 52    :    return -3125.0/24*M[49] - 3125.0/24*M[54] - 3125.0/24*M[55] + 625.0/4*M[29] + 625.0/4*M[33] + 6875.0/24*M[34] - 1375.0/24*M[15] - 1375.0/24*M[18] - 5125.0/24*M[19] + 25.0/4*M[6] + 25.0/4*M[8] + 1525.0/24*M[9] - 25.0/4*M[3];    break;
    	      case 53    :    return 3125.0/24*M[49] - 625.0/4*M[29] + 1375.0/24*M[15] - 25.0/4*M[6];    break;
    	      case 54    :    return 3125.0/24*M[54] - 625.0/4*M[33] + 1375.0/24*M[18] - 25.0/4*M[8];    break;
    	      case 55    :    return 625.0/24*M[55] - 625.0/12*M[34] + 875.0/24*M[19] - 125.0/12*M[9] + M[3];    break;
    	        }    break;
    	}
    }


    /** \brief  Analytic map from local to global coordinates, given explicitly by the polynomial class. Implementation for simplex  */
    PolynomialGlobalCoordinate interpolatoryVectorAnalyticalSimplex() const {
        int dof = dofPerOrder();
        PolynomialGlobalCoordinate rez;

        // Specialization for points - just return a constant polynomial vector
    	if (mydim == 0) {
    		for (int d = 0; d < coorddimension; d++)
    		{
    			rez[d] = LocalPolynomial(Monomial(point_[0][d], std::vector<int>()));
    		}
    		return rez;
    	}

        PolynomialVector monomial_basis;
        PolynomialVector lagrange_basis;

        // Step 0. Construct an interpolatory grid over the element
        IntegerCoordinateVector simplexPoints = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(order_);
        std::vector<LocalCoordinate> localCoordinateSet = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<ctype, mydim>(simplexPoints, order_);

        // Step 1. Construct monomial basis and set of local coordinates
        for (unsigned int i = 0; i < simplexPoints.size(); i++)
        {
                monomial_basis.push_back(LocalPolynomial(Monomial(1.0, simplexPoints[i])));
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

            rez[d] = tmpPoly;  // Write this global coordinate to a vector
            rez[d].cleanUp();  // Remove all terms that are insignificantly small
        }

        return rez;
    }


}; // SimplexInterpolator




} // Namespace Dune

#endif /** DUNE_LAGRANGE_INTERPOLATOR_HH **/
