/*******************************************************************
 * Curvilinear Geometry Helper
 * 
 * author: Aleksejs Fomins
 * date: 01.10.2014 - created
 * 
 * description:
 * 
 * All sort of common functionality
 *
 * 
 * 
 *******************************************************************/



#ifndef DUNE_CURVILINEARGEOMETRYHELPER_HH
#define DUNE_CURVILINEARGEOMETRYHELPER_HH

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


namespace Dune {


class CurvilinearGeometryHelper {
  public:

    // Public typedefs
	// **************************************************************
    typedef  int     InternalIndexType;
    typedef  int     InterpolatoryOrderType;

    typedef std::vector<std::vector<int>> IntegerCoordinateVector;




    /** \brief Empty constructor. This class only contains static members anyway */
    CurvilinearGeometryHelper() {}


    /** \brief Returns the name of the geometry based on the geometry type. Useful for logging messages
     *  \param[in]  gt	    GeometryType of the element of interest
     */
    static std::string geometryName (Dune::GeometryType gt)
    {
    	if (gt.isVertex())         { return "vertex"; }
    	if (gt.isLine())           { return "edge"; }
    	if (gt.isTriangle())       { return "triangle"; }
    	if (gt.isQuadrilateral())  { return "quadrilateral"; }
    	if (gt.isCube())           { return "cube"; }
    	if (gt.isHexahedron())     { return "hexahedron"; }
    	if (gt.isTetrahedron())    { return "tetrahedron"; }
    	if (gt.isPrism())          { return "prism"; }
    	if (gt.isPyramid())        { return "pyramid"; }
    }


    /** \brief Number of degrees of freedom of an element (number of interpolation points)
     *  \param[in]  geomType	GeometryType of the element of interest
     *  \param[in]  order	Interpolation Order of the element of interest
     */
    static int dofPerOrder(Dune::GeometryType geomType, int order)
    {
    	int rez;

        if (!geomType.isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: dofPerOrder() only implemented for Simplex geometries at the moment"); }

        const int curvilinearDofPerOrderEdge_[5]        = {2,3,4,5,6};
        const int curvilinearDofPerOrderTriangle_[5]    = {3,6,10,15,21};
        const int curvilinearDofPerOrderTetrahedron_[5] = {4,10,20,35,56};


        switch (geomType.dim()) {
        case 1:  rez = curvilinearDofPerOrderEdge_[order - 1];         break;
        case 2:  rez = curvilinearDofPerOrderTriangle_[order - 1];     break;
        case 3:  rez = curvilinearDofPerOrderTetrahedron_[order - 1];  break;
        }

        return rez;
    }


    /** \brief Maps from internal corner index to internal vertex index of an element
     *  \param[in]  geomType	GeometryType of the element of interest
     *  \param[in]  order	Interpolation Order of the element of interest
     *  \param[in]  i		Internal corner index of the corner of interest
     */
    static InternalIndexType cornerIndex(Dune::GeometryType geomType, int order, InternalIndexType i)
    {
        if (!geomType.isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: cornerIndex() only implemented for Simplex geometries at the moment"); }

        int ind = 0;
        int dofperorder = dofPerOrder(geomType, order);

        switch(geomType.dim()) {
        case 1 :  // EDGE
        {
            switch (i)
            {
            case 0 : ind = 0;                  break;
            case 1 : ind = order;              break;
            }
        } break;
        case 2 : // TRIANGLE
        {
            switch (i)
            {
            case 0 : ind = 0;                break;
            case 1 : ind = order;            break;
            case 2 : ind = dofperorder - 1;  break;
            }
        } break;

        case 3 : // TETRAHEDRON
        {
            switch (i)
            {
            case 0 : ind = 0;                       break;
            case 1 : ind = order;                   break;
            case 2 : ind = order*(order + 3) / 2;   break;
            case 3 : ind = dofperorder - 1;         break;
            }
        } break;
        }

        return ind;
    }


    /** \brief Generates an integer vector corresponding to the internal coordinate of one of the corners of an entity
     *
     *  \param[in]  gt 		geometry type of the entity
     *  \param[in]  subInd	corner subentity index
     */
    template <typename ctype, int cdim>
    static Dune::FieldVector<ctype, cdim> cornerInternalCoordinate(GeometryType gt, InternalIndexType subInd)
    {
    	int mydim = gt.dim();
    	assert(mydim <= cdim);    // Embedded geometry can not exceed the geometry of the world
    	assert(gt.isSimplex());   // So far only simplex geometries allowed
    	assert(subInd <= mydim);  // In simplex geometries nCorners = mydim + 1

    	Dune::FieldVector<ctype, cdim> rez;         // By convention initialized as 0-vector
    	if (subInd > 0) { rez[subInd - 1] = 1; }  // Thus generating coordinates {[0,0,0], [1,0,0], [0,1,0], [0,0,1]}

    	return rez;
    }


    /** \brief Generates a vector of integer vectors which label the positions of a regular grid over a simplex
     *  \param[in]  n 		the number of linear intervals
     */
    template <int mydim>
    static IntegerCoordinateVector simplexGridEnumerate(int n)
    {
        IntegerCoordinateVector rez;

        switch (mydim)
        {
        case 1:
            for (int i = 0; i <= n; i++)
            {
                std::vector<int> tmp;
                tmp.push_back(i);
                rez.push_back(tmp);
            } break;
        case 2:
            for (int i = 0; i <= n; i++)
            {
                for (int j = 0; j <= n - i; j++)
                {
                    std::vector<int> tmp;
                    tmp.push_back(j);
                    tmp.push_back(i);
                    rez.push_back(tmp);
                }
            } break;
        case 3:
            for (int i = 0; i <= n; i++)
            {
                for (int j = 0; j <= n - i; j++)
                {
                    for (int k = 0; k <= n - i - j; k++)
                    {
                        std::vector<int> tmp;
                        tmp.push_back(k);
                        tmp.push_back(j);
                        tmp.push_back(i);
                        rez.push_back(tmp);
                    }
                }
            } break;
        }
        return rez;
    }


    /** \brief Generates a vector of coordinates which correspond to a regular grid over a simplex
     *  \param[in]  n 		the number of linear intervals (aka interpolation order)
     * 
     *  note: the coordinates are just the integere labels of simplexGridEnumerate divided by n
     */
    template <class ct, int mydim>
    static std::vector<Dune::FieldVector<ct, mydim> > simplexGridCoordinateSet(int n)
    {
        return simplexGridCoordinateSet<ct, mydim>(simplexGridEnumerate<mydim>(n), n);
    }


    /** \brief Generates a vector of coordinates which correspond to a regular grid over a simplex. Re-uses already calculated simplexGridEnumerate for speedup
     *  \param[in]  integerGrid	output of simplexGridEnumerate(int n)
     *  \param[in]  n 		the number of linear intervals (aka interpolation order)
     * 
     *  note: the coordinates are just the integer labels of simplexGridEnumerate divided by n
     */
    template <class ct, int mydim>
    static std::vector<Dune::FieldVector<ct, mydim> > simplexGridCoordinateSet(IntegerCoordinateVector integerGrid, int n)
    {
        typedef Dune::FieldVector<ct, mydim>  LocalVector;
    	std::vector<LocalVector > rez;

        for (int i = 0; i < integerGrid.size(); i++)
        {
                LocalVector tmpLocalCoordinate;
                for (int j = 0; j < mydim; j++) { tmpLocalCoordinate[j] = double(integerGrid[i][j]) / n; }
                rez.push_back(tmpLocalCoordinate);
        }
        return rez;
    }


    /** \brief List of internal vertex indices of an element corresponding to the interpolatory vertices of a subentity of the element
     *  \param[in]  entityGeometry	GeometryType of the element of interest
     *  \param[in]  order		Interpolation Order of the element of interest
     *  \param[in]  subentityCodim	Codimension of the subentity of interest
     *  \param[in]  subentityIndex		Subentity internal index inside the element
     */
    template <class ct, int cdim>
    static std::vector<InternalIndexType> subentityInternalCoordinateSet(Dune::GeometryType entityGeometry, int order, int subentityCodim, int subentityIndex)
    {
    	typedef Dune::FieldVector<int, cdim> IntFieldVector;

    	std::vector<InternalIndexType> rez;

    	const Dune::ReferenceElement< ct, cdim > & ref = Dune::ReferenceElements< ct, cdim >::general(entityGeometry);

    	int nSubentity = ref.size(subentityCodim);
    	int nSubCorner = ref.size(0, subentityCodim, cdim);

        if (!entityGeometry.isSimplex())                           { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: subentityInternalCoordinateSet() only implemented for Simplex geometries at the moment"); }
        if ((subentityCodim < 0)||(subentityCodim >= cdim))        { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: subentityInternalCoordinateSet() - Unexpected subentity codimension"); }
        if ((subentityIndex < 0)||(subentityIndex >= nSubentity))  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: subentityInternalCoordinateSet() - Unexpected subentity index"); }


        // Special cases
        // ********************************************************

        // If coordinates of the element itself are requested, just return all coordinates
        if (subentityCodim == 0)  {
        	for (int i = 0; i < nSubCorner; i++ ) { rez.push_back(i); }
        	return rez;
        }

        // If coordinates of a corner are requested, then it is simply 1 corner of the element
        if (subentityCodim == cdim)  { return std::vector<InternalIndexType> (1, cornerIndex(entityGeometry, order, subentityIndex)); }


        // General case
        // ********************************************************

        // Get corner internal indices, and associated internal coordinates
        std::vector<InternalIndexType> cornerInd;
        std::vector<IntFieldVector> cornerInternalCoord;
        for (int i = 0; i < nSubCorner; i++)  {
        	cornerInd.push_back(ref.subEntity(subentityIndex, subentityCodim, i, cdim));
        	cornerInternalCoord.push_back(cornerInternalCoordinate<int, cdim>(entityGeometry, cornerInd[i]));
        }

        // Consider each geometry separately
             if ((cdim == 3) && (subentityCodim == 1) )  { rez = tet2TriIndexSet (order, cornerInternalCoord); }
        else if ((cdim == 3) && (subentityCodim == 2) )  { rez = tet2EdgeIndexSet(order, cornerInternalCoord); }
        else if ((cdim == 2) && (subentityCodim == 1) )  { rez = tri2EdgeIndexSet(order, cornerInternalCoord); }
        else { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: subentityInternalCoordinates() - unexpected dimension-codimension pair"); }

        return rez;
    }



    // Returns corner id's of this entity
    template<class ct, int mydim>
    static std::vector<int> entityVertexCornerSubset(
            Dune::GeometryType gt,
            const std::vector<int> & vertexIndexSet,
            InterpolatoryOrderType order)
    {
        std::vector<int> corner;

        // Get corner number
        int cornerNo = Dune::ReferenceElements<ct, mydim>::general(gt).size(gt.dim());
        //int cornerNo = gt.dim() + 1;  // for simplices

        // Get corners
        for (int j = 0; j < cornerNo; j++) {
            InternalIndexType internalId = Dune::CurvilinearGeometryHelper::cornerIndex(gt, order, j );
            corner.push_back(vertexIndexSet[internalId]);
        }

        return corner;
    }


  private:


    /** \brief Computes internal indices of all interpolatory vertices of a face subentity of a tetrahedron
     *
     *  \param[in]  order                  interpolatory order of the entity
     *  \param[in]  cornerInternalCoord    internal integer coordinates of the corners of the face within the reference tetrahedron
     *
     *  [FIXME] - there should exist a nicer algorithm to do this...
     *
     *  Algorithm:
     *    1) Write all coordinate numbers in a ((order + 1) x (order + 1) x (order + 1)) matrix
     *    2) Construct a parametric face connecting provided subentity corners
     *    3) Collect all points on that line in an array and return it
     * */
    template <int cdim>
    static std::vector<InternalIndexType> tet2TriIndexSet(int order, std::vector<Dune::FieldVector<int, cdim> > & cornerInternalCoord)
    {
    	typedef Dune::FieldVector<int, cdim> IntFieldVector;

        std::vector<InternalIndexType> rez;

        // Build a 3D matrix consisting of coordinate numbers in the point_ vector in the shape of reference simplex
        int coord_map [order + 1][order + 1][order + 1];
        IntegerCoordinateVector simplexGrid = simplexGridEnumerate<3>(order);
        for (int i = 0; i < simplexGrid.size(); i++) { coord_map[simplexGrid[i][0]][simplexGrid[i][1]][simplexGrid[i][2]] = i; }

        // Find the direction vectors
        IntFieldVector v1 = cornerInternalCoord[1] - cornerInternalCoord[0];
        IntFieldVector v2 = cornerInternalCoord[2] - cornerInternalCoord[0];


        for (int i = 0; i <= order; i++)
        {
            for (int j = 0; j <= order - i; j++)
            {
            	IntFieldVector mapIndex;
            	mapIndex.axpy(order, cornerInternalCoord[0]);
            	mapIndex.axpy(i, v1);
            	mapIndex.axpy(j, v2);
            	rez.push_back(coord_map[mapIndex[0]][mapIndex[1]][mapIndex[2]]);
            }
        }
        return rez;
    }


    /** \brief Computes internal indices of all interpolatory vertices of an edge subentity of a tetrahedron
     *
     *  \param[in]  order                  interpolatory order of the entity
     *  \param[in]  cornerInternalCoord    internal integer coordinates of the corners of the edge within the reference tetrahedron
     *
     *  \note Can't figure out a nice way to do this.
     *
     *  Algorithm:
     *    1) Write all coordinate numbers in a ((order + 1) x (order + 1) x (order + 1)) matrix
     *    2) Construct a parametric line connecting provided subentity corners
     *    3) Collect all points on that line in an array and return it
     * */
    template <int cdim>
    static std::vector<InternalIndexType> tet2EdgeIndexSet(int order, std::vector<Dune::FieldVector<int, cdim> > & cornerInternalCoord)
    {
    	typedef Dune::FieldVector<int, cdim> IntFieldVector;

        std::vector<InternalIndexType> rez;

        // Build a 3D matrix consisting of coordinate numbers in the point_ vector in the shape of reference simplex
        int coord_map [order + 1][order + 1][order + 1];
        IntegerCoordinateVector simplexGrid = simplexGridEnumerate<3>(order);
        for (int i = 0; i < simplexGrid.size(); i++) { coord_map[simplexGrid[i][0]][simplexGrid[i][1]][simplexGrid[i][2]] = i; }

        // Find the direction vectors
        IntFieldVector v = cornerInternalCoord[1] - cornerInternalCoord[0];

        for (int i = 0; i <= order; i++)
        {
        	IntFieldVector mapIndex;
        	mapIndex.axpy(order, cornerInternalCoord[0]);
        	mapIndex.axpy(i, v);
        	rez.push_back(coord_map[mapIndex[0]][mapIndex[1]][mapIndex[2]]);
        }
        return rez;
    }


    /** \brief Computes internal indices of all interpolatory vertices of an edge subentity of a triangle
     *
     *  \param[in]  order                  interpolatory order of the entity
     *  \param[in]  cornerInternalCoord    internal integer coordinates of the corners of the edge within the reference triangle
     *
     *  \note Can't figure out a nice way to do this.
     *
     *  Algorithm:
     *    1) Write all coordinate numbers in a ((order + 1) x (order + 1)) matrix
     *    2) Construct a parametric line connecting provided subentity corners
     *    3) Collect all points on that line in an array and return it
     * */
    template <int cdim>
    static std::vector<InternalIndexType> tri2EdgeIndexSet(int order, std::vector<Dune::FieldVector<int, cdim> > & cornerInternalCoord)
    {
    	typedef Dune::FieldVector<int, cdim> IntFieldVector;

        std::vector<InternalIndexType> rez;

        // Build a 2D matrix consisting of coordinate numbers in the point_ vector in the shape of reference simplex
        int coord_map [order + 1][order + 1];
        IntegerCoordinateVector simplexGrid = simplexGridEnumerate<2>(order);
        for (int i = 0; i < simplexGrid.size(); i++) { coord_map[simplexGrid[i][0]][simplexGrid[i][1]] = i; }

        // Find the direction vector
        IntFieldVector v = cornerInternalCoord[1] - cornerInternalCoord[0];


        for (int i = 0; i <= order; i++)
        {
        	IntFieldVector mapIndex;
        	mapIndex.axpy(order, cornerInternalCoord[0]);
        	mapIndex.axpy(i, v);
        	rez.push_back(coord_map[mapIndex[0]][mapIndex[1]]);
        }
        return rez;
    }

}; // SimplexInterpolator




} // Namespace Dune

#endif /** DUNE_CURVILINEARGEOMETRYHELPER_HH **/
