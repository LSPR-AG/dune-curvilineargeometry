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
  protected:
    typedef std::vector<std::vector<int>> IntegerCoordinateVector;

  public:


    // Public typedefs
    typedef  int     InternalIndexType;
    typedef  int     InterpolatoryOrderType;




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
    static InternalIndexType cornerID(Dune::GeometryType geomType, int order, InternalIndexType i)
    {
        if (!geomType.isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: cornerID() only implemented for Simplex geometries at the moment"); }

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


    // Returns d-1 subentity corner local indices, sorted
    static std::vector<InternalIndexType> linearElementSubentityCornerInternalIndexSet(GeometryType gt, int subentityCodim, InternalIndexType subentityIndex)
    {
        std::vector<InternalIndexType> rez;

        if (gt.isTriangle() && (subentityCodim == 1))          { rez = linearTriangleSubentityEdgeIndexSet(subentityIndex);  }
        else if (gt.isTetrahedron() && (subentityCodim == 0))  { rez = std::vector<InternalIndexType> {0, 1, 2, 3};          }
        else if (gt.isTetrahedron() && (subentityCodim == 1))  { rez = linearTetrahedronSubentityTriangleIndexSet(subentityIndex);  }
        else if (gt.isTetrahedron() && (subentityCodim == 2))  { rez = linearTetrahedronSubentityEdgeIndexSet(subentityIndex);  }
        else  {  DUNE_THROW(Dune::IOError, "Curvilinear Helper: Not implemented element subentityIndex for this element type " ); }

        return rez;
    }

    /** \brief List of internal vertex indices of an element corresponding to the interpolatory vertices of a subentity of the element
     *  \param[in]  entityGeometry	GeometryType of the element of interest
     *  \param[in]  order		Interpolation Order of the element of interest
     *  \param[in]  subentityCodim	Codimension of the subentity of interest
     *  \param[in]  subentityIndex		Subentity internal index inside the element
     * 
     *  FIXME: The Subentity internal index is not in sync with the Dune convention
     */
    template <class ct, int mydim>
    static std::vector<InternalIndexType> subentityInternalCoordinateSet(Dune::GeometryType entityGeometry, int order, int subentityCodim, int subentityIndex)
    {
    	int nSubentity = Dune::ReferenceElements< ct, mydim >::general(entityGeometry).size(subentityCodim);


        if (!entityGeometry.isSimplex())  { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: subentityInternalCoordinates() only implemented for Simplex geometries at the moment"); }
        if ((subentityIndex < 0)||(subentityIndex >= nSubentity))
        {
        	DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: SubentityInterpolator() - Unexpected subentity index");
        }

        std::vector<InternalIndexType> rez;

             if ((entityGeometry.dim() == 3) && (subentityCodim == 1) ) {  rez = tetrahedronSubentityTriangleIndexSet(order, subentityIndex); }
        else if ((entityGeometry.dim() == 3) && (subentityCodim == 2) ) {  rez = tetrahedronSubentityEdgeIndexSet(order, subentityIndex); }
        else if ((entityGeometry.dim() == 2) && (subentityCodim == 1) ) {  rez = triangleSubentityEdgeIndexSet(order, subentityIndex); }
        else { DUNE_THROW(Dune::IOError, "CURVILINEAR_ELEMENT_INTERPOLATOR: subentityInternalCoordinates() - unexpected dimension-codimension pair"); }

        return rez;
    }



    // Returns corner id's of this entity
    static std::vector<int> entityVertexCornerSubset(
            Dune::GeometryType gt,
            const std::vector<int> & vertexIndexSet,
            InterpolatoryOrderType order) const
    {
        std::vector<int> corner;

        // Get corner number
        int cornerNo = Dune::ReferenceElements::general(gt).size(gt.dim());
        //int cornerNo = gt.dim() + 1;  // for simplices

        // Get corners
        for (int j = 0; j < cornerNo; j++) {
            InternalIndexType internalId = Dune::CurvilinearGeometryHelper::cornerID(gt, order, j );
            corner.push_back(vertexIndexSet[internalId]);
        }

        return corner;
    }


  private:

    /** \brief  internal corner vertex index subset corresponding to subentity corners of a linear element
     *  \param[in]  subentityIndex		Subentity internal index inside the element
     *
     *  FIXME - the orientation convention has not been sync with Dune
     */
    static std::vector<InternalIndexType> linearTriangleSubentityEdgeIndexSet(InternalIndexType subentityIndex)
    {
        std::vector<InternalIndexType> rez;

        switch (subentityIndex)
        {
        case 0 :  rez = std::vector<InternalIndexType> {0, 1};  break;
        case 1 :  rez = std::vector<InternalIndexType> {1, 2};  break;
        case 2 :  rez = std::vector<InternalIndexType> {0, 2};  break;
        default : DUNE_THROW(Dune::IOError, "Curvilinear Helper: Wrong input arguments for SubentityCorners " );
        }

        return rez;
    }

    static std::vector<InternalIndexType> linearTetrahedronSubentityEdgeIndexSet(InternalIndexType subentityIndex)
    {
        std::vector<InternalIndexType> rez;

        switch (subentityIndex)
        {
        case 0 :  rez = std::vector<InternalIndexType> {0, 1};  break;
        case 1 :  rez = std::vector<InternalIndexType> {1, 2};  break;
        case 2 :  rez = std::vector<InternalIndexType> {2, 0};  break;
        case 3 :  rez = std::vector<InternalIndexType> {3, 0};  break;
        case 4 :  rez = std::vector<InternalIndexType> {3, 2};  break;
        case 5 :  rez = std::vector<InternalIndexType> {3, 1};  break;
        default : DUNE_THROW(Dune::IOError, "Curvilinear Helper: Wrong input arguments for SubentityCorners " );
        }

        return rez;
    }

    static std::vector<InternalIndexType> linearTetrahedronSubentityTriangleIndexSet(InternalIndexType subentityIndex)
    {
        std::vector<InternalIndexType> rez;

        switch (subentityIndex)
        {
        case 0 :  rez = std::vector<InternalIndexType> {0, 1, 2};  break;
        case 1 :  rez = std::vector<InternalIndexType> {0, 2, 3};  break;
        case 2 :  rez = std::vector<InternalIndexType> {2, 1, 3};  break;
        case 3 :  rez = std::vector<InternalIndexType> {0, 3, 1};  break;
        default : DUNE_THROW(Dune::IOError, "Curvilinear Helper: Wrong input arguments for SubentityCorners " );
        }

        return rez;
    }


    /** \brief  internal interpolatory vertex index subset corresponding to the subentity interpolatory vertices
     *  \param[in]  gt			GeometryType of the subentity
     *  \param[in]  order		Interpolation Order of the element of interest
     *  \param[in]  subentityIndex		Subentity internal index inside the element
     * 
     *  FIXME - there should exist a nicer algorithm to do this...
     *
     *  FIXME - the orientation convention has not been sync with Dune
     */
    static std::vector<InternalIndexType> tetrahedronSubentityTriangleIndexSet(int order, InternalIndexType subentityIndex)
    {
    	std::vector<InternalIndexType> rez;
        // Can't figure out a nice way to do this
        // For now will write all coordinate numbers in a ((order + 1) x (order + 1) x (order + 1)) matrix and map from there

        // Build a 3D matrix consisting of coordinate numbers in the point_ vector in the shape of reference simplex
        int coord_map [order + 1][order + 1][order + 1];
        IntegerCoordinateVector simplexGrid = simplexGridEnumerate<3>(order);
        for (int i = 0; i < simplexGrid.size(); i++) { coord_map[simplexGrid[i][0]][simplexGrid[i][1]][simplexGrid[i][2]] = i; }

        std::cout << "subFaceReqInd=" << subentityIndex << std::endl;

        for (int i = 0; i <= order; i++)
        {
            for (int j = 0; j <= order - i; j++)
            {
                switch (subentityIndex)
                {
                    case 0 : rez.push_back(coord_map[j][i][0]);               break;  // Face (0 1 2)
                    case 1 : rez.push_back(coord_map[0][j][i]);               break;  // Face (0 2 3)
                    case 2 : rez.push_back(coord_map[j][order - i - j][i]);   break;  // Face (2 1 3)
                    case 3 : rez.push_back(coord_map[i][0][j]);               break;  // Face (0 3 1)
                }
            }
        }
        return rez;
    }

    static std::vector<InternalIndexType> tetrahedronSubentityEdgeIndexSet(int order, InternalIndexType subentityIndex)
    {
        std::vector<InternalIndexType> rez;
        // Can't figure out a nice way to do this
        // For now will write all coordinate numbers in a ((order + 1) x (order + 1) x (order + 1)) matrix and map from there

        // Build a 3D matrix consisting of coordinate numbers in the point_ vector in the shape of reference simplex
        int coord_map [order + 1][order + 1][order + 1];
        IntegerCoordinateVector simplexGrid = simplexGridEnumerate<3>(order);
        for (int i = 0; i < simplexGrid.size(); i++) { coord_map[simplexGrid[i][0]][simplexGrid[i][1]][simplexGrid[i][2]] = i; }

        for (int i = 0; i <= order; i++)
        {
            switch (subentityIndex)
            {
                case 0 : rez.push_back(coord_map[i][0][0]);          break;  // Edge (0 1)
                case 1 : rez.push_back(coord_map[order - i][i][0]);  break;  // Edge (1 2)
                case 2 : rez.push_back(coord_map[0][order - i][0]);  break;  // Edge (2 0)
                case 3 : rez.push_back(coord_map[0][0][order - i]);  break;  // Edge (3 0)
                case 4 : rez.push_back(coord_map[0][i][order - i]);  break;  // Edge (3 2)
                case 5 : rez.push_back(coord_map[i][0][order - i]);  break;  // Edge (3 1)
            }
        }
        return rez;
    }

    static std::vector<InternalIndexType> triangleSubentityEdgeIndexSet(int order, InternalIndexType subentityIndex)
    {
        std::vector<InternalIndexType> rez;

        // Can't figure out a nice way to do this
        // For now will write all coordinate numbers in a ((order + 1) x (order + 1)) matrix and map from there

        // Build a 2D matrix consisting of coordinate numbers in the point_ vector in the shape of reference simplex
        int coord_map [order + 1][order + 1];
        IntegerCoordinateVector simplexGrid = simplexGridEnumerate<2>(order);
        for (int i = 0; i < simplexGrid.size(); i++) { coord_map[simplexGrid[i][0]][simplexGrid[i][1]] = i; }

        for (int i = 0; i <= order; i++)
        {
            switch(subentityIndex)
            {
                case 0: rez.push_back(coord_map[i][0]);          break;  // Edge (0 1)
                case 1: rez.push_back(coord_map[order - i][i]);  break;  // Edge (1 2)
                case 2: rez.push_back(coord_map[0][order - i]);  break;  // Edge (2 0)
            }
        }
        return rez;
    }

}; // SimplexInterpolator




} // Namespace Dune

#endif /** DUNE_CURVILINEARGEOMETRYHELPER_HH **/
