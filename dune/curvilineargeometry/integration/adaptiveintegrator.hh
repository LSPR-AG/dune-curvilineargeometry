/*******************************************************************
 * Numerical Recursive Interpolation Integrator
 * 
 * author: Aleksejs Fomins
 * date: 01.08.2014 - created
 * 
 * description:
 * Integrates arbitrary function over reference simplex (edge or triangle) using piecewise interpolation.
 * - Function is given in terms of a functor object. It need not be continuous
 * - The algorithm uses lagrange interpolation to sample the function at the interpolation points for 2nd and 4th order
 * - The integral of the real function over a simplex is approximated by an analytical integral of the interpolated function over this simplex
 * - If the total expected error is too high, the 4th order simplex with highest expected error is split into 2nd order simplices, which are then refined to 4th order and dealt with recursively. The algorithm re-uses previous sampling points when splitting
 * - The running error of a simplex is given by the difference between the approximated integrals for orders 2 and 4. All sub-simplices are always sorted by expected integration error regardless of their size
 * 
 * 
 * TODO: The algorithm is rather slow. It is possible that a higher interpolation order may be beneficial.
 * 
 *******************************************************************/


#ifndef DUNE_ADAPTIVE_INTEGRATOR_HH
#define DUNE_ADAPTIVE_INTEGRATOR_HH


#include <iostream>
#include <vector>
#include <queue>          // std::queue
#include <algorithm>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>
#include <config.h>

#include <dune/common/function.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/interpolation/lagrangeinterpolator.hh>



namespace Dune {


template<class ct, int mydim, int order1, int order2>
class ElementTwoOrders {

    typedef FieldVector< ct, mydim > GlobalVector;
    typedef FieldVector< ct, mydim + 1 > CompositeVector;
    typedef std::vector< GlobalVector > GlobalVectorVector;
    typedef std::vector< CompositeVector > CompositeVectorVector;
    typedef std::vector<ct> sampleVector;
    typedef Polynomial<ct, mydim> LocalPolynomial;
    typedef std::vector<LocalPolynomial> PolynomialVector;
    typedef Function<GlobalVector, ct> GlobalFunction;
    typedef LagrangeInterpolator<double, mydim, mydim + 1> CompositeInterpolator;

private:
    GlobalVectorVector vertexSet_;
    std::vector<double> value_;
    double integralO1_;
    double integralO2_;
    double error_;

public:
    ElementTwoOrders() {}

    ElementTwoOrders(
              Dune::GeometryType gt,
              const GlobalVectorVector & verticesO1, const std::vector<double> & valuesO1,
              const GlobalVectorVector & verticesO2, const std::vector<double> & valuesO2
            ) {
        vertexSet_ = verticesO2;
        value_ = valuesO2;

        // Merge vertices and values to one single vector of vertices of dimension cdim+1
        CompositeVectorVector composite_verticesO1 = compositeVertices(verticesO1, valuesO1);
        CompositeVectorVector composite_verticesO2 = compositeVertices(verticesO2, valuesO2);

        CompositeInterpolator interp_composite_O1(gt, composite_verticesO1, order1);
        CompositeInterpolator interp_composite_O2(gt, composite_verticesO2, order2);

        PolynomialVector intPoly_composite_O1 = interp_composite_O1.interpolatoryVectorAnalytical();
        PolynomialVector intPoly_composite_O2 = interp_composite_O2.interpolatoryVectorAnalytical();

        PolynomialVector intPoly_global_O1 = intPoly_composite_O1;      intPoly_global_O1.erase(intPoly_global_O1.end());
        PolynomialVector intPoly_global_O2 = intPoly_composite_O2;      intPoly_global_O2.erase(intPoly_global_O2.end());

        LocalPolynomial J_O1 = JacobianDeterminantAnalytical(intPoly_global_O1);
        LocalPolynomial J_O2 = JacobianDeterminantAnalytical(intPoly_global_O2);

        integralO1_    = (intPoly_composite_O1[mydim] * J_O1).integrateRefSimplex();
        integralO2_ = (intPoly_composite_O2[mydim] * J_O2).integrateRefSimplex();

        error_ = fabs(integralO1_ - integralO2_);
    }

    // Returns vertex of the interpolatory grid
    GlobalVector vertex(int i) const{return vertexSet_[i]; }

    // Returns value of the integrand at the interpolatory vertex
    double value(int i) const {return value_[i]; }

    // Returns analytical integral over element using lower order interpolation
    double integralLowerOrder() const { return integralO1_; }

    // Returns analytical integral over element using higher order interpolation
    double integralHigherOrder() const { return integralO2_; }

    // Returns the estimated integration error associated with this element
    double error() const { return error_; }



private:
    CompositeVectorVector compositeVertices(const GlobalVectorVector & vertices, const std::vector<double> & values) const
    {
        CompositeVectorVector rez;
        for (int i = 0; i < vertices.size(); i++)
        {
            CompositeVector tmp;
            for (int j = 0; j < mydim; j++) { tmp[j] = vertices[i][j]; }
            tmp[mydim] = values[i];
            rez.push_back(tmp);
        }
        return rez;
    }

    LocalPolynomial JacobianDeterminantAnalytical(const PolynomialVector & analyticalMap) const
    {
        LocalPolynomial rez;
        GlobalVector mid;
        switch(mydim)
        {
        case 1:
            rez = analyticalMap[0].derivative(0);
            mid[0] = 0.1;
            break;
        case 2:
            rez = analyticalMap[0].derivative(0) * analyticalMap[1].derivative(1) - analyticalMap[0].derivative(1) * analyticalMap[1].derivative(0);
            mid[0] = 0.1;
            mid[1] = 0.1;
            break;
        case 3:
            rez += analyticalMap[0].derivative(0) * ( analyticalMap[1].derivative(1) * analyticalMap[2].derivative(2) - analyticalMap[1].derivative(2) * analyticalMap[2].derivative(1) );
            rez += analyticalMap[0].derivative(1) * ( analyticalMap[1].derivative(2) * analyticalMap[2].derivative(0) - analyticalMap[1].derivative(0) * analyticalMap[2].derivative(2) );
            rez += analyticalMap[0].derivative(2) * ( analyticalMap[1].derivative(0) * analyticalMap[2].derivative(1) - analyticalMap[1].derivative(1) * analyticalMap[2].derivative(0) );
            mid[0] = 0.1;
            mid[1] = 0.1;
            mid[2] = 0.1;
            break;
        }
        // Change sign if determinant is negative
        if (rez.evaluate(mid) < 0) { rez *= -1; }

        return rez;
    }
};



template<class ct, int mydim>
class AdaptiveIntegrator {
    typedef FieldVector< ct, mydim > GlobalVector;
    typedef std::vector< GlobalVector > GlobalVectorVector;

    typedef ElementTwoOrders<ct, mydim, 2, 4> Element;
    typedef std::vector<Element> ElementVector;

    typedef Function<GlobalVector, ct> GlobalFunction;

    typedef Dune::ReferenceElement< ct, mydim > ReferenceElement;
    typedef Dune::ReferenceElements< ct, mydim > ReferenceElements;

protected:

  GlobalVector mid(const GlobalVector & a, const GlobalVector & b) const
  {
      GlobalVector rez;
      for (int i = 0; i < mydim; i++) { rez[i] = (a[i] + b[i]) / 2.0; }
      return rez;
  }

  template<typename Functor>
  ElementVector splitRefine(const Element & elem, const Functor & f) const {
      ElementVector rez;

      switch (mydim)
      {
        case 1 : rez = splitRefineEdgeOrder4(elem, f);      break;
        case 2 : rez = splitRefineTriangleOrder4(elem, f);  break;
        //case 3 : rez = splitRefineTetrahedron();            break;
      }
      return rez;
  }

  template<typename Functor>
  ElementVector splitRefineEdgeOrder2(const Element & edge, const Functor & f) const
  {
      ElementVector rez;

      int indices[2][2] = { {0,1}, {1,2} };

      for (int i = 0; i < 2; i++)
      {
      GlobalVectorVector edge_vertices;
      std::vector<double> edge_values;

      GlobalVectorVector edge_refined_vertices;
          std::vector<double> edge_refined_values;

          edge_vertices.push_back(edge.vertex(indices[i][0]));  edge_values.push_back(edge.value(indices[i][0]));
          edge_vertices.push_back(edge.vertex(indices[i][1]));  edge_values.push_back(edge.value(indices[i][1]));

          edge_refined_vertices.push_back(edge_vertices[0]);
          edge_refined_vertices.push_back(mid(edge_vertices[0], edge_vertices[1]));
          edge_refined_vertices.push_back(edge_vertices[1]);

          edge_refined_values.push_back(edge_values[0]);
          edge_refined_values.push_back(f(edge_refined_vertices[1]));
          edge_refined_values.push_back(edge_values[1]);

          rez.push_back(Element(type(), edge_vertices, edge_values, edge_refined_vertices, edge_refined_values));
      }

      return rez;
  }

  template<typename Functor>
  ElementVector splitRefineEdgeOrder4(const Element & edge, const Functor & f) const
  {
      ElementVector rez;

      int indices[2][3] = { {0,1,2}, {2,3,4} };

      for (int i = 0; i < 2; i++)
      {
      GlobalVectorVector edge_vertices;
      std::vector<double> edge_values;

      GlobalVectorVector edge_refined_vertices;
          std::vector<double> edge_refined_values;

          edge_vertices.push_back(edge.vertex(indices[i][0]));  edge_values.push_back(edge.value(indices[i][0]));
          edge_vertices.push_back(edge.vertex(indices[i][1]));  edge_values.push_back(edge.value(indices[i][1]));
          edge_vertices.push_back(edge.vertex(indices[i][2]));  edge_values.push_back(edge.value(indices[i][2]));

          edge_refined_vertices.push_back(edge_vertices[0]);
          edge_refined_vertices.push_back(mid(edge_vertices[0], edge_vertices[1]));
          edge_refined_vertices.push_back(edge_vertices[1]);
          edge_refined_vertices.push_back(mid(edge_vertices[1], edge_vertices[2]));
          edge_refined_vertices.push_back(edge_vertices[2]);

          edge_refined_values.push_back(edge_values[0]);
          edge_refined_values.push_back(f(edge_refined_vertices[1]));
          edge_refined_values.push_back(edge_values[1]);
          edge_refined_values.push_back(f(edge_refined_vertices[3]));
          edge_refined_values.push_back(edge_values[2]);

          rez.push_back(Element(type(), edge_vertices, edge_values, edge_refined_vertices, edge_refined_values));
      }

      return rez;
  }

  template<typename Functor>
  ElementVector splitRefineTriangleOrder2(const Element & tri, const Functor & f) const
  {
      ElementVector rez;

      int indices[4][3] = {{0,1,3}, {1,2,4}, {4,3,1}, {3,4,5} };

      for (int i = 0; i < 4; i++)
      {
          GlobalVectorVector tri_vertices;
          std::vector<double> tri_values;

          GlobalVectorVector tri_refined_vertices;
          std::vector<double> tri_refined_values;

          tri_vertices.push_back(tri.vertex(indices[i][0]));  tri_values.push_back(tri.value(indices[i][0]));
          tri_vertices.push_back(tri.vertex(indices[i][1]));  tri_values.push_back(tri.value(indices[i][1]));
          tri_vertices.push_back(tri.vertex(indices[i][2]));  tri_values.push_back(tri.value(indices[i][2]));


          tri_refined_vertices.push_back(tri_vertices[0]);
          tri_refined_vertices.push_back(mid(tri_vertices[0], tri_vertices[1]));
          tri_refined_vertices.push_back(tri_vertices[1]);
          tri_refined_vertices.push_back(mid(tri_vertices[0], tri_vertices[2]));
          tri_refined_vertices.push_back(mid(tri_vertices[1], tri_vertices[2]));
          tri_refined_vertices.push_back(tri_vertices[2]);

          tri_refined_values.push_back(tri_values[0]);
          tri_refined_values.push_back(f(tri_refined_vertices[1]));
          tri_refined_values.push_back(tri_values[1]);
          tri_refined_values.push_back(f(tri_refined_vertices[3]));
          tri_refined_values.push_back(f(tri_refined_vertices[4]));
          tri_refined_values.push_back(tri_values[2]);

          rez.push_back(Element(type(), tri_vertices, tri_values, tri_refined_vertices, tri_refined_values));
      }

      return rez;
  }

  template<typename Functor>
  ElementVector splitRefineTriangleOrder4(const Element & tri, const Functor & f) const
  {
      ElementVector rez;

      int indices[4][6] = {{0,1,2,5,6,9}, {2,3,4,7,8,11}, {11,10,9,7,6,2}, {9,10,11,12,13,14} };

      for (int i = 0; i < 4; i++)
      {
          GlobalVectorVector tri_vertices;
          std::vector<double> tri_values;

          GlobalVectorVector tri_refined_vertices;
          std::vector<double> tri_refined_values;

          tri_vertices.push_back(tri.vertex(indices[i][0]));  tri_values.push_back(tri.value(indices[i][0]));
          tri_vertices.push_back(tri.vertex(indices[i][1]));  tri_values.push_back(tri.value(indices[i][1]));
          tri_vertices.push_back(tri.vertex(indices[i][2]));  tri_values.push_back(tri.value(indices[i][2]));
          tri_vertices.push_back(tri.vertex(indices[i][3]));  tri_values.push_back(tri.value(indices[i][3]));
          tri_vertices.push_back(tri.vertex(indices[i][4]));  tri_values.push_back(tri.value(indices[i][4]));
          tri_vertices.push_back(tri.vertex(indices[i][5]));  tri_values.push_back(tri.value(indices[i][5]));


          tri_refined_vertices.push_back(tri_vertices[0]);
          tri_refined_vertices.push_back(mid(tri_vertices[0], tri_vertices[1]));
          tri_refined_vertices.push_back(tri_vertices[1]);
          tri_refined_vertices.push_back(mid(tri_vertices[1], tri_vertices[2]));
          tri_refined_vertices.push_back(tri_vertices[2]);
          tri_refined_vertices.push_back(mid(tri_vertices[0], tri_vertices[3]));
          tri_refined_vertices.push_back(mid(tri_vertices[0], tri_vertices[4]));
          tri_refined_vertices.push_back(mid(tri_vertices[1], tri_vertices[4]));
          tri_refined_vertices.push_back(mid(tri_vertices[2], tri_vertices[4]));
          tri_refined_vertices.push_back(tri_vertices[3]);
          tri_refined_vertices.push_back(mid(tri_vertices[3], tri_vertices[4]));
          tri_refined_vertices.push_back(tri_vertices[4]);
          tri_refined_vertices.push_back(mid(tri_vertices[3], tri_vertices[5]));
          tri_refined_vertices.push_back(mid(tri_vertices[4], tri_vertices[5]));
          tri_refined_vertices.push_back(tri_vertices[5]);

          tri_refined_values.push_back(tri_values[0]);
          tri_refined_values.push_back(f(tri_refined_vertices[1]));
          tri_refined_values.push_back(tri_values[1]);
          tri_refined_values.push_back(f(tri_refined_vertices[3]));
          tri_refined_values.push_back(tri_values[2]);
          tri_refined_values.push_back(f(tri_refined_vertices[5]));
          tri_refined_values.push_back(f(tri_refined_vertices[6]));
          tri_refined_values.push_back(f(tri_refined_vertices[7]));
          tri_refined_values.push_back(f(tri_refined_vertices[8]));
          tri_refined_values.push_back(tri_values[3]);
          tri_refined_values.push_back(f(tri_refined_vertices[10]));
          tri_refined_values.push_back(tri_values[4]);
          tri_refined_values.push_back(f(tri_refined_vertices[12]));
          tri_refined_values.push_back(f(tri_refined_vertices[13]));
          tri_refined_values.push_back(tri_values[5]);

          rez.push_back(Element(type(), tri_vertices, tri_values, tri_refined_vertices, tri_refined_values));
      }

      return rez;
  }

  template<typename Functor>
  Element referenceElementRefined(int order1, int order2, const Functor & f) const {
      GlobalVectorVector vertices_O1 = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<ct, mydim>(order1);
      GlobalVectorVector vertices_O2 = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<ct, mydim>(order2);

      std::vector<double> values_O1;
      std::vector<double> values_O2;

      for (int i = 0; i < vertices_O1.size(); i++) { values_O1.push_back(f(vertices_O1[i])); }
      for (int i = 0; i < vertices_O2.size(); i++) { values_O2.push_back(f(vertices_O2[i])); }

      return Element(type(), vertices_O1, values_O1, vertices_O2, values_O2);
  }


  // Sinks last edge until it encounters edge with error smaller than own, or until it hits the bottom
  void heapSink (ElementVector & elementHeap)
  {
      int i = elementHeap.size() - 1;
      while( (i > 0)&&( elementHeap[i].error() < elementHeap[i - 1].error())   ) { std::swap(elementHeap[i], elementHeap[i - 1]);  i--; }
  }

  // Decides when integration has finished
  bool still_integrating(double rezIntegral, double totalError, double tolerance, int sampleNo) const
  {
      double expected_err = fabs(rezIntegral * tolerance);

      // First check if the expected absolute precision of the error is within reasonable limits
      // If not, because of too small integral, then it is almost 0 anyway.
      // If not, because of too small tolerance, then this is a computationally infeasible request.
      bool feasible_expected_error = ( expected_err > 1e-20 );

      // Check if relative error within tolerance (aka finish if err/val <= tolerance)
      bool feasible_error = (totalError > expected_err);

      // For now, also limit maximum number of iterations such that the computation time does not explode
      bool feasible_sample_number = (sampleNo < 100000);

      return (feasible_expected_error && feasible_error && feasible_sample_number );
  }

public:

  AdaptiveIntegrator(Dune::GeometryType gt) : refElement_( &ReferenceElements::general( gt ) ) {  }

  /** \brief obtain the name of the reference element */
  Dune::GeometryType type () const { return refElement_->type(); }

  // Integrates the provided function over reference tetrahedron, returns result
  template<typename Functor>
  double integrate( const Functor & f, double tolerance) {
    Functor integrand_ = f;
    ElementVector elementHeap;
    elementHeap.push_back(referenceElementRefined(2,4,f));

    double rezIntegral = elementHeap[0].integralHigherOrder();
    double totalError = elementHeap[0].error();

    //int iterPoints = (mydim + 1) * (mydim + 2) / 2;
    int iterPoints = (2 * mydim + 1) * (mydim * mydim + mydim + 3) / 3;


    //std::cout << "initial error estimate: " << totalError << std::endl;
    std::cout << "---------------------------------- started integrating element (dim " << mydim << ") ---------------------------------" << std::endl;

    while (still_integrating(rezIntegral, totalError, tolerance, iterPoints))
    {
      Element thisElement = elementHeap[elementHeap.size() - 1];
      elementHeap.erase(elementHeap.end());

      ElementVector newElements = splitRefine(thisElement, f);
      //iterPoints += mydim * (mydim + 1) * (mydim + 2) * (mydim + 7) / 24;
      iterPoints += 4 * (2 * mydim - 1) * (2 * mydim - 1);

      double this_err = 0;
      double this_integral = 0;
      for (int i = 0; i < newElements.size(); i++) {
          this_integral += newElements[i].integralHigherOrder();
          this_err += newElements[i].error();
          elementHeap.push_back(newElements[i]);
          heapSink(elementHeap);

          //std::cout << "**** Integral-linear: " << newElements[i].integralLowerOrder() << ", Integral-quadratic: "<< newElements[i].integralHigherOrder() << ", error: " << newElements[i].error() << std::endl;
      }

      rezIntegral -= thisElement.integralHigherOrder();
      rezIntegral += this_integral;

      totalError -= thisElement.error();
      totalError += this_err;

/*      if (iterPoints < 200)
      {
          std::cout << "******* Printing element heap: ";
          for (int l = 0; l < elementHeap.size(); l++) {
              std::cout << l << " " << elementHeap[l].error() << ": ";
              for (int r = 0; r < mydim + 1; r++) { std::cout << elementHeap[l].vertex(r) << "; ";  }
              std::cout << std::endl;
          }
          std::cout << std::endl;
      }*/


      double lengthscale = fabs(thisElement.vertex(0)[0] - thisElement.vertex(2)[0]);

     // std::cout << "points used: " << iterPoints << ", tmp result: " << rezIntegral << ", Queue Size: " << elementHeap.size();
     // std::cout << ", lengthscale: " << lengthscale << ", expected total error: " << totalError << ", previous local error: "<< thisElement.error() <<", expected local error " << this_err << std::endl;
    }
    std::cout << "---------------------------------- done integrating element ---------------------------------" << std::endl << std::endl;

    return rezIntegral;

  }


protected:
  const ReferenceElement *refElement_;


};

} // Namespace Dune

#endif /** DUNE_ADAPTIVE_INTEGRATOR_HH **/
