/*******************************************************************
 * Test of Numerical Recursive Interpolation Integrator
 * 
 * author: Aleksejs Fomins
 * date: 01.08.2014 - created
 * 
 * description:
 * Integrates a set of functions in 1D and 2D. Compares result to analytic or known numeric result
 * 
 * TODO: Record integration time
 * TODO: Implement automatic testing by evaluating analytic results, thus providing tests of the form PASS/FAIL
 * 
 *******************************************************************/

#include <config.h>
#include<iostream>
#include <math.h>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/curvilineargeometry/integration/adaptiveintegrator.hh>


using namespace Dune;

typedef FieldVector< double, 1 > GlobalVector1D;
typedef FieldVector< double, 2 > GlobalVector2D;
typedef FieldVector< double, 3 > GlobalVector3D;

typedef Function<GlobalVector1D, double> GlobalFunction1D;
typedef Function<GlobalVector2D, double> GlobalFunction2D;
typedef Function<GlobalVector3D, double> GlobalFunction3D;


struct function1d1 {   double operator()(const GlobalVector1D & in) const  { return 1; }  };
struct function1d2 {   double operator()(const GlobalVector1D & in) const  { return in[0]; }  };
struct function1d3 {   double operator()(const GlobalVector1D & in) const  { return in[0] * in[0] * in[0] - 3 * in[0] + 3; }  };
struct function1d4 {   double operator()(const GlobalVector1D & in) const  { return sqrt(in[0]); }  };

struct function2d1 {   double operator()(const GlobalVector2D & in) const  { return 1; }  };
struct function2d2 {   double operator()(const GlobalVector2D & in) const  { return 1 + in[0]; }  };
struct function2d3 {   double operator()(const GlobalVector2D & in) const  { return 1 + in[0] * in[0] + in[1] * in[1]; }  };
struct function2d4 {   double operator()(const GlobalVector2D & in) const  { return in[0] * in[1] * in[1]; }  };
struct function2d5 {   double operator()(const GlobalVector2D & in) const  { return sqrt(in[0] * in[1]); }  };
struct function2d6 {   double operator()(const GlobalVector2D & in) const  { return 2000 * in[0] * in[0] * in[0] * in[1] * in[1] * in[1]; }  };
struct function2d7 {   double operator()(const GlobalVector2D & in) const  { return 3628800 * pow(in[0], 7) * pow(in[1], 10); }  };
struct function2d8 {   double operator()(const GlobalVector2D & in) const  { return sqrt(pow(in[0], 7) * pow(in[1], 10) + 0.5); }  };


int main ()
{
  std::cout << "initialized" << std::endl;

  Dune::GeometryType edgeGeometry=Dune::GeometryTypes::simplex(1);
  Dune::GeometryType faceGeometry=Dune::GeometryTypes::simplex(2);


  AdaptiveIntegrator<double, 1> funIntegrator1D(edgeGeometry);
  AdaptiveIntegrator<double, 2> funIntegrator2D(faceGeometry);

  std::cout << "started 1" << std::endl;    double rez1D_1 = funIntegrator1D.integrate(function1d1(), 1.0e-5);
  std::cout << "started 2" << std::endl;    double rez1D_2 = funIntegrator1D.integrate(function1d2(), 1.0e-5);
  std::cout << "started 3" << std::endl;    double rez1D_3 = funIntegrator1D.integrate(function1d3(), 1.0e-5);
  std::cout << "started 4" << std::endl;    double rez1D_4 = funIntegrator1D.integrate(function1d4(), 1.0e-5);

  std::cout << "started 5" << std::endl;    double rez2D_1 = funIntegrator2D.integrate(function2d1(), 1.0e-5);
  std::cout << "started 6" << std::endl;    double rez2D_2 = funIntegrator2D.integrate(function2d2(), 1.0e-5);
  std::cout << "started 7" << std::endl;    double rez2D_3 = funIntegrator2D.integrate(function2d3(), 1.0e-5);
  std::cout << "started 8" << std::endl;    double rez2D_4 = funIntegrator2D.integrate(function2d4(), 1.0e-5);
  std::cout << "started 9" << std::endl;    double rez2D_5 = funIntegrator2D.integrate(function2d5(), 1.0e-5);
  std::cout << "started 10" << std::endl;    double rez2D_6 = funIntegrator2D.integrate(function2d6(), 1.0e-5);
  std::cout << "started 11" << std::endl;    double rez2D_7 = funIntegrator2D.integrate(function2d7(), 1.0e-5);
  std::cout << "started 12" << std::endl;    double rez2D_8 = funIntegrator2D.integrate(function2d8(), 1.0e-5);


  std::cout << "---------------------------Summary: -----------------------------------------" << std::endl;
  std::cout << "Integrating f(x) = 1            over unit edge. Expected: 1.00000, result: " << rez1D_1 << ", error: " << fabs(rez1D_1 - 1.0) << std::endl;
  std::cout << "Integrating f(x) = x            over unit edge. Expected: 0.50000, result: " << rez1D_2 << ", error: " << fabs(rez1D_2 - 0.5) << std::endl;
  std::cout << "Integrating f(x) = x^3 - 3x + 3 over unit edge. Expected: 1.75000, result: " << rez1D_3 << ", error: " << fabs(rez1D_3 - 1.75) << std::endl;
  std::cout << "Integrating f(x) = sqrt(x)      over unit edge. Expected: 0.66667, result: " << rez1D_4 << ", error: " << fabs(rez1D_4 - 2.0/3) << std::endl;
  std::cout << "-----------------------------------------------------------------------------" << std::endl;
  std::cout << "Integrating f(x,y) = 1                      over unit triangle. Expected: 0.500000, result: " << rez2D_1 << ", error: " << fabs(rez2D_1 - 0.5) << std::endl;
  std::cout << "Integrating f(x,y) = 1 + x                  over unit triangle. Expected: 0.666667, result: " << rez2D_2 << ", error: " << fabs(rez2D_2 - 2.0/3) << std::endl;
  std::cout << "Integrating f(x,y) = 1 + x^2 + y^2          over unit triangle. Expected: 0.666667, result: " << rez2D_3 << ", error: " << fabs(rez2D_3 - 2.0/3) << std::endl;
  std::cout << "Integrating f(x,y) = x * y^2                over unit triangle. Expected: 0.016667, result: " << rez2D_4 << ", error: " << fabs(rez2D_4 - 1.0/60) << std::endl;
  std::cout << "Integrating f(x,y) = sqrt(x * y)            over unit triangle. Expected: 0.130900, result: " << rez2D_5 << ", error: " << fabs(rez2D_5 - M_PI/24) << std::endl;
  std::cout << "Integrating f(x,y) = 2000 * x^3 * y^3       over unit triangle. Expected: 1.785714, result: " << rez2D_6 << ", error: " << fabs(rez2D_6 - 25.0/14) << std::endl;
  std::cout << "Integrating f(x,y) = 10! * x^7 * y^10       over unit triangle. Expected: 0.545584, result: " << rez2D_7 << ", error: " << fabs(rez2D_7 - 0.545584) << std::endl;
  std::cout << "Integrating f(x,y) = sqrt(x^7 * y^10 + 0.5) over unit triangle. Expected: 0.353554, result: " << rez2D_8 << ", error: " << fabs(rez2D_8 - 0.353554) << std::endl;

  return 0;
}
