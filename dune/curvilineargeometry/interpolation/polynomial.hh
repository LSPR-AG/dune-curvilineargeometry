/*******************************************************************
 * LocalPolynomial Class
 * 
 * author: Aleksejs Fomins
 * date: 01.07.2014 - created.
 * 
 * description:
 * Provides the implementation of a Polynomial represented by the sum $p(x,y,z) = \sum_i A_i x^{m_i} y^{n_i} z^{k_i}$
 * - The number of summands of the Polynomial and its dimension are not restricted
 * - The polynomials can be added, subtracted, multiplied by scalar and by Polynomial, evaluated at a given point, differentiated wrt given dimension, integrated over reference simplex for dimensions 1-3. Can also print Polynomial to screen
 * - Runs compactification routine, which adds up summands with the same power whenever needed, thus saving space
 * - Runs cleanup routine, which deletes terms with low prefactor, thus saving space
 * 
 * 
 * FIXME: cleanup routine too naive, since it is based on absolute prefactor magnitude. A better algorithm would find the magnitude of the leading order. If that magnitude is effectively zero (with very low tolerance), then delete all terms except of a constant zero. Otherwise, only delete terms which are small relatively to the leading term. Introduce constants POLYNOMIAL_ABSOLUTE_TOLERANCE = 10e-25 and POLYNOMIAL_RELATIVE_TOLERANCE = 10e-10
 * 
 * TODO: Replace print to screen to print to string, more scalable. Also, implement printing for higher dimensions by using alphabet.
 * TODO: Implement integration over reference simplex in arbitrary dimension. Not sure what for, but I know the analytic expression.
 * TODO: Implement integration over reference hypercube in arbitrary dimension. It is trivial, and may be necessary
 * 
 *******************************************************************/



#ifndef DUNE_POLYNOMIAL_HH
#define DUNE_POLYNOMIAL_HH


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <config.h>
#include <dune/common/fvector.hh>



namespace Dune {


struct PolynomialTraits
{
	class Monomial {
	public:
	  double pref_;
	  std::vector<int> power_;

	  Monomial(double prefNew, std::vector<int> powerNew)
	  {
	    pref_ = prefNew;
	    power_ = powerNew;
	  }

	  // 1D initializer
	  Monomial(double prefNew, int x)
	  {
	    pref_ = prefNew;
	    power_.push_back(x);
	  }

	  // 2D initializer
	  Monomial(double prefNew, int x, int y)
	  {
	    pref_ = prefNew;
	    power_.push_back(x);
	    power_.push_back(y);
	  }

	  // 3D initializer
	  Monomial(double prefNew, int x, int y, int z)
	  {
	    pref_ = prefNew;
	    power_.push_back(x);
	    power_.push_back(y);
	    power_.push_back(z);
	  }
	};

	// In order to compactify the Polynomial, we want to be able to sort its summands with respect
	// to the powers of each of the parameters
	static bool polySortOrder( Monomial A, Monomial B )
	{
		int i = 0;
		int thisDim = A.power_.size();
		while ( (i + 1 < thisDim) && (A.power_[i] == B.power_[i]) )  { i++; }

		return (i < thisDim) ? A.power_[i] < B.power_[i] : 0;
	}
};



template<class ctype, int dim>
class Polynomial {

  typedef Polynomial<ctype, dim> LocalPolynomial;
  typedef PolynomialTraits::Monomial     Monomial;
  typedef typename std::vector<Monomial> SummandVector;
  typedef FieldVector< ctype, dim > LocalCoordinate;

public:
  SummandVector poly_;

  // Constructor - creates Polynomial with 1 summand
  Polynomial() { }
  Polynomial(SummandVector polyNew) { poly_ = polyNew; }
  Polynomial(Monomial polySM) { poly_.push_back(polySM); }

  /** \brief Add a summand to a Polynomial
   *
   *  \param[in]  polySM  a summand
   */
  void append(Monomial polySM) { poly_.push_back(polySM); }

  /** \brief Adds another Polynomial to this one
   *
   *  \param[in]  polyNew  a Polynomial to add
   */
  void mergeTo(LocalPolynomial polyNew) {
    for (uint i = 0; i < polyNew.poly_.size(); i++) { poly_.push_back(polyNew.poly_[i]); }

    // Addition is likely to generate repeating powers, need to compactify
    compactify();
  }

  /** \brief Multiplies this Polynomial by a scalar
   *
   *  \param[in]  c     scalar to multiply by
   */
  void multScalar(double c) {
    // If we multiply by zero, return zero Polynomial
    if (fabs(c) < 1.0e-25) {
        SummandVector zeroPoly;
        zeroPoly.push_back(Monomial(0, std::vector<int> (dim, 0)));
        poly_ = zeroPoly;
    } else
    {
        for (uint i = 0; i < poly_.size(); i++) { poly_[i].pref_ *= c; }
    }
  }

  /** \brief Multiply incoming Polynomial by scalar, then add to this one
   *
   *  \param[in]  polyNew  Polynomial to scale
   *  \param[in]  c        scalar to multiply the summand by
   */
  void axpy(LocalPolynomial polyNew, double c) {
      LocalPolynomial polyTmp = polyNew;
      polyTmp.multScalar(c);
      mergeTo(polyTmp);
  }

  /** \brief Multiply this Polynomial by another one
   *
   *  \param[in]  polyNew     Polynomial to multiply by
   */
  void multPoly(LocalPolynomial polyNew) {
    SummandVector rez;

    double magnTmp = magnitude() * polyNew.magnitude();

    for (uint i = 0; i < poly_.size(); i++) {
      for (uint j = 0; j < polyNew.poly_.size(); j++) {
          double prefTmp = poly_[i].pref_ * polyNew.poly_[j].pref_;

          // Only add a summand if it is not too small
          if (fabs(prefTmp) > magnTmp * 1.0e-10)
          {
              std::vector<int> powerRez;
              for (int d = 0; d < dim; d++) { powerRez.push_back(poly_[i].power_[d] + polyNew.poly_[j].power_[d]); }
              rez.push_back(Monomial(prefTmp , powerRez));
          }
      }
    }
    poly_ = rez;

    if (poly_.size() == 0)
    {
        poly_.push_back(Monomial(0, std::vector<int> (dim, 0)));
    } else
    {
        // The product is likely to have produced several summands with the same power. Need to compactify
        compactify();
    }
  }

  /** \brief Addition for two polynomials
   *
   *  \param[in]  a     Polynomial to add
   *  \returns    sum of polynomials
   */
  LocalPolynomial operator+(const LocalPolynomial & a) const {
      LocalPolynomial rez = *this;
      rez.mergeTo(a);
      return rez;
  }

  /** \brief Addition of a scalar
   *
   *  \param[in]  a     a scalar
   *  \returns    new Polynomial with added scalar
   */
  LocalPolynomial operator+(const ctype a) const {
      LocalPolynomial rez = *this;
      rez.mergeTo(LocalPolynomial(Monomial(a, std::vector<int> (dim, 0))));
      return rez;
  }

  /** \brief Subtraction for two polynomials
   *
   *  \param[in]  a     Polynomial to add
   *  \returns    difference of polynomials
   */
  LocalPolynomial operator-(const LocalPolynomial & a) const {
      LocalPolynomial rez = *this;
      rez.axpy(a, -1);
      return rez;
  }

  /** \brief Subtraction of a scalar
   *
   *  \param[in]  a     a scalar
   *  \returns    new Polynomial with subtracted scalar
   */
  LocalPolynomial operator-(const ctype a) const {
      LocalPolynomial rez = *this;
      rez.mergeTo(LocalPolynomial(Monomial(-a, std::vector<int> (dim, 0))));
      return rez;
  }

  /** \brief Multiplication by scalar
   *
   *  \param[in]  a     a scalar
   *  \returns    new Polynomial multiplied by a scalar
   */
  LocalPolynomial operator*(const ctype a) const {
      LocalPolynomial rez = *this;
      rez.multScalar(a);
      return rez;
  }

  /** \brief Multiplication of two polynomials
   *
   *  \param[in]  a     Polynomial to multiply
   *  \returns    multiplication of polynomials
   */
  LocalPolynomial operator*(const LocalPolynomial & a) const {
      LocalPolynomial rez = *this;
      rez.multPoly(a);
      return rez;
  }

  // Returns the total summand power (order) maximized over summands

  /** \brief Order of Polynomial
   *
   *  \returns    Highest total poly_ of a summand of a Polynomial
   */
  uint order() const {
    int rez = 0;

    for (uint i = 0; i < poly_.size(); i++) {
      int maxTmp = 0;
      for (int d = 0; d < dim; d++) { maxTmp += poly_[i].power_[d]; }
      rez = std::max(rez, maxTmp);
    }

    return rez;
  }

  /** \brief Magnitude of Polynomial
   *
   *  \returns    Highest absolute value of a prefactor of a summand of a Polynomial
   */
  double magnitude() const {
        double rez = 0;

        for (uint i = 0; i < poly_.size(); i++) { rez = std::max(rez, fabs(poly_[i].pref_)); }

        return rez;
  }


  /** \brief Takes derivative of a Polynomial
   *
   *  \param[in]  paramNo     index to parameter to differentiate with {0,1,2}
   *  \returns    Polynomial - derivative of this Polynomial
   */
  LocalPolynomial derivative(int paramNo) const {
    LocalPolynomial rez (poly_);

    for (uint i = 0; i < rez.poly_.size(); i++) {
      rez.poly_[i].pref_ *= poly_[i].power_[paramNo];
      rez.poly_[i].power_[paramNo] -= 1;
    }

    // Delete all summands that are now 0
    rez.cleanUp();
    return rez;
  }

  /** \brief Evaluate the Polynomial at given coordinates
   *
   *  \param[in]  point     Local coordinate to evaluate at
   *  \returns    value of the polynoimal at the given point
   */
  double evaluate(const LocalCoordinate & point) const {
    double rez = 0;

    for (uint i = 0; i < poly_.size(); i++) {
      double powTmp = 1;
      for (int d = 0; d < dim; d++) { powTmp *= pow(point[d], poly_[i].power_[d]); }
      rez += poly_[i].pref_ * powTmp;
    }
    return rez;
  }

  /** \brief Integrates the Polynomial over a reference simplex (edge, triangle or tetrahedron, depending on dim)
   *
   *  \returns    Returns the value of the analytical integral of the Polynomial
   */
  double integrateRefSimplex() const {
    double rez = 0;

    // Generate factorial list up to the needed order
    uint polyorder = order();
    std::vector<double> factorial(2, 1.0);
    for (uint i = 2; i <= polyorder + dim; i++) { factorial.push_back( i * factorial[i - 1] ); }

    // Compute the analytical integral
    for (uint i = 0; i < poly_.size(); i++) {
      switch (dim)
      {
          case 1 : rez += poly_[i].pref_ / (poly_[i].power_[0] + 1);  break;
          case 2 : rez += poly_[i].pref_ * factorial[poly_[i].power_[0]] * factorial[poly_[i].power_[1]] / factorial[poly_[i].power_[0] + poly_[i].power_[1] + 2];  break;
          case 3 : rez += poly_[i].pref_ * factorial[poly_[i].power_[0]] * factorial[poly_[i].power_[1]] * factorial[poly_[i].power_[2]] / factorial[poly_[i].power_[0] + poly_[i].power_[1] + poly_[i].power_[2] + 3];  break;
      }
    }
    return rez;
  }

  /** \brief Removes all summands that are effectively 0
   *
   */
  void cleanUp() {
    SummandVector poly_cleaned;

    // Everything that is 10^10 less than the principal term is insignificant
    double tolerance = magnitude() * 1.0e-10;

    for (uint i = 0; i < poly_.size(); i++) {
      bool isZero = fabs(poly_[i].pref_) < tolerance;

      for (int d = 0; d < dim; d++) { isZero = isZero || (poly_[i].power_[d] < 0); }

      if (!isZero) { poly_cleaned.push_back(poly_[i]); }
    }

    // If the Polynomial has no non-zero terms, it must have 1 zero-term not to have zero length
    if (poly_cleaned.size() == 0) { poly_cleaned.push_back(Monomial(0, std::vector<int> (dim, 0))); }

    poly_ = poly_cleaned;
  }

  /** \brief Adds up repeating powers.
   *
   * Sorts the summands by increasing 1st -> 2nd -> 3rd parameter power (if available),
   * then merges the neighboring summands if they have the same powers
   *
   */

  void compactify() {
    // This way we make sure that the identical powers are consecutive
    std::sort(poly_.begin(), poly_.end(), PolynomialTraits::polySortOrder);

    SummandVector polyNew;
    polyNew.push_back(poly_[0]);

    uint i = 0;
    uint j = 1;

    while (j < poly_.size())
    {
      // Check if the powers match
      bool samePow = true;
      for (int d = 0; d < dim; d++) { samePow = samePow && ( polyNew[i].power_[d] == poly_[j].power_[d] ); }

      if (samePow) { polyNew[i].pref_ += poly_[j].pref_; j++;} // If they match, merge the summands
      else { polyNew.push_back(poly_[j]); i++; j++; } // Otherwise construct new summand
    }

    poly_ = polyNew;
  }

  /** \brief Convert the Polynomial to a string for future output
   *
   */
  std::string to_string() const {
    std::vector<char> coordNames;
    coordNames.push_back('x');
    coordNames.push_back('y');
    coordNames.push_back('z');

    std::stringstream out_str;

    //std::cout << "Printing Polynomial: ";
    for (uint i = 0; i < poly_.size(); i++) {
      if (poly_[i].pref_ >= 0) { out_str << "+"; }

      out_str << poly_[i].pref_ << " ";
      for (int d = 0; d < dim; d++) { out_str << coordNames[d] << "^" << poly_[i].power_[d] << " "; }
      std::cout << " ";
    }
    return out_str.str();
  }

};

/** \brief Generate an identity Polynomial
 *
 */
template<class ctype, int dim>
Polynomial<ctype, dim> identityPolynomial() {
    return Polynomial<ctype, dim> ( PolynomialTraits::Monomial(1.0, std::vector<int>(dim, 0)) );
}


} // namespace Dune

#endif // #ifndef DUNE_POLYNOMIAL_HH
