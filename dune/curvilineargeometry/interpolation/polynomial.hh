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

#include <config.h>

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <dune/common/fvector.hh>



namespace Dune {

template<class ctype>
struct PolynomialTraits
{
	class Monomial {

	public:
	  ctype pref_;
      std::vector<int> power_;


	  Monomial(ctype prefNew, std::vector<int> powerNew)
	  {
	    pref_ = prefNew;
	    power_ = powerNew;
	  }

	  // 1D initializer
	  Monomial(ctype prefNew, int x)
	  {
	    pref_ = prefNew;
	    power_.push_back(x);
	  }

	  // 2D initializer
	  Monomial(ctype prefNew, int x, int y)
	  {
	    pref_ = prefNew;
	    power_.push_back(x);
	    power_.push_back(y);
	  }

	  // 3D initializer
	  Monomial(ctype prefNew, int x, int y, int z)
	  {
	    pref_ = prefNew;
	    power_.push_back(x);
	    power_.push_back(y);
	    power_.push_back(z);
	  }

	  Monomial(const Monomial& other) :
		  pref_(other.pref_),
		  power_(other.power_)
	  {

	  }



	  int dim() const { return power_.size(); }


	  Monomial operator*(const Monomial & other) const {
		  assert(dim() == other.dim());
		  std::vector<int> newPower;
		  for(int i = 0; i < dim(); i++)  { newPower.push_back(power_[i] + other.power_[i]); }
		  return Monomial(pref_ * other.pref_, newPower);
	  }

	  Monomial derivative(int paramNo) const
	  {
		  assert((paramNo >=0 )&&(paramNo < dim()));
		  if (power_[paramNo] <= 0)  { return Monomial(0.0, std::vector<int>(dim(), 0)); }

		  Monomial monNew = *this;
		  monNew.pref_ *= monNew.power_[paramNo];
		  monNew.power_[paramNo]--;
		  return monNew;
	  }

	};

	// In order to compactify the Polynomial, we want to be able to sort its summands with respect
	// to the powers of each of the parameters
	static bool polySortOrder( Monomial A, Monomial B )
	{
		//assert(A.dim() == B.dim());
		int i = 0;
		int thisDim = A.power_.size();

		//while ( (i + 1 < thisDim) && (A.power_[i] == B.power_[i]) )  { i++; }
		//return (i < thisDim) ? A.power_[i] < B.power_[i] : 0;

		while ( (i < thisDim) && (A.power_[i] == B.power_[i]) )  { i++; }
		return (i < thisDim) ? A.power_[i] < B.power_[i] : 0;
	}
};



template<class ctype, int dim>
class Polynomial {

  typedef Polynomial<ctype, dim> LocalPolynomial;
  typedef typename PolynomialTraits<ctype>::Monomial Monomial;
  typedef typename std::vector<Monomial> SummandVector;
  typedef FieldVector< ctype, dim > LocalCoordinate;

  typedef unsigned int uint;

public:
  SummandVector poly_;

  Polynomial()                      { poly_.push_back(zeroMonomial()); }  // Empty polynomial

  // Polynomial from a summand
  Polynomial(Monomial polySM)  {
	  if (dim != polySM.dim())  { std::cout << "attempting to push a " << polySM.dim() << " monomial to a poly " << dim << std::endl;  }
	  assert(dim == polySM.dim());  // Must push polynomial of correct dimension
	  poly_.push_back(polySM);
  }
  //Polynomial(SummandVector polyNew) { poly_ = polyNew; }                  // Polynomial from a vector of summands

  Polynomial(const Polynomial & other) :
    poly_(other.poly_)
  {

  }

  Monomial zeroMonomial()  { return Monomial(0.0, std::vector<int> (dim, 0)); }


  /** \brief Add a summand to a Polynomial
   *  \param[in]  polySM  a summand
   */
  LocalPolynomial & operator+=(const Monomial & polySM)
  {
	  assert(dim == polySM.dim());  // Must push monomial of correct dimension
	  poly_.push_back(polySM);
	  return *this;
  }


  /** \brief Adds another Polynomial to this one
   *  \param[in]  polyNew  a Polynomial to add
   */
  LocalPolynomial & operator+=(const LocalPolynomial & other)
  {
	assert(other.size() > 0);   // Non-initialized polynomials should not exist
	for (uint i = 0; i < other.size(); i++)  { poly_.push_back(other.poly_[i]); }

	// Addition is likely to generate repeating powers, need to compactify
	compactify();
	return *this;
  }


  /** \brief Adds another Polynomial to this one
   *  \param[in]  polyNew  a Polynomial to add
   */
  LocalPolynomial & operator+=(const ctype c)
  {
	  poly_.push_back(Monomial(c, std::vector<int> (dim, 0)));
	  return *this;
  }


  /** \brief Multiplies this Polynomial by a scalar
   *
   *  \param[in]  c     scalar to multiply by
   */
  LocalPolynomial & operator*=(const ctype c)
  {
    // If we multiply by zero, return zero Polynomial
    if (fabs(c) < 1.0e-25)  { poly_ = SummandVector(1, zeroMonomial()); }
    else
    {
        for (uint i = 0; i < poly_.size(); i++) { poly_[i].pref_ *= c; }
    }
    return *this;
  }


  /** \brief Multiply incoming Polynomial by scalar, then add to this one
   *
   *  \param[in]  polyNew  Polynomial to scale
   *  \param[in]  c        scalar to multiply the summand by
   */
  void axpy(const LocalPolynomial & other, ctype c) {
	  assert(other.size() > 0);   // Non-initialized polynomials should not exist
      *this += other * c;
  }


  /** \brief Multiply this Polynomial by another one
   *
   *  \param[in]  polyNew     Polynomial to multiply by
   */
  LocalPolynomial & operator*=(const LocalPolynomial & other)
  {
	assert(other.size() > 0);   // Non-initialized polynomials should not exist
    SummandVector rez;

    ctype magnTmp = magnitude() * other.magnitude();

    for (uint i = 0; i < poly_.size(); i++) {
      for (uint j = 0; j < other.size(); j++) {
    	Monomial mult = poly_[i] * other.poly_[j];

        // Only add a summand if it is not too small
        if (fabs(mult.pref_) > magnTmp * 1.0e-10)  { rez.push_back(mult); }
      }
    }

    // If by some magic polynomial has no entries, make it zero
    if (rez.size() == 0)  {
    	//std::cout << "--- Warining: Polynomial Multiplication resulted in empty polynomial when multiplying " << to_string() << " and " << other.to_string() << std::endl;
    	poly_ = SummandVector(1, zeroMonomial());
    }
    // The product is likely to have produced several summands with the same power. Need to compactify
    else  {
    	poly_ = rez;
    	compactify();
    }

    return *this;
  }


  /** \brief Addition for two polynomials
   *
   *  \param[in]  a     Polynomial to add
   *  \returns    sum of polynomials
   */
  LocalPolynomial operator+(const LocalPolynomial & other) const {
      LocalPolynomial rez = *this;
      rez += other;
      return rez;
  }


  /** \brief Addition of a scalar
   *
   *  \param[in]  a     a scalar
   *  \returns    new Polynomial with added scalar
   */
  LocalPolynomial operator+(const ctype a) const {
      LocalPolynomial rez = *this;
      rez += Monomial(a, std::vector<int> (dim, 0));
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
      rez += Monomial(-a, std::vector<int> (dim, 0));
      return rez;
  }

  /** \brief Multiplication by scalar
   *
   *  \param[in]  a     a scalar
   *  \returns    new Polynomial multiplied by a scalar
   */
  LocalPolynomial operator*(const ctype a) const {
      LocalPolynomial rez = *this;
      rez *= a;
      return rez;
  }

  /** \brief Multiplication of two polynomials
   *
   *  \param[in]  a     Polynomial to multiply
   *  \returns    multiplication of polynomials
   */
  LocalPolynomial operator*(const LocalPolynomial & a) const {
      LocalPolynomial rez = *this;
      rez *= a;
      return rez;
  }


  /** \brief Number of monomials */
  uint size() const { return poly_.size(); }


  /** \brief Order of Polynomial
   *
   *  \returns    Returns the total summand power (order) maximized over summands
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
  ctype magnitude() const {
	ctype rez = 0;

    for (uint i = 0; i < poly_.size(); i++) { rez = std::max(rez, fabs(poly_[i].pref_)); }

    return rez;
  }


  /** \brief Takes derivative of a Polynomial
   *
   *  \param[in]  paramNo     index to parameter to differentiate with {0,1,2}
   *  \returns    Polynomial - derivative of this Polynomial
   */
  LocalPolynomial derivative(int paramNo) const {
	assert((paramNo>=0)&&(paramNo < dim));

	LocalPolynomial newPoly;

    for (uint i = 0; i < poly_.size(); i++) {
      // Only add this monomial if it does not differentiate to 0
      if (poly_[i].power_[paramNo] > 0)  { newPoly += poly_[i].derivative(paramNo); }
    }
    return newPoly;
  }

  /** \brief Evaluate the Polynomial at given coordinates
   *
   *  \param[in]  point     Local coordinate to evaluate at
   *  \returns    value of the polynoimal at the given point
   */
  ctype evaluate(const LocalCoordinate & point) const {
	ctype rez = 0;

    for (uint i = 0; i < poly_.size(); i++) {
      ctype powTmp = 1.0;
      for (int d = 0; d < dim; d++) { powTmp *= pow(point[d], poly_[i].power_[d]); }
      rez += poly_[i].pref_ * powTmp;
    }
    return rez;
  }

  /** \brief Integrates the Polynomial over a reference simplex (edge, triangle or tetrahedron, depending on dim)
   *
   *  \returns    Returns the value of the analytical integral of the Polynomial
   */
  ctype integrateRefSimplex() const {
	assert((dim > 0) && (dim <= 3));

	ctype rez = 0.0;

    // Generate factorial list up to the needed order
    uint polyorder = order();
    std::vector<double> factorial(2, 1.0);  // NOTE: !!! factorial should be a Real variable, do not use ctype here
    for (uint i = 2; i <= polyorder + dim; i++) { factorial.push_back( i * factorial[i - 1] ); }

    // Compute the analytical integral
    for (uint i = 0; i < poly_.size(); i++) {
      switch (dim)
      {
          case 1 : rez += poly_[i].pref_ / double(poly_[i].power_[0] + 1);  break;
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
    ctype tolerance = magnitude() * 1.0e-10;

    for (uint i = 0; i < poly_.size(); i++) {
      if (fabs(poly_[i].pref_) > tolerance) { poly_cleaned.push_back(poly_[i]); }
    }

    // If the Polynomial has no non-zero terms, it must have 1 zero-term not to have zero length
    if (poly_cleaned.size() == 0) { poly_cleaned.push_back(zeroMonomial()); }

    poly_ = poly_cleaned;
  }

  /** \brief Adds up repeating powers.
   *
   * Sorts the summands by increasing 1st -> 2nd -> 3rd parameter power (if available),
   * then merges the neighboring summands if they have the same powers
   *
   */

  void compactify() {
	assert(size() > 0);

    // This way we make sure that the identical powers are consecutive
    std::sort(poly_.begin(), poly_.end(), PolynomialTraits<ctype>::polySortOrder);

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

    for (uint i = 0; i < poly_.size(); i++) {
      if (poly_[i].pref_ >= 0) { out_str << "+"; }

      out_str << poly_[i].pref_ << " ";
      for (int d = 0; d < dim; d++) { out_str << coordNames[d] << "^" << poly_[i].power_[d] << " "; }
    }
    return out_str.str();
  }

};


/** \brief Generate an identity Polynomial*/
template<class ctype, int dim>
Polynomial<ctype, dim> zeroPolynomial() {
    return Polynomial<ctype, dim> ( typename PolynomialTraits<ctype>::Monomial(0.0, std::vector<int>(dim, 0)) );
}


/** \brief Generate an identity Polynomial*/
template<class ctype, int dim>
Polynomial<ctype, dim> identityPolynomial() {
    return Polynomial<ctype, dim> ( typename PolynomialTraits<ctype>::Monomial(1.0, std::vector<int>(dim, 0)) );
}


} // namespace Dune

#endif // #ifndef DUNE_POLYNOMIAL_HH
