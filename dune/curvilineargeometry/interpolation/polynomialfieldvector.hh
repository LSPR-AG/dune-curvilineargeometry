/*******************************************************************
 * Polynomial Field Vector Class
 * 
 * author: Aleksejs Fomins
 * date: 18.12.2014 - created.
 * 
 * description:
 * Provides the implementation of a polynomial vector,
 * 
 * 
 * [FIXME] Throw errors if the dimensions or lengths of polynomials/monomials mismatch
 * 
 *******************************************************************/



#ifndef DUNE_POLYNOMIAL_FIELD_VECTOR_HH
#define DUNE_POLYNOMIAL_FIELD_VECTOR_HH


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <config.h>
#include <dune/common/fvector.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>



namespace Dune {



// [TODO]  Merge with DifferentialHelper at some point
template<class ctype, int mydim, int cdim>
class PolynomialFieldVector {
  typedef PolynomialFieldVector<ctype, mydim, cdim>  This;

  typedef Polynomial<ctype, mydim>              LocalPolynomial;
  typedef typename LocalPolynomial::Monomial    Monomial;
  typedef FieldVector< ctype, mydim >           LocalCoordinate;
  typedef FieldVector< ctype, cdim >            GlobalCoordinate;
  typedef FieldVector< LocalPolynomial, cdim >  PolyVector;

  typedef unsigned int   uint;

public:

  PolynomialFieldVector()  {  }

  // Brackets operator for reading
  LocalPolynomial operator [](uint i) const  {return polyvec_[i];}

  // Brackets operator for assign
  LocalPolynomial & operator [](uint i)      {return polyvec_[i];}


  /** \brief Adds another Polynomial vector to this one
   *  \param[in]  other  a Polynomial vector to add
   */
  This & operator+=(const This & other)
  {
	  for (uint i = 0; i < cdim; i++)  { polyvec_[i] += other[i]; }
	  return *this;
  }


  /** \brief Subtracts another Polynomial vector to this one
   *  \param[in]  other  a Polynomial vector to subtract
   */
  This & operator-=(const This & other)
  {
	  for (uint i = 0; i < cdim; i++)  { polyvec_[i] -= other[i]; }
	  return *this;
  }


  /** \brief Multiplies this Polynomial vector by a scalar
   *
   *  \param[in]  c     scalar to multiply by
   */
  This & operator*=(const ctype c)
  {
	  for (uint i = 0; i < cdim; i++) { polyvec_[i] *= c; }
	  return *this;
  }


  /** \brief Multiplies this Polynomial vector by a polynomial
   *
   *  \param[in]  poly     a polynomial to multiply by
   */
  This & operator*=(const LocalPolynomial & poly)
  {
	  for (uint i = 0; i < cdim; i++) { polyvec_[i] *= poly; }
	  return *this;
  }


  /** \brief Multiply another Polynomial vector by scalar, then add to this one
   *
   *  \param[in]  other    Polynomial vector to scale
   *  \param[in]  c        scalar to multiply the summand by
   */
  void axpy(const LocalPolynomial & other, ctype c) {
	  LocalPolynomial tmp = other;
	  tmp *= c;
      *this += tmp;
  }


  /** \brief Evaluate a polynomial vector at a local cooridnate */
  GlobalCoordinate evaluate(LocalCoordinate x) const
  {
	  GlobalCoordinate rez;
	  for (int i = 0; i < cdim; i++)  { rez[i] = polyvec_[i].evaluate(x); }
	  return rez;
  }


  /** \brief Compute the gradient of a scalar polynomial */
  static This grad(const LocalPolynomial & other) {
	  assert(mydim == cdim);  // Gradient only defined if dimension of result vector is equal to the number of variables
	  This rez;
	  for (uint i = 0; i < cdim; i++) { rez[i] = other.derivative(i); }
	  return rez;
  }


  /** \brief Divergence differential operator. Only defined for 3D in 3D at the moment */
  LocalPolynomial div() const
  {
	  assert(cdim == 3);
	  assert(mydim == 3);

	  LocalPolynomial rez;
	  for (int i = 0; i < cdim; i++)  { rez += polyvec_[i].derivative(i); }

	  rez.compactify();   // Add-up all monomials that are the same
	  rez.cleanUp();      // Clean-up in case some summands summed up to 0

	  return rez;
  }


  /** \brief Rotor (AKA Curl) differential operator. Only defined for 3D in 3D at the moment */
  This rot() const {
	  assert(cdim == 3);
	  assert(mydim == 3);

	  This rez;
	  rez[0] = polyvec_[2].derivative(1) - polyvec_[1].derivative(2);
	  rez[1] = polyvec_[0].derivative(2) - polyvec_[2].derivative(0);
	  rez[2] = polyvec_[1].derivative(0) - polyvec_[0].derivative(1);

	  for (int i = 0; i < cdim; i++) {
		  rez[i].compactify();   // Add-up all monomials that are the same
		  rez[i].cleanUp();      // Clean-up in case some summands summed up to 0
	  }

	  return rez;
  }


  /** \brief Polynomial map composition. If the argument is a map between new and old local coordinates,
   * then the result is this polynomial vector written in new local coordinates. This is achieved by simply
   * substituting polynomials in place of variables, and then converting to canonical polynomial form
   *
   * Algoritm-polySubstitution
   *
   *  1) Get P1, P2
   *  2) For P1, for each variable xyz, evaluate max power
   *  3) For P2, for each coordinate, obtain all its powers up to the one evaluated for P1
   *  4) Loop over P1 monomials, obtain products of polynomials in question
   *
   *  Perhaps there is a more efficient algorithm by considering 3d monomials instead of tensor product of 1D
   *  [TODO] Store all computed monomials, as they may be reused again in next coordinate
   *  [TODO] Store every computed product of first two coordinates. Then, after computing xyz, computing xyz^2 will be 1 multiplication instead of 2
   *
   */
  template <int mydim1>
  PolynomialFieldVector<ctype, mydim1, cdim> composite(const PolynomialFieldVector<ctype, mydim1, mydim> & map) const
  {
	  typedef Polynomial<ctype, mydim1> LocalPolyNew;
	  std::vector<std::vector<LocalPolyNew > > powerExpansion(mydim, std::vector<LocalPolyNew>());


	  // For each variable of original polynomial, find highest power in which it appears in monomials
	  // Then, pre-calculate all polynomials x, x*x, x*x*x, ... until that power and store them in a vector
	  for (uint iVar = 0; iVar < mydim; iVar++) {
		  // Find the maximal power of this variable among all monomials of all coordinates of this polyvector
		  uint orderThis = 0;
		  for (uint iDim = 0; iDim < cdim; iDim++)  {
			  int ordNew = polyvec_[iDim].order(iVar);

			  assert(ordNew >= 0);
			  orderThis = std::max(orderThis, (uint)ordNew);
		  }

		  // Compute all powers that contribute to this variable
		  if (orderThis > 0) {
			  powerExpansion[iVar] = std::vector<LocalPolyNew>(orderThis);
			  powerExpansion[iVar][0] = map[iVar];

			  //std::cout << "Power expansion var=" << iVar << ", ordExt=" << 0 << ", size=" << powerExpansion[iVar][0].size() << ", mag=" << powerExpansion[iVar][0].magnitude() << ", ordInt=" << powerExpansion[iVar][0].order() << std::endl;

			  for (uint iOrd = 1; iOrd < orderThis; iOrd++) {
				  powerExpansion[iVar][iOrd] = powerExpansion[iVar][iOrd-1] * map[iVar];
				  //std::cout << "Power expansion var=" << iVar << ", ordExt=" << iOrd << ", size=" << powerExpansion[iVar][iOrd].size() << ", mag=" << powerExpansion[iVar][iOrd].magnitude() << ", ordInt=" << powerExpansion[iVar][iOrd].order() << std::endl;
			  }
		  }


	  }

	  // For each coordinate, for each monomial, calculate it as a product of expanded coordinates to the corresponding power, given by the original polynomial
	  PolynomialFieldVector<ctype, mydim1, cdim> rez;
	  for (uint iDim = 0; iDim < cdim; iDim++)  {
		  for (uint iMon = 0; iMon < polyvec_[iDim].size(); iMon++)  {

			  const Monomial & mon = polyvec_[iDim].poly_[iMon];

			  // 1) Loop over mydim until find coordinate with power > 0
			  // 2) Multiply by expansion of that coordinate
			  // 3) Loop over remaining mydim, and multiply if power > 0
			  // 4) Multiply by prefactor of original monomial
			  // VERY IMPORTANT! powerExpansion stores powers with index=power-1, because storing power 0 is pointless. DO NOT MESS UP THE INDEX!!!
			  int iVarTmp = 0;
			  while ((mon.power_[iVarTmp] == 0)&&(iVarTmp < mydim)) {  iVarTmp++; }

			  LocalPolyNew p;

			  if (iVarTmp == mydim)  {  // If this is a constant monomial, there is nothing to expand
				  p = LocalPolyNew( LocalPolyNew::constMonomial(mon.pref_)); }
			  else {
				  p = powerExpansion[iVarTmp][mon.power_[iVarTmp] - 1];

				  // If there is more than one variable contributing to this monomial, all variable expansions should be multiplied together
				  for (uint iVar = iVarTmp + 1; iVar < mydim; iVar++)  {
					  if (mon.power_[iVar] != 0)  { p *= powerExpansion[iVar][mon.power_[iVar] - 1]; }
				  }

				  // Note: Instead, one could have created a constant polynomial with pref index, and simply multiplied it by all expansions
				  // However, multiplying a polynomial by a constant is cheaper than multiplying a polynomial by a constant polynomial, because
				  // in polynomial multiplication, a few special cases need to be checked for which do not apply here, and that slows down the code.
				  p *= mon.pref_;
			  }

			  rez[iDim] += p;
			  //std::cout << "Dim=" << iDim << ", mag=" << p.magnitude() << ", ord=" << p.order() << std::endl;
		  }
	  }

	  return rez;
  }

  template<class T>
  std::string vec2str(const std::vector<T> & vec) const {
	  std::stringstream s;
	  s << "{";
	  for (int i = 0; i < vec.size(); i++)  { s << vec[i] << " "; }
	  s << "}";
	  return s.str();
  }


  /** \brief Convert the polynomial vector to a string for text output */
  std::string to_string() const
  {
	  std::stringstream rez;

	  rez << "{" << std::endl;
	  for (uint i = 0; i < cdim; i++)  { rez << "  {" << polyvec_[i].to_string() << "}" << std::endl;}
	  rez << "}" << std::endl;

	  return rez.str();
  }



private:

  PolyVector polyvec_;


};

} // namespace Dune

#endif // #ifndef DUNE_POLYNOMIAL_FIELD_VECTOR_HH
