#ifndef DUNE_CURVILINEAR_SINGULARITY_HELPER
#define DUNE_CURVILINEAR_SINGULARITY_HELPER



namespace Dune {



/** This class transforms the domain of a singular
 integrand

             1 / |p(x) - p(x0)|

 from the curvilinear reference triangle to the
 reference square, eliminating the singularity

 Algorithm:
 0) This algorithm ONLY works if x0 = (1,0), thus
 the singularity has to be at a particular corner
 of the reference triangle. While this is definitely
 not the general case, in general it is possible
 to split any triangle into 2-3 triangles, such that
 the singularity is at a corner of each of them

  |\           |\
  | \          ||\
  |  \         |\ \
  |  /\        | \ \
  | /  \       | /\ \
  |/    \      |/  \ \
  -------      -------

*/
template <class ctype, int cdim>
class CurvilinearSingularityHelper {

	typedef LagrangeInterpolator<ctype, 2, cdim>     Interpolator;

	typedef typename Interpolator::LocalCoordinate   LocalCoordinate;
	typedef typename Interpolator::GlobalCoordinate  GlobalCoordinate;
	typedef typename Interpolator::Monomial          Monomial;
	typedef typename Interpolator::LocalPolynomial   LocalPolynomial;
	typedef typename Interpolator::PolynomialVector  PolynomialVector;
	typedef typename Interpolator::PolynomialGlobalCoordinate  PolyGlobalCoord;

	typedef std::pair<LocalPolynomial, LocalPolynomial>  polyPair;

	const unsigned int X_COORD = 0;
	const unsigned int Y_COORD = 1;

public:

	/** \brief Constructor
	 * \param[in]  map    The mapping from local to global coordinates */
	CurvilinearSingularityHelper(PolyGlobalCoord & map, unsigned int order) {
		// A good polynomial map should be of the right size
		assert(map.size() == cdim);

		// Generate set of polynomials for (x^n - 1) / (x - 1)
		// Namely {1, 1+x, 1+x+x^2, 1+x+x^2+x^3, ...}
		PolynomialVector powPolyTmp(order);
		for (unsigned int i = 0; i < order; i++)  { powPolyTmp[i] = LocalPolynomial(Monomial(1.0, i)); }
		for (unsigned int i = 1; i < order; i++)  { powPolyTmp[i] += powPolyTmp[i-1]; }

		// Split into y-independent and y-proportional part
		// Namely p(x,y) = P0 + y*R = p(x,0) + y*R(x,y)
		PolyGlobalCoord P0(cdim);
		D_ = PolyGlobalCoord(cdim);
		R_ = PolyGlobalCoord(cdim);
		for (unsigned int i = 0; i < cdim; i++) {
			// Initialize the polynomial classes
			P0[i] = LocalPolynomial::zeroMonomial();
			D_[i] = P0[i];
			R_[i] = P0[i];

			for (unsigned int j = 0; j < map[i].poly_.size(); j++) {

				Monomial thisM = map[i].poly_[j];
				if (thisM.power_[Y_COORD] == 0) {  P0[i] += thisM; }
				else                            {
					thisM.power_[Y_COORD]--;
					R_[i] += thisM;
				}
			}
		}

		// Perform analytical division of D = (p(x,0)-p(1,0))/(1-x)
		// This is effectively a sum of terms pref*((x^n - 1) / (x - 1))
		for (unsigned int i = 0; i < cdim; i++) {
			for (unsigned int j = 0; j < P0[i].poly_.size(); j++) {

				unsigned int xpow = P0[i].poly_[j].power_[X_COORD];
				ctype        pref = P0[i].poly_[j].pref_;

				if (xpow > 0) {
					LocalPolynomial summand = powPolyTmp[xpow - 1];
					summand *= pref * (-1);
					D_ += summand;
				}
			}
		}
	}


	/** Evaluate the resulting fraction in transformed coordinates
	 *
	 * VERY IMPORTANT: Transformed coordinates (x,t) are to be integrated
	 * over the reference square [0,1]^2, not the reference triangle,
	 * as were the original coordinates (x,y)
	 *
	 * Note that evaluating D and R separately is much cheaper than combining them into
	 * a single polynomial and taking the norm. At the moment we just need to evaluate
	 * a polynomial of order n six times and do a few small arithmetic operations. Otherwise,
	 * a polynomial of order 2 * (2n + 1) needs to be evaluated once, which is definitely worse.
	 * */
	ctype evalDenominator(ctype x, ctype t) const {

		LocalCoordinate r;
		r[0] = x;
		r[1] = t * (1-x);

		GlobalCoordinate drez = evalPolyVec(D_, r);
		GlobalCoordinate rrez = evalPolyVec(R_, r);
		rrez *= t;

		return (drez + rrez).two_norm();
	}



private:

	GlobalCoordinate evalPolyVec(const PolyGlobalCoord & pv, const LocalCoordinate & r) {
		GlobalCoordinate rez;
		for (unsigned int i = 0; i < cdim; i++)  { rez[i] = PolyGlobalCoord[i].evaluate(r); }
		return rez;
	}

	PolyGlobalCoord D_;
	PolyGlobalCoord R_;


};


}





#endif // DUNE_CURVILINEAR_SINGULARITY_HELPER
