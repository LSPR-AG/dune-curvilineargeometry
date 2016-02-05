#ifndef DUNE_CURVILINEAR_SINGULARITY_HELPER
#define DUNE_CURVILINEAR_SINGULARITY_HELPER


#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>


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

	static const int mydim = 2;

	typedef LagrangeInterpolator<ctype, mydim, mydim>   InterpolatorLocal;
	typedef LagrangeInterpolator<ctype, mydim, cdim>    InterpolatorGlobal;

	typedef typename InterpolatorGlobal::LocalCoordinate   LocalCoordinate;
	typedef typename InterpolatorGlobal::GlobalCoordinate  GlobalCoordinate;
	typedef typename InterpolatorGlobal::Monomial          Monomial;
	typedef typename InterpolatorGlobal::LocalPolynomial   LocalPolynomial;
	typedef typename InterpolatorGlobal::PolynomialVector  PolynomialVector;
	typedef typename InterpolatorGlobal::PolynomialLocalCoordinate   PolyLocalCoord;
	typedef typename InterpolatorGlobal::PolynomialGlobalCoordinate  PolyGlobalCoord;

	typedef std::pair<LocalPolynomial, LocalPolynomial>  polyPair;

	typedef std::vector<LocalCoordinate> LocalCoordVec;
	typedef std::vector<LocalCoordVec>   LocalCoordVecVec;

	const unsigned int X_COORD = 0;
	const unsigned int Y_COORD = 1;


public:

	typedef std::pair<LocalCoordinate, GlobalCoordinate> ResultType;
	typedef std::vector<ResultType>                      ResultVector;


	/* \brief Constructor for singularity elimination. Pre-computes power expansion of x-coordinate
	 * \param[in]  curvOrder           Order of curvature of this entity
	 */
	CurvilinearSingularityHelper(unsigned int curvOrder) : curvOrder_(curvOrder) {
		gt_.makeTriangle();

		// Extract the corners of reference symplex
		pTriangle_ = LocalCoordVec {
			CurvilinearGeometryHelper::cornerInternalCoordinate<ctype, 2, cdim>(gt_, 0),
			CurvilinearGeometryHelper::cornerInternalCoordinate<ctype, 2, cdim>(gt_, 1),
			CurvilinearGeometryHelper::cornerInternalCoordinate<ctype, 2, cdim>(gt_, 2)
		};

		// Generate set of polynomials for (x^n - 1) / (x - 1)
		// Namely {1, 1+x, 1+x+x^2, 1+x+x^2+x^3, ...}
		// NOTE: Local coordinate is 2D, however, we only expand terms that have y=0
		// [TODO] Optimization: this operation can be done once for all singularities.
		//                      It makes sense to have a more generic constructor
		powerDiffExpansionX_ = PolynomialVector(curvOrder_);
		for (unsigned int i = 0; i < curvOrder_; i++)  { powerDiffExpansionX_[i] = LocalPolynomial(Monomial(1.0, i, 0)); }
		for (unsigned int i = 1; i < curvOrder_; i++)  { powerDiffExpansionX_[i] += powerDiffExpansionX_[i-1]; }
	}


	/** Algorithm-constructor:
	 *   1) Get PolyMap p(r), polyOrder, and singularity r0
	 *   2) Find if singularity is in the corner, side, or middle of ref. triangle
	 *   3) For each effective sub-triangle, generate map r = psi(r1), such that singularity is at r1=(1,0) for each sub-triangle.
	 *   4) For each effective sub-triangle, generate map p'(r1) = p(psi(r1)) by polynomial vector substitution
	 *   5) Calculate zeta(r1) = (p'(r1) - p'(1,0)) / (1 - r1x) = D(qx, 0) + qy * R(qx, qy(1-qx))
	 *        Where r1 = (qx, qy(1-qx)), and q are the square coordinates, mapped to the triangle coordinates r1
	 */

	/** \brief Constructor
	 * \param[in]  map                 The mapping from local to global coordinates
	 * \param[in]  singularityOrder    Order of singularity of the kind 1 / |r - r0|^n
	 * \param[in]  r0                  Local coordinate of the singularity
	 *
	 * */
	void init(
			PolyGlobalCoord & map,
			unsigned int singularityOrder,
			LocalCoordinate r0)
	{
		// A good polynomial map should be of the right size
		assert(map.size() == cdim);
		assert((singularityOrder == 1) || (singularityOrder == 2));

		// Generate point sets for each sub-triangle. Take care to
		// 1) Put the singularity @ 2nd corner
		// 2) New coordinate axis need to have correct orientation (non-negative cross-product)
		LocalCoordVecVec point {
			LocalCoordVec {pTriangle_[0], r0, pTriangle_[1]},
			LocalCoordVec {pTriangle_[1], r0, pTriangle_[2]},
			LocalCoordVec {pTriangle_[2], r0, pTriangle_[0]}
		};

		// If singularity is at side or corner, there will be less than 3 sub-triangles
		std::vector<unsigned int> triIndex;
		for (unsigned int i = 0; i < point.size(); i++)  { if(!isDegenerate(point[i]))  { triIndex.push_back(i); } }

		// Reserve space for resulting singularity-cancelled function, for each sub-triangle
		subTriInterp_ = std::vector<PolyLocalCoord>(triIndex.size());
		zetaD_        = std::vector<PolyGlobalCoord>(triIndex.size());
		zetaR_        = std::vector<PolyGlobalCoord>(triIndex.size());



		// [FIXME] Consider det(J) in 2nd coord transform




		for (unsigned int iTriangle = 0; iTriangle < triIndex.size(); iTriangle++)
		{
			subTriInterp_[iTriangle] = InterpolatorLocal(gt_, point[triIndex[iTriangle]], curvOrder_).interpolatoryVectorAnalytical();

			// generate composite map p'(r1) = p(psi(r1))
			PolyGlobalCoord compositeSubTriMap = map.composite(subTriInterp_[iTriangle]);

			localTriangleSingularityElimination(compositeSubTriMap, zetaD_[iTriangle], zetaR_[iTriangle]);


			// If have order 2, here more stuff
			if (singularityOrder == 2) {

				// Split cube into 2 triangles, one with singularity
				// Get

			}


		}
	}



	/**
	 * Algorithm-evaluate-O1-scalar
	 *   1) Get quadrature point q1
	 *   2) Calculate r1 = [q1x, q1y*(1-q1x) ]
	 *   3) Evaluate vector zeta(q) and return it
	 */

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
	ResultVector evalOrder1(LocalCoordinate q) const {

		LocalCoordinate r1;
		r1[0] = q[0];
		r1[1] = q[1] * (1 - q[0]);

		ResultVector rez(subTriInterp_.size());

		for (unsigned int i = 0; i < subTriInterp_.size(); i++)
		{
			GlobalCoordinate zeta  = zetaD_[i].evaluate(r1);
			GlobalCoordinate zetaR = zetaR_[i].evaluate(r1);
			zetaR *= q[1];
			zeta += zetaR;

			rez[i] = ResultType(subTriInterp_[i].evaluate(r1), zeta);
		}

		return rez;
	}



private:

	void localTriangleSingularityElimination(
			const PolyGlobalCoord & compositeSubTriMap,
			PolyGlobalCoord & D,
			PolyGlobalCoord & R
	) const {
		// Split into y-independent and y-proportional part
		// Namely p(x,y) = P0 + y*R = p(x,0) + y*R(x,y)
		PolyGlobalCoord P0(cdim);
		for (unsigned int iDim = 0; iDim < cdim; iDim++) {
			for (unsigned int iMon = 0; iMon < compositeSubTriMap[iDim].poly_.size(); iMon++) {

				// If monomial proportional to y, take y in front of the brackets, and add the remainder to zetaR
				// If monomial !proportional to y, add it to P0
				Monomial thisM = compositeSubTriMap[iDim].poly_[iMon];
				if (thisM.power_[Y_COORD] == 0) {  P0[iDim] += thisM; }
				else                            {
					thisM.power_[Y_COORD]--;
					R[iDim] += thisM;
				}
			}
		}

		// Perform analytical division of D = (p(x,0)-p(1,0))/(1-x)
		// This is effectively a sum of terms pref*((x^n - 1) / (x - 1))
		for (unsigned int iDim = 0; iDim < cdim; iDim++) {
			for (unsigned int iMon = 0; iMon < P0[iDim].poly_.size(); iMon++) {

				unsigned int xpow = P0[iDim].poly_[iMon].power_[X_COORD];
				ctype        pref = P0[iDim].poly_[iMon].pref_;

				// The constant term cancels out, for all other terms the power is decreased by 1 according to the (x^n - 1) / (x - 1) expansion
				if (xpow > 0) {
					LocalPolynomial summand = powerDiffExpansionX_[xpow - 1];
					summand *= pref * (-1);
					D += summand;
				}
			}
		}
	}

	// Determine if a 2D triangle is degenerate by its 3 coordinates
	bool isDegenerate(const std::vector<LocalCoordinate> & r) const {
		LocalCoordinate px = r[1];   px-= r[0];
		LocalCoordinate py = r[2];   py-= r[0];

		// Compute 2D cross product and check if it is non-zero
		ctype residual = px[0]*py[1] - py[0]*px[1];
		return (residual < 1.0e-10);
	}


	GeometryType  gt_;
	unsigned int curvOrder_;
	LocalCoordVec pTriangle_;
	PolynomialVector powerDiffExpansionX_;

	std::vector<PolyLocalCoord>   subTriInterp_;  // SubTriangle -> SubLocal(r1) -> Local(r).   AKA psi-map
	std::vector<PolyGlobalCoord>  zetaD_;         // SubTriangle -> SubLocal(r1) -> Global(x).  y-independent part of zeta
	std::vector<PolyGlobalCoord>  zetaR_;         // SubTriangle -> SubLocal(r1) -> Global(x).  y-dependent part of zeta


};


}





#endif // DUNE_CURVILINEAR_SINGULARITY_HELPER
