#ifndef DUNE_OPTIMIZATION_
#define DUNE_OPTIMIZATION_


template <class CurvGrid>
struct EnergyFunctionL2Squared
{
	static const int coorddimension = CurvGrid::coorddimension;

	typedef typename CurvGrid::LocalPolynomial             LocalPolynomial;
    typedef typename CurvGrid::GlobalCoordinate            GlobalCoordinate;
	typedef typename CurvGrid::PolynomialGlobalCoordinate  PolyMap;

	static LocalPolynomial eval(const GlobalCoordinate & global, const PolyMap & polymap)
	{
		LocalPolynomial rez;

		for (int i = 0; i < CurvGrid::coorddimension; i++)
		{
			LocalPolynomial tmp = polymap[i] - global[i];
			rez += tmp * tmp;
		}
	}
};


template <class CurvGrid>
struct EnergyFunctionCoordinate
{
	static const int coorddimension = CurvGrid::coorddimension;

	typedef typename CurvGrid::LocalPolynomial             LocalPolynomial;
    typedef typename CurvGrid::GlobalCoordinate            GlobalCoordinate;
	typedef typename CurvGrid::PolynomialGlobalCoordinate  PolyMap;

	static LocalPolynomial eval(int coord, const PolyMap & polymap)
	{
		assert((coord > 0)&&(coord < coorddimension));
		return polymap[coord];
	}
};


// Implement method to newton-minimize constrained for element for functor


#endif // DUNE_OPTIMIZATION_
