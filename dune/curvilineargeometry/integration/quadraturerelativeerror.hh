#ifndef DUNE_QUADRATURE_RELATIVE_ERROR
#define DUNE_QUADRATURE_RELATIVE_ERROR

namespace Dune
{

enum QuadratureNorm {
	QUADRATURE_NORM_L1,
	QUADRATURE_NORM_L2,			// Ratio of two-norms of error / value
	QUADRATURE_NORM_RELATIVE_L2,   // two-norm of vector of ratios (error_i / value_i)
	QUADRATURE_NORM_INDIVIDUAL     // max of all ratios error_i / value_i
};



// ***************************************************
// Calculates relative 1-norm of a scalar
// ***************************************************


template<class ctype, class ValueType, int NormType>
struct QuadratureRelativeError {
	static ctype eval(
		ValueType & delta,     // The difference between two values used as absolute error estimate
		ValueType & mag,       // A value used as a magnitude estimate
		ctype ACCURACY_GOAL)   // This number determines how small an absolute number should be to be considered effectively zero
	{
		ctype abs_d = std::abs(delta);
		ctype abs_m = std::abs(mag);

		// Only calculate relative error if the absolute value of the actual quantity is large enough
		// Otherwise just return absolute error
		return (abs_m > ACCURACY_GOAL) ? abs_d / abs_m : abs_d;
	}
};


// ***************************************************
// Calculates relative norm of a DynamicVector
// ***************************************************


template<class ctype, class ValueType, int NormType>
struct QuadratureRelativeError<ctype, Dune::DynamicVector<ValueType>, NormType>
{
	typedef Dune::DynamicVector<ValueType>  ResultType;

	static ctype eval( ResultType & delta, ResultType & mag, ctype ACCURACY_GOAL) {  }
};


template<class ctype, class ValueType>
struct QuadratureRelativeError<ctype, Dune::DynamicVector<ValueType>, QUADRATURE_NORM_INDIVIDUAL>
{
	typedef Dune::DynamicVector<ValueType>  ResultType;

	static ctype eval( ResultType & delta, ResultType & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;
		for (int i = 0; i < delta.N(); i++) { rez = std::max(rez, QuadratureRelativeError<ctype, ValueType, 0>::eval(delta[i], mag[i], ACCURACY_GOAL)); }
		return rez;
	}
};


template<class ctype, class ValueType>
struct QuadratureRelativeError<ctype, Dune::DynamicVector<ValueType>, QUADRATURE_NORM_RELATIVE_L2>
{
	typedef Dune::DynamicVector<ValueType>  ResultType;

	static ctype eval( ResultType & delta, ResultType & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;
		for (int i = 0; i < delta.N(); i++) {
			ctype tmp = QuadratureRelativeError<ctype, ValueType, 0>::eval(delta[i], mag[i], ACCURACY_GOAL);
			rez += tmp * tmp;
		}
		return sqrt(rez / delta.N());
	}
};


template<class ctype, class ValueType>
struct QuadratureRelativeError<ctype, Dune::DynamicVector<ValueType>, QUADRATURE_NORM_L2>
{
	typedef Dune::DynamicVector<ValueType>  ResultType;

	static ctype eval( ResultType & delta, ResultType & mag, ctype ACCURACY_GOAL) {
		return delta.two_norm() / mag.two_norm();
	}
};


// ***************************************************
// Calculates relative one-norm of a DynamicMatrix
// ***************************************************

template<class ctype, class ValueType, int NormType>
struct QuadratureRelativeError<ctype, Dune::DynamicMatrix<ValueType>, NormType>
{
	typedef Dune::DynamicMatrix<ValueType>  ResultType;

	static ctype eval( ResultType & delta, ResultType & mag, ctype ACCURACY_GOAL) {  }
};


template<class ctype, class ValueType>
struct QuadratureRelativeError<ctype, Dune::DynamicMatrix<ValueType>, QUADRATURE_NORM_INDIVIDUAL>
{
	typedef Dune::DynamicMatrix<ValueType>  ResultType;

	static ctype eval( ResultType & delta, ResultType & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;
		for (int i = 0; i < delta.N(); i++) {
			for (int j = 0; j < delta.M(); j++) {
				rez = std::max(rez, QuadratureRelativeError<ctype, ValueType, 0>::eval(delta[i][j], mag[i][j], ACCURACY_GOAL));
			}
		}
		return rez;
	}
};


template<class ctype, class ValueType>
struct QuadratureRelativeError<ctype, Dune::DynamicMatrix<ValueType>, QUADRATURE_NORM_RELATIVE_L2>
{
	typedef Dune::DynamicMatrix<ValueType>  ResultType;

	static ctype eval( ResultType & delta, ResultType & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;
		for (int i = 0; i < delta.N(); i++) {
			for (int j = 0; j < delta.M(); j++) {
				ctype tmp = QuadratureRelativeError<ctype, ValueType, 0>::eval(delta[i][j], mag[i][j], ACCURACY_GOAL);
				rez = tmp * tmp;
			}
		}
		return sqrt(rez / (delta.N() * delta.M()));
	}
};


template<class ctype, class ValueType>
struct QuadratureRelativeError<ctype, Dune::DynamicMatrix<ValueType>, QUADRATURE_NORM_L2>
{
	typedef Dune::DynamicMatrix<ValueType>  ResultType;

	static ctype eval( ResultType & delta, ResultType & mag, ctype ACCURACY_GOAL) {
		return delta.frobenius_norm() / mag.frobenius_norm();
	}
};




}





#endif  //DUNE_QUADRATURE_RELATIVE_ERROR
