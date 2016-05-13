#ifndef DUNE_QUADRATURE_RELATIVE_ERROR
#define DUNE_QUADRATURE_RELATIVE_ERROR

#include<dune/curvilineargeometry/integration/quadraturetypedef.hh>

namespace Dune
{

// ***************************************************
// Calculates relative 1-norm of a scalar
// ***************************************************


template<class ctype, int NormType>
struct QuadratureRelativeError {

	template <class ValueType>
	static ctype eval( ValueType & delta, ValueType & mag, ctype ACCURACY_GOAL) {	DUNE_THROW(NotImplemented, "Trying to call generic method"); }

	// FieldVector
	template<class ValueType, int mydim>
	static ctype eval( Dune::FieldVector<ValueType, mydim> & delta, Dune::FieldVector<ValueType, mydim> & mag, ctype ACCURACY_GOAL) {
		DUNE_THROW(NotImplemented, "Trying to call generic method");
	}

	// DynamicVector
	template<class ValueType>
	static ctype eval( Dune::DynamicVector<ValueType> & delta, Dune::DynamicVector<ValueType> & mag, ctype ACCURACY_GOAL) {
		DUNE_THROW(NotImplemented, "Trying to call generic method");
	}

	// DynamicMatrix
	template<class ValueType>
	static ctype eval( Dune::DynamicMatrix<ValueType> & delta, Dune::DynamicMatrix<ValueType> & mag, ctype ACCURACY_GOAL) {
		DUNE_THROW(NotImplemented, "Trying to call generic method");
	}

};


// ***************************************************
// Calculates relative 1-norm for scalars
// ***************************************************


template<class ctype>
struct QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT> {

	template <class ValueType>
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
// Calculates relative infinity-norm
// ***************************************************


template<class ctype>
struct QuadratureRelativeError<ctype, QUADRATURE_NORM_INDIVIDUAL> {

	// Scalar
	template <class ValueType>
	static ctype eval( ValueType & delta, ValueType & mag, ctype ACCURACY_GOAL) {
		return QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(delta, mag, ACCURACY_GOAL);
	}

	// FieldVector
	template<class ValueType, int mydim>
	static ctype eval( Dune::FieldVector<ValueType, mydim> & delta, Dune::FieldVector<ValueType, mydim> & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;
		for (int i = 0; i < mydim; i++) {
			ctype tmp = QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(delta[i], mag[i], ACCURACY_GOAL);
			rez = std::max(rez, tmp);
		}
		return rez;
	}

	// DynamicVector
	template<class ValueType>
	static ctype eval( Dune::DynamicVector<ValueType> & delta, Dune::DynamicVector<ValueType> & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;

		for (int i = 0; i < delta.N(); i++) {
			ctype tmp = QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(delta[i], mag[i], ACCURACY_GOAL);
			rez = std::max(rez, tmp);
		}
		return rez;
	}

	// DynamicMatrix
	template<class ValueType>
	static ctype eval( Dune::DynamicMatrix<ValueType> & delta, Dune::DynamicMatrix<ValueType> & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;
		for (int i = 0; i < delta.N(); i++) {
			for (int j = 0; j < delta.M(); j++) {
				ctype tmp = QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(delta[i][j], mag[i][j], ACCURACY_GOAL);
				rez = std::max(rez, tmp);
			}
		}
		return rez;
	}

};





// ***************************************************
// Calculates 2-norm of relative vectors
// ***************************************************


template<class ctype>
struct QuadratureRelativeError<ctype, QUADRATURE_NORM_RELATIVE_L2> {

	// Scalar
	template <class ValueType>
	static ctype eval( ValueType & delta, ValueType & mag, ctype ACCURACY_GOAL) {
		return QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(delta, mag, ACCURACY_GOAL);
	}

	// FieldVector
	template<class ValueType, int mydim>
	static ctype eval( Dune::FieldVector<ValueType, mydim> & delta, Dune::FieldVector<ValueType, mydim> & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;
		for (int i = 0; i < mydim; i++) {
			ctype tmp = QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(delta[i], mag[i], ACCURACY_GOAL);
			rez += tmp * tmp;
		}
		return sqrt(rez / mydim);
	}

	// DynamicVector
	template<class ValueType>
	static ctype eval( Dune::DynamicVector<ValueType> & delta, Dune::DynamicVector<ValueType> & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;
		for (int i = 0; i < delta.N(); i++) {
			ctype tmp = QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(delta[i], mag[i], ACCURACY_GOAL);
			rez += tmp * tmp;
		}
		return sqrt(rez / delta.N());
	}

	// DynamicMatrix
	template<class ValueType>
	static ctype eval( Dune::DynamicMatrix<ValueType> & delta, Dune::DynamicMatrix<ValueType> & mag, ctype ACCURACY_GOAL) {
		ctype rez = 0;
		for (int i = 0; i < delta.N(); i++) {
			for (int j = 0; j < delta.M(); j++) {
				ctype tmp = QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(delta[i][j], mag[i][j], ACCURACY_GOAL);
				rez = tmp * tmp;
			}
		}
		return sqrt(rez / (delta.N() * delta.M()));
	}

};



// ***************************************************
// Calculates relative 2-norm
// ***************************************************


template<class ctype>
struct QuadratureRelativeError<ctype, QUADRATURE_NORM_L2> {

	// Scalar
	template <class ValueType>
	static ctype eval( ValueType & delta, ValueType & mag, ctype ACCURACY_GOAL) {
		return QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(delta, mag, ACCURACY_GOAL);
	}

	// FieldVector
	template<class ValueType, int mydim>
	static ctype eval( Dune::FieldVector<ValueType, mydim> & delta, Dune::FieldVector<ValueType, mydim> & mag, ctype ACCURACY_GOAL) {
		ctype normDelta = delta.two_norm();
		ctype normMag = mag.two_norm();
		return QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(normDelta, normMag, ACCURACY_GOAL);
	}

	// DynamicVector
	template<class ValueType>
	static ctype eval( Dune::DynamicVector<ValueType> & delta, Dune::DynamicVector<ValueType> & mag, ctype ACCURACY_GOAL) {
		ctype normDelta = delta.two_norm();
		ctype normMag = mag.two_norm();
		return QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(normDelta, normMag, ACCURACY_GOAL);
	}

	// DynamicMatrix
	template<class ValueType>
	static ctype eval( Dune::DynamicMatrix<ValueType> & delta, Dune::DynamicMatrix<ValueType> & mag, ctype ACCURACY_GOAL) {
		ctype normDelta = delta.frobenius_norm();
		ctype normMag = mag.frobenius_norm();
		return QuadratureRelativeError<ctype, QUADRATURE_NORM_DEFAULT>::eval(normDelta, normMag, ACCURACY_GOAL);
	}

};






}





#endif  //DUNE_QUADRATURE_RELATIVE_ERROR
