#ifndef DUNE_QUADRATURE_ABSOLUTE_VALUE
#define DUNE_QUADRATURE_ABSOLUTE_VALUE

#include<dune/curvilineargeometry/integration/quadraturetypedef.hh>

namespace Dune
{

// ***************************************************
// Calculates relative 1-norm of a scalar
// ***************************************************


template<class ctype, int NormType>
struct QuadratureAbsoluteValue {

	template <class ValueType>
	static ctype eval(ValueType & v) {	DUNE_THROW(NotImplemented, "Trying to call generic method"); }

	// FieldVector
	template<class ValueType, int mydim>
	static ctype eval( Dune::FieldVector<ValueType, mydim> & v) {
		DUNE_THROW(NotImplemented, "Trying to call generic method");
	}

	// DynamicVector
	template<class ValueType>
	static ctype eval( Dune::DynamicVector<ValueType> & v) {
		DUNE_THROW(NotImplemented, "Trying to call generic method");
	}

	// DynamicMatrix
	template<class ValueType>
	static ctype eval( Dune::DynamicMatrix<ValueType> & v) {
		DUNE_THROW(NotImplemented, "Trying to call generic method");
	}

};


// ***************************************************
// Calculates relative 1-norm for scalars
// ***************************************************


template<class ctype>
struct QuadratureAbsoluteValue<ctype, QUADRATURE_NORM_DEFAULT> {
	template <class ValueType>
	static ctype eval(ValueType & v) { return std::abs(v); }
};


// ***************************************************
// Calculates relative infinity-norm
// ***************************************************


template<class ctype>
struct QuadratureAbsoluteValue<ctype, QUADRATURE_NORM_INDIVIDUAL> {

	// Scalar
	template <class ValueType>
	static ctype eval( ValueType & v) {
		return QuadratureAbsoluteValue<ctype, QUADRATURE_NORM_DEFAULT>::eval(v);
	}

	// FieldVector
	template<class ValueType, int mydim>
	static ctype eval( Dune::FieldVector<ValueType, mydim> & v) {
		ctype rez = 0;
		for (int i = 0; i < mydim; i++) {
			ctype tmp = QuadratureAbsoluteValue<ctype, QUADRATURE_NORM_DEFAULT>::eval(v[i]);
			rez = std::max(rez, tmp);
		}
		return rez;
	}

	// DynamicVector
	template<class ValueType>
	static ctype eval( Dune::DynamicVector<ValueType> & v) {
		ctype rez = 0;

		for (int i = 0; i < v.N(); i++) {
			ctype tmp = QuadratureAbsoluteValue<ctype, QUADRATURE_NORM_DEFAULT>::eval(v[i]);
			rez = std::max(rez, tmp);
		}
		return rez;
	}

	// DynamicMatrix
	template<class ValueType>
	static ctype eval( Dune::DynamicMatrix<ValueType> & v) {
		ctype rez = 0;
		for (int i = 0; i < v.N(); i++) {
			for (int j = 0; j < v.M(); j++) {
				ctype tmp = QuadratureAbsoluteValue<ctype, QUADRATURE_NORM_DEFAULT>::eval(v[i][j]);
				rez = std::max(rez, tmp);
			}
		}
		return rez;
	}

};


// ***************************************************
// Calculates relative 2-norm
// ***************************************************


template<class ctype>
struct QuadratureAbsoluteValue<ctype, QUADRATURE_NORM_L2> {

	// Scalar
	template <class ValueType>
	static ctype eval( ValueType & v) {
		return QuadratureAbsoluteValue<ctype, QUADRATURE_NORM_DEFAULT>::eval(v);
	}

	// FieldVector
	template<class ValueType, int mydim>
	static ctype eval( Dune::FieldVector<ValueType, mydim> & v) { return v.two_norm(); }

	// DynamicVector
	template<class ValueType>
	static ctype eval( Dune::DynamicVector<ValueType> & v) { return v.two_norm(); }

	// DynamicMatrix
	template<class ValueType>
	static ctype eval( Dune::DynamicMatrix<ValueType> & v) { return v.frobenius_norm(); }

};






}





#endif  //DUNE_QUADRATURE_ABSOLUTE_VALUE
