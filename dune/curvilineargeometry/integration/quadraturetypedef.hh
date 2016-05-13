#ifndef DUNE_QUADRATURE_TYPEDEF
#define DUNE_QUADRATURE_TYPEDEF

namespace Dune
{

enum QuadratureNorm {
	QUADRATURE_NORM_DEFAULT = 0,       // Default norm to be used for scalars
	QUADRATURE_NORM_L1 = 1,
	QUADRATURE_NORM_L2 = 2,			   // Ratio of two-norms of error / value
	QUADRATURE_NORM_RELATIVE_L2 = 3,   // two-norm of vector of ratios (error_i / value_i)
	QUADRATURE_NORM_INDIVIDUAL = 4     // max of all ratios error_i / value_i
};

}

#endif  //DUNE_QUADRATURE_TYPEDEF
