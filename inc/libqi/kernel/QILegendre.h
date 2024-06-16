// Things concerning parameterizations
// curve_param's, surface_param's and quad_inter's

#ifndef _qi_legendre_h_
#define _qi_legendre_h_

/** rpl */
#include <lidia/bigint_matrix.h>


using namespace rpl;

// Enter namespace QI
namespace QI {

// Computes an integer solution on a conic with integer coefficients
// Returns 0 if one such point has been found
lidia_size_t FindRationalPointOnConic(const bigint_matrix &q,
				      math_vector <bigint> &sol, std::ostream &s);

} // end of namespace QI

#endif
