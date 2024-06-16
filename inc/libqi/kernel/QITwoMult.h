// Intersection when the determinantal equation has two double roots

#ifndef _qi_two_mult_h_
#define _qi_two_mult_h_

/** QI */
#include <libqi/kernel/QIHompoly.h>
#include <libqi/kernel/QIQsicStruct.h>

using namespace std;
using namespace rpl;

// Enter namespace QI
namespace QI {

// Main procedure when the determinantal equation has two double roots
quad_inter <bigint> inter_two_mult(const bigint_matrix &q1, const bigint_matrix &q2, 
				   const hom_polynomial <bigint> &det_p, 
				   const hom_polynomial <bigint> &det_p_orig,
				   const hom_polynomial <bigint> &gcd_p, 
				   const int opt_level, std::ostream &s);

} // end of namespace QI

#endif
