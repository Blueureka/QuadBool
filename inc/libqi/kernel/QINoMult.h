// Intersection when the determinantal equation has no multiple root

#ifndef _qi_no_mult_h_
#define _qi_no_mult_h_

/** QI */
#include <libqi/kernel/QIQsicStruct.h>

using namespace std;
using namespace rpl;

// Enter namespace QI
namespace QI {

quad_inter <bigint> inter_no_mult(const bigint_matrix &q1, const bigint_matrix &q2, 
				  const hom_polynomial <bigint> &det_p, 
				  const hom_polynomial <bigint> &det_p_orig, 
				  const int opt_level, ostream &s);

} // end of namespace QI

#endif
