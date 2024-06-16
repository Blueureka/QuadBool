// Intersection when the determinantal equation vanishes identically

#ifndef _qi_vanish_det_h_
#define _qi_vanish_det_h_

/** rpl */
#include <libqi/rpl/rational_factorization.h>

/** QI */
#include <libqi/kernel/QIHompoly.h>
#include <libqi/kernel/QIQsicStruct.h>

using namespace std;
using namespace rpl;

// Enter namespace QI
namespace QI {

// The main intersection procedure when the determinantal equation vanishes
quad_inter <bigint> inter_vanish_det(const bigint_matrix &q1, const bigint_matrix &q2, 
				     const hom_polynomial <bigint> &det_p, 
				     const int opt_level, 
				     ostream &s);

} // end of namespace QI

#endif
