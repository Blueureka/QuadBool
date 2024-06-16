// Main intersection loop

#ifndef _qi_inter_h_ 
#define _qi_inter_h_

/** rpl */
#include <libqi/rpl/bigint_matrix.h>

/** QI */
#include <libqi/kernel/QIQsicStruct.h>

using namespace std;
using namespace rpl;

// Enter namespace QI
namespace QI {

// The main intersection procedure
quad_inter <bigint> intersection(const bigint_matrix &q1, const bigint_matrix &q2, 
				 const int opt_level, ostream &s);

} // end of namespace QI

#endif
