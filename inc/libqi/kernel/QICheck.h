// Things for checking whether computed params are alright

#ifndef _qi_check_h_
#define _qi_check_h_

/** rpl   */
#include <libqi/rpl/bigint.h>
#include <libqi/rpl/bigint_matrix.h>

/** QI */
#include "QIParamStruct.h"
#include "QIParam.h"

using namespace std;
using namespace rpl;

// Enter namespace QI
namespace QI {

// Checks if the computed param is ok by replugging in initial quadrics
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
		 const surface_param <bigint> &p, ostream &s);

// Checks if the computed param is ok by replugging in initial quadrics
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
		 const curve_param <bigint> &p, ostream &s);

// Checks if p1+sqrt(D)*p2 is ok
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
		 const curve_param <bigint> &p1, const curve_param <bigint> &p2, 
		 const bigint &D, ostream &s);

// Checks if p1+sqrt(D1)*p2+sqrt(D2)*p3 is ok
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
		 const curve_param <bigint> &p1, const curve_param <bigint> &p2, 
		 const curve_param <bigint> &p3, 
		 const bigint &D1, const bigint &D2, ostream &s);

// Checks if param c1+sqrt(a)*c2+sqrt(b)*c3+sqrt(a*b)*c4 is ok
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
		 const curve_param <bigint> &c1, const curve_param <bigint> &c2, 
		 const curve_param <bigint> &c3, const curve_param <bigint> &c4,
 		 const bigint &a, const bigint &b, ostream &s);

// Checks if param c1+sqrt(a)*c2+(c3+sqrt(a)*c4)*sqrt(b+c*sqrt(a)) is ok
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
		 const curve_param <bigint> &c1, const curve_param <bigint> &c2, 
		 const curve_param <bigint> &c3, const curve_param <bigint> &c4,
		 const bigint &a, const bigint &b, const bigint &c, ostream &s);

// Checks if the computed param is ok by replugging in initial quadrics
// Second case, when instead we have a surface param and a polynomial
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
		 const surface_param <bigint> &p, 
		 const hom_polynomial <bigint> &pol, ostream &s);

// Checks if the parameterization s1 + sqrt(det_q[1]) s2 of the quadric 22 q lies on q
void check_param(const bigint_matrix &q, const surface_param <bigint> &s1, 
		 const surface_param <bigint> &s2, const bigint &delta, ostream &s);

// Checks if the computed quartic param c is ok by replugging in initial quadrics
//   c = c1 + sqrt(d). c2 + eps. sqrt(Delta). (c3 + sqrt(d). c4)
// Delta = Delta1 + sqrt(d). Delta2
// Note if check_param is ok for c then it is also ok after replacing 
// sqrt(d) by -sqrt(d) and sqrt(Delta) by -sqrt(Delta)
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
		 const curve_param <bigint> &c1, const curve_param <bigint> &c2, 
		 const curve_param <bigint> &c3, const curve_param <bigint> &c4, 
		 const bigint &delta, const hom_polynomial <bigint> &Delta1, 
		 const hom_polynomial <bigint> &Delta2, ostream &s); 

} // end of namespace QI

#endif
