// Special procedures for the ``generic'' algorithm, i.e. going through the param
// of a quadric of inertia (2,2)

#ifndef _qi_inter_with_22_
#define _qi_inter_with_22_

/** QI */
#include <libqi/kernel/QIHompoly.h>
#include <libqi/kernel/QIQsicStruct.h>

using namespace std;
using namespace rpl;

// Enter namespace QI
namespace QI {

// Input : initial quadrics q1, q2
//         determinental equation of the pencil : det_p
//         a quadric 22 : q
//         case_flag  = 0 if the intersection is a smooth quartic, 2 real affinely finite
//                    = 1           "                              1 real affinely finite
//                    = 2           "                              2 real affinely infinite
//                    = 3 if intersection is a cubic and a line secant
//                    = 4           "                           non-secant
//                    = 5 if intersection is two real skew lines
// Output : Compute the Parameterization s1 + sqrt(xi) s2 of a quadric 22 in the pencil 
//           and the Associated biquadratic equation poly1 + sqrt(xi) poly2 = 0
//    Call smooth_quartic(...) or cubic_line(...) for computing the intersection curve
quad_inter <bigint> intersection_with_quadric22(const bigint_matrix &q1, 
						const bigint_matrix &q2, 
						const hom_polynomial <bigint> &det_p, 
						const bigint_matrix &q, 
						const unsigned int case_flag, 
						const int opt_level, 
						ostream &s);

// Input : initial quadrics q1, q2
//         determinental equation of the pencil : det_p
//         a quadric 22  q
// Output : a quadric 22  quadric_sol and a point on it point_sol
//          the determinant a^2.b of the quadric 22 in the form [a,b]
void find_quadric_22_and_point(const bigint_matrix &q1, const bigint_matrix &q2, 
			       const hom_polynomial <bigint> &det_p, 
			       const bigint_matrix &q_cp, 
			       bigint_matrix &quadric_sol, 
			       math_vector <bigint> &point_sol, 
			       bigint &det_R,
			       math_vector <bigint> &det_sol, 
			       const int opt_level, ostream &s);

// Input:  a matrix q of inertia 2,2 
//   a point  on it and the determinant of q in the form [a,b] such that det(q)=a^2. b
// Output : Parameterization of the quadric :   s1 + sqrt(det_q[1]) s2
// If det_q[1]=1 then  s2=[0,0,0,0]
void find_parameterization_quadric_22(const math_matrix <bigint> &q_cp, 
				      const math_vector <bigint> &point_on_q, 
				      const math_vector <bigint> &det_q_cp,
				      const bigint &det_R,
				      surface_param <bigint> &s1, 
				      surface_param <bigint> &s2,
				      const int opt_level, ostream &s);

} // end of namespace QI

#endif
