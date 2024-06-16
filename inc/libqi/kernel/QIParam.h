// Things concerning parameterizations

#ifndef _qi_param_h_
#define _qi_param_h_

/** QI */
#include <libqi/kernel/QIHompoly.h>
#include <libqi/kernel/QIParamStruct.h>

using namespace std;

// Enter namespace QI
namespace QI {

// Output the parameterization of the quadric of inertia (2,2) of implicit
// equation xy - zw = 0
void surface_22_param(surface_param <bigint> &s);

// Print the bi-quadratic homogeneous polynomial 
void print_bi_quadratic_hom_poly(const hom_hom_polynomial <bigint> &poly1, 
				 const hom_hom_polynomial <bigint> &poly2, 
				 const bigint &delta, ostream &s);

// Plug a param in the equation of a quadric
hom_hom_polynomial <bigint> plug_param_in_quadric(const surface_param <bigint> &s1,
						  const bigint_matrix &q, 
						  const surface_param <bigint> &s2);

// Plug a curve_param in the equation of a quadric
hom_polynomial <bigint> plug_param_in_quadric(const curve_param <bigint> &s1,
					      const bigint_matrix &q, 
					      const curve_param <bigint> &s2);

// Compute the parameterization of a non-singular conic through a rational point
void conic_param_through_ratpoint(const bigint_matrix &q, 
				  const math_vector <bigint> &rat_point,
				  bigint_matrix &p_trans, curve_param <bigint> &p,
				  math_vector <bigint> &l, ostream &s);

// Compute the parameterization of a non-singular conic through the point [0,0,1]
void conic_param_through_origin(const bigint_matrix &q, bigint_matrix &p_trans, 
				curve_param <bigint> &p, math_vector <bigint> &l);

// Compute the parameterization of a cone through a rational point
void cone_param_through_ratpoint(const bigint_matrix &q, 
				 const math_vector <bigint> &sing,
				 const math_vector <bigint> &rat_point, 
				 bigint_matrix &p_trans, surface_param <bigint> &p, 
				 math_vector <bigint> &l, ostream &s);

// Compute the parameterization of 2 x 2 matrix, used for lines and planes
void two_by_two_param(const bigint_matrix &q, bigint &D, bigint_matrix &m1, 
		      bigint_matrix &m2, const int opt_level, ostream &s);

// Compute the parameterization of a pair of lines in the plane
void pair_of_lines_param(const bigint_matrix &q, const math_vector <bigint> &sing, 
			 bigint &D, bigint_matrix &m1, bigint_matrix &m2, 
			 const int opt_level, ostream &s);

// Compute the parameterization of a double line
void double_line_param(const bigint_matrix &q, const bigint_matrix &sing,
		       bigint_matrix &m1, ostream &s);

// Compute the parameterization of a double plane
void double_plane_param(const bigint_matrix &q, const bigint_matrix &sing,
			bigint_matrix &m1, ostream &s);

// Compute the parameterization of a pair of planes
void pair_of_planes_param(const bigint_matrix &q, const bigint_matrix &sing, bigint &D, 
			  bigint_matrix &m1, bigint_matrix &m2, 
			  const int opt_level, ostream &s);

// The procedure for parameterizing conics
// The param is (p_trans + sqrt(D)*p_trans2) * [u^2 v^2 uv]
// par = [u^2 v^2 uv]
void conic_param(const bigint_matrix &q, bigint &D, 
		 bigint_matrix &p_trans, bigint_matrix &p_trans2, 
		 curve_param <bigint> &par, const int opt_level, ostream &s);

// The procedure for parameterizing a conics whose corresponding 3x3 matrix is 
// C1 +sqrt(delta).C2 when no rational point is known
// The param is 
//      (p_trans0 + sqrt(delta)*p_trans1 + sqrt(D)*p_trans2 
//                + sqrt(delta)*sqrt(D)*p_trans3) * [u^2 v^2 uv]
// par = [u^2 v^2 uv]
// D = D[0] + sqrt(delta)*D[1]
void non_rational_conic_param(const bigint_matrix &C1, const bigint_matrix &C2, 
			      const bigint &delta, math_vector <bigint> &D, 
			      bigint_matrix &p_trans0, bigint_matrix &p_trans1, 
			      bigint_matrix &p_trans2, bigint_matrix &p_trans3,
			      curve_param <bigint> &par, const int opt_level,
			      ostream &s);

// Compute the parameterization of a cone
void cone_param(const bigint_matrix &q, const math_vector <bigint> &sing, bigint &D,
		bigint_matrix &p_trans, bigint_matrix &p_trans2,
		surface_param <bigint> &p, const int opt_level, ostream &s);

// Try to reparameterize a rational line by taking as "endpoints" (ie. the points
// (u,v) = (1,0) and (0,1)) points on the planes x = 0, y = 0, z = 0, w = 0. 
// In other words, one of the coordinates is 0 when u=0 and another is zero when v=0
// Returns the value of the new parameters corresponding to the first endpoint of
// the old parameterization
math_vector <bigint> improved_line_param(curve_param <bigint> &line_par);

// Try to reparameterize a line c = c1+sqrt(xi). c2 such that 
// one of the coordinates is 0 when u=0 and another is zero when v=0
math_vector <bigint> improved_line_param(curve_param <bigint> &c1, 
					 curve_param <bigint> &c2, const bigint &xi);

// Try to reparameterize a line c = c1+sqrt(xi). c2 + sqrt(D). c3 + sqrt(D.xi). c4
// such that one of the coordinates is 0 when u=0 and another is zero when v=0
// (D = D1 + sqrt(xi). D2)
void improved_line_param(curve_param <bigint> &c1, curve_param <bigint> &c2, 
			 curve_param <bigint> &c3, curve_param <bigint> &c4, 
			 const bigint &xi, const bigint &D1, const bigint &D2);

} // end of namespace QI

#endif

