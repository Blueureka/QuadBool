// Number theory related procedures

#ifndef _qi_number_h_
#define _qi_number_h_

/** QI */
/** $$ Redondance */
#include "QIHompoly.h"
#include "QIParam.h"

using namespace std;

// Enter namespace QI
namespace QI {

// Optimize coefficients of polynomial using the pp function
void optimize(hom_polynomial <bigint> &pol);

// Optimize coefficients of a pair of hom_polynomials by taking the gcd of their coeff
void optimize(hom_polynomial <bigint> &pol1, hom_polynomial <bigint> &pol2); 

// Optimize coefficients of a pair of hom_hom_polynomials by taking the gcd of the
// contents of the coefficients (we only eliminate constants, not polynomial factors)
void optimize(hom_hom_polynomial <bigint> &pol1, hom_hom_polynomial <bigint> &pol2);

// Output the gcd of the content of the coordinates of a curve_param
bigint cont(const curve_param <bigint> &par);

// Compute the content of a hom_hom_polynomial, ie the gcd of the coefficients in (u,v)
hom_polynomial <bigint> cont(const hom_hom_polynomial <bigint> &a);

// Compute the primitive part  of a hom_hom_polynomial a, i.e. a divided by its content
hom_hom_polynomial <bigint> pp(const hom_hom_polynomial <bigint> &a);

// Compute the primitive part of a hom_hom_polynomial poly knowing its content
hom_hom_polynomial <bigint> pp(const hom_hom_polynomial <bigint> &a, 
			       const hom_polynomial <bigint> &content);

// Optimize a curve_param
void optimize(curve_param <bigint> &par);

// Optimize a pair of curve_params by taking the gcd of their content
void optimize(curve_param <bigint> &par1, curve_param <bigint> &par2);

// Message for factor extraction
void extract_message(const int opt_level, ostream &s, 
		     const string &f = "extracting square factors");

// Take the integer nearest to the cubic root of an integer
bigint cubic_root(const bigint &b);

// Compute the square factors of a bigint
math_vector <bigint> extract_square_factors(const bigint &b, const int opt_level, 
					    ostream &s);

// Content of the coefficients of a vector
bigint content(const math_vector <bigint> &vec);

// Optimize a vector by dividing by its content
void optimize(math_vector <bigint> &v1);

// Optimize a pair of vectors by dividing by the gcd of the contents
void optimize(math_vector <bigint> &v1, math_vector <bigint> &v2);

// Optimize first two, then last two coordinates
void optimize_by_half(math_vector <bigint> &v1);

// Optimize a triple of vectors by dividing by the gcd of the contents
void optimize(math_vector <bigint> &v1, math_vector <bigint> &v2, 
	      math_vector <bigint> &v3);

// Optimize a quadruple of vectors by dividing by the gcd of the contents
void optimize(math_vector <bigint> &v1, math_vector <bigint> &v2, 
	      math_vector <bigint> &v3, math_vector <bigint> &v4);

// Optimize by dividing each column by its content
void optimize_trans1(bigint_matrix &q);

// Optimize: each column should be divided by the gcd of the contents of each
// matrix column
void optimize_trans2(bigint_matrix &q1, bigint_matrix &q2);

// Optimize by dividing columns 1, 2, 3 by a, b, c such that a b = c^2
void optimize_trans3(bigint_matrix &q, const int opt_level, ostream &s);

// Same, but with additional line values to be updated
void optimize_trans3(bigint_matrix &q, math_vector <bigint> &l, const int opt_level, 
		     ostream &s);

// Load balancing of columns of a matrix such that a*b = c*d
void load_balancing(bigint_matrix &q, const bool &whatcase, ostream &s);

// Load balancing of vectors of a matrix such that a*b = c*d
bool load_balancing(const math_vector<bigint> &v1, const math_vector<bigint> &v2,
			bigint &mult, ostream &s);


// Optimize by dividing columns 1, 2, 3 by a, b, c such that a b = c^2
// But here the two matrices should be divided simultaneously
void optimize_trans4(bigint_matrix &q1, bigint_matrix &q2, const int opt_level, 
		     ostream &s);

// Optimize the parameterization of a smooth quartic c1 + sqrt(xi). c2 + 
//     eps. sqrt(Delta). (c3 + sqrt(xi). c4)
// Simplification by the common factor of the parameterization
void optimize(curve_param <bigint> &c1, curve_param <bigint> &c2, 
	      curve_param <bigint> &c3, curve_param <bigint> &c4, 
	      hom_polynomial <bigint> &Delta1, hom_polynomial <bigint> &Delta2, 
	      const int opt_level, ostream &s);

// Optimize the parameterization of a conic c1 + sqrt(xi). c2 + eps. sqrt(D). (c3
// + sqrt(xi). c4) where D and xi are already optimized
// Simplification by the common factor of the parameterization
void optimize(curve_param <bigint> &c1, curve_param <bigint> &c2, 
	      curve_param <bigint> &c3, curve_param <bigint> &c4, ostream &s);

// Compute the pseudo-discriminant of a quadratic homogeneous polynomial
bigint discriminant2(const hom_polynomial <bigint> &pol);

} // end of namespace QI

#endif
