// Number theory related procedures
#include <libqi/rpl/rational_factorization.h>
#include <libqi/rpl/bigfloat.h>
#include <libqi/kernel/QINumber.h>
#include <libqi/kernel/QIGlobal.h>
#include <libqi/rpl/polynomial.h>
#include <libqi/rpl/rpl.h>

// Enter namespace QI
namespace QI {

// Optimize coefficients of polynomial
// I got rid of the pp (primitive part) function because it does not like
// zero coefficients...
void optimize(hom_polynomial <bigint> &pol)
{
	pol = pol/cont(pol);
}

// Optimize coefficients of a pair of hom_polynomials by taking the gcd of their coefficients
void optimize(hom_polynomial <bigint> &pol1, hom_polynomial <bigint> &pol2)
{
	bigint tmp = 0;
	for (rpl_size_t i = 0; i <= pol1.degree(); i++)
		tmp = gcd(tmp, pol1[i]);
	for (rpl_size_t i = 0; i <= pol2.degree(); i++)
		tmp = gcd(tmp, pol2[i]);
	if (tmp > 1)
		{
			for (rpl_size_t i = 0; i <= pol1.degree(); i++)
				div_exact(pol1[i], pol1[i], tmp);
			for (rpl_size_t i = 0; i <= pol2.degree(); i++)
				div_exact(pol2[i], pol2[i], tmp);
		}
}

// Optimize coefficients of a pair of hom_hom_polynomials by taking the gcd of the
// contents of the coefficients (we only eliminate constants, not polynomial factors)
void optimize(hom_hom_polynomial <bigint> &pol1, hom_hom_polynomial <bigint> &pol2)
{
	bigint coef = 0;
	for (rpl_size_t i = 0; i <= pol1.degree(); i++)
		coef = gcd(coef,cont(pol1[i]));
	for (rpl_size_t i = 0; i <= pol2.degree(); i++)
		coef = gcd(coef,cont(pol2[i]));
	if (coef > 1) 
		{
			for (rpl_size_t i = 0; i <= pol1.degree(); i++)
				divide(pol1[i],pol1[i],coef);
			for (rpl_size_t i = 0; i <= pol2.degree(); i++)
				divide(pol2[i],pol2[i],coef); 
		}
}

// Output the gcd of the content of the coordinates of a curve_param
bigint cont(const curve_param <bigint> &par)
{
	// Compare the contents of all polynomial coordinates of par
	bigint c0 = cont(par[0]);
	bigint c1 = cont(par[1]);
	bigint c2 = cont(par[2]);
	bigint c3 = cont(par[3]);
	
	// Now take gcd's
	c0.assign(gcd(c0,c1));
	c1.assign(gcd(c2,c3));
	c0.assign(gcd(c0,c1));
	
	return c0;
}

// Compute the content of a hom_hom_polynomial, ie the gcd of the coefficients in (u,v)
hom_polynomial <bigint> cont(const hom_hom_polynomial <bigint> &a)
{
	hom_polynomial <bigint> tmp; // initialized to 0

	for (rpl_size_t i = 0; i <= a.degree(); i++)
		tmp = gcd(tmp, a[i]);

	return tmp; // if a=0 return 0 
}

// Compute the primitive part	 of a hom_hom_polynomial a, i.e. a divided by its content
hom_hom_polynomial <bigint> pp(const hom_hom_polynomial <bigint> &a)
{
	hom_hom_polynomial <bigint> tmp;
	tmp.set_degree(a.degree());
	hom_polynomial <bigint> content, new_coeff;
	content = cont(a);
	
	for (rpl_size_t i = 0; i <= a.degree(); i++)
		{
			divide(new_coeff, a[i], content);
			tmp[i].assign(new_coeff);
		}

	return tmp; 
}

// Compute the primitive part of a hom_hom_polynomial poly knowing its content
hom_hom_polynomial <bigint> pp(const hom_hom_polynomial <bigint> &a, 
																const hom_polynomial <bigint> &content)
{
	hom_hom_polynomial <bigint> tmp; 
	tmp.set_degree(a.degree());
	hom_polynomial <bigint> new_coeff;

	for (rpl_size_t i = 0; i <= a.degree(); i++)
		{
			divide(new_coeff, a[i], content);
			tmp[i].assign(new_coeff);
		}

	return tmp; 
}

// Optimize a curve_param
void optimize(curve_param <bigint> &par)
{
	bigint c0 = cont(par);

	if (!c0.is_one())
		{
			par[0].assign(par[0]/c0);
			par[1].assign(par[1]/c0);
			par[2].assign(par[2]/c0);
			par[3].assign(par[3]/c0);
		}
}

// Optimize a pair of curve_params by taking the gcd of their content
void optimize(curve_param <bigint> &par1, curve_param <bigint> &par2)
{
	bigint c0 = gcd(cont(par1),cont(par2));

	if (!c0.is_one())
		{
			par1[0].assign(par1[0]/c0);
			par1[1].assign(par1[1]/c0);
			par1[2].assign(par1[2]/c0);
			par1[3].assign(par1[3]/c0);
			par2[0].assign(par2[0]/c0);
			par2[1].assign(par2[1]/c0);
			par2[2].assign(par2[2]/c0);
			par2[3].assign(par2[3]/c0);
		}
}

// Message for factor extraction
void extract_message(const int opt_level, ostream &s, const string &f)
{
	s << ">> " << f << " using ";

	if (opt_level == 0)
		s << "rpl's trialdiv up to MAXFACTOR" << endl;
	else
		s << "rpl's factor" << endl;
}

// Take the integer nearest to the cubic root of an integer
bigint cubic_root(const bigint &b)
{
	assert(b >= 0);
	bigint tmp,tmp2 = b;
	mpz_root(tmp.val.get_mpz_t(), tmp2.val.get_mpz_t(), 3);

	/* to get the nearest integer, 
	 * let x = floor(b^(1/3)), then the nearest
	 integer is x when b <= (x+1/2)^3 = x^3 + 3/2*x^2 + 3/4*x + 1/8,
	 and x+1 otherwise (assuming b > 0). */
	bigfloat x = tmp;
	bigfloat bf = b;
	if ( bf <= ((x + 0.5) * (x + 0.5) * (x+ 0.5)))
	    return tmp;
	else   
	    return tmp+1;
}

// Compute the square factors of a bigint - Output is [b c], a = b^2 c
// NB: I discovered after writing this that there is something similar in
// lidia/nmbrthry_functions.h
math_vector <bigint> extract_square_factors(const bigint &a, const int opt_level, 
							 															ostream &s)
{
	math_vector <bigint> out(2,2);
	out[1].assign(a.sign());
	
	// If a is a perfect square, no need to do the big thing
	if (is_square(out[0],abs(a)))
		return out;
	else // a is not a perfect square
		{
			out[0].assign_one();

			rational_factorization f;
			f.assign(abs(a));

			int maxdiv;

			bigint cub = cubic_root(abs(a));

			if (cub.is_uint()) // cub can be converted to an unsigned int
				{
					if (MAXFACTOR < cub)
						maxdiv = MAXFACTOR;
					else
						maxdiv = static_cast<unsigned int>(mpz_get_si(cub.val.get_mpz_t()));
				}
			else // cub can't be converted: use MAXFACTOR
				maxdiv = MAXFACTOR;

 			if (opt_level == 1)
		 		{
			 		// Find factors with trialdiv up to MAXFACTOR and then refine if needed
			 		if (MAXFACTOR >= cub)
				 		f.trialdiv(maxdiv);
			 		else
				 		f.factor();
		 		}
 			else
		 		f.trialdiv(maxdiv);

			/*catch(const basic_error &ex) {
			s << ">>>>>> " << "Problem factoring " << a << endl;
			s << "	 ---> " << ex.what() << endl;
			s << ">>>>>> " << "We continue without factoring" << endl;
			out[1].assign(a);
			return out;
			}*/

			int tmp_e;
			bigint h,tmp_b;
		 
			for (size_t i = 0; i < f.no_of_comp(); i++)
				{
		 			tmp_e = f.exponent(i);
		 			tmp_b = f.base(i);
		 			
			 		// Odd exponent
			 		if (tmp_e % 2)
				 		{
					 		// Special business if we look at the last base, this base is
					 		// larger than maxdiv and opt_level = 1
					 
					 		if ((i == f.no_of_comp()-1) && (!opt_level) && (tmp_b > maxdiv))
								{
									bigint sq_tmp;
				
									if (is_square(sq_tmp,tmp_b))
										h.assign(sq_tmp);
									else
										{
											multiply(out[1],out[1],tmp_b);
											h.assign_one();
										}
								}
							else
								{
									if (tmp_e == 1)
										h.assign_one();
									else
										power(h,tmp_b,(tmp_e-1)/2);
						
									multiply(out[1],out[1],tmp_b);
								}
						 }
					 else // Even exponent
						 if (tmp_e == 2)
							 h.assign(tmp_b);
						 else
							 power(h,tmp_b,tmp_e/2);
			
					 multiply(out[0],out[0],h);
				}
			
			return out;
		} // a is not a perfect square
}

// Content of the coefficients of a vector
bigint content(const math_vector <bigint> &vec)
{
	if (vec.size() == 1) {
		return vec[0];
	}
	
	bigint gcd_c = gcd(vec[0],vec[1]); 
	
	for (size_t i = 2; i < vec.size(); i++) {
		gcd_c = gcd(gcd_c, vec[i]); 
	}
	
	return gcd_c;
}

// Optimize a vector by dividing by its content
void optimize(math_vector <bigint> &v1)
{
	bigint g = content(v1);

	if (g != 1)
		divide(v1,v1,g);
}

// Optimize a pair of vectors by dividing by the gcd of the contents
void optimize(math_vector <bigint> &v1, math_vector <bigint> &v2)
{
	bigint g = gcd(content(v1),content(v2));

	if (g != 1)
		{
			divide(v1,v1,g);
			divide(v2,v2,g);
		}
}

// Optimize first half, then last half (capacity is even)
void optimize_by_half(math_vector <bigint> &v1)
{
	int l = v1.get_capacity()/2;
	math_vector <bigint> w1(l,l),w2(l,l);

	w1.assign(0, v1, 0, l - 1);
	w2.assign(0, v1, l, 2 * l - 1);
	
	optimize(w1);
	optimize(w2);

	v1.concat(w1,w2);
}

// Optimize a triple of vectors by dividing by the gcd of the contents
void optimize(math_vector <bigint> &v1, math_vector <bigint> &v2, 
							math_vector <bigint> &v3)
{
	bigint g = gcd(gcd(content(v1),content(v2)),content(v3));

	if (g != 1)
		{
			divide(v1,v1,g);
			divide(v2,v2,g);
			divide(v3,v3,g);
		}
}

// Optimize a quadruple of vectors by dividing by the gcd of the contents
void optimize(math_vector <bigint> &v1, math_vector <bigint> &v2, 
							math_vector <bigint> &v3, math_vector <bigint> &v4)
{
	bigint g = gcd(gcd(content(v1),content(v2)),gcd(content(v3),content(v4)));

	if (g != 1)
		{
			divide(v1, v1, g);
			divide(v2, v2, g);
			divide(v3, v3, g);
			divide(v4, v4, g);
		}
}

// Optimize by dividing each column by its content
void optimize_trans1(bigint_matrix &q)
{
	bigint ct;

	for (size_t i = 0; i < q.get_no_of_columns(); i++)
		{
			ct = content(column(q, i));

			if (ct != 1)
				for (size_t j = 0; j < q.get_no_of_rows(); j++)
		 			q.sto(j, i, q(j, i) / ct);
		}
}

// Optimize: each column should be divided by the gcd of the contents of each
// matrix column
void optimize_trans2(bigint_matrix &q1, bigint_matrix &q2)
{
	bigint ct1,ct2,ct;

	for (size_t i = 0; i < q1.get_no_of_columns(); i++)
		{
			ct1 = content(column(q1, i));
			ct2 = content(column(q2, i));
			ct = gcd(ct1,ct2);

			if (ct != 1)
				for (size_t j = 0; j < q1.get_no_of_rows(); j++)
					{
			 			q1.sto(j, i, q1(j,i) / ct);
		 				q2.sto(j, i, q2(j,i) / ct);
					}
		}
}

// Optimize by dividing columns 1, 2, 3 by a, b, c such that a b = c^2
// Divide remaining columns by their content
void optimize_trans3(bigint_matrix &q, const int opt_level, ostream &s)
{
	bigint A = content(column(q, 0));
	bigint B = content(column(q, 1));
	bigint C = content(column(q, 2));

	// A = a0^2*a1*gcdAB and B = b0^2*b1*gcdAB
	bigint gcdAB = gcd(A,B);

	#ifndef NDEBUG
	extract_message(opt_level,s);
	#endif

	bigint a0 = extract_square_factors(A/gcdAB,opt_level,s)[0];
	bigint b0 = extract_square_factors(B/gcdAB,opt_level,s)[0];

	bigint gamma = gcd(gcdAB,C);
	bigint alpha = gcd(a0,C/gamma);
	bigint beta = gcd(b0,C/gamma);
	bigint a = alpha*alpha*gamma;
	bigint b = beta*beta*gamma;
	bigint c = alpha*beta*gamma;

	// There is a simplification only if c != 1
	if (!c.is_one())
		{
			for (size_t i = 0; i < q.get_no_of_rows(); i++)
				{
					// Column 0
					q.sto(i,0,q.member(i,0)/a);
					// Column 1
					q.sto(i,1,q.member(i,1)/b);
					// Column 2
					q.sto(i,2,q.member(i,2)/c);
				}
		}

	// Divide remaining columns by content
	bigint ct;

	for (size_t i = 3; i < q.get_no_of_columns(); i++)
		{
			ct = content(column(q, i));

			if (ct != 1)
				for (size_t j = 0; j < q.get_no_of_rows(); j++)
		 			q.sto(j,i,q.member(j,i)/ct);
		}
}

// Same, but with additional line values to be updated
void optimize_trans3(bigint_matrix &q, math_vector <bigint> &l, const int opt_level, ostream &s)
{
	bigint A = content(column(q, 0));
	bigint B = content(column(q, 1));
	bigint C = content(column(q, 2));

	// A = a0^2*a1*gcdAB and B = b0^2*b1*gcdAB
	bigint gcdAB = gcd(A,B);

	#ifndef NDEBUG
	extract_message(opt_level,s);
	#endif

	bigint a0 = extract_square_factors(A/gcdAB,opt_level,s)[0];
	bigint b0 = extract_square_factors(B/gcdAB,opt_level,s)[0];

	bigint gamma = gcd(gcdAB,C);
	bigint alpha = gcd(a0,C/gamma);
	bigint beta = gcd(b0,C/gamma);
	bigint a = alpha*alpha*gamma;
	bigint b = beta*beta*gamma;
	bigint c = alpha*beta*gamma;

	// There is a simplification only if c != 1
	if (!c.is_one())
		{
			for (size_t i = 0; i < q.get_no_of_rows(); i++)
				{
					// Column 0
					q.sto(i,0,q.member(i,0)/a);
					// Column 1
					q.sto(i,1,q.member(i,1)/b);
					// Column 2
					q.sto(i,2,q.member(i,2)/c);
				}
		}

	// Divide remaining columns by content
	bigint ct;

	for (size_t i = 3; i < q.get_no_of_columns(); i++)
		{
			ct = content(column(q, i));

			if (ct != 1)
				for (size_t j = 0; j < q.get_no_of_rows(); j++)
		 			q.sto(j, i, q.member(j, i) / ct);
		}

	l[0] = l[0]*c;
	l[1] = l[1]*a;
}

// Load balancing of two vectors
bool load_balancing(const math_vector <bigint> &v1, const math_vector <bigint> &v2, 
										bigint &mult, ostream &s)
{
	bigint A = content(v1), B = content(v2);

	bigint maxA = abs(v1[0]), maxB = abs(v2[0]);

	for (int i = 1; i < v1.get_capacity(); i++)
		{
			if (abs(v1[i]) > maxA)
				maxA = abs(v1[i]);
			if (abs(v2[i]) > maxB)
				maxB = abs(v2[i]);
		}

	bigint tmp_e, tmp_b;
	bool finished = 0, maxwhat = 0;
	mult = 1;

	rational_factorization f;

	if (maxB > maxA)
		f = B;
	else
		{
			f = A;
			swap(maxA,maxB);
			maxwhat = 1;
		}

	f.trialdiv(MAXFACTOR);

	for (size_t i = 0; i < f.no_of_comp(); i++)
		{
			tmp_e = f.exponent(i);
			tmp_b = f.base(i);

			while (tmp_e > 0)
				{
		 			if (abs(tmp_b*tmp_b*maxA-maxB) < tmp_b*abs(maxA-maxB))
						{
							mult = mult*tmp_b;
							maxA = maxA*tmp_b;
							maxB = maxB/tmp_b;
			
							tmp_e--;
						}
					else
						{
							finished = 1;
							break;
						}
				}
			
			if (finished)
				break;
		}

	// Returns true if first column has to be divided and second multiplied, false otherwise
	return maxwhat;
}

// Load balancing of columns of a matrix such that a*b = c*d
void load_balancing(bigint_matrix &q, const bool &whatcase, ostream &s)
{
	#ifndef NDEBUG
	extract_message(1, s, "load balancing");
	#endif

	bigint mult;

	if (load_balancing(column(q, 0), column(q, 1), mult, s))
		{
			column(q, 0) = column(q, 0) / mult;
			column(q, 1) = column(q, 1) * mult;
		}
	else
		{
			column(q, 0) = column(q,0) * mult;
			column(q, 1) = column(q,1) / mult;
		}

	if (whatcase)
		{
			if (load_balancing(column(q, 2),column(q, 3), mult, s))
				{
		 			column(q, 2) = column(q, 2) / mult;
		 			column(q, 3) = column(q, 3) * mult;
				}
			else
				{
		 			column(q, 2) = column(q, 2) * mult;
		 			column(q, 3) = column(q, 3) / mult;
				}
		}
}

// Optimize by dividing columns 1, 2, 3 by a, b, c such that a b = c^2
// But here the two matrices should be divided simultaneously
// Dividing remaining columns by gcd of content
void optimize_trans4(bigint_matrix &q1, bigint_matrix &q2, const int opt_level, 
				 							ostream &s)
{
	bigint A = gcd(content(column(q1, 0)), content(column(q2, 0)));
	bigint B = gcd(content(column(q1, 1)), content(column(q2, 1)));
	bigint C = gcd(content(column(q1, 2)), content(column(q2, 2)));

	// A = a0^2*a1*gcdAB and B = b0^2*b1*gcdAB
	bigint gcdAB = gcd(A,B);

	#ifndef NDEBUG
	extract_message(opt_level,s);
	#endif

	bigint a0 = extract_square_factors(A / gcdAB, opt_level, s)[0];
	bigint b0 = extract_square_factors(B / gcdAB, opt_level, s)[0];

	bigint gamma = gcd(gcdAB, C);
	bigint alpha = gcd(a0, C / gamma);
	bigint beta	 = gcd(b0, C / gamma);
	bigint a = alpha * alpha * gamma;
	bigint b = beta	 * beta	 * gamma;
	bigint c = alpha * beta	 * gamma;

	// There is a simplification only if c != 1
	if (!c.is_one())
		{
			for (size_t i = 0; i < q1.get_no_of_rows(); i++)
				{
		 			// Column 0
		 			q1.sto(i, 0, q1(i, 0) / a); q2.sto(i, 0, q2(i, 0) / a);
		 			// Column 1
		 			q1.sto(i, 1, q1(i, 1) / b); q2.sto(i, 1, q2(i, 1) / b);
		 			// Column 2
		 			q1.sto(i, 2, q1(i, 2) / c); q2.sto(i, 2, q2(i, 2) / c);
				}
		}

	// Dividing remaining columns by gcd of content
	bigint ct1,ct2,ct;

	for (size_t i = 3; i < q1.get_no_of_columns(); i++)
		{
			ct1 = content(column(q1, i));
			ct2 = content(column(q2, i));
			ct = gcd(ct1, ct2);

			for (size_t j = 0; j < q1.get_no_of_rows(); j++)
				{
		 			q1.sto(j, i, q1(j, i) / ct);
		 			q2.sto(j, i, q2(j, i) / ct);
				}
		}
}

// Compute the pseudo-discriminant of a quadratic homogeneous polynomial
bigint discriminant2(const hom_polynomial <bigint> &pol)
{
	return pol[1]*pol[1]-4*pol[0]*pol[2];
}

// Optimize the parameterization of a smooth quartic c1 + sqrt(xi). c2 +
// eps. sqrt(Delta). (c3 + sqrt(xi). c4) 
// Delta = Delta1 + sqrt(xi). Delta2
// If xi = 1 then c2, c4	and Delta2 are 0
// Simplification by the common factor of the parameterization
void optimize(curve_param <bigint> &c1, curve_param <bigint> &c2, 
							curve_param <bigint> &c3, curve_param <bigint> &c4, 
							hom_polynomial <bigint> &Delta1, hom_polynomial <bigint> &Delta2, 
							const int opt_level, ostream &s)
{
	bigint factor_Delta = extract_square_factors(gcd(cont(Delta1),
							 					cont(Delta2)),opt_level,s)[0]; 

	bigint fD2 = factor_Delta*factor_Delta;

	divide(Delta1, Delta1, fD2);
	divide(Delta2, Delta2, fD2);

	bigint factor_c12 = gcd(cont(c1), cont(c2));

	if (factor_c12 == 1) // Multiply c3, c4 by factor_Delta and leave
		{
			multiply(c3, c3, factor_Delta);
			multiply(c4, c4, factor_Delta);
		}
	else
		{
			bigint factor_c34 = gcd(cont(c3), cont(c4));

			multiply(c3, c3, factor_Delta);
			multiply(c4, c4, factor_Delta);
			
			bigint factor = gcd(factor_c12, factor_Delta*factor_c34);

			if (factor != 1)
 				{
 					#ifndef NDEBUG
 					s << ">> simplification factor of param: " << factor << endl;
 					#endif

 					divide(c1, c1, factor); 
 					divide(c2, c2, factor); 
 					divide(c3, c3, factor); 
 					divide(c4, c4, factor); 
				}
		}
}

// Optimize the parameterization of a conic c1 + sqrt(xi). c2 + eps. sqrt(D). (c3
// + sqrt(xi). c4) 
// where D and xi are already optimized
// Simplification by the common factor of the parameterization
void optimize(curve_param <bigint> &c1, curve_param <bigint> &c2, 
				 			curve_param <bigint> &c3, curve_param <bigint> &c4, ostream &s)
{
	bigint factor = gcd(cont(c1), cont(c2));
	factor = gcd(factor, cont(c3));
	factor = gcd(factor, cont(c4));

	if (factor != 1)
		{
			#ifndef NDEBUG
			s << ">> simplification factor of param: " << factor << endl;
			#endif

			divide(c1, c1, factor); 
			divide(c2, c2, factor); 
			divide(c3, c3, factor); 
			divide(c4, c4, factor); 
		}
}

} // end of namespace QI

