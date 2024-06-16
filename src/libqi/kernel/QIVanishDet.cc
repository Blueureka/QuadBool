// Intersection when the determinantal equation vanishes identically

#include <libqi/kernel/QIVanishDet.h>
#include <libqi/kernel/QIElem.h>
#include <libqi/kernel/QINumber.h>

using namespace rpl;

// Enter namespace QI
namespace QI {

//////////////////////////////////////////////////////////////////////////////////
// The four concurrent lines case, when we have found one (or two) rational pair
// of planes
quad_inter <bigint> inter_vanish_four_concurrent_lines_special(const bigint_matrix &q1, 
																					  										const bigint_matrix &q2,
																					  	  								const bigint_matrix &q1_r, 
																					  	  								const bigint_matrix &q2_r,
																					  	  								const math_vector <bigint> &sing_p0,
																					  	  								const bigint_matrix &proj_mat,
																					  	  								const int opt_level, 
																					  	  								ostream &s)
{
	// One, two or four components
	quad_inter <bigint> ic;

	// Here q1_r is necessarily a rational pair of lines -- q2_r is a rational pair
	// of lines or a conic with rational coefficients

	// Parameterize the pair of lines q1_r
	bigint D;
	bigint_matrix m1(3,2),m2(3,1);

	pair_of_lines_param(q1_r,column(QI::singular(q1_r), 0),D,m1,m2,opt_level,s);

	// If D is negative, the solution is reduced to one point

	if (D.is_zero())
		{
			// The two lines are rational
			bigint_matrix l1p = m1;
			column(l1p, 0) = column(m1, 0) + column(m2, 0);
			bigint_matrix l1m = m1;
			column(l1m, 0) = column(m1, 0) - column(m2, 0);

			#ifndef NDEBUG
			s << ">> optimization of transformation" << endl;
			#endif

			optimize_trans1(l1p);
			optimize_trans1(l1m);

			// Plug in the other quadric
			bigint_matrix p2p = prod(base_matrix<bigint> (prod(trans(l1p), q2_r)), l1p);
			bigint_matrix p2m = prod(base_matrix<bigint> (prod(trans(l1m), q2_r)), l1m);

			// Parameterize the points
			bigint_matrix m1_1d(2,1),m2_1d(2,1),m1_1d2(2,1),m2_1d2(2,1);

			two_by_two_param(p2p,D,m1_1d,m2_1d,opt_level,s);

			bigint D2;
			two_by_two_param(p2m,D2,m1_1d2,m2_1d2,opt_level,s);

			if ((D < 0) && (D2 < 0))
				{
		 			ic = (quad_inter <bigint>)(1);
		 			ic.set_type(16,1);
				}
		 	else if (((D < 0) && (D2 >= 0)) || ((D >= 0) && (D2 < 0)))
				{
		 			ic = (quad_inter <bigint>)(2);
		 			ic.set_type(16,2);		 
				}
			else
				{
		 			ic = (quad_inter <bigint>)(4);
		 			ic.set_type(16,3);
				}

			#ifndef NDEBUG
			//print_type(ic,s);
			#endif

			// Stack of transformation
			l1p.resize(4,2);
			bigint_matrix tr_stack = prod(proj_mat, l1p);

			// If D is negative and D2 > 0, this means that two of the lines are complex

			if (D.is_zero())
				{
		 			// The two points are rational
		 			bigint_matrix pt1 = prod(tr_stack, m1_1d+m2_1d);
		 			bigint_matrix pt2 = prod(tr_stack, m1_1d-m2_1d);

					#ifndef NDEBUG
		 			s << ">> optimization of parameterization" << endl;
		 			#endif

		 			optimize_trans1(pt1);
		 			optimize_trans1(pt2);

		 			curve_param <bigint> line1(sing_p0,column(pt1, 0)), line2(sing_p0,column(pt2, 0));

		 			// Cut parameter of the lines
		 			math_vector <bigint> cut1(2,2),cut2(2,2);

				 	#ifndef NDEBUG
		 			s << ">> reparameterization of lines" << endl;
					#endif

		 			cut1 = improved_line_param(line1);
		 			cut2 = improved_line_param(line2);

		 			optimize_by_half(cut1);
		 			optimize_by_half(cut2);

		 			ic.cc[0].create_component(INTER_TYPE_LINE,line1);
		 			ic.cc[0].create_cut_parameter(-1,cut1);

					ic.cc[1].create_component(INTER_TYPE_LINE,line2);
					ic.cc[1].create_cut_parameter(-1,cut2);
				}
			else if (D > 0)
				{
		 			// The points are not rational

					// Stack the transformations
		 			bigint_matrix pt1 = prod(tr_stack,m1_1d);
					bigint_matrix pt2 = prod(tr_stack,m2_1d);

					#ifndef NDEBUG
		 			s << ">> optimization of parameterization" << endl;
					#endif

		 			optimize_trans2(pt1,pt2);

		 			math_vector <bigint> zero(4,4);

					curve_param <bigint> line_rat(sing_p0,column(pt1, 0)), line_D(zero,column(pt2, 0));

					math_vector <bigint> cut(6,6);

					// Reparameterize the lines

				 	#ifndef NDEBUG
		 			s << ">> reparameterization of lines" << endl;
					#endif

					cut = improved_line_param(line_rat,line_D,D);
		 
		 			optimize_by_half(cut);

		 			ic.create_two_components(0,INTER_TYPE_LINE,line_rat,line_D,D);
		 			ic.cc[0].create_cut_parameter(-1,cut,D);
					ic.cc[1].create_cut_parameter(-1,cut,D);
				}

	 		// Stack of transformation
		 	l1m.resize(4,2);
			tr_stack = prod(proj_mat, l1m);

			if (D2.is_zero())
				{
		 			// The two points are rational
					bigint_matrix pt1 = prod(tr_stack, m1_1d2+m2_1d2);
					bigint_matrix pt2 = prod(tr_stack, m1_1d2-m2_1d2);

				 	#ifndef NDEBUG
		 			s << ">> optimization of parameterization" << endl;
		 			#endif

					optimize_trans1(pt1);
					optimize_trans1(pt2);

		 			curve_param <bigint> line3(sing_p0,column(pt1, 0)), line4(sing_p0,column(pt2, 0));

		 			// Cut parameter of the lines
		 			math_vector <bigint> cut1(2,2), cut2(2,2);

					#ifndef NDEBUG
		 			s << ">> reparameterization of lines" << endl;
					#endif

		 			cut1 = improved_line_param(line3);
		 			cut2 = improved_line_param(line4);

		 			optimize_by_half(cut1);
					optimize_by_half(cut2);

					unsigned int shift = 0;

		 			if (D >= 0)
		 				shift = 2;

					ic.cc[shift].create_component(INTER_TYPE_LINE,line3);
					ic.cc[shift].create_cut_parameter(-1,cut1);

		 			ic.cc[shift+1].create_component(INTER_TYPE_LINE,line4);
		 			ic.cc[shift+1].create_cut_parameter(-1,cut2);
				}
			else if (D2 > 0)
				{
		 			// The points are not rational

					// Stack the transformations
		 			bigint_matrix pt1 = prod(tr_stack,m1_1d2);
		 			bigint_matrix pt2 = prod(tr_stack,m2_1d2);

				 	#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif

					optimize_trans2(pt1,pt2);

					math_vector <bigint> zero(4,4);

		 			curve_param <bigint> line_rat(sing_p0,column(pt1, 0)),line_D(zero,column(pt2, 0));

		 			math_vector <bigint> cut(6,6);

					// Reparameterize the lines

				 	#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
				 	#endif

					cut = improved_line_param(line_rat,line_D,D2);

		 			optimize_by_half(cut);

					unsigned int shift = 0;

					if (D >= 0)
		 				shift = 2;

					ic.create_two_components(shift,INTER_TYPE_LINE,line_rat,line_D,D2);
					ic.cc[shift].create_cut_parameter(-1,cut,D2);
		 			ic.cc[shift+1].create_cut_parameter(-1,cut,D2);
				}

			if ((D < 0) && (D2 < 0))
				{
		 			// One component

					// Output singular point of all quadrics
					curve_param <bigint> c(sing_p0);

		 			ic.cc[0].create_component(INTER_TYPE_POINT,c);
				}
			
			ic.set_optiflag(true);

			return ic;
		}
	else if (D > 0)
		{
			// The lines are not rational

			// The line transformations are m1 +/- sqrt(D)*m2.

			bigint_matrix m_1 = prod(base_matrix<bigint> (prod(trans(m1), q2_r)), m1); // 2 x 2
			bigint_matrix m_2 = prod(base_matrix<bigint> (prod(D*trans(m2), q2_r)), m2); // 1 x 1
			bigint_matrix m_3 = bigint(2)*prod(base_matrix<bigint> (prod(trans(m1), q2_r)), m2); // 2 x 1

			bigint del_rat = 4*m_1(1,0)*m_1(1,0)+D*m_3(1,0)*m_3(1,0)
							-4*m_1(1,1)*m_1(0,0)-4*m_1(1,1)*m_2(0,0);
			bigint del_D = 4*m_1(1,0)*m_3(1,0)-4*m_1(1,1)*m_3(0,0);

			bigint e = 1;

			#ifndef NDEBUG
			extract_message(opt_level,s,"optimization of square root");
			#endif

			math_vector <bigint> D_fact1 = extract_square_factors(del_rat,opt_level,s);
			math_vector <bigint> D_fact2 = extract_square_factors(del_D,opt_level,s);
		 
			e = gcd(D_fact1[0],D_fact2[0]);

			bigint e2 = e*e;

			del_rat = del_rat/e2;
			del_D = del_D/e2;

			// If delta < 0, only two real lines
			bigint delta = del_rat*del_rat-D*del_D*del_D;

			if (delta < 0)
				ic = (quad_inter <bigint>)(2);
			else
				ic = (quad_inter <bigint>)(4);

			if (delta < 0)
				ic.set_type(16,2);
			else
				ic.set_type(16,3);

			#ifndef NDEBUG
			//print_type(ic,s);
			#endif

			m1.resize(4,2);
			m2.resize(4,1);

			bigint_matrix t1 = prod(proj_mat,m1);
			bigint_matrix t2 = prod(proj_mat,m2);

			math_vector <bigint> k1 = column(t1, 0);
			math_vector <bigint> k2 = column(t1, 1);			
			math_vector <bigint> l1 = column(t2, 0);

			// Solution of quadratic equation is:
			//		(2 m_1(1,1), -2 m_1(1,0) - eps1*m_3(1,0)*sqrt(D) + eps2*e*sqrt(D'))
			// with D' = del_rat + eps1*del_D

			math_vector <bigint> po_1 = 2*m_1(1,1)*k1-2*m_1(1,0)*k2;
			math_vector <bigint> po_2 = 2*m_1(1,1)*l1-m_3(1,0)*k2;
			math_vector <bigint> po_3 = e*k2;

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize(po_1,po_2,po_3);

			math_vector <bigint> zero(4,4);

			// If the Galois group is C4 or V4, the solution can be written in a
			// simpler way 
			bigint sq;

			if (del_D.is_zero()) // V4, first case
				{
					curve_param <bigint> c1(sing_p0,po_1),c2(zero,po_2),c3(zero,po_3),c4(4);

					// For presentation
		 			if (D > del_rat)
			 			{
							swap(D,del_rat);
				 			swap(c2,c3);
						}

					// Reparameterize the lines
					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif

					improved_line_param(c1,c2,c3,c4,D,del_rat,0);

					ic.create_four_components(INTER_TYPE_LINE,c1,c2,c3,c4,D,del_rat);

					// Temp
					ic.cc[0].create_cut_parameter(-1,bigint(0),bigint(0));
					ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));
					ic.cc[2].create_cut_parameter(-1,bigint(0),bigint(0));
					ic.cc[3].create_cut_parameter(-1,bigint(0),bigint(0));
				}
			else if ((delta > 0) && (is_square(sq,delta))) // V4, second case
				{
					// Need two (non-nested) square roots to define the lines
					bigint tmp1 = 2*del_rat;
					bigint tmp2 = 2*sq;

					bigint dp = tmp1+tmp2;
					bigint ep = 1;
					bigint dm = tmp1-tmp2;
					bigint em = 1;

					// The discriminants of the quadratic equation are
					// sqrt(D') = 1/2*(ep*sqrt(dp) +/- em*sqrt(dm))
					bigint Da,Db;
					math_vector <bigint> r1,r2,r3,r4;

					#ifndef NDEBUG		 
					extract_message(opt_level,s,"optimization of second square root");
					#endif

					math_vector <bigint> D_fact = extract_square_factors(dp,opt_level,s);
					ep = D_fact[0];
					dp = D_fact[1];

					D_fact = extract_square_factors(dm,opt_level,s);
					em = D_fact[0];
					dm = D_fact[1];

					if ((dm == 0) || ((dm >= D) && (dm >= dp))) 
						// Smallest square roots are D and dp
						{ 
			 				r1 = ep*ep*dp*po_3;
							r2 = 2*del_D*po_3;
							r3 = 2*ep*po_1;
							r4 = 2*ep*po_2;
				 
							Da = D;
							Db = dp;
						}
					else if ((dp >= D) && (dp >= dm)) // Smallest roots are D and dm
			 			{ 
							r1 = em*em*dm*po_3;
							r2 = 2*del_D*po_3;
							r3 = 2*em*po_1;
							r4 = 2*em*po_2;
				 
							Da = D;
							Db = dm;
						}
					else // D is largest
						{ 
			 				r1 = 2*del_D*po_1;
							r2 = em*del_D*po_3;
							r3 = ep*del_D*po_3;
							r4 = ep*em*po_2;
				 
							Da = dm;
							Db = dp;
						}

					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif

					optimize(r1,r2,r3,r4);

					curve_param <bigint> c1(sing_p0,r1),c2(zero,r2),c3(zero,r3),c4(zero,r4);
					
					// For better presentation		
					if (Da > Db)
			 			{
							swap(Da,Db);
			 				swap(c2,c3);
						}

					// Reparameterize the lines
					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif

					improved_line_param(c1,c2,c3,c4,Da,Db,0);
				 
					ic.create_four_components(INTER_TYPE_LINE,c1,c2,c3,c4,Da,Db);

					// Temp
					ic.cc[0].create_cut_parameter(-1,bigint(0),bigint(0));
					ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));
					ic.cc[2].create_cut_parameter(-1,bigint(0),bigint(0));
					ic.cc[3].create_cut_parameter(-1,bigint(0),bigint(0));
				}
			else if ((delta > 0) && (is_square(sq,D*delta))) // C4
				{
					curve_param <bigint> c1(sing_p0,po_1),c2(zero,po_2),c3(zero,po_3);
					curve_param <bigint> c4(4);

					// Reparameterize the lines

					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif
		
					curve_param <bigint> temp(4);
					improved_line_param(c1,c2,c3,temp,D,del_rat,del_D);

					ic.create_two_components(0,INTER_TYPE_LINE,c1,c2,c3,c4,D,del_rat,del_D);

					math_vector <bigint> po_1_2 = D*(del_D*po_1-del_rat*po_2);
					math_vector <bigint> po_2_2 = del_rat*po_1-D*del_D*po_2;
					math_vector <bigint> po_3_2 = sq*po_3;

					optimize(po_1_2,po_2_2,po_3_2);
	
					curve_param <bigint> c1_2(sing_p0,po_1_2),c2_2(zero,po_2_2),c3_2(zero,po_3_2); 

					// Reparameterize the lines

					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif

					improved_line_param(c1_2,c2_2,c3_2,temp,D,del_rat,del_D);

					ic.create_two_components(2,INTER_TYPE_LINE,c1_2,c2_2,c3_2,c4,D,del_rat,del_D);

					// Temp
					ic.cc[0].create_cut_parameter(-1,bigint(0),bigint(0));
					ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));
					ic.cc[2].create_cut_parameter(-1,bigint(0),bigint(0));
					ic.cc[3].create_cut_parameter(-1,bigint(0),bigint(0));
				}
			else // Full case, i.e. D4
				{
					curve_param <bigint> c1(sing_p0,po_1),c2p(zero,po_2),c2m(zero,-po_2),c3(zero,po_3),c4(4);

					// If delta < 0, only two lines are real

					if (delta < 0)
			 			{
							if (del_D > 0)
								 ic.create_two_components(0,INTER_TYPE_LINE,c1,c2p,c3,c4,D,del_rat,del_D);
			 				else
								 ic.create_two_components(0,INTER_TYPE_LINE,c1,c2m,c3,c4,D,del_rat,-del_D);

							// Temp
			 				ic.cc[0].create_cut_parameter(-1,bigint(0),bigint(0));
							ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));
			 			}
					else
						{
							 ic.create_two_components(0,INTER_TYPE_LINE,c1,c2p,c3,c4,D,del_rat,del_D);
							 ic.create_two_components(2,INTER_TYPE_LINE,c1,c2m,c3,c4,D,del_rat,-del_D);

							// Temp
							ic.cc[0].create_cut_parameter(-1,bigint(0),bigint(0));
							ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));
							ic.cc[2].create_cut_parameter(-1,bigint(0),bigint(0));
							ic.cc[3].create_cut_parameter(-1,bigint(0),bigint(0));
						}
				}
		}
	else // One point
		{
			// One component
			ic = (quad_inter <bigint>)(1);

			ic.set_type(16,1);

			#ifndef NDEBUG
			//print_type(ic,s);
			#endif

			// Output singular point of all quadrics
			curve_param <bigint> c(sing_p0);

			ic.cc[0].create_component(INTER_TYPE_POINT,c);
		}

	ic.set_optiflag(true);

	return ic;
}

// Compute all positive divisors of a given rational factorization (trialdiv)
// maxdiv is to limit the height of the divisors considered
inline void all_positive_divisors1(math_vector <bigint> &div, rational_factorization &f,
					 													const bigint &maxdiv)
{
	rpl_size_t e,i,j,k,l,fno,fexp;
	rpl_size_t cap = 1;
	bigint h;
	bigint base;

	fno = f.no_of_comp();

	for (i = 0; i < fno; i++)
		cap *= (1 + f.exponent(i));

	div.set_capacity(cap);
	div.set_size(cap);

	k = 0;
	l = 1;
	div[0] = bigint (1);

	for (i = 0; i < fno; i++) 
		{
			h.assign_one();
			fexp = f.exponent(i);
			base = f.base(i);

			if (base < maxdiv)
				{
					for (e = 1; e <= fexp; e++) 
						{
			 				multiply(h,h,base);

							for (j = 0; j <= k; j++)
								multiply(div[l++],div[j], h);
						}
			
					k = l - 1;
				}
		}
	
	div.set_size(l);
}

inline math_vector <bigint> all_positive_divisors(rational_factorization &f, 
																									const bigint &maxdiv)
{
	math_vector <bigint> div;
	all_positive_divisors1(div,f,maxdiv);

	return div;
}

// All negative divisors
inline math_vector <bigint> all_negative_divisors(rational_factorization &f, 
																									const bigint &maxdiv)
{
	math_vector <bigint> div;
	all_positive_divisors1(div,f,maxdiv);

	for (size_t i = 0; i < div.size(); i++)
		div[i].negate();

	return div;
}

// Add the negative divisors
inline void all_divisors1(math_vector <bigint> &div, rational_factorization &f,
				 									const bigint &maxdiv)
{
	rpl_size_t i,div_sz,div_cap;

	all_positive_divisors1(div,f,maxdiv);

	div_sz = div.size();
	div_cap = (static_cast<rpl_size_t>(2))*div_sz;
	div.set_capacity(div_cap);
	div.set_size(div_cap);

	div_cap--;
	for (i = 0; i < div_sz; i++) 
		{
			div[div_cap-i] = div[i];
			div[i].negate();
		}
}

inline math_vector <bigint> all_divisors(rational_factorization &f, bigint maxdiv)
{
	math_vector <bigint> div;
	all_divisors1(div,f,maxdiv);

	return div;
}

//////////////////////////////////////////////////////////////////////////////////
// The four concurrent lines case (over the complexes)
quad_inter <bigint> inter_vanish_four_concurrent_lines(const hom_polynomial <bigint> &det_p3, 
																											 const bigint_matrix &q1, 
																											 const bigint_matrix &q2,
																											 const bigint_matrix &q1_r,
																											 const bigint_matrix &q2_r, 
																											 const math_vector <bigint> &sing_p0,
																											 const bigint_matrix &proj_mat, 
																											 const int opt_level, 
																											 ostream &s)
{
	quad_inter <bigint> ic;

	// If one (or both) of the initial quadrics is a pair of planes, special
	// procedure
	if (det_p3[3].is_zero()) // q1 is a pair of planes
		{
			ic = inter_vanish_four_concurrent_lines_special(q1,q2,q1_r,q2_r,sing_p0,proj_mat,opt_level,s);
		}
	else if (det_p3[0].is_zero()) // q2 is a pair of planes
		{
			ic = inter_vanish_four_concurrent_lines_special(q1,q2,q2_r,q1_r,sing_p0,proj_mat,opt_level,s);
		}
	else
		{
			bool rat_found = 0;

			// Try to find rational roots if opt_flag was passed
			if (opt_level)
				{
					rational_factorization f = abs(det_p3[0]), g = abs(det_p3[3]);

					int maxdiv = 100;

					f.trialdiv(maxdiv);
					g.trialdiv(maxdiv);

					rpl_size_t desc = descartes(det_p3);
					math_vector <bigint> div0,div3;

					if (desc == 0) // No positive root
						{
			 				div0 = all_positive_divisors(f,maxdiv);
							div3 = all_negative_divisors(g,maxdiv);
						}
					else if (desc == 3) // No negative root
						{
			 				div0 = all_positive_divisors(f,maxdiv);
							div3 = all_positive_divisors(g,maxdiv);
						}
					else // Negative and positive roots
						{
							div0 = all_positive_divisors(f,maxdiv);
							div3 = all_divisors(g,maxdiv);
						}

					size_t i,j;
					for (i = 0; i < div0.size(); i++)
			 			{
							for (j = 0; j < div3.size(); j++)
								{
									if ((det_p3.eval(div0[i],div3[j])).is_zero())
										rat_found = 1;
									if (rat_found)
										break;
								}
							if (rat_found)
								break;
						}

					if (rat_found)
						{
							#ifndef NDEBUG
			 				s << ">> rational root found" << endl;
							#endif

							bigint_matrix q = div0[i]*q1_r+div3[j]*q2_r;

							ic = inter_vanish_four_concurrent_lines_special(q1,q2,q,q1_r,
																					sing_p0,proj_mat,opt_level,s);
						}
				}

			if (!rat_found)
				{
					#ifndef NDEBUG
		 			if (opt_level)
						s << ">> no rational root found" << endl;
					#endif

					#ifndef NDEBUG
					s << ">> finding separators for the roots of 3 x 3 det equation" << endl;
					#endif

					// Use of inflexion and intersection of inflexion tangent with axis to separate
					// the roots. Inflexion is at -b/(3a). The other point is at (9ad-bc)/2(b^2-3ac)
		
					math_vector <bigint> p1(2,2),p2(2,2);
					p1[0] = -det_p3[2];
					p1[1] = 3*det_p3[3];

					optimize(p1);

					if ((det_p3.eval(p1[0],p1[1])).is_zero())
						{
			 				// We have found a rational root, after all!!!!

							bigint_matrix q = p1[0]*q1_r+p1[1]*q2_r;

							ic = inter_vanish_four_concurrent_lines_special(q1,q2,q,q1_r,
																	sing_p0,proj_mat,opt_level,s);
						}
					else
			 			{
							p2[0] = 9*det_p3[3]*det_p3[0]-det_p3[2]*det_p3[1];
			 				p2[1] = 2*(det_p3[2]*det_p3[2]-3*det_p3[3]*det_p3[1]);

							p1[0] = p1[0]*p1[1].sign();
							p1[1] = abs(p1[1]);
							p2[0] = p2[0]*p2[1].sign();
							p2[1] = abs(p2[1]);	 

			 				if (p2[1] != 0)
								optimize(p2);

							// Re-order the points if needed
			 				if (p1[0]*p2[1]-p2[0]*p1[1] > 0)
								swap(p1,p2);

							// Fake polynomial just to be able to use Descartes' rule of signs
			 				hom_polynomial <bigint> fake;
							fake.set_degree(3);
			 				fake[0] = det_p3.eval(-1,0);
							fake[1] = det_p3.eval(p1[0],p1[1]);
							fake[2] = det_p3.eval(p2[0],p2[1]);
							fake[3] = det_p3.eval(1,0);

			 				// Number of real roots
							unsigned int nbroot = descartes(fake);

							#ifndef NDEBUG
			 				s << ">> number of real roots: " << nbroot << endl;
							#endif

							bool degen_inter = 0;

							// Three real roots
			 				if (nbroot == 3)
								{
									// Decide whether a [3 0] is at p1 and p2

									// The matrix associated to p1
									bigint_matrix q = p1[0]*q1_r+p1[1]*q2_r;

									// Its signed inertia
									math_vector <int> in_q = signed_inertia(q);

									#ifndef NDEBUG
									s << ">> test point 1 at " << p1 << ", sign = " << fake[1].sign() << 
												" -- " << "signed inertia " << in_q << " found" << endl;
									#endif

									// If signed inertia is [1 2] or [2 1]: proceed with next point
									if ((in_q[0] == 1) || (in_q[0] == 2))
										{
											// The matrix associated to p2
											q = p2[0]*q1_r+p2[1]*q2_r;

											// Its signed inertia
											in_q = signed_inertia(q);

											#ifndef NDEBUG
											s << ">> test point 2 at " << p2 << ", sign = " << fake[2].sign() << 
														" -- " << "signed inertia " << in_q << " found" << endl;
											#endif
										}

									// If signed inertia is [0 3] or [3 0]: point
									if ((in_q[0] == 3) || (in_q[0] == 0)) // Point
										{
											degen_inter = 1;

											ic = (quad_inter <bigint>)(1);

											ic.set_type(16,1);

											#ifndef NDEBUG
											//print_type(ic,s);
											#endif

											// Output singular point of all quadrics
											curve_param <bigint> c(sing_p0);

											ic.cc[0].create_component(INTER_TYPE_POINT,c);
										}
								}

				 			if (!degen_inter)
								{
									// Find a conic of the pencil through the point (0,0,1)
									math_vector <bigint> l(2,2), pt(3,3);
									pt[2] = 1;

									bigint_matrix q_rat = pencil_quadric_through_ratpoint(q1_r,q2_r,pt,l);

									if (!(det_p3.eval(l[0],l[1])).is_zero())
										{
											// The point (0,0,1) is on a non-singular conic
											curve_param <bigint> par(3);
											bigint_matrix p_trans(3,3);
											math_vector <bigint> l(2,2);

											conic_param_through_ratpoint(q_rat,pt,p_trans,par,l,s);

											// Now the param is u*sing_p0 + v*proj_mat*p_trans*par, subject to the 
											// constraint quart = 0. We put proj_mat and p_trans in a single
											// transformation matrix and optimize 

											// p_trans is a 3 x 3 matrix so resize first
											p_trans.resize(4,3);

											multiply(p_trans,proj_mat,p_trans);

											#ifndef NDEBUG
											s << ">> optimization of transformation" << endl;
											#endif
						
											optimize_trans3(p_trans,opt_level,s);

											// The param in 3-space
											curve_param <bigint> par_3d(4);
											multiply(par_3d,p_trans,par);

											// Compute the constraint polynomial quart
											hom_polynomial <bigint> quart;

											quart = plug_param_in_quadric(par_3d,q1,par_3d);

											if (quart.is_zero())
												quart = plug_param_in_quadric(par_3d,q2,par_3d);

											#ifndef NDEBUG
											s << ">> optimization of constraint polynomial" << endl;
											#endif

											optimize(quart);

											// The singularity of all cones
											curve_param <bigint> c0(sing_p0);

											if ((quart[0] == 0) || (quart[4] == 0)) // Bit of a simplification here
												{
													ic = (quad_inter <bigint>)(2);

													math_vector <bigint> p_3d;
													if (quart[4] == 0)
														p_3d = par_3d.eval((bigint)1,(bigint)0);
													else
														p_3d = par_3d.eval((bigint)0,(bigint)1);

													optimize(p_3d);

													if (nbroot == 1)
														ic.set_type(16,2);
													else
														ic.set_type(16,3);

													#ifndef NDEBUG
													//print_type(ic,s);
													#endif

													curve_param <bigint> line(sing_p0,p_3d);

													// Cut parameter of the line
													math_vector <bigint> cut(2,2);

													#ifndef NDEBUG
													s << ">> reparameterization of line" << endl;
													#endif
													
													cut = improved_line_param(line);

													optimize_by_half(cut);

													ic.cc[0].create_component(INTER_TYPE_LINE,line);
													ic.cc[0].create_cut_parameter(-1,cut);
									
													hom_polynomial <bigint> div;
													if (quart[4] == 0)
														div.assign_y();
													else
														div.assign_x();
									
													hom_polynomial <bigint> quart2;
													divide(quart2,quart,div);
													
													ic.cc[1].create_component(INTER_TYPE_LINES_WITH_CONSTRAINT,c0,par_3d,quart2);
													ic.cc[1].create_cut_parameter(-1,bigint(1),bigint(0));
												}
											else
												{
										 			ic = (quad_inter <bigint>)(1);
										
										 			if (nbroot == 1)
										 				ic.set_type(16,2);
										 			else
										 				ic.set_type(16,3);
										
													#ifndef NDEBUG
													//print_type(ic,s);
													#endif
										
										 			ic.cc[0].create_component(INTER_TYPE_LINES_WITH_CONSTRAINT,c0,par_3d,quart);
										 			ic.cc[0].create_cut_parameter(-1,bigint(1),bigint(0));
										 		}
										}
									else
										{
											// The point (0,0,1) is on a pair of planes: use these planes to parameterize

											if (l[0].is_zero())
												ic = inter_vanish_four_concurrent_lines_special(q1,q2,q_rat,q1_r,
																		sing_p0,proj_mat,opt_level,s);
											else
												ic = inter_vanish_four_concurrent_lines_special(q1,q2,q_rat,q2_r,
									 									 sing_p0,proj_mat,opt_level,s);
										}
								}
						}
				}
		}

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the two concurrent lines and one double line case
quad_inter <bigint> inter_vanish_two_conc_lines(const bigint_matrix &q1, 
																								const bigint_matrix &q2,
																								const bigint_matrix &q, 
																								const bigint_matrix &q_other,
																								const bigint_matrix &q_sing, 
																								const bigint_matrix &proj_mat, 
																								const math_vector <bigint> &sing_p0,
																								const int in_q, 
																								const int opt_level, ostream &s)
{
	// One or three components
	quad_inter <bigint> ic;

	if (in_q == 1) // Inertia is [1 1]: real pair of planes
		{
			ic = (quad_inter <bigint>)(3);
			ic.set_type(17,2);
		}
	else // Imaginary pair of planes: inertia [2 0]
		{
			ic = (quad_inter <bigint>)(1);
			ic.set_type(17,1);
		}

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// The double line is the singular line of the pair of planes q
	// We already now one point of this line (sing_p0). We compute a second point
	// by taking the singular point of the pair of lines q

	bigint_matrix q_sing3 = q_sing;
	q_sing3.resize(4,1);
	q_sing3 = prod(proj_mat, q_sing3);

	#ifndef NDEBUG
	s << ">> optimization of parameterization" << endl;
	#endif

	optimize_trans1(q_sing3);

	// Let us build a transformation matrix
	bigint_matrix tr(4,2);
	column(tr, 0) = sing_p0;
	column(tr, 1) = column(q_sing3, 0);

	curve_param <bigint> lin(column(tr, 0),column(tr, 1));

	#ifndef NDEBUG
	s << ">> reparameterization of line" << endl;
	#endif

	// Cut parameter of the lines ?????
	math_vector <bigint> cut1(4,4), cut2(4,4);

	cut1 = improved_line_param(lin);			

	ic.cc[0].create_component(INTER_TYPE_LINE,lin);
	ic.cc[0].mult = 2;

	if (in_q == 1)
		{
			// q_other is a pair of lines: parameterize it

			bigint D;
			bigint_matrix m1(3,2),m2(3,1);

			pair_of_lines_param(q_other,column(QI::singular(q_other), 0),D,m1,m2,opt_level,s);

			// The two lines are necessarily rational
			bigint_matrix l1p = m1;
			column(l1p, 0)	= column(m1, 0) + column(m2, 0);
			bigint_matrix l1m = m1;
			column(l1m, 0)	= column(m1, 0) - column(m2, 0);

			// Now, one of trans(l1p)*q_other*l1p and trans(l1m)*q_other*l1m has
			// determinant zero: one line of q_other is in q

			bigint_matrix it_pts = prod(base_matrix<bigint> (prod(trans(l1p), q)), l1p);
			bigint_matrix tr_mat = l1p;

			if (det(it_pts).is_zero())
				{
		 			it_pts = prod(base_matrix<bigint> (prod(trans(l1m), q)), l1m);
					tr_mat = l1m;
				}

			// Parameterize it_pts
			bigint_matrix m1_1d(2,1),m2_1d(2,1);

			two_by_two_param(it_pts,D,m1_1d,m2_1d,opt_level,s);

			if (D.is_zero())
				{
					// The two lines are rational
					bigint_matrix m1_1dp = m1_1d+m2_1d;
					bigint_matrix m1_1dm = m1_1d-m2_1d;
			
					// Return to 2D
					m1_1dp = prod(tr_mat,m1_1dp);
					m1_1dm = prod(tr_mat,m1_1dm);
			
					optimize_trans1(m1_1dp);
					optimize_trans1(m1_1dm);
			
					// And go to 3D: we now have two points on the pair of lines
					m1_1dp.resize(4,1);
					m1_1dm.resize(4,1);
					m1_1dp = prod(proj_mat,m1_1dp);
					m1_1dm = prod(proj_mat,m1_1dm);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
							
					optimize_trans1(m1_1dp);
					optimize_trans1(m1_1dm);
			
					curve_param <bigint> linep(sing_p0,column(m1_1dp, 0)),
															 linem(sing_p0,column(m1_1dm, 0));

					#ifndef NDEBUG
		 			s << ">> reparameterization of lines" << endl;
					#endif

					// Cut parameter of the lines ?????
					math_vector <bigint> cut1(4,4), cut2(4,4);
			
					cut1 = improved_line_param(linep);			
					cut2 = improved_line_param(linem);			
			
					ic.cc[1].create_component(INTER_TYPE_LINE,linep);
					ic.cc[2].create_component(INTER_TYPE_LINE,linem);
				}
			else
				{
					// Lines are not rational
			
					// Stack the transformations
					m1_1d = prod(tr_mat,m1_1d);
					m2_1d = prod(tr_mat,m2_1d);
			
					// And go to 3D: we now have two points on the pair of lines
					m1_1d.resize(4,1);
					m2_1d.resize(4,1);
					m1_1d = prod(proj_mat,m1_1d);
					m2_1d = prod(proj_mat,m2_1d);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
							
					optimize_trans2(m1_1d,m2_1d);
			
					// The components of the lines
					math_vector <bigint> null(4,4);
					curve_param <bigint> line_rat(sing_p0,column(m1_1d, 0)),
															 line_D(null,column(m2_1d, 0));
			
					// Reparameterize the lines
			
					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif
			
					improved_line_param(line_rat,line_D,D);
			
					ic.create_two_components(1,INTER_TYPE_LINE,line_rat,line_D,D);
				}

			#ifndef NDEBUG
			s << ">> the three lines meet at " << sing_p0 << endl;
			#endif
		}

	if (in_q == 1)
		{
			ic.cc[0].create_cut_parameter(-1,bigint(0),bigint(0));
			ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));
			ic.cc[2].create_cut_parameter(-1,bigint(0),bigint(0));
		}

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the two double lines case
quad_inter <bigint> inter_vanish_two_double_lines(const bigint_matrix &q1, 
																									const bigint_matrix &q2,
																									const bigint_matrix &q, 
																									const bigint_matrix &q_other,
																									const bigint_matrix &q_sing, 
																									const bigint_matrix &proj_mat,
																									const math_vector <bigint> &sing_p0, 
																									const int opt_level, ostream &s)
{
	int in_qo = inertia(q_other)[0];

	// One or two components
	quad_inter <bigint> ic;

	if (in_qo == 1) // Inertia of other singular quadric is [1 1]
		{
			ic = (quad_inter <bigint>)(2);
			ic.set_type(19,2);
		}
	else // Inertia of other singular quadric is [2 0]
		{
			ic = (quad_inter <bigint>)(1);
			ic.set_type(19,1);
		}

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// q is a double plane, q_other a pair of planes. When q_other is imaginary, the
	// intersection is reduced to the common singular point sing_p0

	if (in_qo == 2)
		{
			curve_param <bigint> point(sing_p0);

			ic.cc[0].create_component(INTER_TYPE_POINT,point);
		}
	else
		{
			// First parameterize q
			bigint_matrix m1(3,2);
			double_line_param(q,q_sing,m1,s);

			// Now parameterize the pair of lines
			bigint_matrix q_other_1d = prod(base_matrix<bigint> (prod(trans(m1), q_other)), m1);
			bigint D;
			bigint_matrix m1_1d(2,1),m2_1d(2,1);

			two_by_two_param(q_other_1d,D,m1_1d,m2_1d,opt_level,s);

			if (D.is_zero())
				{
					// The two lines are rational
					bigint_matrix m1_1dp = m1_1d+m2_1d;
					bigint_matrix m1_1dm = m1_1d-m2_1d;
			
					// Return to 2D
					m1_1dp = prod(m1,m1_1dp);
					m1_1dm = prod(m1,m1_1dm);
			
					optimize_trans1(m1_1dp);
					optimize_trans1(m1_1dm);
			
					// And go to 3D: we now have two points on the pair of lines
					m1_1dp.resize(4,1);
					m1_1dm.resize(4,1);
					m1_1dp = prod(proj_mat,m1_1dp);
					m1_1dm = prod(proj_mat,m1_1dm);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
							
					optimize_trans1(m1_1dp);
					optimize_trans1(m1_1dm);
			
					curve_param <bigint> linep(sing_p0,column(m1_1dp, 0)),
															 linem(sing_p0,column(m1_1dm, 0));
			
					// Reparameterize the lines
			
					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif
			
					// Cut parameter of the lines ?????
					math_vector <bigint> cut1(4,4), cut2(4,4);
			
					cut1 = improved_line_param(linep);
					cut2 = improved_line_param(linem);
			
					ic.cc[0].create_component(INTER_TYPE_LINE,linep);
					ic.cc[1].create_component(INTER_TYPE_LINE,linem);
				}
			else
				{
					// The lines are not rational
			
					// Stack the transformations
					m1_1d = prod(m1, m1_1d);
					m2_1d = prod(m1, m2_1d);
			
					// And go to 3D: we now have two points on the pair of lines
					m1_1d.resize(4,1);
					m2_1d.resize(4,1);
					m1_1d = prod(proj_mat, m1_1d);
					m2_1d = prod(proj_mat, m2_1d);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
							
					optimize_trans2(m1_1d,m2_1d);
			
					// The components of the lines
					math_vector <bigint> null(4,4);
					curve_param <bigint> line_rat(sing_p0,column(m1_1d, 0)),
															 line_D(null,column(m2_1d, 0));
			
					// Reparameterize the lines
			
					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif
			
					improved_line_param(line_rat,line_D,D);
			
					ic.create_two_components(0,INTER_TYPE_LINE,line_rat,line_D,D);
				}

			ic.cc[0].mult = 2;
			ic.cc[1].mult = 2;
			
			#ifndef NDEBUG
			s << ">> the two lines meet at " << sing_p0 << endl;
			#endif
		}

	// Temp
	if (in_qo != 2)
		{
			ic.cc[0].create_cut_parameter(-1, bigint(0), bigint(0));
			ic.cc[1].create_cut_parameter(-1, bigint(0), bigint(0));
		}

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the quadruple line case
quad_inter <bigint> inter_vanish_quadruple_line(const bigint_matrix &q1, 
																								const bigint_matrix &q2,
																								const bigint_matrix &q, 
																								const bigint_matrix &q_other,
																								const bigint_matrix &q_sing, 
																								const bigint_matrix &proj_mat,
																								const math_vector <bigint> &sing_p0, 
																								const int opt_level, ostream &s)
{
	// One component here
	quad_inter <bigint> ic(1);

	ic.set_type(20,1);

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// We already know one point of the multiple line: sing_p0. The other point if
	// found by parameterizing the double plane q and computing its intersection
	// with q_other

	bigint_matrix m1(3,2);
	double_line_param(q,q_sing,m1,s);			 

	bigint_matrix two_two = prod(base_matrix<bigint> (prod(trans(m1), q_other)), m1);

	// This 2 x 2 matrix is singular and represents a double point
	bigint_matrix point1d(2,1);
	if (!two_two(0,0).is_zero())
		{
			point1d.sto(0,0,-two_two(0,1));
			point1d.sto(1,0,two_two(0,0));
		}
	else
		{
			point1d.sto(0,0,two_two(1,1));
			point1d.sto(1,0,-two_two(0,1));
		}

	// Stack the transformations
	m1.resize(4,2);
	bigint_matrix point3d = prod(base_matrix<bigint> (prod(proj_mat, m1)), point1d);

	#ifndef NDEBUG
	s << ">> optimization of parameterization" << endl;
	#endif
				 
	optimize_trans1(point3d);

	curve_param <bigint> lin(sing_p0, column(point3d, 0));

	// Reparameterize the line

	#ifndef NDEBUG
	s << ">> reparameterization of line" << endl;
	#endif

	// Cut parameter of the lines ?????
	math_vector <bigint> cut1(4,4), cut2(4,4);

	cut1 = improved_line_param(lin);

	ic.cc[0].create_component(INTER_TYPE_LINE,lin);
	ic.cc[0].mult = 4;

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the line and triple line case
quad_inter <bigint> inter_vanish_triple_and_line(const bigint_matrix &q1, 
																								 const bigint_matrix &q2,
																								 const bigint_matrix &q, 
																								 const bigint_matrix &q_other,
																								 const bigint_matrix &q_sing, 
																								 const bigint_matrix &proj_mat,
																								 const math_vector <bigint> &sing_p0, 
																								 const int opt_level, ostream &s)
{
	// Two components
	quad_inter <bigint> ic(2);

	ic.set_type(18,1);

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// The multiple line is the singular line of the 3D pair of planes q
	// We already now one point of this line (sing_p0). We compute a second point
	// by taking the singular point of the pair of lines q.

	bigint_matrix q_sing3 = q_sing;
	q_sing3.resize(4,1);
	q_sing3 = prod(proj_mat, q_sing3);

	#ifndef NDEBUG
	s << ">> optimization of parameterization" << endl;
	#endif

	optimize_trans1(q_sing3);

	// Let us build a transformation matrix
	bigint_matrix tr(4,2);
	column(tr, 0) = sing_p0;
	column(tr, 1) = column(q_sing3, 0);

	curve_param <bigint> lin(column(tr, 0),column(tr, 1));

	#ifndef NDEBUG
	s << ">> reparameterization of line" << endl;
	#endif 

	// Cut parameter of the lines ?????
	math_vector <bigint> cut1(4,4), cut2(4,4);

	cut1 = improved_line_param(lin);			

	// The triple line
	ic.cc[0].create_component(INTER_TYPE_LINE,lin);
	ic.cc[0].mult = 3;

	// Parameterize the pair of lines q
	bigint D;
	bigint_matrix m1(3,2),m2(3,1);

	pair_of_lines_param(q,column(q_sing, 0),D,m1,m2,opt_level,s);

	// The two lines are necessarily rational
	bigint_matrix l1p = m1;
	column(l1p, 0) = column(m1, 0) + column(m2, 0);
	bigint_matrix l1m = m1;
	column(l1m, 0) = column(m1, 0) - column(m2, 0);

	// Now, one of trans(l1p)*q_other*l1p and trans(l1m)*q_other*l1m represents
	// a double point, and the other a pair of points

	bigint_matrix it_pts = prod(base_matrix<bigint> (prod(trans(l1p), q_other)), l1p);
	bigint_matrix tr_mat = l1p;

	if (it_pts(0,1).is_zero())
		{
			// The pair of points is on the other
			it_pts = prod(base_matrix<bigint> (prod(trans(l1m), q_other)), l1m);
			tr_mat = l1m;
		}

	// The two solutions are (0,1) and (2*it_pts(0,1),-it_pts(0,0)): the first
	// point is already known
	bigint_matrix p1d(2,1);
	p1d.sto(0,0,2*it_pts(0,1));
	p1d.sto(1,0,-it_pts(0,0));

	optimize_trans1(p1d);

	// The point in 2D
	p1d = prod(tr_mat, p1d);

	// The point in 3D
	p1d.resize(4,1);
	p1d = prod(proj_mat, p1d);

	#ifndef NDEBUG
	s << ">> optimization of parameterization" << endl;
	#endif

	optimize_trans1(p1d);

	// Let us build a transformation matrix
	bigint_matrix tr2(4,2);
	column(tr2, 0) = sing_p0;
	column(tr2, 1) = column(p1d, 0);

	curve_param <bigint> lin2(column(tr2, 0),column(tr2, 1));

	#ifndef NDEBUG
	s << ">> reparameterization of line" << endl;
	#endif 

	// Cut parameter of the lines ?????
	cut1 = improved_line_param(lin2);			 

	// The simple line
	ic.cc[1].create_component(INTER_TYPE_LINE,lin2);

	#ifndef NDEBUG
	s << ">> the two lines meet at " << sing_p0 << endl;
	#endif	

	// Temp
	ic.cc[0].create_cut_parameter(-1,bigint(0),bigint(0));
	ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the conic and double line case
quad_inter <bigint> inter_vanish_conic_db_line(const bigint_matrix &q1, 
																								const bigint_matrix &q2,
																								const bigint_matrix &sing1, 
																								const bigint_matrix &sing2,
																								const int opt_level, ostream &s)
{
	// Two components
	quad_inter <bigint> ic(2);

	ic.set_type(15,1);

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// There are two cases: 1. when one of q1 or q2 is a rational pair of planes,
	// 2. when the two given quadrics are cones

	if ((sing1.get_no_of_columns() == 2) || (sing2.get_no_of_columns() == 2))
		{
			// q1 or q2 is a pair of planes
			bigint_matrix q,q_other,sing;
			if (sing1.get_no_of_columns() == 2)
				{
					q = q1;
					sing = sing1;
					q_other = q2;
				}
						 else
				{
					q = q2;
					sing = sing2;
					q_other = q1;
				}

			// Parameterize the pair of planes q
			bigint D;
			bigint_matrix m1(4,3),m2(4,1);

			pair_of_planes_param(q,sing,D,m1,m2,opt_level,s);

			// The two planes are rational so we can directly replace the first column
			// of m1 by m1(0) +/- m2(0) 
			bigint_matrix m1p = m1;
			column(m1p, 0) = column(m1, 0) + column(m2, 0);
			bigint_matrix m1m = m1;
			column(m1m, 0) = column(m1, 0) - column(m2, 0);

			#ifndef NDEBUG
			s << ">> optimization of transformation" << endl;
			#endif

			optimize_trans1(m1p);
			optimize_trans1(m1m);

			// Now apply the transformations to q_other
			bigint_matrix qop = prod(base_matrix<bigint> (prod(trans(m1p), q_other)), m1p);
			bigint_matrix qom = prod(base_matrix<bigint> (prod(trans(m1m), q_other)), m1m);

			// One of two matrices represents a double line, the other a conic
			if (!det(qop).is_zero())
				{
					swap(qop,qom);
					swap(m1p,m1m);
				}

			// Parameterize the double line
			bigint_matrix m1db(3,2);
			double_line_param(qop,QI::singular(qop),m1db,s);			

			// Stack the transformations
			bigint_matrix tr = prod(m1p, m1db);

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize_trans1(tr);

			// The double line is this line
			curve_param <bigint> line(column(tr, 0),column(tr, 1));

			math_vector <bigint> cut(4,4);

			#ifndef NDEBUG
			s << ">> reparameterization of line" << endl;
			#endif

			cut = improved_line_param(line);			

			//		 optimize_by_half(cut);

			ic.cc[0].create_component(INTER_TYPE_LINE,line);
			ic.cc[0].mult = 2;
			//			ic.cc[0].create_cut_parameter(cut);

			// Parameterize the conic: it goes through a rational point which is the
			// intersection of the singular line of q with q_other
			curve_param <bigint> tmp(4),sgline4(column(sing, 0),column(sing, 1));

			hom_polynomial <bigint> res;
			multiply(tmp,q_other,sgline4);
			multiply(res,sgline4,tmp);

			math_vector <bigint> dbpoint(4,4);

			if (!res[2].is_zero())
				dbpoint = sgline4.eval(-res[1],2*res[2]);
			else
				dbpoint = sgline4.eval(2*res[0],-res[1]);

			#ifndef NDEBUG
			s << ">> optimization of coordinates of rational point" << endl;
			#endif

			optimize(dbpoint);

			#ifndef NDEBUG
			s << ">> rational point on conic: " << dbpoint << endl;
			#endif

			math_vector <bigint> dbpoint2d(3,3);
			dbpoint2d = prod(base_matrix<bigint>(adj(prod(trans(m1m), m1m))) , math_vector<bigint> (prod(trans(m1m), dbpoint)));

			curve_param <bigint> par(3);
			bigint_matrix m1co(3,3);
			math_vector <bigint> l(2,2);

			conic_param_through_ratpoint(qom,dbpoint2d,m1co,par,l,s);
			
			// Stack the transformations
			tr = prod(m1m, m1co);

			// Rescale so that the double line and the conic meet at the point
			// corresponding to the parameters (1,0)

			math_vector <bigint> m(2,2);
			if (l[1].is_zero())
				{
					m[0] = 0;
					m[1] = 1;
				}
			else
				{
					m[0] = 1;
					m[1] = 0;
				}
			
			// Now rewrite the par accordingly, i.e. replace l[0]*u+l[1]*v by v and 
			// m[0]*u+m[1]*v by u
			
			bigint_matrix par_trans(3,3);
			par_trans.sto(0,0,l[1]*l[1]);
			par_trans.sto(0,1,m[1]*m[1]);
			par_trans.sto(0,2,-2*m[1]*l[1]);
			par_trans.sto(1,0,l[0]*l[0]);
			par_trans.sto(1,1,m[0]*m[0]);
			par_trans.sto(1,2,-2*m[0]*l[0]);
			par_trans.sto(2,0,-l[0]*l[1]);
			par_trans.sto(2,1,-m[0]*m[1]);
			par_trans.sto(2,2,l[0]*m[1]+l[1]*m[0]);

			tr = prod(tr, par_trans);

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize_trans3(tr,l,opt_level,s);

			curve_param <bigint> conic(4);
			multiply(conic,tr,par);

			ic.cc[1].create_component(INTER_TYPE_CONIC,conic);

			#ifndef NDEBUG
			s << ">> the line and the conic meet at " << dbpoint << endl;
			#endif
		}
	else // The two initial quadrics are cones
		{
			// The apex of each cone is a rational point on the other cone - We find a
			// simple rational point on the line linking the two vertices

			// The line param
			curve_param <bigint> line(column(sing1, 0), column(sing2, 0));

			// The rational point
			math_vector <bigint> rat_point = column(sing2, 0);

			// Try to find a point of small height on the line by looking at the
			// parameter values such that one of the components is zero

			#ifndef NDEBUG
			s << ">> reparameterization of line" << endl;
			#endif

			// Cut parameter of the lines ?????
			math_vector <bigint> cut1(4,4), cut2(4,4);
			
			cut1 = improved_line_param(line);
		 
			if (are_equal(column(sing1, 0),line.eval(1,0)))
				rat_point.assign(line.eval(0,1));
			else
				rat_point.assign(line.eval(1,0));

			ic.cc[0].create_component(INTER_TYPE_LINE,line);
			ic.cc[0].mult = 2;

			#ifndef NDEBUG
			s << ">> singular point of cone: " << column(sing1, 0) << endl;
			s << ">> rational point on cone: " << rat_point << endl; 
			#endif

			// Parameterize one of the cones
			surface_param <bigint> par(4);
			bigint_matrix p_trans(4,4);
			math_vector <bigint> l(2,2);

			cone_param_through_ratpoint(q1,column(sing1, 0),rat_point,p_trans,par,l,s);

			#ifndef NDEBUG
			s << ">> optimization of transformation" << endl;
			#endif
						
			optimize_trans3(p_trans,l,opt_level,s);

			// The quadric in the canonical frame of the cone
			bigint_matrix q_tmp = prod(base_matrix<bigint> (prod(trans(p_trans), q2)), p_trans);

			// Plug in the other quadric
			hom_hom_polynomial <bigint> res = plug_param_in_quadric(par,q_tmp,par);

			// The common factor of res (see below)
			hom_polynomial <bigint> tmp;

			// res has the form:
			//			 ---> (l[0]*s+l[1]*t)^2*(u^2*p(s,t) + u*v)

			// The solution (-l[1],l[0]) corresponds to the double line. The
			// solution v = -u*p(s,t) corresponds to the conic. Let us rescale so
			// that (1,0) corresponds to the point of intersection of the line and
			// the conic
			
			math_vector <bigint> m(2,2);
			if (l[1].is_zero())
				{
					m[0] = 0;
					m[1] = 1;
				}
			else
				{
					m[0] = 1;
					m[1] = 0;
				}
			
			// Now rewrite the par accordingly, i.e. replace l[0]*s+l[1]*t by t and 
			// m[0]*s+m[1]*t by s
			bigint_matrix par_trans(4,4);
			par_trans.sto(0,0,l[1]*l[1]);
			par_trans.sto(0,1,m[1]*m[1]);
			par_trans.sto(0,2,-2*m[1]*l[1]);
			par_trans.sto(1,0,l[0]*l[0]);
			par_trans.sto(1,1,m[0]*m[0]);
			par_trans.sto(1,2,-2*m[0]*l[0]);
			par_trans.sto(2,0,-l[0]*l[1]);
			par_trans.sto(2,1,-m[0]*m[1]);
			par_trans.sto(2,2,l[0]*m[1]+l[1]*m[0]);
			par_trans.sto(3,3,1);
			
			// And do the plug again
			res = plug_param_in_quadric(par,prod(base_matrix<bigint> (prod(trans(par_trans), q_tmp)), par_trans),par);
			p_trans = prod(p_trans, par_trans);
			
			tmp.assign_y2();
			
			l[0] = 0;
			l[1] = -1;

			hom_polynomial <bigint> p1,p2;
			divide(p1,res[2],tmp);
			divide(p2,res[1],tmp);
			negate(p2,p2);

			// For the conic, build the param [u^2 v^2 u*v] and the transformation
			// matrix as follows
			curve_param <bigint> conic2d(3);
			conic2d[0].assign_x2();
			conic2d[1].assign_y2();
			conic2d[2].assign_xy();

			bigint_matrix tr(4,3);

			tr.sto(0,0,p2[0]);
			tr.sto(1,1,p2[0]);
			tr.sto(2,2,p2[0]);

			if (!p1.is_zero())
				{
					tr.sto(3,0,p1[2]);
					tr.sto(3,1,p1[0]);
					tr.sto(3,2,p1[1]);
				}

			// Stack the transformations
			tr = prod(p_trans, tr);

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif
						
			optimize_trans3(tr,opt_level,s);

			curve_param <bigint> conic(4);
			multiply(conic,tr,conic2d);

			ic.cc[1].create_component(INTER_TYPE_CONIC,conic);

			#ifndef NDEBUG
			math_vector <bigint> meet_pt = conic.eval(-l[1],l[0]);
			optimize(meet_pt);

			s << ">> the line and the conic meet at " << meet_pt << endl;
			#endif
		}

	// Temp
	ic.cc[0].create_cut_parameter(-1,bigint(0),bigint(0));
	ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The main procedure when the determinantal equation vanishes identically
quad_inter <bigint> inter_vanish_det(const bigint_matrix &q1, const bigint_matrix &q2, 
																		 const hom_polynomial <bigint> &det_p, 
																		 const int opt_level, 
																		 ostream &s)
{
	quad_inter <bigint> ic;

	///////////////////////////////////////////////////////////////////////////////
	// Preliminaries: compute the common singular locus of the two quadrics and its
	// dimension
	///////////////////////////////////////////////////////////////////////////////

	#ifndef NDEBUG
	s << ">> " << "vanishing 4 x 4 determinantal equation" << endl;
	#endif

	// The singular loci of the two quadrics
	bigint_matrix sing1 = QI::singular(q1);
	bigint_matrix sing2 = QI::singular(q2);

	// The common singular locus of the two quadrics
	bigint_matrix sing = linear_intersect(sing1,sing2);

	// The dimension of this singular locus
	rpl_size_t d = sing.get_no_of_columns()-1;

	#ifndef NDEBUG
	s << ">> dimension of common singular locus: " << d << endl;
	#endif

	///////////////////////////////////////////////////////////////////////////////
	// The main switch
	///////////////////////////////////////////////////////////////////////////////
	// d = -1: no singular point in common
	///////////////////////////////////////////////////////////////////////////////
	if (d == -1) 
		ic = inter_vanish_conic_db_line(q1,q2,sing1,sing2,opt_level,s);
	else
		{
			///////////////////////////////////////////////////////////////////////////
			// d >= 0: at least one common singular point

			// The first common singular point
			math_vector <bigint> sing_p0 = column(sing, 0);

			#ifndef NDEBUG
			s << ">> optimization of coordinates of common singular point" << endl;
			#endif
	
			optimize(sing_p0);

			#ifndef NDEBUG
			s << ">> common singular point of quadrics: " << sing_p0 << endl;
			s << ">> computing matrix sending singular point to [ 0 0 0 1 ]" << endl;
			#endif

			// Compute a 4x4 projective transformation sending infinity to the first
			// point of sing
			bigint_matrix proj_mat = send_to_infinity(sing_p0);

			// Its transpose (to avoid recomputing)
			bigint_matrix proj_mat_t = trans(proj_mat);

			// ``Rotate'' q1 and q2, and remove last line and column
			// After that: we enumerate the possible intersections of conics
			bigint_matrix q1_r = prod(base_matrix<bigint> (prod(proj_mat_t, q1)), proj_mat);
			q1_r = mat_minor(q1_r,3,3);
			bigint_matrix q2_r = prod(base_matrix<bigint> (prod(proj_mat_t, q2)), proj_mat);
			q2_r = mat_minor(q2_r,3,3);

			hom_polynomial <bigint> det_p3 = det_pencil3(q1_r,q2_r);

			///////////////////////////////////////////////////////////////////////////
			// New 3 x 3 determinantal equation vanishes
			if (det_p3.is_zero())
				{
					#ifndef NDEBUG
					s << ">> " << "vanishing 3 x 3 determinantal equation" << endl;
					#endif
			
					// No other singular point in common
					if (d == 0)
						{
							#ifndef NDEBUG
							s << ">> no other singular point in common" << endl;
							#endif
			
							// Two components in the intersection
							ic = (quad_inter <bigint>)(2);
			
							ic.set_type(22,1);
			
							#ifndef NDEBUG
							//print_type(ic,s);
							#endif
			
							// All the quadrics of the pencil are pairs of planes
							// Parameterize the pair of planes q1
							bigint D;
							bigint_matrix m1(4,3),m2(4,1);
			
							pair_of_planes_param(q1,sing1,D,m1,m2,opt_level,s);
			
							// The two planes are rational so we can directly replace the first column
							// of m1 by m1(0) +/- m2(0) 
							bigint_matrix m1p = m1;
							column(m1p, 0) = column(m1, 0) + column(m2, 0);
							bigint_matrix m1m = m1;
							column(m1m, 0) = column(m1, 0) - column(m2, 0);
			
							#ifndef NDEBUG
							s << ">> optimization of transformation" << endl;
							#endif
			
							optimize_trans1(m1p);
							optimize_trans1(m1m);
			
							// Now apply the transformations to q2
							bigint_matrix qop = prod(base_matrix<bigint> (prod(trans(m1p), q2)), m1p);
							bigint_matrix qom = prod(base_matrix<bigint> (prod(trans(m1m), q2)), m1m);
			
							// One matrix has rank 0, the other has rank 2
							if (rpl::rank(qop) != 0)
								{
									swap(qop,qom);
									swap(m1p,m1m);
								}

							// Here is the plane
							ic.cc[0].create_component(INTER_TYPE_PLANE,m1p);
							
							bigint_matrix m1d(3,2),m2d(3,1);
			
							pair_of_lines_param(qom,column(QI::singular(qom), 0),D,m1d,m2d,opt_level,s);
			
							// The two lines are necessarily rational
							bigint_matrix l1p = m1d;
							column(l1p, 0) = column(m1d, 0) + column(m2d, 0);
							bigint_matrix l1m = m1d;
							column(l1m, 0) = column(m1d, 0) - column(m2d, 0);
			
							// The lines in 3D
							l1p = prod(m1m, l1p);
							l1m = prod(m1m, l1m);
							
							// The line we want is the one that is not on the plane
							bigint_matrix l1;
			
							if (rpl::rank(prod(trans(kernel(trans(m1p))),l1m)) == 0)
								l1 = l1p;
							else
								l1 = l1m;

							#ifndef NDEBUG
							s << ">> optimization of parameterization" << endl;
							#endif

							optimize_trans1(l1);
			
							curve_param <bigint> lin(column(l1, 0),column(l1, 1));
			
							#ifndef NDEBUG
							s << ">> reparameterization of line" << endl;
							#endif
			
							// Cut parameter of the lines ?????
							math_vector <bigint> cut1(4,4), cut2(4,4);
			
							cut1 = improved_line_param(lin);
			
							ic.cc[1].create_component(INTER_TYPE_LINE,lin);
			
							#ifndef NDEBUG
							s << ">> the line and the plane meet at " << sing_p0 << endl;
							#endif				
			
							ic.set_optiflag(true);
						}
					else
						{
							// More than one singular point in common: send the singular line to
							// z = w = 0 and ``rotate''
			
							// The second common singular point
							math_vector <bigint> sing_p1 = column(sing, 1);
			
							#ifndef NDEBUG
							s << ">> optimization of coordinates of second common singular " 
								<< "point" << endl;
							#endif
			
							optimize(sing_p1);
			
							#ifndef NDEBUG
							s << ">> second common singular point of quadrics: " << sing_p1 
								<< endl;
							s << ">> computing matrix sending singular line to z = w = 0" 
								<< endl;
							#endif

							// Compute a 4x4 projective transformation sending the two common
							// singular points to the line z = w = 0
							bigint_matrix proj_mat = send_to_zw(sing_p0,sing_p1);
			
							// Its transpose (to avoid recomputing)
							bigint_matrix proj_mat_t = trans(proj_mat);
			
							// Note we could have started from q1_r and have rotated them... it
							// would have made aggregation of transformation matrices harder
							// ``Rotate'' q1 and q2, and remove the two last lines and columns
							bigint_matrix q1_r = prod(base_matrix<bigint> (prod(proj_mat_t, q1)), proj_mat);
							q1_r = mat_minor(q1_r,3,3);
							q1_r = mat_minor(q1_r,2,2);
					
							bigint_matrix q2_r = prod(base_matrix<bigint> (prod(proj_mat_t, q2)), proj_mat);
							q2_r = mat_minor(q2_r,3,3);
							q2_r = mat_minor(q2_r,2,2);
			
							hom_polynomial <bigint> det_p2 = det_pencil2(q1_r,q2_r);
							
							// The reduced 2 x 2 determinantal equation vanishes
							if (det_p2.is_zero())
								{
									#ifndef NDEBUG
									s << ">> " << "vanishing 2 x 2 determinantal equation" 
										<< endl;
									#endif
						
									bigint_matrix q,q_r,sing;
						
									// If q1_r is zero, consider q2_r
									if (((q1_r(0,0)).is_zero()) && ((q1_r(1,1).is_zero())))
										{
											q_r = q2_r;
											q = q2;
											sing = sing2;
										}
									else
										{
											q_r = q1_r;
											q = q1;
											sing = sing1;
										}
						
									ic = (quad_inter <bigint>)(1);
						
									// All quadrics of pencil are zero
									if (((q_r(0,0)).is_zero()) && ((q_r(1,1).is_zero())))
										{
											ic.set_type(27,1);					
						
											#ifndef NDEBUG
											//print_type(ic,s);
											#endif
									
											ic.cc[0].create_universe();
										}
									else
										{
											// One component in the intersection
											ic.set_type(26,1);
						
											#ifndef NDEBUG
											//print_type(ic,s);
											#endif

											bigint_matrix m1(4,3);
						
											// Parameterization of double plane
											double_plane_param(q,sing,m1,s);
						
											#ifndef NDEBUG
											s << ">> optimization of parameterization" << endl;
											#endif
						
											optimize_trans1(m1);
						
											ic.cc[0].create_component(INTER_TYPE_PLANE,m1);
										}
								}
							// The reduced 2 x 2 determinantal equation does not vanish
							 else
								{
									#ifndef NDEBUG
									s << ">> optimization of coefficients of 2 x 2 determinantal"
										<< " equation" << endl;
									#endif
										
									optimize(det_p2);
										
									#ifndef NDEBUG
									s << ">> 2 x 2 determinantal equation: ";
									det_p2.print_verbose(s,'l','m');
									s << endl; 
									#endif
						
									bigint delta_p2 = discriminant2(det_p2);
						
									// The discriminant of the 2 x 2 equation is not zero
									if (delta_p2 != 0)
										{
											// One component in the intersection
											ic = (quad_inter <bigint>)(1);

											if (delta_p2 > 0)
												ic.set_type(23,1);
											else
												ic.set_type(23,2);			
						
											#ifndef NDEBUG
											//print_type(ic,s);
											#endif		
						
											curve_param <bigint> lin(sing_p0,sing_p1);
						
											#ifndef NDEBUG
											s << ">> reparameterization of line" << endl;
											#endif
						
											// Cut parameter of the lines ?????
											math_vector <bigint> cut1(4,4), cut2(4,4);
						
											cut1 = improved_line_param(lin);
						
											ic.cc[0].create_component(INTER_TYPE_LINE,lin);
											ic.cc[0].mult = 4;
						
											ic.set_optiflag(true);
										}
									// The discriminant of the 2 x 2 equation vanishes
									else
										{
											// The (lambda,mu) of the multiple root
											math_vector <bigint> root(2,2);					 
						
											if (det_p2[0] == 0)
												{
													root[0].assign(det_p2[1]);
													multiply(root[1],det_p2[2],2);
													root[1].negate();
												}
											else
												{
													multiply(root[0],det_p2[0],2);
													root[0].negate();
													root[1].assign(det_p2[1]);
												}

											#ifndef NDEBUG
											s << ">> optimization of coordinates of root" << endl;
											#endif
						
											optimize(root);
									
											#ifndef NDEBUG
											s << ">> " << "double real root" << ": " << root 
												<< endl;
											#endif
						
											// The associated matrix
											bigint_matrix q = root[0]*q1_r+root[1]*q2_r;
						
											// Its inertia
											math_vector <int> in_q = inertia(q);
						
											#ifndef NDEBUG
											s << ">> inertia: " << in_q << endl;
											#endif
						
											// The singular quadric of the 2 x 2 pencil corresponds to a plane
											if (in_q[0] == 1)
												{
													// One component in the intersection
													ic = (quad_inter <bigint>)(1);
									
													ic.set_type(24,1);
									
													#ifndef NDEBUG
													//print_type(ic,s);
													#endif
									
													// The two quadrics share a common (rational) plane
													// Let us parameterize q
									
													bigint_matrix m1(4,3);
													if (!q(0,0).is_zero())
														{
															m1.sto(0,0,-q(0,1));
															m1.sto(1,0,q(0,0));
														}
													else
														{
															m1.sto(0,0,q(1,1));
															m1.sto(1,0,-q(0,1));
														}
													m1.sto(2,1,1);
													m1.sto(3,2,1);
									
													// Stack the transformations
													m1 = prod(proj_mat, m1);
									
													#ifndef NDEBUG
													s << ">> optimization of parameterization" << endl;
													#endif
									
													optimize_trans1(m1);
									
													ic.cc[0].create_component(INTER_TYPE_PLANE,m1);
												}
											else
												{
													// To avoid picking the zero matrix
													bigint_matrix sing;
													if (root[0].is_zero())
														{
															q = q1;
															sing = sing1;
														}
													else
														{
															q = q2;
															sing = sing2;
														}
									
													// Imaginary pair of planes
													if (inertia_known_rank(q,2)[0] == 2) 
														{
															// One component in the intersection
															ic = (quad_inter <bigint>)(1);
									
															ic.set_type(25,1);
									
															#ifndef NDEBUG
															//print_type(ic,s);
															#endif
									
															curve_param <bigint> lin(sing_p0,sing_p1);
									
															#ifndef NDEBUG
															s << ">> reparameterization of line" << endl;
															#endif
									
															// Cut parameter of the lines ?????
															math_vector <bigint> cut1(4,4), cut2(4,4);
															
															cut1 = improved_line_param(lin);
									
															ic.cc[0].create_component(INTER_TYPE_LINE,lin);
									
															ic.set_optiflag(true);
														}
													else
														{
															// We separate the planes if the discrim is a square
															// Parameterize the pair of planes q
															bigint D;
															bigint_matrix m1(4,3),m2(4,1);
									
															pair_of_planes_param(q,sing,D,m1,m2,opt_level,s);
									
															if (D.is_zero())
																{
																	// The planes are rational
																	bigint_matrix m1p = m1;
																	column(m1p, 0) = column(m1, 0)+ column(m2, 0);
																	bigint_matrix m1m = m1;
																	column(m1m, 0) = column(m1, 0)- column(m2, 0);
												
																	#ifndef NDEBUG
																	s << ">> optimization of transformation" << endl;
																	#endif
												
																	optimize_trans1(m1p);
																	optimize_trans1(m1m);
																	
																	// Two components in the intersection
																	ic = (quad_inter <bigint>)(2);
												
																	ic.set_type(25,2);
																					
																	#ifndef NDEBUG
																	//print_type(ic,s);
																	#endif
												
																	ic.cc[0].create_component(INTER_TYPE_PLANE,m1p);
																	ic.cc[1].create_component(INTER_TYPE_PLANE,m1m);
																}
															else
																{
																	// D is not a square: do not separate
												
																	// One component in the intersection
																	ic = (quad_inter <bigint>)(1);
												
																	ic.set_type(25,2);
												
																	#ifndef NDEBUG
																	//print_type(ic,s);
																	#endif
												
																	ic.cc[0].create_component(INTER_TYPE_PAIR_OF_PLANES,q);
																}
														}
												}
										}
								}
			 			}
				}
			///////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////
			// New 3 x 3 determinantal equation does not identically vanish
			else
				{
					#ifndef NDEBUG
		 			s << ">> optimization of coefficients of 3 x 3 determinantal equation" 
						<< endl; 
					#endif
			 
					optimize(det_p3);

					#ifndef NDEBUG
					s << ">> 3 x 3 determinantal equation: ";
					det_p3.print_verbose(s,'l','m');
					s << endl;
					#endif

					// The gcd of the derivatives of the 3 x 3 det equation
					hom_polynomial <bigint> gcd_p3 = gcd(derivative(det_p3,'x'),
																	 derivative(det_p3,'y'));

					// Its degree
					rpl_size_t d = gcd_p3.degree();

					// Always optimize, otherwise may run into trouble in division
					#ifndef NDEBUG
					s << ">> optimization of coefficients of gcd" << endl;
					#endif
			 
					optimize(gcd_p3);

					#ifndef NDEBUG
					s << ">> gcd of derivatives of 3 x 3 determinantal equation: ";
					gcd_p3.print_verbose(s,'l','m');
					s << endl;
					#endif		
		 
					///////////////////////////////////////////////////////////////////////
					// The reduced 3 x 3 det equation has no multiple root
					if (d == 0) // Only simple roots
						ic = inter_vanish_four_concurrent_lines(det_p3,q1,q2,q1_r,q2_r,sing_p0,
											 proj_mat,opt_level,s); 
					// The reduced 3 x 3 det equation has a multiple root
					else // Double or triple real root
						{
							// The (lambda,mu) of the multiple root
							math_vector <bigint> root(2,2);
			
							if (d == 1) // Double real root
								{
									negate(root[0],gcd_p3[0]);
									root[1].assign(gcd_p3[1]);
								}
							else // Triple (d = 2) root
								if (!gcd_p3[0].is_zero())
									{
										multiply(root[0],gcd_p3[0],2);
										root[0].negate();
										root[1].assign(gcd_p3[1]);
									}
								else
									{
										negate(root[0],gcd_p3[1]);
										multiply(root[1],gcd_p3[2],2);
									}
			
							// Simplify the coordinates of the root if optimization flag was passed
							#ifndef NDEBUG
							s << ">> optimization of coordinates of root" << endl;
							#endif			 
			
							optimize(root);
			
							// The associated matrix
							bigint_matrix q = root[0]*q1_r+root[1]*q2_r;
			
							// Its inertia
							math_vector <int> in_q = inertia(q);
			
							// Its rank
							rpl_size_t rank_q = in_q[0]+in_q[1];
			
							// The singular locus of q (not needed in all cases, but eases
							// readability of code) 
							bigint_matrix q_sing = QI::singular(q);
			
							///////////////////////////////////////////////////////////////////
							// The reduced 3 x 3 det equation has a double real root: the rank
							// of the associated quadric is either 2 or 1
							if (d == 1)
								{
									#ifndef NDEBUG
									s << ">> " << "double real root" << ": " << root << endl;
									#endif
						
									// Divide det_p3 by the gcd squared
									hom_polynomial <bigint> lin3;
									multiply(lin3,gcd_p3,gcd_p3);
									hom_polynomial <bigint> det_e3;
									divide(det_e3,det_p3,lin3);
						
									// Compute the other simple root of the pencil
									math_vector <bigint> root2(2,2);
						
									negate(root2[0],det_e3[0]);
									root2[1].assign(det_e3[1]);
						
									#ifndef NDEBUG
									s << ">> optimization of coordinates of second root" << endl;
									#endif
										 
									optimize(root2);
						
									#ifndef NDEBUG
									s << ">> second root: " << root2 << endl;
									#endif
						
									bigint_matrix q_other = root2[0]*q1_r+root2[1]*q2_r;
									
									// Rank of singular quadric is 2
									if (rank_q == 2)
										ic = inter_vanish_two_conc_lines(q1,q2,q,q_other,q_sing,proj_mat,
														 sing_p0,in_q[0],opt_level,s);
									// Rank of singular quadric is 1
									else
										ic = inter_vanish_two_double_lines(q1,q2,q,q_other,q_sing,proj_mat,
															 sing_p0,opt_level,s);
								}
							///////////////////////////////////////////////////////////////////
							// The reduced 3 x 3 det equation has a triple real root: the rank
							// of the associated quadric is either 2, 1 or 0
							else
								{
									#ifndef NDEBUG
									s << ">> " << "triple real root" << ": " << root << endl;
									#endif
									
									bigint_matrix q_other;
									if (root[0].is_zero())
										q_other = q1_r;
									else
										q_other = q2_r;
						
									if (rank_q == 2)
										ic = inter_vanish_triple_and_line(q1,q2,q,q_other,q_sing,proj_mat,
															sing_p0,opt_level,s);
									else if (rank_q == 1)
										ic = inter_vanish_quadruple_line(q1,q2,q,q_other,q_sing,proj_mat,
														 sing_p0,opt_level,s);
									// Rank of singular quadric is 0
									else
										{
											// To avoid picking the zero matrix
											if (root[0].is_zero()) // wrong : root is for 3x3 det eq,
												// not 4x4
												q = q1;
											else
												q = q2;
						
								 			ic = (quad_inter <bigint>)(1);
			
								 			// The pencil is only made of imaginary cones
								 			if (inertia_known_rank(q,3)[0] == 3)
												{
													// One component in the intersection
													ic.set_type(21,1);
									
													#ifndef NDEBUG
													//print_type(ic,s);
													#endif
									
													// Output singular point of cone
													 curve_param <bigint> cone_p(sing_p0);
									
													ic.cc[0].create_component(INTER_TYPE_POINT,cone_p);
									
													ic.set_optiflag(true);
												}
											else // The pencil is only made of real cones
												{
													// One component in the intersection
													ic.set_type(21,2);
									
													#ifndef NDEBUG
													//print_type(ic,s);
													#endif
									
													ic.cc[0].create_component(INTER_TYPE_CONE,q);
												}
							 			}
								}
						}
				}
		}

	return ic;
}

} // end of namespace QI
