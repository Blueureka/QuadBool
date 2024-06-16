// Special procedures for the ``generic'' algorithm, i.e. going through the param
// of a quadric of inertia (2,2)

#include <libqi/kernel/QIInterWith22.h>

#include <libqi/kernel/QIElem.h>
#include <libqi/kernel/QINumber.h>
#include <libqi/rpl/bigfloat.h>
#include <libqi/rpl/math_matrix.h>

#ifndef NDEBUG
#include <libqi/kernel/QICheck.h>
#endif

using namespace rpl;

// Enter namespace QI
namespace QI {

//// Returns max coeff of discriminant of polynomial p+sqrt(d)*q, after simplification
//bigint h_of_discr4(const hom_polynomial <bigint> &p, const hom_polynomial <bigint> &q, const bigint &d)
//{
//	bigint p0 = p[0], p1 = p[1], p2 = p[2], p3 = p[3], p4 = p[4];
//	
//	bigint q0, q1, q2, q3, q4;
//
//	if (q.is_zero())
//		{
//			q0 = 0;
//			q1 = 0;
//			q2 = 0;
//			q3 = 0;
//			q4 = 0;
//		}
//	else
//		{
//			q0 = q[0];
//			q1 = q[1];
//			q2 = q[2];
//			q3 = q[3];
//			q4 = q[4];
//		}
//
//	bigint ir = 12*p0*p4+12*d*q0*q4-3*p1*p3-3*d*q1*q3+p2*p2+d*q2*q2;
//	bigint id = 12*q0*p4+12*p0*q4-3*p1*q3-3*q1*p3+2*p2*q2;
//
//	bigint jr = 72*p0*p2*p4+9*p1*p2*p3+72*d*p0*q2*q4-54*d*p1*q1*q4
//							 +9*d*p1*q2*q3+72*d*q0*p2*q4-54*d*q0*p3*q3+72*d*q0*q2*p4
//							 +9*d*p2*q1*q3+9*d*q1*p3*q2-2*p2*p2*p2-27*p0*p3*p3-27*p1*p1*p4
//							 -27*d*p0*q3*q3-6*d*p2*q2*q2-27*d*q1*q1*p4;
//	bigint jd = 72*p0*p2*q4-54*p0*p3*q3+72*p0*q2*p4+9*p1*p2*q3
//							 -54*p1*q1*p4+9*p1*p3*q2+72*q0*p2*p4+9*p2*q1*p3+72*d*q0*q2*q4
//							 +9*d*q1*q2*q3-2*d*q2*q2*q2-27*q0*p3*p3-6*p2*p2*q2-27*p1*p1*q4
//							 -27*d*q0*q3*q3-27*d*q1*q1*q4;
//
////////////////////		cout << "I := " << ir << " + " << id << "*sqrt(" << d << ")" << endl;
///////		cout << "J := " << jr << " + " << jd << "*sqrt(" << d << ")" << endl;
//
//	bigint Dr = 4*ir*ir*ir+12*d*ir*id*id-jr*jr-d*jd*jd;
//	bigint Dd = 12*id*ir*ir+4*d*id*id*id-2*jr*jd;
//
////////////////////		cout << "Delta := " << Dr << " + " << Dd << "*sqrt(" << d << ")" << endl;
//
//	bigint pr = 3*p3*p3+3*d*q3*q3-8*p4*p2-8*d*q4*q2;
//	bigint pd = 6*p3*q3-8*p2*q4-8*p4*q2;
//
////////////////////		cout << "p : " << pr << "	 " << pd << endl;
//
//	bigint rr = p3*p3*p3+8*p4*p4*p1-4*p4*p3*p2;
//
////////////////////		cout << "r = " << rr << endl;
//
//	return Dr;
//}
//
//int disc_co(const hom_hom_polynomial <bigint> &c1, const hom_hom_polynomial <bigint> &c2, 
//			 const bigint &xi, const rpl_size_t &i)
//{
//	bigint c10,c11,c12,c20,c21,c22;
//
//	if ((!c1.is_zero()) && (!c1[i].is_zero()))
//		{
//			c10 = c1[i][0];
//			c11 = c1[i][1];
//			c12 = c1[i][2];
//		}
//	if ((!c2.is_zero()) && (!c2[i].is_zero()))
//		{
//			c20 = c2[i][0];
//			c21 = c2[i][1];
//			c22 = c2[i][2];
//		}
//
//	bigint p1 = c11*c11+xi*c21*c21-4*c12*c10-4*xi*c22*c10;
//	bigint p2 = 2*c11*c21-4*c12*c20-4*c22*c10;
//
//	// disc is p1+sqrt(xi)*p2
//
//	rpl_size_t sg = sign(p1,p2,xi);
//
//	if (sg != 0)
//		return(sg);
//	else
//		{
//			//			std::cout << "disc 0 ";
//			
//			// Double root
//
//
//			return(3);
//		}	 
//}

// Compute the parameterization of a smooth quartic
// Input : initial quadrics q1, q2 (only for checking)
//				 Parameterization s1 + sqrt(xi) s2 of a quadric 22 in the pencil
//				 Associated biquadratic equation poly1 + sqrt(xi) poly2 = 0
//					 (if xi = 1 then s2 and poly2 are 0)
// Output : compute of the parameterization of the intersection curve
//		 c1 + sqrt(xi). c2 + eps. sqrt(Delta). (c3 + sqrt(xi). c4)
// where Delta = Delta1 + sqrt(xi). Delta2
// (if xi = 1 then c2, c4, and Delta2 are 0)
quad_inter <bigint> smooth_quartic(const bigint_matrix &q1, const bigint_matrix &q2, 
																	 const surface_param <bigint> &s1, 
																	 const surface_param <bigint> &s2, 
																	 const hom_hom_polynomial <bigint> &poly1,
																	 const hom_hom_polynomial <bigint> &poly2, 
																	 const bigint &xi, const unsigned int case_flag,
																	 const int opt_level, ostream &s)
{
	#ifndef NDEBUG
	s << ">> entering smooth_quartic" << endl;
	#endif
	
	// One component here, but two branches
	quad_inter <bigint> ic(2);

	if (case_flag == 0)
		ic.set_type(1,2);
	else if (case_flag == 1)
		ic.set_type(1,3);
	else
		ic.set_type(1,4);

	#ifndef NDEBUG
	// print_type(ic,s);
	#endif

	// Store the (2,2) quadric surface param
	ic.s1 = s1;
	ic.s2 = s2;

	// rpl_size_t co = 0;

	// s << "A1: " << disc_co(poly1,poly2,xi,2) << endl;
	// if (disc_co(poly1,poly2,xi,2) == 1)
	//	 co++;
	// s << "B1: " << disc_co(poly1,poly2,xi,1) << endl;
	// if (disc_co(poly1,poly2,xi,1) == 1)
	//	 co++;
	// s << "C1: " << disc_co(poly1,poly2,xi,0) << endl;
	// if (disc_co(poly1,poly2,xi,0) == 1)
	//	 co++;

	hom_hom_polynomial <bigint> poly1sw = exchange_uv_st(poly1);
	hom_hom_polynomial <bigint> poly2sw = exchange_uv_st(poly2);

	// s << "A2: " << disc_co(poly1sw,poly2sw,xi,2) << endl;
	// if (disc_co(poly1sw,poly2sw,xi,2) == 1)
	//	 co++;
	// s << "B2: " << disc_co(poly1sw,poly2sw,xi,1) << endl;
	// if (disc_co(poly1sw,poly2sw,xi,1) == 1)
	//	 co++;
	// s << "C2: " << disc_co(poly1sw,poly2sw,xi,0) << endl;
	// if (disc_co(poly1sw,poly2sw,xi,0) == 1)
	//	 co++;

	// s << "-----------> count = " << co << endl;
	//	if (co == 6)
	//		exit(0);

	// >> quadric 1: 4*x^2 - 4*x*y - 3*y^2 + 8*y*z - 2*y*w + z^2 - 3*z*w + w^2
	// >> quadric 2: 9*x^2 - 6*x*y - 4*x*w - y^2 + 16*y*z - 10*y*w - 3*z^2 - z*w + 5*w^2


	
	// Compute the discriminant Delta of the bi-quadratic	 equation poly
	// seen as a degree 2 equation in one of the (projective) variables
	// poly = poly1 + sqrt(xi) poly2
	// =	(poly1[2] + sqrt(xi). poly2[2]). u^2 + (poly1[1] + sqrt(xi). poly2[1]). u.v
	//																		 + (poly1[0] + sqrt(xi). poly2[0]). v^2
	// Delta = Delta1 + sqrt(xi). Delta2 where 
	// Delta1 = poly1[1]^2 + d.poly2[1]^2 - 4.poly1[0].poly1[2] - 4.d.poly2[0].poly2[2]
	// Delta2 = 2.poly1[1].poly2[1] - 4.poly1[2].poly2[0] - 4.poly1[0].poly2[2]
	hom_polynomial <bigint> Delta1, Delta2, htmp; 
	if (!poly1.is_zero())
		{
			multiply(Delta1, poly1[1], poly1[1]);
			multiply(htmp, poly1[0], poly1[2]);
			multiply(htmp, htmp, (bigint) 4);
			htmp.negate(htmp);
			add(Delta1, Delta1, htmp);
		}
	if (!poly2.is_zero())
		{
			multiply(htmp, poly2[1], poly2[1]);
			multiply(htmp, htmp, xi);
			add(Delta1, Delta1, htmp);
			multiply(htmp, poly2[0], poly2[2]);
			multiply(htmp, htmp, xi);
			multiply(htmp, htmp, (bigint) 4);
			htmp.negate(htmp);
			add(Delta1, Delta1, htmp);
		}
	if ((!poly1.is_zero()) && (!poly2.is_zero()))
		{
			// Delta2 = 2.poly1[1].poly2[1] - 4.poly1[2].poly2[0] - 4.poly1[0].poly2[2]
			multiply(Delta2, poly1[0], poly2[2]);
			multiply(htmp, poly1[2], poly2[0]);
			add(Delta2, Delta2, htmp);
			Delta2.negate(Delta2);
			multiply(Delta2, Delta2, (bigint) 2);
			multiply(htmp, poly1[1], poly2[1]);
			add(Delta2, Delta2, htmp);
			multiply(Delta2, Delta2, (bigint) 2);
		}
	// End of the computation of the discriminant Delta of the bi-quadratic equation poly

	// Recall the parameterization of the quadric 22 : s = s1 + sqrt(xi) s2		 
	// and the bi-quadratic homogeneous polynomial : poly = poly1 + sqrt(xi) poly2
	// =	(poly1[2] + sqrt(xi). poly2[2]). u^2 + (poly1[1] + sqrt(xi). poly2[1]). u.v
	//																			 + (poly1[0] + sqrt(xi). poly2[0]). v^2
	// Note that neither the coefficient in u^2 nor the coeff in v^2 is zero because 
	// poly is irreducible since we are in the smooth quartic case. 

	// Compute the output	 parameterized curve : 
	//										 c = c1 + sqrt(xi). c2 + eps. sqrt(Delta). (c3 + sqrt(xi). c4)
	// c1 = s1(u = -poly1[1], v = 2.poly1[2]) + xi. s2(u = -poly2[1], v = 2.poly2[2])
	// c2 = s1(u = -poly2[1], v = 2.poly2[2]) + s2(u = -poly1[1], v = 2.poly1[2])
	// c3 = s1(u = 1, v = 0) 
	// c4 = s2(u = 1, v = 0) 
	curve_param <bigint> c1(4), c2(4), c3(4), c4(4), c_tmp(4); 
	hom_polynomial <bigint> htmp1, htmp2;
	c3 = s1.eval((bigint) 1,(bigint) 0);
	c4 = s2.eval((bigint) 1,(bigint) 0);
	if (!poly1.is_zero())
		{
			QI::negate(htmp1, poly1[1]);	
			multiply(htmp2, poly1[2], (bigint) 2);	
			c1 = s1.eval(htmp1, htmp2);	 
			c2 = s2.eval(htmp1, htmp2);	 
		}
	if (!poly2.is_zero())
		{
			QI::negate(htmp1, poly2[1]);	
			multiply(htmp2, poly2[2], (bigint) 2);	
			c_tmp = s2.eval(htmp1, htmp2);	
			multiply(c_tmp, c_tmp,	xi);	
			add(c1, c1, c_tmp);	 
			c_tmp = s1.eval(htmp1, htmp2);	
			add(c2, c2, c_tmp);
		}

	// The parameterization of the smooth quartic is 
	//			c1 + sqrt(xi). c2 + eps. sqrt(Delta). (c3 + sqrt(xi). c4) 
	// where Delta = Delta1 + sqrt(xi). Delta2 

	#ifndef NDEBUG
	s << ">> optimization of parameterization of quartic" << endl;
	#endif

	// square factors are extracted from Delta
	optimize(c1, c2, c3, c4, Delta1, Delta2, opt_level,s);

	// Output parameterized curve : c = c1 + sqrt(xi). c2 + eps. sqrt(Delta). (c3 +
	// sqrt(xi). c4) 
	// if xi=1 c2=[0,0,0,0] and c4 =[0,0,0,0]

	// cout << height_of_polynomial((hom_polynomial <bigint>)h_of_discr4(Delta1,Delta2,xi),		 size_of_input(q1,q2)) << endl;

	if (xi == 1)
		{
			ic.cc[0].create_component(INTER_TYPE_SMOOTH_QUARTIC_BRANCH_1,c1,c3,Delta1);
			ic.cc[1].create_component(INTER_TYPE_SMOOTH_QUARTIC_BRANCH_2,c1,-c3,Delta1);
			
			ic.set_optiflag(true);
		}
	else
		{
			ic.cc[0].create_component(INTER_TYPE_SMOOTH_QUARTIC_BRANCH_1,c1,c2,c3,c4,xi,Delta1,Delta2);
			ic.cc[1].create_component(INTER_TYPE_SMOOTH_QUARTIC_BRANCH_2,c1,c2,-c3,-c4,xi,Delta1,Delta2);

			ic.set_optiflag(false);
		}
		
	#ifndef NDEBUG
	s << ">> exiting smooth_quartic" << endl;
	#endif

	ic.numberOfSqrt=xi;
	return ic;
}

// Compute the parameterization of two real skew	lines
// Input : initial quadrics q1, q2 (only for checking)
//				 Parameterization s1_cp + sqrt(xi) s2_cp of a quadric 22 in the pencil
//				 Associated biquadratic equation poly1 + sqrt(xi) poly2 = 0
//					(if xi = 1 then s2_cp and poly2 are 0)
// Output : print of the parameterization of the intersection curves 
// (one curve_param for each of the two lines)
// c1 + sqrt(xi). c2 + eps. sqrt(Delta). (c3 + sqrt(xi). c4)
// where the coeff of the ci are linear and Delta = (Delta1 + sqrt(xi). Delta2) in R
// (if xi = 1 then c2, c4, and Delta2 are 0)
quad_inter <bigint> two_real_skew_lines(const bigint_matrix &q1, 
																				const bigint_matrix &q2, 
																				const surface_param <bigint> &s1_cp, 
																				const surface_param <bigint> &s2_cp, 
																				const hom_hom_polynomial <bigint> &poly1, 
																				const hom_hom_polynomial <bigint> &poly2, 
																				const bigint &xi, const int opt_level, 
																				ostream &s) 
{
	// Two components here
	quad_inter <bigint> ic(2);

	ic.set_type(14,4);

	#ifndef NDEBUG
	////print_type(ic,s);
	#endif

	// The bi-variate homogeneous degree (2,2) polynomial poly can be factored into 
	// two univariate polynomials of degree 2. One of these has two real roots and 
	// can be factored into two linear terms (the two real lines). The other has complex 
	// roots and correspond to two complex lines. 

	// The 2 factors are (up to a complex factor) 
	//					coeff(poly, u^2)	 and coeff(poly, s^2).
	// Indeed, poly = F1.F2 = (a1.u^2 + b1.uv + c1.v^2) . (a2.s^2 + b2.st + c2.t^2) 
	// coeff(poly, u^2) = a1.F2 
	// coeff(poly, s^2) = a2.F1 
	// coeff(poly, u^2) * coeff(poly, s^2) = a1.a2.F1.F2 = coeff(poly, u^2s^2). poly
	// More precesely the 2 factors are, up to a complex factor, 
	// (1) any of the coeff of poly of u^2, uv or v^2 that is non-zero and
	// (2) any of the coeff of poly of s^2, st or t^2 that is non-zero

	// Let content = coeff(poly, u^2 or uv or v^2) and primpart = coeff(poly, s^2 or
	// st or t^2) such that they are non-zero and the sum_of_squares of the coeff is
	// minimized 
	// content = content1 + sqrt(xi) content2
	// content_i = s^2.content_i[2] + st.content_i[1] + t^2.content_i[0]
	// where the content_i[*] are of degree 0 in (s,t)
	// Primpart = primpart1 + sqrt(xi) primpart2
	// primpart_i = u^2.primpart_i[2] + uv.primpart_i[1] + v^2.primpart_i[0]

	// s1_cp and s2_cp are const: make copies
	surface_param <bigint> s1 = s1_cp, s2 = s2_cp;

	hom_polynomial <bigint> content1, content2, primpart1, primpart2, htmp1, htmp2;
	bigint norm, norm_htmp;
	// Choose content1 = poly1[i] and content2 = poly2[i] such that there are not both 0 
	// and sum_of_squares of the coeff is minimized
	norm = 0;
	for (unsigned int	 i = 0; i < 3; i++)
		{
			htmp1 = 0;
			htmp2 = 0;
			norm_htmp = 0;
			if (!poly1.is_zero())
				{
		 			htmp1 = poly1[i];
		 			if (!htmp1.is_zero())
			 			norm_htmp = htmp1[0]*htmp1[0]+htmp1[1]*htmp1[1]+htmp1[2]*htmp1[2];
				}
			if (!poly2.is_zero())
				{
					htmp2 = poly2[i];
					if (!htmp2.is_zero())
						norm_htmp = norm_htmp + htmp2[0]*htmp2[0]+htmp2[1]*htmp2[1]
							+htmp2[2]*htmp2[2];
				}
			if (norm == 0) // content = 0 
				{
					content1 = htmp1;
					content2 = htmp2;
					norm = norm_htmp;
				}
			else if (norm_htmp != 0) // content !=0 and htmp !=0
				{
					if (((htmp1.is_zero()) || (htmp2.is_zero())) && 
					 	 ( ((!content1.is_zero()) && (!content2.is_zero())) || (norm_htmp<norm)))
						// if htmpi =0 and 
						// ((content1 !=0 and content2 !=0) or ||htmp||<||content||)
						{
					 	 	content1 = htmp1;
					 	 	content2 = htmp2;
					 	 	norm = norm_htmp;
						}
					else if ((!content1.is_zero()) && (!content2.is_zero()) && (norm_htmp<norm))
						{
					 	 	content1 = htmp1;
					 	 	content2 = htmp2;
					 	 	norm = norm_htmp;
						}
				}
		}

	// poly1 and poly2 are const
	hom_hom_polynomial <bigint> poly1_cp,poly2_cp;

	// We now do exactly the same thing for primpart
	poly1_cp = exchange_uv_st(poly1);
	poly2_cp = exchange_uv_st(poly2);
	norm = 0;
	for (unsigned int	 i = 0; i < 3; i++)
		{
			htmp1 = 0;
			htmp2 = 0;
			norm_htmp = 0;
			if (!poly1_cp.is_zero())
				{
					htmp1 = poly1_cp[i];
					if (!htmp1.is_zero())		 
						norm_htmp = htmp1[0]*htmp1[0]+htmp1[1]*htmp1[1]+htmp1[2]*htmp1[2];
				}
			if (!poly2_cp.is_zero())
				{
					htmp2 = poly2_cp[i];
					if (!htmp2.is_zero())
						norm_htmp = norm_htmp + htmp2[0]*htmp2[0]+htmp2[1]*htmp2[1]
							+htmp2[2]*htmp2[2];
				}
			if (norm == 0) // primpart = 0 
				{
					primpart1 = htmp1;
					primpart2 = htmp2;
					norm = norm_htmp;
				}
			else if (norm_htmp != 0) // primpart !=0 and htmp !=0
				{
					if (((htmp1.is_zero()) || (htmp2.is_zero())) && 
					 	 (((!primpart1.is_zero()) && (!primpart2.is_zero())) || (norm_htmp<norm)))
						// if htmpi =0 and 
						// ((primpart1 !=0 and primpart2 !=0) or ||htmp||<||primpart||)
						{
					 	 	primpart1 = htmp1;
					 	 	primpart2 = htmp2;
					 	 	norm = norm_htmp;
						}
					else if ((!primpart1.is_zero()) && (!primpart2.is_zero()) && (norm_htmp<norm))
						{
					 	 	primpart1 = htmp1;
					 	 	primpart2 = htmp2;
					 		norm = norm_htmp;
						}
				}
		}

	// Some optimization
	optimize(content1, content2);
	optimize(primpart1, primpart2);
	
	// Delta = b^2 -4ac, discriminant	 of the polynomial content
	//	 = (b1 + sqrt(xi) b2)^2 - 4 (a1 + sqrt(xi) a2) (c1 + sqrt(xi) c2) 
	//	 = b1^2 + xi b2^2 - 4 a1 c1 - 4 xi a2 c2 + sqrt(xi) [2 b1 b2 - 4 a1 c2 -4 a2 c1]
	// Delta = Delta1 + sqrt(xi) Delta2
	bigint	Delta1, Delta2;		
	if (!content1.is_zero())
		Delta1 = content1[1]*content1[1] - 4*content1[2]*content1[0] ; 
	if (!content2.is_zero())
		Delta1 = Delta1 + xi*content2[1]*content2[1]- 4*xi*content2[2]*content2[0];
	if ((!content1.is_zero()) && (!content2.is_zero()))
		Delta2 = 2*content1[1]*content2[1] - 4*content1[2]*content2[0] -
							4*content2[2]*content1[0]; 

	if (sign(Delta1, Delta2, xi) < 0) // if	 Delta < 0 
		{ // The content correspond to the two complex lines
			// we consider instead the primitive part of poly (ie poly =
			// coeff*prim_part*content) which we recall content for simplicity
			content1 = primpart1;
			content2 = primpart2;
			// We redefine as before	Delta
			Delta1 =0; Delta2 =0; 
			if (!content1.is_zero())
				Delta1 = content1[1]*content1[1] - 4*content1[2]*content1[0] ; 
			if (!content2.is_zero())
				Delta1 = Delta1 + xi*content2[1]*content2[1]- 4*xi*content2[2]*content2[0];
			if ((!content1.is_zero()) && (!content2.is_zero()))
				Delta2 = 2*content1[1]*content2[1] - 4*content1[2]*content2[0] -
		 							4*content2[2]*content1[0]; 
		}
	else
		{
			// If (u,v) and (s,t) haven't been exchanged in poly (ie. content and
			// primpart haven't been exchanged) then the solutions of content are in (s,t)
			// and we have to exchange (u,v) and (s,t) in the parameterization of the
			// quadric (2,2) 
			s1 = exchange_uv_st(s1);
			s2 = exchange_uv_st(s2);
		}
 
	// Take out the square factor from Delta1+sqrt(xi). Delta2 
	bigint delta = extract_square_factors(gcd(Delta1, Delta2), opt_level, s)[0];

	if (delta > 1) 
		{
			Delta1 = Delta1/(delta*delta);
			Delta2 = Delta2/(delta*delta);
		}

	// Delta1, Delta2 are constants but we consider them as constant hom_polynomial
	// for conveniency 
	hom_polynomial <bigint>	Delta1h, Delta2h;	 
	Delta1h = Delta1;
	Delta2h = Delta2;

	// Recall the parameterization of the quadric 22 : s = s1 + sqrt(xi) s2		 
	// and the quadratic homogeneous polynomial : content	 = content1 + sqrt(xi) content2
	// =	(content1[2] + sqrt(xi). content2[2]). u^2 + (content1[1] +
	// sqrt(xi). content2[1]). u.v + (content1[0] + sqrt(xi). content2[0]). v^2

	// Compute the output	 parameterized curve : 
	//										 c = c1 + sqrt(xi). c2 + eps. sqrt(Delta). (c3 + sqrt(xi). c4)
	// If the coeff of u^2 is not zero then : 
	// c1 = s1(u = -content1[1], v = 2.content1[2]) + xi. s2(u = -content2[1], v =
	//		 2.content2[2]) 
	// c2 = s1(u = -content2[1], v = 2.content2[2]) + s2(u = -content1[1], v =
	//		 2.content1[2]) 
	// c3 = s1(u = 1, v = 0) 
	// c4 = s2(u = 1, v = 0) 
	// If the coeff of v^2 is not zero then : 
	// c1 = s1(u = 2.content1[0], v = -content1[1]) + xi. s2(u = 2.content2[0], v =
	//		 -content2[1]) 
	// c2 = s1(u = 2.content2[0], v = -content2[1]) +			 s2(u = 2.content1[0], v =
	//		 -content1[1]) 
	// c3 = s1(u = 0, v = 1) 
	// c4 = s2(u = 0, v = 1) 
	// else is the coeffs of u^2 and v^2 are 0 then the two lines are 
	//			 c= s1(u = 0, v = 1) + sqrt(xi) s2(u = 0, v = 1)
	// and c= s1(u = 1, v = 0) + sqrt(xi) s2(u = 1, v = 0)
	curve_param <bigint> c1(4), c2(4), c3(4), c4(4), c_tmp(4);

	if (((!content1.is_zero()) && (content1[2] != 0)) 
			|| ((!content2.is_zero()) && (content2[2] != 0))) // if coeff of u^2 is not 0
		{ // solutions of content are (u,v) = (-b+-sqrt(Delta) , 2a)
			if (!content1.is_zero()) // otherwise crash when is zero and call content1[i]
				{
					c1 = s1.eval(-content1[1], 2*content1[2]);
					c2 = s2.eval(-content1[1], 2*content1[2]);
				}
			if (!content2.is_zero())
				{
					c_tmp = s2.eval(-content2[1], 2*content2[2]);	
			
					multiply(c_tmp, c_tmp,	 xi);	 
			
					add(c1, c1, c_tmp);	
					c_tmp = s1.eval(-content2[1], 2*content2[2]);	
			
					add(c2, c2, c_tmp);
				}

			c3 = s1.eval((bigint) delta, (bigint) 0);
			c4 = s2.eval((bigint) delta, (bigint) 0);
		}
	else if (((!content1.is_zero()) && (content1[0] != 0)) 
			|| ((!content2.is_zero()) && (content2[0] != 0))) // if coeff of v^2 is not 0
		{ // solutions of content are (u,v) = (2c , -b+-sqrt(Delta))
			if (!content1.is_zero()) 
				{
					c1 = s1.eval(2*content1[0], -content1[1]);
					c2 = s2.eval(2*content1[0], -content1[1]);
				}
			if (!content2.is_zero())
				{
					c_tmp = s2.eval(2*content2[0], -content2[1]);	
					multiply(c_tmp, c_tmp,	 xi);	 
					add(c1, c1, c_tmp);	
					c_tmp = s1.eval(2*content2[0], -content2[1]);
					add(c2, c2, c_tmp);
				}
			c3 = s1.eval((bigint) 0, (bigint) delta);
			c4 = s2.eval((bigint) 0, (bigint) delta);
		}
	else // if coeff of u^2 and v^2 are	 0
		{ // solutions are (u,v) = (0, 1) and (1, 0)
			// The two lines are c= s1(u = 0, v = 1) + sqrt(xi) s2(u = 0, v = 1)
			// and c= s1(u = 1, v = 0) + sqrt(xi) s2(u = 1, v = 0)
			// That is	c1 + sqrt(xi) c2 
			// and			 c1b + sqrt(xi) c2b 
			curve_param <bigint> c1b(4), c2b(4);
			c1 = s1.eval((bigint) 0, (bigint) 1);
			c2 = s2.eval((bigint) 0, (bigint) 1);
			c1b = s1.eval((bigint) 1, (bigint) 0);
			c2b = s2.eval((bigint) 1, (bigint) 0);
			
			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize(c1,c2);
			optimize(c1b,c2b);

			ic.cc[0].create_component(INTER_TYPE_LINE,c1,c2,xi);
			ic.cc[1].create_component(INTER_TYPE_LINE,c1b,c2b,xi);

			ic.set_optiflag(true);

			return ic;
		}

	if ((Delta1h.is_one()) && (Delta2h.is_zero())) // Delta = 1 and xi = 1, c2 and c4 are 0
		{
			curve_param <bigint> ca(4), cb(4);
			add(ca,c1,c3);	// ca = c1 + c3
			QI::negate(cb,c3);
			add(cb,c1,cb); // cb = c1 - c3

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize(ca);
			optimize(cb);

			#ifndef NDEBUG
			s << ">> reparameterization of lines" << endl;
			#endif

			// Cut parameter of the lines ?????
			math_vector <bigint> cut1(4,4), cut2(4,4);

			cut1 = improved_line_param(ca);
			cut2 = improved_line_param(cb);

			ic.cc[0].create_component(INTER_TYPE_LINE,ca);
			ic.cc[1].create_component(INTER_TYPE_LINE,cb);
		}
	else if (xi == 1)
		{ 
			// the 2 lines are c1 + sqrt(Delta1h[0]) c3 and c1 - sqrt(Delta1h[0]) c3

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize(c1,c3);

			#ifndef NDEBUG
			s << ">> reparameterization of lines" << endl;
			#endif

			improved_line_param(c1,c3,Delta1h[0]);

			ic.create_two_components(0,INTER_TYPE_LINE,c1,c3,Delta1h[0]);
		}
	else
		{
			// Output parameterized curve : c = c1 + sqrt(xi). c2 +
			// eps. sqrt(Delta). (c3 + sqrt(xi). c4) 
			// Delta = Delta1 + sqrt(xi). Delta2

			#ifndef NDEBUG
			s << ">> reparameterization of lines" << endl;
			#endif

			improved_line_param(c1,c2,c3,c4,xi,Delta1h[0],Delta2h[0]);

			ic.create_two_components(0,INTER_TYPE_LINE,c1,c2,c3,c4,xi,Delta1h[0],Delta2h[0]);
		}

	ic.set_optiflag(true);

	return ic;
} // void two_real_skew_lines


// Compute the parameterization of a cubic and a line
// Input : initial quadrics q1, q2 (only for checking)
//				 Parameterization par of a quadric 22 in the pencil
//				 Associated biquadratic equation poly = 0
// Output : print of the parameterization of the intersection curves 
// (one curve_param for the cubic and one for the line)
quad_inter <bigint> cubic_line(const bigint_matrix &q1, const bigint_matrix &q2, 
																const surface_param <bigint> &par_cp,
																const hom_hom_polynomial <bigint> &poly_cp, 
																const unsigned int case_flag, ostream &s) 
{
	// Two components here
	quad_inter <bigint> ic(2);

	if (case_flag == 3)
		ic.set_type(12,1);
	else
		ic.set_type(12,2);

	#ifndef NDEBUG
	// print_type(ic,s);
	#endif

	// poly_cp is const: make copy
	hom_hom_polynomial <bigint> poly = poly_cp;
	// par_cp is const: make copy
	surface_param <bigint> par = par_cp;

	// A linear term `content' can be factored in poly : 
	hom_polynomial <bigint> content;
	
	content = cont(poly);
	
	bool exchange_flag = 0;

	if (content.degree() < 1)
		{
			exchange_flag = 1;
			poly = exchange_uv_st(poly);
			content = cont(poly);
		}

	// The primitive part of poly is (ie poly = prim_part*content) :
	hom_hom_polynomial <bigint> prim_part;
	prim_part = pp(poly, content); 

	// prim_part = u^2.prim_part[2](s,t) + uv.prim_part[1](s,t) + v^2.prim_part[0](s,t) 
	// where the prim_part[i](s,t) are of degree 1
	// We exchange (u,v) and (s,t) in	 prim_part :
	prim_part = exchange_uv_st(prim_part);
	// prim_part = u.prim_part[1](s,t) + v.prim_part[0](s,t) 
	// where the prim_part[i](s,t) are of degree 2

	// We also	exchange (u,v) and (s,t) in	 the surface param par if 
	// the exhange of (u,v) and (s,t) has been done only once in poly/prim_part
	if (!exchange_flag)
		par = exchange_uv_st(par);

	// We now replace u=prim_part[0] and v=-prim_part[1] in the surface parameterization
	hom_polynomial <bigint> uu, vv;
	uu = prim_part[0];
	QI::negate(vv, prim_part[1]);
	curve_param <bigint> cubic = par.eval(uu, vv);

	// The content part of poly is content = s.content[1] + t.content[0]
	// We replace s=content[0] and t=-content[1] in the surface parameterization
	// where (u,v) and (s,t) have to be exchanged if content is the content of poly 
	// (vs. of the content of exchange_uv_st(poly)); this has already been done at 
	// if (!exchange_flag)	par = exchange_uv_st(par);
	curve_param <bigint> line = par.eval(content[0], -content[1]);

	#ifndef NDEBUG
	s << ">> optimization of parameterization" << endl;
	#endif

	optimize(cubic);
	optimize(line);

	#ifndef NDEBUG
	s << ">> reparameterization of lines" << endl;
	#endif

	// Cut parameter of the lines ?????
	math_vector <bigint> cut1(4,4), cut2(4,4);

	cut1 = improved_line_param(line);

	ic.cc[0].create_component(INTER_TYPE_CUBIC,cubic);
	ic.cc[1].create_component(INTER_TYPE_LINE,line);

	ic.set_optiflag(true);

	if (case_flag == 3)
		{
			ic.cc[0].create_cut_parameters(-1,0,0,-1,0,0);
			ic.cc[1].create_cut_parameters(-1,0,0,-1,0,0);
		}

	return ic;
}

// Input : initial quadrics q1, q2
//				 determinental equation of the pencil : det_p
//				 a quadric 22 : q
//				 case_flag	= 0 if the intersection is a smooth quartic, 2 real affinely finite
//										= 1						"															 1 real affinely finite
//										= 2						"															 2 real affinely infinite
//										= 3 if intersection is a cubic and a line secant
//										= 4						"														non-secant
//										= 5 if intersection is two real skew lines
// Output : Compute the Parameterization s1 + sqrt(xi) s2 of a quadric 22 in the pencil 
//					and the associated biquadratic equation poly1 + sqrt(xi) poly2 = 0
// At the end, call smooth_quartic(...) or cubic_line(...) for computing the
// intersection curve 
quad_inter <bigint> intersection_with_quadric22(const bigint_matrix &q1, 
																								const bigint_matrix &q2, 
																								const hom_polynomial <bigint> &det_p, 
																								const bigint_matrix &q_cp, 
																								const unsigned int case_flag, 
																								const int opt_level, 
																								ostream &s)
{
	#ifndef NDEBUG
	s << ">> entering intersection_with_quadric22" << endl;
	#endif
	quad_inter <bigint> ic;

	// Compute the new quadric q	of inertia 2,2; 
	// point_on_q is a point on it; det_q its determinant in the form [a,b] with det
	// = a^2. b	 
	math_vector <bigint> point_on_q(4,4); 
	math_vector <bigint> det_q(2,2);

	bigint xi = 1;
	hom_polynomial <bigint> tmp_p;
	//////////////////	cout << "det_p" << endl;
	//////////////////	cout << height_of_polynomial((hom_polynomial <bigint>)h_of_discr4(det_p,tmp_p,xi),	 size_of_input(q1,q2)) << endl;

	// q_cp is const: create new variable
	bigint_matrix q;
	bigint det_R;
	find_quadric_22_and_point(q1,q2,det_p,q_cp,q,point_on_q,det_R,det_q,opt_level,s);

	#ifndef NDEBUG
	s << ">> quadric (2,2) found: ";
	print_quadric(q,s);
	s << endl;
	s << ">> decomposition of its determinant [a,b] (det = a^2*b): " << det_q << endl;
	s << ">> a point on the quadric: " << point_on_q << endl;
	#endif
	
	surface_param <bigint> s1(4), s2(4); // s1 and s2 are initialized to [0,0,0,0]

	find_parameterization_quadric_22(q,point_on_q,det_q,det_R,s1,s2,opt_level,s);
	// Parameterization of the quadric is	 s1 + sqrt(det_q[1]) s2
	// If det_q[1]=1 then	 s2=[0,0,0,0]

	// Compute the bi-quadratic homogeneous polynomial poly = poly1 + sqrt(xi) poly2
	// poly = (s1 + sqrt(det_q[1]) s2)^t. q1. (s1 + sqrt(det_q[1]) s2)
	// (or with q2 if poly = 0)
	// poly = s1.q1.s1 + det_q[1]. s2.q1.s2 + 2. sqrt(det_q[1]). s1.q1.s2
	// We set poly1 = s1.q1.s1 + det_q[1]. s2.q1.s2 
	// and poly2 = 2. s1.q1.s2 

	hom_hom_polynomial <bigint> poly1, poly2, hhtmp; 
	poly1 = plug_param_in_quadric(s1, q1, s1);
	hhtmp = plug_param_in_quadric(s2, q1, s2);
	multiply(hhtmp, hhtmp, det_q[1]);
	add(poly1, poly1, hhtmp);
	poly2 = plug_param_in_quadric(s1, q1, s2);
	multiply(poly2, poly2, (bigint) 2);

	if ((poly1.is_zero()) && (poly2.is_zero()))
		{
			// We recompute poly1 and poly2 after replacing q1 by q2
			poly1 = plug_param_in_quadric(s1, q2, s1);
			hhtmp = plug_param_in_quadric(s2, q2, s2);
			multiply(hhtmp, hhtmp, det_q[1]);
			add(poly1, poly1, hhtmp);
			poly2 = plug_param_in_quadric(s1, q2, s2);
			multiply(poly2, poly2, (bigint) 2);
		}

	optimize(poly1, poly2); // remove the gcd of all (integer) coefficients 
	// End of the computation of the bi-quadratic homogeneous polynomial poly

	// If xi = 1 then s2 and poly2 are 0

	if ((case_flag >= 0) && (case_flag <= 2))// the intersection is a smooth quartic
		ic = smooth_quartic(q1,q2,s1,s2,poly1,poly2,det_q[1],case_flag,opt_level,s);
	else if ((case_flag == 3) || (case_flag == 4)) // intersection is a cubic and a line 
		{
			#ifndef NDEBUG
			if ((!poly2.is_zero()) || (!s2[0].is_zero()) || (!s2[1].is_zero()) ||
		 			(!s2[2].is_zero()) || (!s2[3].is_zero())) 
				s << ">>	 Error :	poly2 = " << poly2 << ", s2 = " << s2 
		 			<< " should both be zero" << endl; 
			#endif

			ic = cubic_line(q1, q2, s1, poly1, case_flag, s);
		}
	else // (case_flag == 5)	the intersection is two real skew	 lines	
		ic = two_real_skew_lines(q1, q2, s1, s2, poly1,	 poly2, det_q[1], opt_level, s);

	#ifndef NDEBUG
	s << ">> exiting intersection_with_quadric22" << endl;
	#endif

	return ic;
}

// Input : initial quadrics q1, q2
//				 determinental equation of the pencil : det_p
//				 a quadric 22	 q
// Output : a quadric 22	quadric_sol and a point on it point_sol
//						the determinant a^2.b of the quadric 22 in the form [a,b]
void find_quadric_22_and_point(const bigint_matrix &q1, const bigint_matrix &q2, 
																const hom_polynomial <bigint> &det_p, 
																const bigint_matrix &q_cp, 
																bigint_matrix &quadric_sol, 
																math_vector <bigint> &point_sol, 
																bigint &det_R,
																math_vector <bigint> &det_sol, const int opt_level, 
																ostream &s) 
{	 
	#ifndef NDEBUG
	s << ">> entering find_quadric_22_and_point" << endl;
	#endif

	// q_cp is const: make a copy
	bigint_matrix q = q_cp;

	// bigfloat_precision is the length of the mantissa (in base 10)
	unsigned int bigfloat_precision = 5;
	unsigned int increment = 0;
	bool not_found = 1;
	math_vector <bigint> srf (2,2);
	math_vector <bigint> coeffs_sol (2,2);
	bigint sol_a2k; 
	// axis is a 4x2 matrix such that axis.(m,n) is a parametrization of one of the
	// 6 axis in P^3, with (m,n) in P^1. axis is initialized to 0
	bigint_matrix tm; 
	math_matrix <bigfloat> tm_float (4,4);
	math_matrix <bigfloat> axis (4,2);
	axis.sto(0,0,0);	axis.sto(0,1,0);	
	axis.sto(1,0,0);	axis.sto(1,1,0);	
	axis.sto(2,0,0);	axis.sto(2,1,0);	
	axis.sto(3,0,0);	axis.sto(3,1,0);	

	bigfloat a, b, c, Delta, real_point_max_abs;
	math_vector <bigfloat> m_sols (2,2), n_sols (2,2);
	math_vector <bigfloat> real_point (4,4), mn_sol (2,2);
	math_vector <bigint> approximated_point (4,4), tmp3 (4,4), tmp4 (1,1);
	unsigned int k = 0, k_current_precision; 
	bigint det_new_R,tmp_bi, cx, cy, gcd_c, multiply_factor=1, multiply_factor_current_precision; 

	while (not_found)
		{
			bigfloat::set_precision(bigfloat_precision);
			k_current_precision = k;
			multiply_factor_current_precision = multiply_factor;

			if (increment == 1) 
				{	
					gauss(q, tm, q); // replace q by its canonical form 
					k_current_precision =0; // reinitialize the current_precision to 0
					multiply_factor_current_precision=1;	 // and the corresponding
					// multiply_factor_current_precision 
				}
			if (increment) // if increment>0
				for (unsigned int	i = 0; i < 4; i++)
					for (unsigned int j = 0; j < 4; j++)
						tm_float.sto(i,j,bigfloat(tm.member(i,j)));

			for (unsigned int i = 0; i < 4; i++)
				{
		 			axis.sto(i,0,1);
		 			for (unsigned int j = i+1; j < 4; j++)
						{
							axis.sto(j,1,1); 
							// axis.(m,n) = (m,n,0,0), (m,0,n,0), (m,0,0,n), (0,m,n,0),
							// (0,m,0,n) or (0,0,m,n) 
							
							// (axis.(m,n))^t . q . axis.(m,n) = a.m^2 + 2.b.m.n +c.n^2
							a = bigfloat (q.member(i,i));
							b = bigfloat (q.member(i,j));
							c = bigfloat (q.member(j,j));
							Delta = b*b-a*c;
			
							if ((Delta>= 0) && !(a.is_zero() && b.is_zero() && c.is_zero()))
								{
									// Compute the solutions (m,n) of the 2nd degree equation
									if (!(a.is_zero()))
										{
											m_sols[0] = -b + sqrt(Delta);
											m_sols[1] = -b - sqrt(Delta);
											n_sols[0] = a;
											n_sols[1] = a;
										}
									else 
										{
											m_sols[0] = 1;
											m_sols[1] = -c;
											n_sols[0] = 0;
											n_sols[1] = 2*b;
										}
						
									for (unsigned int sol_index = 0; sol_index < 2; sol_index++) 
										{
											// real_point =
											// (axis.(m_sols[sol_index],n_sols[sol_index])) +
											// normalization ou	tm_float * axis * (m_sols[sol_index],
											// n_sols[sol_index]) + normalization 
											mn_sol[0] = m_sols[sol_index];
											mn_sol[1] = n_sols[sol_index];
											multiply(real_point, axis,mn_sol); 
											if (increment) // if increment>0
									 			multiply(real_point, tm_float, real_point); 
											real_point_max_abs = abs(real_point[0]) ;
											if (real_point_max_abs < abs(real_point[1])) 
									 			real_point_max_abs = abs(real_point[1]);
											if (real_point_max_abs < abs(real_point[2])) 
									 			real_point_max_abs = abs(real_point[2]);
											if (real_point_max_abs < abs(real_point[3])) 
									 			real_point_max_abs = abs(real_point[3]);
											real_point = real_point/real_point_max_abs;
							
											// We now approzimate real_point
											for (unsigned int floor_ceil = 0; floor_ceil < 2; floor_ceil++) 
									 			{
													k = k_current_precision;
													multiply_factor = multiply_factor_current_precision;
													det_new_R =0;
													while ((det_new_R <= 0) && (k<=bigfloat_precision))
														{
															if (floor_ceil == 0)
										 						{
																	floor(tmp_bi, real_point[0]*multiply_factor);
																	approximated_point[0] = tmp_bi;
																	floor(tmp_bi, real_point[1]*multiply_factor);
																	approximated_point[1] = tmp_bi;
																	floor(tmp_bi, real_point[2]*multiply_factor);
																	approximated_point[2] = tmp_bi;
																	floor(tmp_bi, real_point[3]*multiply_factor);
																	approximated_point[3] = tmp_bi;
										 						}
															else
																{
																	ceil(tmp_bi, real_point[0]*multiply_factor);
																	approximated_point[0] = tmp_bi;
																	ceil(tmp_bi, real_point[1]*multiply_factor);
																	approximated_point[1] = tmp_bi;
																	ceil(tmp_bi, real_point[2]*multiply_factor);
																	approximated_point[2] = tmp_bi;
																	ceil(tmp_bi, real_point[3]*multiply_factor);
																	approximated_point[3] = tmp_bi;
																}
															// The approximated point is	computed
															// We now compute cx, cy such that 
															// (approx_point)^t. (x.q1 + y.q2). approx_point =
															// cx.x +cy.y 
															multiply(tmp3,q1,approximated_point); 
															// original version
															//multiply(tmp4,trans(approximated_point), tmp3);
															//cx = tmp4[0];
															cx = inner_prod(approximated_point, tmp3);
															multiply(tmp3,q2,approximated_point);
															//multiply(tmp4,trans(approximated_point), tmp3);
															//cy = tmp4[0];
															cy = inner_prod(approximated_point, tmp3);
								
															if ((cx != 0) && (cy != 0))
																{
																	gcd_c = gcd(cx,cy);
																	divide(cx, cx, gcd_c);
																	divide(cy, cy, gcd_c);
																}
												
															// The quadric through approx_point is -cy.q1 + cx.q2
															// unless cx=cy=0 in which case q1 and q2 contain
															// approx_point; 
															// Determinant of the quadric -cy.q1 + cx.q2
															// through approx_point is: 
															if ((cx!=0) || (cy!=0))
												 				det_new_R = det_p.eval(-cy,cx);
															else // cx=cy=0
																{
																	cx = 1; 
																	det_new_R = det_p.eval(-cy,cx);
																	if (det_new_R <=0)
																 	 {
																 		 cx = 0; 
																 		 cy = 1;
																 		 det_new_R = det_p.eval(-cy,cx);
																 	 }
																}
							
															if (det_new_R>0) 
																{
													 				srf = extract_square_factors(det_new_R,opt_level,s); 
													 				tmp_bi = srf[0] * multiply_factor; // tmp_bi =
													 				// srf[0] * 10^k 
								 									if (not_found==1)
												 						{
																			not_found = 0;
																			det_sol.assign(srf);
																			det_R = det_new_R;
																			sol_a2k = tmp_bi;
																			coeffs_sol[0] = -cy;
																			coeffs_sol[1] = cx;
																			point_sol.assign(approximated_point);
																		}
																	else
																		{
																			if ((srf[1]<det_sol[1]) ||
																				(srf[1]==det_sol[1] && tmp_bi<sol_a2k))
																	 			{
																					det_sol.assign(srf);	
																					det_R = det_new_R;
																					sol_a2k = tmp_bi;
																					coeffs_sol[0] = -cy;
																					coeffs_sol[1] = cx;
																					point_sol.assign(approximated_point);
																	 			}
						 												}
																} // end if (det_new_R>0) 
															k = k+1;
															multiply_factor = 10*multiply_factor;
														} // end while (det_new_R <= 0 && k<bigfloat_precision)

			 									} // end for	(unsigned int floor_ceil = 0; floor_ceil <
					 								// 2; floor_ceil++) 
				 						} // for	 (unsigned int sol_index = 0; sol_index < 2; sol_index++)
								} // end if (Delta>= 0 && !(a.is_zero() && b.is_zero() && c.is_zero()))
				 			axis.sto(j,1,0); // re-initialize the coeff axis[j,1] to 0
			 			} // end of the loop for (unsigned int j = i+1; j < 4; j++)
		 			axis.sto(i,0,0); // re-initialize the coeff axis[i,0] to 0
				} // end of the loop for (unsigned int i = 0; i < 4; i++)

			if (not_found) // no quadric 22 has been found
				{
					increment++;
					if (increment > 1)
						 bigfloat_precision = bigfloat_precision + 10;
					
					#ifndef NDEBUG
					s << "INCREMENT, bigfloat_precision	 " << increment << "		"
						 << bigfloat_precision << endl;
					// I have found an example where this had to be set to 50 and it worked!
					// Well of course then, the size of the rational point had a big
					// influence on the height of the output
					#endif
				}
			else	// a quadric 22 has been found
				{
					// The matrix associated to p
					bigint_matrix q1_tmp;
					q1_tmp = q1 * coeffs_sol[0];
					bigint_matrix q2_tmp;
					q2_tmp = q2 * coeffs_sol[1];
					quadric_sol = q1_tmp + q2_tmp;
			
					#ifndef NDEBUG
					optimize(coeffs_sol);
			
					s << ">> ending up using lambda " << coeffs_sol << endl;
					#endif
				}
		} // end while (not_found)

	optimize(point_sol);
	
	#ifndef NDEBUG
	s << ">> exiting find_quadric_22_and_point" << endl;
	#endif
} // find_quadric_22_and_point

// Input:	 a matrix q of inertia (2,2)
//		a point on it and the determinant of q in the form [a,b] such that det(q) = a^2.b
// Output : Parameterization of the quadric : s1 + sqrt(det_q[1]) s2
// If det_q[1] = 1 then s2 = [0,0,0,0]
void find_parameterization_quadric_22(const math_matrix <bigint> &q_cp, 
																			const math_vector <bigint> &point_on_q, 
																			const math_vector <bigint> &det_q_cp,
																			const bigint &det_R,
																			surface_param <bigint> &s1, 
																			surface_param <bigint> &s2,
																			const int opt_level, ostream &s)	
{
	#ifndef NDEBUG
	s << ">> entering find_parameterization_quadric_22" << endl;
	#endif
	
	// q_cp is const: make a copy
	math_matrix <bigint> q = q_cp;	

	// det_q is const: make a copy
	math_vector <bigint> det_q = det_q_cp;
	math_vector <bigint> p0 = point_on_q;

	// The transformation matrix
	// First find a point of the canonical basis that is not on the tangent plane to
	// q at p0
	math_vector <bigint> Rp0 = prod(q, p0), v(4,4);
	
	// Try (1,0,0,0), (0,1,0,0), (0,0,1,0) and (0,0,0,1) in that order -- At least
	// one of the four points is good
	v[0] = 1;

	bigint m0 = inner_prod(Rp0, v);
			
	if (m0.is_zero()) // (1,0,0,0) is no good
		{
			v[0] = 0;
			v[1] = 1;

			m0 = inner_prod(Rp0, v);

			if (m0.is_zero()) // (0,1,0,0) is no good
				{
					v[1] = 0;
					v[2] = 1;
			
					m0 = inner_prod(Rp0, v);
					 	 
					if (m0.is_zero()) // (0,0,1,0) is no good
						{
					 	 	v[2] = 0;
					 	 	v[3] = 1;
			
					 	 	m0 = inner_prod(Rp0, v); // (0,0,0,1) is good!
						}
				}
		}

	bigint_matrix P = send_to_zw(v,p0);

	// The determinant of the transformation matrix
	bigint theta = det(P);
	
	// The transformed matrix
	bigint_matrix Rp = prod(base_matrix<bigint>(prod(trans(P), q)), P);

	// Entries of Rp
	bigint a1 = Rp.member(0,0), a2 = Rp.member(0,1), a3 = Rp.member(1,1), a4 = Rp.member(0,2), 
		a5 = Rp.member(1,2), a6 = Rp.member(2,2), a7 = Rp.member(0,3), a8 = Rp.member(1,3),
		a9 = Rp.member(2,3);

	// Minors
	bigint l1 = a1*a9-a4*a7, l2 = a6*a7-a4*a9, l3 = a2*a9-a5*a7;

	// The coefficients alpha and beta
	bigint alpha = a9*l1+a7*l2, beta = a9*l3+a8*l2;

	// The matrices of the result
	bigint_matrix m1(4,4), m2(4,4);
	
	// We have two cases according to whether alpha is 0 or not
	if (alpha == 0)
		{
			bigint tp = abs(theta)*det_q[0]*beta.sign()*a9.sign();

			bigint gamma = a3*a9*a9+a6*a8*a8-2*a5*a8*a9;

			bigint sig1 = 2*beta*a8-gamma*a7;
			bigint sig2 = (2*beta*(a5*a9-a6*a8)+l2*gamma)/a9;

			m1.sto(0,2,gamma*a9);
			m1.sto(0,3,a9*beta);
			m1.sto(1,2,-2*m1(0,3));
			m1.sto(2,1,-m1(1,2));
			m1.sto(2,2,sig1);
			m1.sto(2,3,-a7*beta);
			m1.sto(3,0,theta*theta*det_R);
			m1.sto(3,1,-beta*a6);
			m1.sto(3,2,sig2);
			m1.sto(3,3,l2*tp);
		}
	else
		{
			// Coefficients tau1 and tau2
			bigint tau1 = a8*l1-a7*l3;
			bigint tau2 = a5*l1+a2*l2+a8*(a4*a4-a1*a6);

			// Param matrices
			m1.sto(0,2,a9*beta);
			m1.sto(0,3,m1(0,2));
			m1.sto(1,2,-a9*alpha);
			m1.sto(1,3,m1(1,2));
			m1.sto(2,1,-2*m1(1,2));
			m1.sto(2,2,a9*tau1);
			m1.sto(2,3,m1(2,2));
			m1.sto(3,0,theta*theta*det_R);
			m1.sto(3,1,-a6*alpha);
			m1.sto(3,2,a9*tau2);
			m1.sto(3,3,m1(3,2));

			m2.sto(0,2,-a9*a9);
			m2.sto(0,3,-m2(0,2));
			m2.sto(2,2,a7*a9);
			m2.sto(2,3,-m2(2,2));
			m2.sto(3,2,-l2);
			m2.sto(3,3,l2);
		}

	// Put back in global frame
	m1 = prod(P, m1);
	m2 = prod(P, m2);

	// m2 should be multiplied by this
	bigint cm2 = theta * det_q[0];

	// sp = [u*t, v*s, u*s, v*t]
	surface_param <bigint> sp(4);
	surface_22_param(sp);

	bigint A = content(column(m1, 0)), B = content(column(m1,1));

	// Now simplification
	if (det_q[1] == 1)
		{
			m1 = m1 + cm2 * m2;

			// Find the factors a, b, c, d of the 4 columns of M such that a*b = c*d

			// First, the gcds of the columns
			bigint C = content(column(m1, 2)), D = content(column(m1, 3));

			bigint a = gcd(A,D), d = a, b = gcd(B,C), c = b;

			A = A/a;
			B = B/b;
			C = C/c;
			D = D/d;

			bigint gcdBD = gcd(B,D); 
			b = b*gcdBD;
			d = d*gcdBD;
			bigint gcdAC = gcd(A,C); 
			a = a*gcdAC;
			c = c*gcdAC;

			column(m1, 0) = column(m1, 0) / a;
			column(m1, 1) = column(m1, 1) / b;
			column(m1, 2) = column(m1, 2) / c;
			column(m1, 3) = column(m1, 3) / d;

			load_balancing(m1, 1, s);
			multiply(s1, m1, sp);
		}
	else
		{
			// Find the factors a, b, c, d of the 4 columns of M such that a*b = c*d
			// This time with two matrices!

			bigint C = gcd(content(column(m1, 2)),cm2 * content(column(m2, 2)));

			bigint gcdAB = gcd(A,B);
			bigint a0 = extract_square_factors(A/gcdAB,opt_level,s)[0]; 
			bigint b0 = extract_square_factors(B/gcdAB,opt_level,s)[0]; 
			// A=a0^2.a1.gcdAB and B=b0^2.b1.gcdAB
		 
			bigint gamma = gcd(gcdAB,C);
			bigint alpha = gcd(a0, C/gamma);
			bigint beta = gcd(b0, C/gamma);
			bigint a = alpha*alpha*gamma;
			bigint b = beta*beta*gamma;
			bigint c = alpha*beta*gamma;
			// a|A, b|B, sqrt(ab)|C because gamma | A, B and C,	 alpha^2 | A, 
			// alpha | C/gamma, beta^2 | B, beta | C/gamma, and alpha and beta
			// are relatively prime (since a0 and b0 are)

			column(m1, 0) = column(m1, 0) / a;
			column(m1, 1) = column(m1, 1) / b;
			column(m1, 2) = column(m1, 2) / c;
			column(m1, 3) = column(m1, 3) / c;

			bigint_matrix m2tmp = cm2*m2;

			column(m2, 2) = column(m2tmp, 2) / c;
			column(m2, 3) = column(m2tmp, 3) / c;

			// We can still do load_balancing on the first two columns of m1
			load_balancing(m1,0,s);

			QI::multiply(s1, m1, sp);
			QI::multiply(s2, m2, sp);
		}

	#ifndef NDEBUG
	s << ">> param of quadric (2,2): " << s1;
	if (det_q[1] != 1)
		s << " + sqrt(" << det_q[1] << ")*" << s2;
	s << endl;

	check_param(q, s1, s2, det_q[1],s);
	#endif
	
	#ifndef NDEBUG
	s << ">> exiting find_parameterization_quadric_22" << endl;
	#endif
}

} // end of namespace QI
