// Parsing of cases: what is the type of the intersection?

#include <libqi/kernel/QIParam.h>
#include <libqi/kernel/QINumber.h>
#include <libqi/kernel/QIElem.h>
//#include <libqi/kernel/QILegendre.h>

#include <string>

using namespace rpl;

// Enter namespace QI
namespace QI {

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////// Parameterization of (2,2) quadrics

// Output the parameterization of the quadric of inertia (2,2) of implicit
// equation xy - zw = 0	 :	s = [u.t, v.s, u.s, v.t]	
void surface_22_param(surface_param <bigint> &s)
{
	// Polynomials u and v
	hom_polynomial <bigint> u_pol,v_pol;
	u_pol.assign_x();
	v_pol.assign_y();
	
	// cu = [t, 0, s, 0] and cv = [0, s, 0, t]
	curve_param <bigint> cu(4),cv(4);
	cu[0].assign_y();
	cu[1].assign_zero();
	cu[2].assign_x();
	cu[3].assign_zero();
	cv[0].assign_zero();
	cv[1].assign_x();
	cv[2].assign_zero();
	cv[3].assign_y();

	// Needed for temporary storage
	surface_param <bigint> tmp(4);

	QI::multiply(s, cu, u_pol);
	QI::multiply(tmp, cv, v_pol);
	add(s,s,tmp);
}

// Print the bi-quadratic homogeneous polynomial 
void print_bi_quadratic_hom_poly(const hom_hom_polynomial <bigint> &poly1, 
																	const hom_hom_polynomial <bigint> &poly2, 
																	const bigint &delta, ostream &s)
{
	s << ">> bi-quadratic homogeneous polynomial:" << endl;

	if (!poly1.is_zero())
		s << poly1;
	if (!poly2.is_zero())
		{
			if	(!poly1.is_zero())
				s << " + "; 
			s << "sqrt(" << delta << ") * (" << poly2 << ")";
		}
	s << endl;
}

// Plug a surface_param in the equation of a quadric
hom_hom_polynomial <bigint> plug_param_in_quadric(const surface_param <bigint> &s1,
																									const bigint_matrix &q, 
																									const surface_param <bigint> &s2)
{
	surface_param <bigint> tmp(4);
	multiply(tmp,q,s2);
	
	hom_hom_polynomial <bigint> res;
	multiply(res,s1,tmp);
	
	return res;
}

// Plug a curve_param in the equation of a quadric
hom_polynomial <bigint> plug_param_in_quadric(const curve_param <bigint> &s1,
																							const bigint_matrix &q, 
																							const curve_param <bigint> &s2)
{
	curve_param <bigint> tmp(q.get_no_of_columns());

	multiply(tmp,q,s2);
	hom_polynomial <bigint> res;
	multiply(res,s1,tmp);

	return res;
}

///////////// Parameterization of conics

// The structure is this:
// conic_param is the high-level function. If opt_level = 1, look for a rational
// point (FindRationalPointOnConic). If such a point is found, launch
// conic_param_through_ratpoint, otherwise launch conic_param_no_ratpoint.
// If opt_level = 0, launch conic_param_no_opt, were some heuristics try to find a
// rational point anyway and branch accordingly.

// Compute the parameterization of a non-singular conic through a rational point
// Conic should be in the plane (x,y,z)
void conic_param_through_ratpoint(const bigint_matrix &q, 
																	const math_vector <bigint> &rat_point,
																	bigint_matrix &p_trans, curve_param <bigint> &p,
																	math_vector <bigint> &l, ostream &s)
{
	#ifndef NDEBUG
	s << ">> parameterization of conic with rational point" << endl;
	#endif

	if ((rat_point[0].is_zero()) && (rat_point[1].is_zero()))
		{
			// No need to rotate
			conic_param_through_origin(q,p_trans,p,l);
		}
	else if ((rat_point[0].is_zero()) && (rat_point[2].is_zero()))
		{
			// The conic goes through [0,1,0]: switch y and z
			bigint_matrix p_trans2(3,3);
			p_trans2.sto(0, 0,1);
			p_trans2.sto(1, 2, 1);
			p_trans2.sto(2, 1, 1);
			conic_param_through_origin(prod(trans(p_trans2), base_matrix<bigint>( prod(q,p_trans2))),
																 p_trans, p, l);
			p_trans = prod(p_trans2, p_trans);
		}
	else if ((rat_point[1].is_zero()) && (rat_point[2].is_zero())) 
		{
			// The conic goes through [1,0,0]: switch x and z
			bigint_matrix p_trans2(3,3);
			p_trans2.sto(1,1,1);
			p_trans2.sto(0,2,1);
			p_trans2.sto(2,0,1);
			conic_param_through_origin(prod(trans(p_trans2), base_matrix<bigint> (prod(q,p_trans2))), 
																 p_trans, p, l);
			p_trans = prod(p_trans2, p_trans);
		}
	else
		{
			// Rotate and param
			bigint_matrix proj_mat = send_to_infinity(rat_point);

			bigint_matrix q_r;
			multiply(q_r, q, proj_mat);
			multiply(q_r, (math_matrix<bigint>)trans(proj_mat), q_r);

			conic_param_through_origin(q_r, p_trans, p, l);

			multiply(p_trans, proj_mat, p_trans);
		}
}

// Compute the parameterization of a non-singular conic going through the point [0,0,1]
void conic_param_through_origin(const bigint_matrix &q, bigint_matrix &p_trans, 
																curve_param <bigint> &p, math_vector <bigint> &l)
{
	// Conic is a x^2 + b xy + c y^2 + d yz + e xz
	bigint a = q(0,0);
	bigint b = 2*q(0,1);
	bigint c = q(1,1);
	bigint d = 2*q(1,2);
	bigint e = 2*q(0,2);

	// The curve_param (u^2, v^2, uv)
	p[0].assign_x2();
	p[1].assign_y2();
	p[2].assign_xy();

	// Param is u^2 [e,0,-a] + v^2 [0,d,-c] + uv [d,e,-b]
	p_trans.sto(0,0,e);
	p_trans.sto(2,0,-a);
	p_trans.sto(1,1,d);
	p_trans.sto(2,1,-c);
	p_trans.sto(0,2,d);
	p_trans.sto(1,2,e);
	p_trans.sto(2,2,-b);

	// The ``line'' values (parameters for the ratpoint)
	l[0] = q(0,2);
	l[1] = q(1,2);

	optimize(l);
}

// Parameterization by discriminant D = e^2-a*f > 0
// The param is (p_trans + sqrt(D)*p_trans2) * [u^2 v^2 uv]
// par = [u^2 v^2 uv]
void conic_param1(const bigint_matrix &q, bigint &D, const math_vector <bigint> &nd, 
									bigint_matrix &p_trans, bigint_matrix &p_trans2,
									curve_param <bigint> &par) 
{
	// Conic is ax^2 + 2bxy + cy^2 + 2dyz + 2exz + fz^2
	// Inertia is supposed to be [2 1]
	bigint a = q(0,0);
	bigint b = q(0,1);
	bigint c = q(1,1);
	bigint d = q(1,2);
	bigint e = q(0,2);
	bigint f = q(2,2);

	p_trans.sto(0,0,a*e);
	p_trans.sto(0,1,c*e);
	p_trans.sto(0,2,2*a*d);
	p_trans.sto(1,1,2*(a*d-b*e));
	p_trans.sto(2,0,-a*a);
	p_trans.sto(2,1,-a*c);
	p_trans.sto(2,2,-2*a*b);

	p_trans2.sto(0,0,-a);
	p_trans2.sto(0,1,c);
	p_trans2.sto(1,1,-2*b);
	p_trans2.sto(1,2,-2*a);
	p_trans2 = p_trans2*nd[0];

	// The param
	par[0].assign_x2();
	par[1].assign_y2();
	par[2].assign_xy();

	if (nd[1] == 1)
		{
			// Rational param
			D = 0;

			p_trans = p_trans+p_trans2;
		}
	else
		D = nd[1];
}

// Parameterization by -f*p > 0
// Equation is (a*f-e^2)*(f*z+e*x+d*y)^2 + (x*(a*f-e^2)+y*(b*f-d*e))^2 + f*p*y^2 = 0
// The param is (p_trans + sqrt(D)*p_trans2) * [u^2 v^2 uv]
// par = [u^2 v^2 uv]
void conic_param2(const bigint_matrix &q, bigint &D, const bigint &de, const bigint &p3, 
									const math_vector <bigint> &nd, const math_vector <bigint> &np, 
									bigint_matrix &p_trans, bigint_matrix &p_trans2, 
									curve_param <bigint> &par, const int opt_level) 
{
	// Conic is ax^2 + 2bxy + cy^2 + 2dyz + 2exz + fz^2
	// Inertia is supposed to be [2 1]
	bigint a = q(0,0);
	bigint b = q(0,1);
	bigint c = q(1,1);
	bigint d = q(1,2);
	bigint e = q(0,2);
	bigint f = q(2,2);

	bigint t1 = b*f-d*e;
	bigint t2 = a*d-b*e;

	// e^2-a*f < 0
	if (nd[1] < 0)
		{
			p_trans.sto(0,0,f*t1);
			p_trans.sto(0,1,-p_trans(0,0)*nd[1]);
			p_trans.sto(1,0,de*f);
			p_trans.sto(1,1,-p_trans(1,0)*nd[1]);
			p_trans.sto(2,0,f*t2);
			p_trans.sto(2,1,-p_trans(2,0)*nd[1]);

			p_trans2.sto(0,0,-f);
			p_trans2.sto(0,1,p_trans2(0,0)*nd[1]);
			p_trans2.sto(2,0,e);
			p_trans2.sto(2,1,p_trans2(2,0)*nd[1]);
			p_trans2.sto(2,2,2*nd[0]*nd[1]);
			p_trans2 = p_trans2*np[0];
		}
	else
		{
			math_vector <bigint> vf = extract_square_factors(de*np[1],opt_level,cout);
			
			p_trans.sto(0,0,f*p3);
			p_trans.sto(0,1,p_trans(0,0)*vf[1]);
			p_trans.sto(2,0,-e*p3);
			p_trans.sto(2,1,p_trans(2,0)*vf[1]);

			p_trans2.sto(0,0,f*t1);
			p_trans2.sto(0,1,-p_trans2(0,0)*vf[1]);
			p_trans2.sto(1,0,f*de);
			p_trans2.sto(1,1,-p_trans2(1,0)*vf[1]);
			p_trans2.sto(2,0,f*t2);
			p_trans2.sto(2,1,-p_trans2(2,0)*vf[1]);
			p_trans2.sto(2,2,-2*np[0]*vf[0]*vf[1]);
			p_trans2 = p_trans2*np[0];
		}

	// The param
	par[0].assign_x2();
	par[1].assign_y2();
	par[2].assign_xy();

	if (np[1] == 1)
		{
			// Rational param
			D = 0;

			p_trans = p_trans+p_trans2;
		}
	else
		D = np[1];
}

// Compute the parameterization of a non-singular conic when there is no rational
// point or we were unable to find a simple one
// The optional parameters come from conic_param_no_opt.
void conic_param_no_ratpoint(const bigint_matrix &q, bigint &D,	 
															bigint_matrix &p_trans, bigint_matrix &p_trans2,	
															curve_param <bigint> &par, const int opt_level, ostream &s,
															const bigint &dxyv = 0, const bigint &dxzv = 0, 
															const bigint &dyzv = 0,
															const bigint &p1v = 0, const bigint &p2v = 0, const bigint &p3v = 0) 
{
	// If opt_level = 1, we know there is no rational point. Otherwise, there still
	// is a chance that we find one such point.

	bigint dxy,dxz,dyz,p1,p2,p3;

	if (opt_level)
		{
			// The optional parameters were not initialized
			bigint a = q(0,0); bigint b = q(0,1); bigint c = q(1,1);
			bigint d = q(1,2); bigint e = q(0,2); bigint f = q(2,2);

			dxy = b*b-a*c; dxz = e*e-a*f; dyz = d*d-c*f;

			bigint p = det(q);
			p1 = -a*p; p2 = -c*p; p3 = -f*p;
		}
	else
		{
			dxy = dxyv; dxz = dxzv; dyz = dyzv; p1 = p1v; p2 = p2v; p3 = p3v;
		}

	// Try to find the smallest square root

	math_vector <bigint> v1(2,2),v2(2,2),v3(2,2),v4(2,2),v5(2,2),v6(2,2);
	v1[0] = 1; v2[0] = 1; v3[0] = 1; v4[0] = 1; v5[0] = 1; v6[0] = 1;
	v1[1] = dxy;
	v2[1] = dxz;
	v3[1] = dyz;
	v4[1] = p1;
	v5[1] = p2;
	v6[1] = p3;
		 
	// Optimize	 dxy, dxz, dyz, p1, p2, p3
	if (v1[1] > 0)
		v1 = extract_square_factors(v1[1],opt_level,s);
	if (v2[1] > 0)
		v2 = extract_square_factors(v2[1],opt_level,s);
	if (v3[1] > 0)
		v3 = extract_square_factors(v3[1],opt_level,s);
	if (v4[1] > 0)
		v4 = extract_square_factors(v4[1],opt_level,s);
	if (v5[1] > 0)
		v5 = extract_square_factors(v5[1],opt_level,s);
	if (v6[1] > 0)
		v6 = extract_square_factors(v6[1],opt_level,s);

	int id;
	bigint sq;
	// Find the smallest square root among dxy, dxz, dyz, p1, p2, p3
	// If dxy < 0, then p1 > 0 and conversely (initialization)
	if (v1[1] < 0)
		{
			id = 4;
			sq = v4[1];
		}
	else
		{
			id = 1;
			sq = v1[1];
		}

	if ((v2[1] < sq) && (v2[1] > 0))
		{
			id = 2;
			sq = v2[1];
		}
	if ((v3[1] < sq) && (v3[1] > 0))
		{
			id = 3;
			sq = v3[1];
		}
	if ((id != 4) && (v4[1] < sq) && (v4[1] > 0))
		{
			id = 4;
			sq = v4[1];
		}
	if ((v5[1] < sq) && (v5[1] > 0))
		{
			id = 5;
			sq = v5[1];
		}
	if ((v6[1] < sq) && (v6[1] > 0))
		id = 6;

	if (id == 1)
		{
			// Switch y and z
			bigint_matrix tr(3,3);
			tr.sto(0,0,1);
			tr.sto(1,2,1);
			tr.sto(2,1,1);

			conic_param1(prod(tr, prod<base_matrix<bigint> >(q, tr)) ,D,v1,p_trans,p_trans2,par);

			p_trans = prod(tr, p_trans);
			p_trans2 = prod(tr, p_trans2);
		}
	else if (id == 2)
		conic_param1(q,D,v2,p_trans,p_trans2,par);
	else if (id == 3)
		{
			// Switch x and y
			bigint_matrix tr(3,3);
			tr.sto(0,1,1);
			tr.sto(1,0,1);
			tr.sto(2,2,1);
			
			conic_param1(prod(tr, prod<base_matrix<bigint> >(q, tr)), D, v3, p_trans,p_trans2,par);
			
			p_trans = prod(tr, p_trans);
			p_trans2 = prod(tr, p_trans2);			
		}
	else if (id == 4)
		{
			// Switch x and z
			bigint_matrix tr(3,3);
			tr.sto(0,2,1);
			tr.sto(1,1,1);
			tr.sto(2,0,1);
			
			conic_param2(prod(tr,prod<base_matrix<bigint> >(q,tr)) ,D,dxz,p1,v2,v4,p_trans,p_trans2,par,opt_level);
			
			p_trans = prod(tr, p_trans);
			p_trans2 = prod(tr, p_trans2);
		}
	else if (id == 5)
		{
			// Switch y and z
			bigint_matrix tr(3,3);
			tr.sto(0,0,1);
			tr.sto(1,2,1);
			tr.sto(2,1,1);
			
			conic_param2(prod(tr,prod<base_matrix<bigint> >(q,tr)), D, dxy, p2, v1, v5, p_trans, p_trans2, par, opt_level);
			
			p_trans = prod(tr, p_trans);
			p_trans2 = prod(tr, p_trans2);
		}
	else // id = 6
		conic_param2(q,D,dxz,p3,v2,v6,p_trans,p_trans2,par,opt_level);
}

// The procedure for parameterizing conics when opt_level = 0
void conic_param_no_opt(const bigint_matrix &q, bigint &D,	
												bigint_matrix &p_trans, bigint_matrix &p_trans2,	 
												curve_param <bigint> &par, ostream &s) 
{
	// If one of the diagonal elements is zero, we know how to get a rational param
	if (q(0,0) == 0)
		{
			#ifndef NDEBUG
			s << ">> found conic through (1 0 0)!" << endl;
			#endif

			math_vector <bigint> rat_point(3,3);
			rat_point[0].assign_one();

			math_vector <bigint> l(2,2);
			conic_param_through_ratpoint(q,rat_point,p_trans,par,l,s);
		}
	else if (q(1,1) == 0)
		{
			#ifndef NDEBUG
			s << ">> found conic through (0 1 0)!" << endl;
			#endif

			math_vector <bigint> rat_point(3,3);
			rat_point[1].assign_one();

			math_vector <bigint> l(2,2);
			conic_param_through_ratpoint(q,rat_point,p_trans,par,l,s);
		}
	else if (q(2,2) == 0)
		{
			#ifndef NDEBUG
			s << ">> found conic through (0 0 1)!" << endl;
			#endif

			math_vector <bigint> rat_point(3,3);
			rat_point[2].assign_one();

			math_vector <bigint> l(2,2);
			conic_param_through_ratpoint(q,rat_point,p_trans,par,l,s);
		}
	else
		{
			// Now, if one of the following 2x2 discriminants is a square, we also have
			// a rational param. 
			bigint a = q(0,0); bigint b = q(0,1); bigint c = q(1,1);
			bigint d = q(1,2); bigint e = q(0,2); bigint f = q(2,2);

			bigint dxy = b*b-a*c; bigint dxz = e*e-a*f; bigint dyz = d*d-c*f;

			bigint p = det(q);
			bigint p1 = -a*p; bigint p2 = -c*p; bigint p3 = -f*p;

			bigint sq;

			if ((dxy > 0) && (is_square(sq,dxy)))
				{
					#ifndef NDEBUG
					s << ">> dxy is a square!" << endl;
					#endif
			
					// Switch y and z
					bigint_matrix tr(3,3);
					tr.sto(0,0,1);
					tr.sto(1,2,1);
					tr.sto(2,1,1);
			
					math_vector <bigint> v(2,2);
					v[0] = sq;
					v[1] = 1;
					conic_param1(prod(tr, prod<base_matrix<bigint> >(q,tr)) ,D,v,p_trans,p_trans2,par);
			
					p_trans = prod(tr, p_trans);
					p_trans2 = prod(tr, p_trans2);
				}
			else if (is_square(sq,dxz))
				{
					#ifndef NDEBUG
					s << ">> dxz is a square!" << endl;
					#endif
			
					math_vector <bigint> v(2,2);
					v[0] = sq;
					v[1] = 1;
					conic_param1(q,D,v,p_trans,p_trans2,par);
				}
			else if (is_square(sq,dyz))
				{
					#ifndef NDEBUG
					s << ">> dyz is a square!" << endl;
					#endif
			
					// Switch x and y
					bigint_matrix tr(3,3);
					tr.sto(0,1,1);
					tr.sto(1,0,1);
					tr.sto(2,2,1);
			
					math_vector <bigint> v(2,2);
					v[0] = sq;
					v[1] = 1;
					conic_param1(prod(tr,prod<base_matrix<bigint> >(q,tr)),D,v,p_trans,p_trans2,par);
			
					p_trans = prod(tr, p_trans);
					p_trans2 = prod(tr, p_trans2);
				}
			else if (is_square(sq,p1))
				{
					#ifndef NDEBUG
					s << ">> p1 is a square!" << endl;
					#endif
			
					// Switch x and z
					bigint_matrix tr(3,3);
					tr.sto(0,2,1);
					tr.sto(1,1,1);
					tr.sto(2,0,1);
			
					math_vector <bigint> v(2,2);
					v[0] = sq;
					v[1] = 1;
			
					math_vector <bigint> v2 = extract_square_factors(dxz,0,s);
			
					conic_param2(prod(tr,prod<base_matrix<bigint> >(q,tr)),D,dxz,p1,v2,v,p_trans,p_trans2,par,0);
			
					p_trans = prod(tr, p_trans);
					p_trans2 = prod(tr, p_trans2);
				}
			else if (is_square(sq,p2))
				{
					#ifndef NDEBUG
					s << ">> p2 is a square!" << endl;
					#endif
			
					// Switch y and z
					bigint_matrix tr(3,3);
					tr.sto(0,0,1);
					tr.sto(1,2,1);
					tr.sto(2,1,1);
			
					math_vector <bigint> v(2,2);
					v[0] = sq;
					v[1] = 1;
			
					math_vector <bigint> v1 = extract_square_factors(dxy,0,s);
			
					conic_param2(prod(tr,prod<base_matrix<bigint> >(q,tr)) ,D,dxy,p2,v1,v,p_trans,p_trans2,par,0);
			
					p_trans = prod(tr, p_trans);
					p_trans2 = prod(tr, p_trans2);
				}
			else if (is_square(sq,p3))
				{
					#ifndef NDEBUG
					s << ">> p3 is a square!" << endl;
					#endif
			
					math_vector <bigint> v(2,2);
					v[0] = sq;
					v[1] = 1;
			
					math_vector <bigint> v2 = extract_square_factors(dxz,0,s);
			
					conic_param2(q,D,dxz,p3,v2,v,p_trans,p_trans2,par,0);
				}
			else
				{
					// No simple rational point found
			
					conic_param_no_ratpoint(q,D,p_trans,p_trans2,par,0,s,dxy,dxz,dyz,p1,p2,p3);
				}
		}

	#ifndef NDEBUG
	if (D.is_zero())
		s << ">> conic with rational point found" << endl;
	#endif
}

// The procedure for parameterizing conics when no rational point is known a
// priori, but there may be some.
// The param is (p_trans + sqrt(D)*p_trans2) * [u^2 v^2 uv], par = [u^2 v^2 uv] 
void conic_param(const bigint_matrix &q, bigint &D,	 
									bigint_matrix &p_trans, bigint_matrix &p_trans2,	 
									curve_param <bigint> &par, const int opt_level, ostream &s) 
{
	#ifndef NDEBUG
	s << ">> parameterization of conic" << endl;
	#endif

	math_vector <bigint> sol(3,3);

	// OPTIMIZATION	 :	a faire plus tard
	//if ((opt_level) && (FindRationalPointOnConic(q,sol,s)))
	//	{
	//		// Rational point found

	//		#ifndef NDEBUG
	//		s << ">> found rational point on conic: " << sol << endl;
	//		#endif
	//			
	//		// Parameterization of the conic
	//		math_vector <bigint> l(2,2);

	//		conic_param_through_ratpoint(q,sol,p_trans,par,l,s);
	//	}
	// else
		{
			if (opt_level) // There is no rational point
				conic_param_no_ratpoint(q,D,p_trans,p_trans2,par,opt_level,s);
			else // Let us try some simple heuristics to find a rational point
				conic_param_no_opt(q,D,p_trans,p_trans2,par,s);
		}
}

///////////// Parameterization of cones

// Compute the parameterization of a cone through a rational point
// Output is p_trans*[u*s^2 u*t^2 u*s*t v]
// s*l[0]+t*l[1] gives the line (sing,rat_point) on the cone
void cone_param_through_ratpoint(const bigint_matrix &q, 
																	const math_vector <bigint> &sing,
																	const math_vector <bigint> &rat_point, 
																	bigint_matrix &p_trans, surface_param <bigint> &p, 
																	math_vector <bigint> &l, ostream &s)
{
	#ifndef NDEBUG
	s << ">> parameterization of cone with rational point" << endl;
	#endif

	p_trans = send_to_zw(rat_point,sing);
	bigint_matrix q_r = prod( base_matrix<bigint> (prod(trans(p_trans), q)), p_trans);

	// Resize to 3x3 matrix
	q_r.resize(3,3);

	// Parameterize the conic
	curve_param <bigint> par(3);
	bigint_matrix p_trans_2d(3,3);
	conic_param_through_origin(q_r,p_trans_2d,par,l);

	// Stack the transformations
	p_trans_2d.resize(4,4);
	p_trans_2d.sto(3,3,1);
	multiply(p_trans,p_trans,p_trans_2d);

	// Complete curve_param to make a surface_param
	par.set_capacity(4);
	par.set_size(4);
	hom_polynomial <bigint> pol_s;
	pol_s.assign_x();
	multiply(p,par,pol_s);
	p[3].assign_y();
}

// Compute the parameterization of a cone
// Output is p_trans*[u*s^2 u*t^2 u*s*t v] + sqrt(D)*p_trans2*[u*s^2 u*t^2 u*s*t v]
void cone_param(const bigint_matrix &q, const math_vector <bigint> &sing, bigint &D,
								bigint_matrix &p_trans, bigint_matrix &p_trans2,
								surface_param <bigint> &p, const int opt_level, ostream &s)
{
	#ifndef NDEBUG
	s << ">> parameterization of cone" << endl;
	#endif

	// Rotate q
	bigint_matrix p_inf = send_to_infinity(sing);
	bigint_matrix q_r = prod(trans(p_inf), prod<base_matrix<bigint> >(q,p_inf));

	// Resize
	q_r.resize(3,3);

	curve_param <bigint> par(3);
	bigint_matrix p_trans2d(3,3),p_trans2d2(3,3);

	conic_param(q_r,D,p_trans2d,p_trans2d2,par,opt_level,s);

	#ifndef NDEBUG
	if (D.is_zero())
		s << ">> cone with rational point found" << endl;
	#endif

	// Stack the transformations
	p_trans2d.resize(4,4);
	p_trans2d.sto(3,3,1);
	p_trans = prod(p_inf, p_trans2d);

	if (!D.is_zero())
		{
			p_trans2d2.resize(4,4);
			p_trans2d2.sto(3,3,1);
			p_trans2 = prod(p_inf, p_trans2d2);
		}

	// Complete curve_param to make a surface_param
	par.set_capacity(4);
	par.set_size(4);

	hom_polynomial <bigint> pol_s;
	pol_s.assign_x();
	multiply(p,par,pol_s);

	p[3].assign_y();
}

///////////// Parameterization of pairs of planes, pairs of lines, and related

// Parameterize a 2 x 2 matrix when the det is 0
// Output is m1*[u]
void two_by_two_param_sing(const bigint_matrix &q, bigint_matrix &m1)
{
	// Thing is a x^2 + 2b xy + c y^2, b^2-ac = 0
	bigint a = q(0,0);
	bigint b = q(0,1);
	bigint c = q(1,1);

	if (!a.is_zero())
		{
			// a is not zero
			// m1 = (-b a)
			m1.sto(0,0,-b);
			m1.sto(1,0,a);
		}
	else
		{
			// c is not zero
			// m1 = (c -b)
			m1.sto(0,0,c);
			m1.sto(1,0,-b);
		}
}

// Compute the parameterization of 2 x 2 matrix, used for lines and planes
// Output is m1*[u] +/- sqrt(D)*m2*[u] if D is not a square
// Output is m1*[u] +/- m2*[u] if D is a square
// When b^2-ac is a square, D is set to 0
void two_by_two_param(const bigint_matrix &q, bigint &D, 
											bigint_matrix &m1, bigint_matrix &m2, 
											const int opt_level, ostream &s)
{
	// Thing is a x^2 + 2b xy + c y^2
	bigint a = q(0,0);
	bigint b = q(0,1);
	bigint c = q(1,1);

	bigint g = gcd(gcd(a,b),c);

	a = a/g;
	b = b/g;
	c = c/g;

	// D = b^2-ac
	bigint tmp;
	power(D,b,2);
	multiply(tmp,a,c);
	subtract(D,D,tmp);

	// Square root of D
	bigint sqrt_D;

	if (is_square(sqrt_D,D))
		{
			// b^2-ac is a square
			D.assign_zero();

			if (!a.is_zero())
				{
					// a is not zero
					// m1 = (-b a), m2 = (sqrt_D 0)
					m1.sto(0,0,-b);
					m1.sto(1,0,a);
					m2.sto(0,0,sqrt_D);
				}
			else if (!c.is_zero())
				{
					// c is not zero
					// m1 = (c -b), m2 = (0 sqrt_D)
					m1.sto(0,0,c);
					m1.sto(1,0,-b);
					m2.sto(1,0,sqrt_D);
				}
			else
				{
					// a = c = 0, b is not zero
					m1.sto(0,0,1);
					m1.sto(1,0,1);
					m2.sto(0,0,1);
					m2.sto(1,0,-1);
				}
		}
	else 
		{
			// D = b^2-ac is not a square
			math_vector <bigint> DD = extract_square_factors(D,opt_level,s);
			D = DD[1];

			if (!a.is_zero())
				{
					// a is not zero
					// m1 = (-b a), m2 = (1 0)
					m1.sto(0,0,-b);
					m1.sto(1,0,a);
					m2.sto(0,0,DD[0]);
				}
			else
				{
					// c is not zero
					// m1 = (c -b), m2 = (0 1)
					m1.sto(0,0,c);
					m1.sto(1,0,-b);
					m2.sto(1,0,DD[0]);
				}
		}
}

// Compute the parameterization of a pair of lines in the plane
// Output is m1*[u v] +/- sqrt(D)*m2*[u] if D is non zero
// Output is m1*[u v] +/- m2*[u] otherwise
void pair_of_lines_param(const bigint_matrix &q, const math_vector <bigint> &sing, 
													bigint &D, bigint_matrix &m1, bigint_matrix &m2, 
													const int opt_level, ostream &s)
{
	#ifndef NDEBUG
	s << ">> parameterization of pair of lines" << endl;
	#endif

	// Rotate the pair of lines by sending the singular point to infinity
	bigint_matrix p_trans = send_to_infinity(sing);
	bigint_matrix q_r;
	multiply(q_r,q,p_trans);
	multiply(q_r,math_matrix<bigint>(trans(p_trans)),q_r);

	// Resize to a 2D matrix
	q_r.resize(2,2);

	// Parameterize the reduced 2D matrix
	two_by_two_param(q_r,D,m1,m2,opt_level,s);

	// Augment the result and stack the transformation matrices
	// Add a 1 at the right place in m1 for the new variable
	m1.sto(2,1,1);

	multiply(m1,p_trans,m1);
	multiply(m2,p_trans,m2);
}

// Compute the parameterization of a double line
// Output is m1*[u u]
void double_line_param(const bigint_matrix &q, const bigint_matrix &sing,
												bigint_matrix &m1, ostream &s)
{
	#ifndef NDEBUG
	s << ">> parameterization of double line" << endl;
	#endif

	// Rotate the line by sending one of its singular points
	bigint_matrix p_trans = send_to_infinity(column(sing,0));
	bigint_matrix q_r;
	multiply(q_r,q,p_trans);
	multiply(q_r,math_matrix<bigint>(trans(p_trans)),q_r);

	// Resize to a 2D matrix
	q_r.resize(2,2);

	// Parameterize the reduced 2D matrix
	two_by_two_param_sing(q_r,m1);

	// Augment the result and stack the transformation matrices
	// Add a 1 at the right place in m1 for the new variable
	m1.sto(2,1,1);

	multiply(m1,p_trans,m1);
}

// Compute the parameterization of a double plane
// Output is m1*[s*u t*u v]
void double_plane_param(const bigint_matrix &q, const bigint_matrix &sing,
												bigint_matrix &m1, ostream &s)
{
	#ifndef NDEBUG
	s << ">> parameterization of double plane" << endl;
	#endif

	// Rotate the plane by sending one of its singular line to infinity
	bigint_matrix p_trans = send_to_zw(column(sing, 0), column(sing,1));
	bigint_matrix q_r;
	multiply(q_r,q,p_trans);
	multiply(q_r, math_matrix<bigint>(trans(p_trans)), q_r);

	// Resize to a 2D matrix
	q_r.resize(2,2);

	// Parameterize the reduced 2D matrix
	two_by_two_param_sing(q_r,m1);

	// Augment the result and stack the transformation matrices
	// Add a 1 at the right places in m1 for the two new variables of par
	m1.sto(2,1,1);
	m1.sto(3,2,1);

	multiply(m1,p_trans,m1);
}

// Compute the parameterization of a pair of planes
// Output is m1*[s*u t*u v] +/- sqrt(D)*m2*[s*u] if D is non zero
// Output is m1*[s*u t*u v] +/- m2*[s*u] otherwise
void pair_of_planes_param(const bigint_matrix &q, const bigint_matrix &sing, bigint &D, 
													bigint_matrix &m1, bigint_matrix &m2, 
													const int opt_level, ostream &s)
{
	#ifndef NDEBUG
	s << ">> parameterization of pair of planes" << endl;
	#endif

	// Rotate the pair of planes by sending its singular line to infinity
	bigint_matrix p_trans = send_to_zw(column(sing, 0), column(sing, 1));
	bigint_matrix q_r;
	multiply(q_r,q,p_trans);
	multiply(q_r,math_matrix<bigint>(trans(p_trans)),q_r);

	// Resize to a 2D matrix
	q_r.resize(2,2);

	// Parameterize the reduced 2D matrix
	two_by_two_param(q_r,D,m1,m2,opt_level,s);

	// Augment the result and stack the transformation matrices
	// Add a 1 at the right places in m1 for the two new variables of par
	m1.sto(2,1,1);
	m1.sto(3,2,1);

	multiply(m1,p_trans,m1);
	multiply(m2,p_trans,m2);
}

///////////// Parameterization of conics with non-rational coefficients

// Parameterization of the conic C1+sqrt(delta)*C2 with	 discriminant D = (eg)
// e^2-a*f > 0 
// The param is :
// (p_trans0 + sqrt(delta)*p_trans1 + sqrt(D)*p_trans2 +
//								sqrt(delta)*sqrt(D)*p_trans3) * [u^2 v^2 uv] 
// par = [u^2 v^2 uv]
// D = D[0] + sqrt(delta)*D[1]
void non_rational_conic_param1(const bigint_matrix &C1, const bigint_matrix &C2, 
																const bigint &delta, math_vector <bigint> &D, 
																const math_vector <bigint> &nd, 
																bigint_matrix &p_trans0, bigint_matrix &p_trans1, 
																bigint_matrix &p_trans2, bigint_matrix &p_trans3,
																curve_param <bigint> &par) 
{
	// Here we simply recopy the procedure conic_param1 after replacing conic q by
	// C1+sqrt(delta)*C2 
	// Conic is ax^2 + 2bxy + cy^2 + 2dyz + 2exz + fz^2
	// a = a[0] + a[1]*sqrt(delta) b=...
	// Inertia of the conic is supposed to be [2 1]
	math_vector <bigint> a(2,2), b(2,2), c(2,2), d(2,2), e(2,2), f(2,2);
	a[0] = C1(0,0);
	b[0] = C1(0,1); 
	c[0] = C1(1,1);
	d[0] = C1(1,2);
	e[0] = C1(0,2);
	f[0] = C1(2,2);

	a[1] = C2(0,0);
	b[1] = C2(0,1); 
	c[1] = C2(1,1);
	d[1] = C2(1,2);
	e[1] = C2(0,2);
	f[1] = C2(2,2);
	
	// Computed with Maple from the parameterization of conic_param1() the
	// parameterization matrix is	 
	// p(0,0) = a[0]*e[0] +delta*a[1]*e[1] + sqrtdelta*(a[0]*e[1]+a[1]*e[0])
	//			-a[0]*sqrtD	 -a[1]*sqrtdelta*sqrtD 
	// p(0,1) = +c[0]*e[0]+delta*c[1]*e[1] + sqrtdelta*(c[0]*e[1]+c[1]*e[0])
	//			+c[0]*sqrtD +c[1]*sqrtdelta*sqrtD 
	// p(0,2) = 2*a[0]*d[0]+2*a[1]*d[1]*delta + sqrtdelta*(2*a[0]*d[1]+2*a[1]*d[0])
	// p(1,0) = 0
	// p(1,1) = +2*a[0]*d[0]+2*a[1]*d[1]*delta-2*b[0]*e[0]-2*delta*b[1]*e[1] 
	//							+ sqrtdelta*(2*a[0]*d[1]+2*a[1]*d[0]-2*b[0]*e[1]-2*b[1]*e[0])
	//							-2*b[0]*sqrtD -2*b[1]*sqrtdelta*sqrtD
	// p(1,2) = -2*a[0]*sqrtD	 -2*a[1]*sqrtdelta*sqrtD
	// p(2,0) = -a[0]^2-a[1]^2*delta	-2*a[0]*a[1]*sqrtdelta
	// p(2,1) = -a[0]*c[0]-a[1]*c[1]*delta + sqrtdelta*(-a[0]*c[1]-a[1]*c[0])
	// p(2,2) = -2*a[0]*b[0]-2*a[1]*b[1]*delta + sqrtdelta*(-2*a[0]*b[1]-2*a[1]*b[0])

	// Lets store the values
	// p(0,0) = a[0]*e[0] +delta*a[1]*e[1] + sqrtdelta*(a[0]*e[1]+a[1]*e[0])
	// -a[0]*sqrtD	-a[1]*sqrtdelta*sqrtD 
	p_trans0.sto(0,0, a[0]*e[0] +delta*a[1]*e[1]);
	p_trans1.sto(0,0, a[0]*e[1]+a[1]*e[0]);
	p_trans2.sto(0,0, -a[0]);
	p_trans3.sto(0,0, -a[1]);
	// p(0,1) = +c[0]*e[0]+delta*c[1]*e[1] + sqrtdelta*(c[0]*e[1]+c[1]*e[0])
	// +c[0]*sqrtD +c[1]*sqrtdelta*sqrtD 
	p_trans0.sto(0,1, c[0]*e[0]+delta*c[1]*e[1]);
	p_trans1.sto(0,1, c[0]*e[1]+c[1]*e[0]);
	p_trans2.sto(0,1, c[0]);
	p_trans3.sto(0,1, c[1]);
	// p(0,2) = 2*a[0]*d[0]+2*a[1]*d[1]*delta + sqrtdelta*(2*a[0]*d[1]+2*a[1]*d[0])
	p_trans0.sto(0,2, 2*a[0]*d[0]+2*a[1]*d[1]*delta);
	p_trans1.sto(0,2, 2*a[0]*d[1]+2*a[1]*d[0]);
	p_trans2.sto(0,2, 0);
	p_trans3.sto(0,2, 0);
	// p(1,0) = 0
	p_trans0.sto(1,0, 0);
	p_trans1.sto(1,0, 0);
	p_trans2.sto(1,0, 0);
	p_trans3.sto(1,0, 0);
	// p(1,1) = +2*a[0]*d[0]+2*a[1]*d[1]*delta-2*b[0]*e[0]-2*delta*b[1]*e[1] 
	//							+ sqrtdelta*(2*a[0]*d[1]+2*a[1]*d[0]-2*b[0]*e[1]-2*b[1]*e[0])
	//							-2*b[0]*sqrtD -2*b[1]*sqrtdelta*sqrtD
	p_trans0.sto(1,1, 2*(a[0]*d[0]+a[1]*d[1]*delta-b[0]*e[0]-delta*b[1]*e[1]));
	p_trans1.sto(1,1, 2*(a[0]*d[1]+a[1]*d[0]-b[0]*e[1]-b[1]*e[0]));
	p_trans2.sto(1,1, -2*b[0]);
	p_trans3.sto(1,1, -2*b[1]);
	// p(1,2) = -2*a[0]*sqrtD	 -2*a[1]*sqrtdelta*sqrtD
	p_trans0.sto(1,2, 0);
	p_trans1.sto(1,2, 0);
	p_trans2.sto(1,2, -2*a[0]);
	p_trans3.sto(1,2, -2*a[1]);
	// p(2,0) = -a[0]^2-a[1]^2*delta	-2*a[0]*a[1]*sqrtdelta
	p_trans0.sto(2,0, -a[0]*a[0]-a[1]*a[1]*delta);
	p_trans1.sto(2,0, -2*a[0]*a[1]);
	p_trans2.sto(2,0, 0);
	p_trans3.sto(2,0, 0);
	// p(2,1) = -a[0]*c[0]-a[1]*c[1]*delta + sqrtdelta*(-a[0]*c[1]-a[1]*c[0])
	p_trans0.sto(2,1, -a[0]*c[0]-a[1]*c[1]*delta);
	p_trans1.sto(2,1, -a[0]*c[1]-a[1]*c[0]);
	p_trans2.sto(2,1, 0);
	p_trans3.sto(2,1, 0);
	// p(2,2) = -2*a[0]*b[0]-2*a[1]*b[1]*delta + sqrtdelta*(-2*a[0]*b[1]-2*a[1]*b[0])
	p_trans0.sto(2,2, -2*a[0]*b[0]-2*a[1]*b[1]*delta);
	p_trans1.sto(2,2, -2*a[0]*b[1]-2*a[1]*b[0]);
	p_trans2.sto(2,2, 0);
	p_trans3.sto(2,2, 0);

	// The param
	par[0].assign_x2();
	par[1].assign_y2();
	par[2].assign_xy();

	// nd is = to eg v1 such that the subdiscriminant dxy =
	// nd[2]^2*(nd[0]+sqrt(delta)*nd[1]) 
	D[0] = nd[0];
	D[1] = nd[1];
	if (nd[2] != 1) 
		{
			// Recall that unlike in the rational case,	 here the square factor is
			// stored in vi[2] (instead vi[0]) 
			p_trans2 = p_trans2*nd[2];
			p_trans3 = p_trans3*nd[2];
		}
}

// Parameterization by -f*p > 0
// Equation is (a*f-e^2)*(f*z+e*x+d*y)^2+(x*(a*f-e^2)+y*(b*f-d*e))^2+f*p*y^2 = 0
// The param is 
// (p_trans0 + sqrt(delta)*p_trans1 + sqrt(D)*p_trans2 +
//									sqrt(delta)*sqrt(D)*p_trans3) * [u^2 v^2 uv] 
// par = [u^2 v^2 uv]
// D = D[0] + sqrt(delta)*D[1]
void non_rational_conic_param2(const bigint_matrix &C1, const bigint_matrix &C2, 
																const bigint &delta, math_vector <bigint> &D, 
																const math_vector <bigint> &de,
																const math_vector <bigint> &p3, 
																const math_vector <bigint> &nd, 
																const math_vector <bigint> &np, 
																bigint_matrix &p_trans0, bigint_matrix &p_trans1, 
																bigint_matrix &p_trans2, bigint_matrix &p_trans3,
																curve_param <bigint> &par, const int opt_level) 
{
	// Here we simply recopy the procedure conic_param2() after replacing conic q by
	// C1+sqrt(delta)*C2 
	// Conic is ax^2 + 2bxy + cy^2 + 2dyz + 2exz + fz^2
	// a = a[0] + a[1]*sqrt(delta) b=...
	// Inertia of the conic is supposed to be [2 1]
	math_vector <bigint> a(2,2), b(2,2), c(2,2), d(2,2), e(2,2), f(2,2), t1(2,2), t2(2,2);
	a[0] = C1(0,0);
	b[0] = C1(0,1); 
	c[0] = C1(1,1);
	d[0] = C1(1,2);
	e[0] = C1(0,2);
	f[0] = C1(2,2);

	a[1] = C2(0,0);
	b[1] = C2(0,1); 
	c[1] = C2(1,1);
	d[1] = C2(1,2);
	e[1] = C2(0,2);
	f[1] = C2(2,2);
	
	// Computed from the parameterization of conic_param2() 
	// the	parameterization matrix is as follows
	// Let	 t1 = b*f-d*e;	 t2 = a*d-b*e;
	// t1 = (b[0]+sqrt(delta)*b[1])*(f[0]+sqrt(delta)*f[1])
	//					-(d[0]+sqrt(delta)*d[1])*(e[0]+sqrt(delta)*e[1])
	// t2 = (a[0]+sqrt(delta)*a[1])*(d[0]+sqrt(delta)*d[1])
	//					-(b[0]+sqrt(delta)*b[1])*(e[0]+sqrt(delta)*e[1])
	t1[0] = b[0]*f[0]-d[0]*e[0] + (b[1]*f[1]-d[1]*e[1])*delta;
	t1[1] = b[0]*f[1]+b[1]*f[0]-d[0]*e[1]-d[1]*e[0];
	t2[0] = a[0]*d[0]-b[0]*e[0] + (a[1]*d[1]-b[1]*e[1])*delta;
	t2[1] = a[0]*d[1]+a[1]*d[0]-b[0]*e[1]-b[1]*e[0];

	// p = (p_trans0 + sqrt(delta)*p_trans1 + sqrt(D)*p_trans2 +
	// sqrt(delta)*sqrt(D)*p_trans3) * [u^2 v^2 uv] 

	// e^2-a*f < 0
	// Commented are the values from conic_param2() 
	if (sign(nd[0], nd[1], delta) < 0)
		{ 
			// p_trans.sto(0,0,f*t1);		
			// f*t1 = (f[0]+sqrt(delta)*f[1])*(t[0]+sqrt(delta)*t[1]) 
			p_trans0.sto(0,0, f[0]*t1[0] +delta*f[1]*t1[1]);
			p_trans1.sto(0,0, f[0]*t1[1] + f[1]*t1[0]);
			// p_trans.sto(0,1,-p_trans(0,0)*nd[1]); 
			// =-(p_trans0(0,0)+sqrt(delta)*p_trans1(0,0))*(nd[0]+sqrt(delta)*nd[1])
			p_trans0.sto(0,1, -p_trans0(0,0)*nd[0] - delta*p_trans1(0,0)*nd[1]);
			p_trans1.sto(0,1, -p_trans1(0,0)*nd[0] - p_trans0(0,0)*nd[1]);
			// p_trans.sto(1,0,de*f);
			p_trans0.sto(1,0, f[0]*de[0] +delta*f[1]*de[1]);
			p_trans1.sto(1,0, f[0]*de[1] + f[1]*de[0]);
			// p_trans.sto(1,1,-p_trans(1,0)*nd[1]);
			p_trans0.sto(1,1, -p_trans0(1,0)*nd[0] - delta*p_trans1(1,0)*nd[1]);
			p_trans1.sto(1,1, -p_trans1(1,0)*nd[0] - p_trans0(1,0)*nd[1]);
			// p_trans.sto(2,0,f*t2);
			p_trans0.sto(2,0, f[0]*t2[0] +delta*f[1]*t2[1]);
			p_trans1.sto(2,0, f[0]*t2[1] + f[1]*t2[0]);
			// p_trans.sto(2,1,-p_trans(2,0)*nd[1]);
			p_trans0.sto(2,1, -p_trans0(2,0)*nd[0] - delta*p_trans1(2,0)*nd[1]);
			p_trans1.sto(2,1, -p_trans1(2,0)*nd[0] - p_trans0(2,0)*nd[1]);
			
			// p_trans2.sto(0,0,-f);
			p_trans2.sto(0,0,-f[0]);
			p_trans3.sto(0,0,-f[1]);
			// p_trans2.sto(0,1,p_trans2(0,0)*nd[1]);
			p_trans2.sto(0,1,p_trans2(0,0)*nd[0] + delta*p_trans3(0,0)*nd[1]);
			p_trans3.sto(0,1,p_trans3(0,0)*nd[0] + p_trans2(0,0)*nd[1]);
			// p_trans2.sto(2,0,e);
			p_trans2.sto(2,0,e[0]);
			p_trans3.sto(2,0,e[1]);
			// p_trans2.sto(2,1,p_trans2(2,0)*nd[1]);
			p_trans2.sto(2,1,p_trans2(2,0)*nd[0] + delta*p_trans3(2,0)*nd[1]);
			p_trans3.sto(2,1,p_trans3(2,0)*nd[0] + p_trans2(2,0)*nd[1]);
			// p_trans2.sto(2,2,2*nd[0]*nd[1]);	 = 2*nd[2]*(nd[0] + sqrt(delta)*nd[1])
			// Recall that unlike in the rational case,	 here the square factor is
			// stored in vi[2] (instead vi[0]) 
			p_trans2.sto(2,2, 2*nd[2]*nd[0]);
			p_trans3.sto(2,2, 2*nd[2]*nd[1]);
		}
	else
		{
			// math_vector <bigint> vf = extract_square_factors(de*np[1],1,cout);
			// de*np[1] is here = (de[0]+sqrt(delta)*de[1]) * (np[0]+sqrt(delta)*np[1])
			// = (de[0]*np[0] + de[1]*np[1]*delta) +sqrt(delta)*(de[0]*np[1]+de[1]*np[0])
			math_vector <bigint> vf(3,3);
			vf[0] = de[0]*np[0] + de[1]*np[1]*delta;
			vf[1] = de[0]*np[1] + de[1]*np[0];
			vf[2] = 1;

			bigint tmp = gcd(vf[0], vf[1]);
			tmp = extract_square_factors(tmp,opt_level,cout)[0];
			if (tmp != 1)
				{
					bigint tmp2 = tmp*tmp;
					vf[0] = vf[0]/tmp2;
					vf[1] = vf[1]/tmp2;
					vf[2] = tmp;
				}

			// p_trans.sto(0,0,f*p3);
			p_trans0.sto(0,0, f[0]*p3[0] +delta*f[1]*p3[1]);
			p_trans1.sto(0,0, f[0]*p3[1] + f[1]*p3[0]);
			// p_trans.sto(0,1,p_trans(0,0)*vf[1]);
			p_trans0.sto(0,1, p_trans0(0,0)*vf[0] + delta*p_trans1(0,0)*vf[1]);
			p_trans1.sto(0,1, p_trans1(0,0)*vf[0] + p_trans0(0,0)*vf[1]);
			// p_trans.sto(2,0,-e*p3);
			p_trans0.sto(2,0, -e[0]*p3[0] -delta*e[1]*p3[1]);
			p_trans1.sto(2,0, -e[0]*p3[1] - e[1]*p3[0]); 
			// p_trans.sto(2,1,p_trans(2,0)*vf[1]);
			p_trans0.sto(2,1, p_trans0(2,0)*vf[0] + delta*p_trans1(2,0)*vf[1]);
			p_trans1.sto(2,1, p_trans1(2,0)*vf[0] + p_trans0(2,0)*vf[1]);

			// p_trans2.sto(0,0,f*t1);
			p_trans2.sto(0,0, f[0]*t1[0] + delta*f[1]*t1[1]);
			p_trans3.sto(0,0, f[0]*t1[1] + f[1]*t1[0]);
			// p_trans2.sto(0,1,-p_trans2(0,0)*vf[1]);
			p_trans2.sto(0,1, -p_trans2(0,0)*vf[0] - delta*p_trans3(0,0)*vf[1]);
			p_trans3.sto(0,1, -p_trans3(0,0)*vf[0] - p_trans2(0,0)*vf[1]);
			// p_trans2.sto(1,0,f*de);
			p_trans2.sto(1,0, f[0]*de[0] +delta*f[1]*de[1]);
			p_trans3.sto(1,0, f[0]*de[1] + f[1]*de[0]);
			// p_trans2.sto(1,1,-p_trans2(1,0)*vf[1]);
			p_trans2.sto(1,1, -p_trans2(1,0)*vf[0] - delta*p_trans3(1,0)*vf[1]);
			p_trans3.sto(1,1, -p_trans3(1,0)*vf[0] - p_trans2(1,0)*vf[1]);
			// p_trans2.sto(2,0,f*t2);
			p_trans2.sto(2,0, f[0]*t2[0] +delta*f[1]*t2[1]);
			p_trans3.sto(2,0, f[0]*t2[1] + f[1]*t2[0]);
			// p_trans2.sto(2,1,-p_trans2(2,0)*vf[1]);
			p_trans2.sto(2,1, -p_trans2(2,0)*vf[0] - delta*p_trans3(2,0)*vf[1]);
			p_trans3.sto(2,1, -p_trans3(2,0)*vf[0] - p_trans2(2,0)*vf[1]);
			// p_trans2.sto(2,2,-2*np[0]*vf[0]*vf[1]); =
			// -2*np[2]*vf[2]*(vf[0]+sqrt(delta)*vf[1]) 
			// Recall that unlike in the rational case,	 here the square factor is
			// stored in vi[2] (instead vi[0]) 
			p_trans2.sto(2,2, -2*np[2]*vf[2]*vf[0]);
			p_trans3.sto(2,2, -2*np[2]*vf[2]*vf[1]);
		}

	// The param
	par[0].assign_x2();
	par[1].assign_y2();
	par[2].assign_xy();

	// nd is = to eg v1 such that the subdiscriminant dxy =
	// nd[2]^2*(nd[0]+sqrt(delta)*nd[1]) 
	D[0] = np[0];
	D[1] = np[1];

	if (np[2] != 1) 
		{
			p_trans2 = p_trans2*np[2];
			p_trans3 = p_trans3*np[2];
		}
}

// The procedure for parameterizing a conics whose corresponding 3x3 matrix is 
// C1 +sqrt(delta).C2 when no rational point is known
// The param is
//					 (p_trans0 + sqrt(delta)*p_trans1 + sqrt(D)*p_trans2 
//											+ sqrt(delta)*sqrt(D)*p_trans3) * [u^2 v^2 uv]
// par = [u^2 v^2 uv]
// D = D[0] + sqrt(delta)*D[1]
void non_rational_conic_param(const bigint_matrix &C1, const bigint_matrix &C2, 
															const bigint &delta, math_vector <bigint> &D, 
															bigint_matrix &p_trans0, bigint_matrix &p_trans1, 
															bigint_matrix &p_trans2, bigint_matrix &p_trans3,
															curve_param <bigint> &par, 
															const int opt_level, ostream &s)
{
	#ifndef NDEBUG
	s << ">> parameterization of conic with non-rational coefficients" << endl;
	#endif

	// The conic has no rational point because they would necessarily be on the 
	// rational singular line of the pair of planes which is impossible since the
	// intersection curve on the other plane is entirely imaginary

	// Now, if one of the following 2x2 discriminants is a square, we also have
	// a rational param. 
	// dxy, dxz, dyz are minus the 2x2 discriminants 
	// (ie -det(minor(C1+sqrt(delta)*C2, 2,2)), -det(minor(C1+sqrt(delta)*C2, 1,1)),
	// -det(minor(C1+sqrt(delta)*C2, 0,0))) dxy = dxy[0] +sqrt(delta).dxy[1]
	math_vector <bigint> dxy(2,2),dxz(2,2),dyz(2,2),detq(2,2),p0(2,2),p1(2,2),p2(2,2);

	dxy[0] = (-C2.member(0,0)*C2.member(1,1)+C2.member(0,1)*C2.member(0,1))*delta 
		-C1.member(0,0)*C1.member(1,1)+C1.member(0,1)*C1.member(0,1);
	dxy[1] = -C1.member(0,0)*C2.member(1,1)-C2.member(0,0)*C1.member(1,1)
		+2*C1.member(0,1)*C2.member(0,1);

	dxz[0] = (-C2.member(0,0)*C2.member(2,2)+C2.member(0,2)*C2.member(0,2))*delta 
		-C1.member(0,0)*C1.member(2,2)+C1.member(0,2)*C1.member(0,2);
	dxz[1] = -C1.member(0,0)*C2.member(2,2)-C2.member(0,0)*C1.member(2,2)
		+2*C1.member(0,2)*C2.member(0,2);

	dyz[0] = (-C2.member(1,1)*C2.member(2,2)+C2.member(1,2)*C2.member(1,2))*delta 
		-C1.member(1,1)*C1.member(2,2)+C1.member(1,2)*C1.member(1,2);
	dyz[1] = -C1.member(1,1)*C2.member(2,2)-C2.member(1,1)*C1.member(2,2)
		+2*C1.member(1,2)*C2.member(1,2);

	detq[0] = (-C2.member(0,1)*C2.member(0,1)*C1.member(2,2)
				-2*C2.member(0,0)*C1.member(1,2)*C2.member(1,2)
				-2*C1.member(0,2)*C2.member(0,2)*C2.member(1,1)
				-C2.member(0,2)*C2.member(0,2)*C1.member(1,1)
				+C2.member(0,0)*C1.member(1,1)*C2.member(2,2)
				+2*C1.member(0,1)*C2.member(0,2)*C2.member(1,2)
				-C1.member(0,0)*C2.member(1,2)*C2.member(1,2)
				+2*C2.member(0,1)*C2.member(0,2)*C1.member(1,2)
				+C2.member(0,0)*C2.member(1,1)*C1.member(2,2)
				-2*C1.member(0,1)*C2.member(0,1)*C2.member(2,2)
				+2*C2.member(0,1)*C1.member(0,2)*C2.member(1,2)
				+C1.member(0,0)*C2.member(1,1)*C2.member(2,2))*delta
		-C1.member(0,2)*C1.member(0,2)*C1.member(1,1)
		-C1.member(0,1)*C1.member(0,1)*C1.member(2,2)
		-C1.member(0,0)*C1.member(1,2)*C1.member(1,2)
		+C1.member(0,0)*C1.member(1,1)*C1.member(2,2)
		+2*C1.member(0,1)*C1.member(0,2)*C1.member(1,2);
	detq[1] = (-C2.member(0,1)*C2.member(0,1)*C2.member(2,2)
				+C2.member(0,0)*C2.member(1,1)*C2.member(2,2)
				-C2.member(0,0)*C2.member(1,2)*C2.member(1,2)
				-C2.member(0,2)*C2.member(0,2)*C2.member(1,1)
				+2*C2.member(0,1)*C2.member(0,2)*C2.member(1,2))*delta
		+(-C2.member(0,0)*C1.member(1,2)*C1.member(1,2)
			+C2.member(0,0)*C1.member(1,1)*C1.member(2,2)
			+C1.member(0,0)*C1.member(1,1)*C2.member(2,2)
			+2*C1.member(0,1)*C2.member(0,2)*C1.member(1,2)
			+2*C2.member(0,1)*C1.member(0,2)*C1.member(1,2)
			+C1.member(0,0)*C2.member(1,1)*C1.member(2,2)
			+2*C1.member(0,1)*C1.member(0,2)*C2.member(1,2)
			-2*C1.member(0,0)*C1.member(1,2)*C2.member(1,2)
			-C1.member(0,1)*C1.member(0,1)*C2.member(2,2)
			-2*C1.member(0,2)*C2.member(0,2)*C1.member(1,1)
			-2*C1.member(0,1)*C2.member(0,1)*C1.member(2,2)
			-C1.member(0,2)*C1.member(0,2)*C2.member(1,1));
			
	// pi = -q(0,0)* det(q) = (-C1(i,i) - sqrt(delta). C2(i,i)) * (detq[0] +
	// sqrt(delta). detq[1]) = -C1(i,i)*detq[0] -delta*C2(i,i)*detq[1]
	// +sqrt(delta)*[-C2(i,i)*detq[0] -C1(i,i)*detq[1]] 
	p0[0] = -C1.member(0,0)*detq[0] -C2.member(0,0)*detq[1]*delta;
	p0[1] = -C2.member(0,0)*detq[0] -C1.member(0,0)*detq[1];

	p1[0] = -C1.member(1,1)*detq[0] -C2.member(1,1)*detq[1]*delta;
	p1[1] = -C2.member(1,1)*detq[0] -C1.member(1,1)*detq[1];

	p2[0] = -C1.member(2,2)*detq[0] -C2.member(2,2)*detq[1]*delta;
	p2[1] = -C2.member(2,2)*detq[0] -C1.member(2,2)*detq[1];

	// if ((dxy > 0) && (is_square(sq,dxy)))
	// else if (is_square(sq,dxz))
	// else if (is_square(sq,dyz))
	// else if (is_square(sq,p1))
	// else if (is_square(sq,p2))
	// else if (is_square(sq,p3))
	// else
	{
		// No rational param found. If opt_level was passed, try to find the
		// smallest square root. If not, take the smallest value among
		// dxy, dxz, dyz, p1, p2, p3...

		math_vector <bigint> v1(3,3),v2(3,3),v3(3,3),v4(3,3),v5(3,3),v6(3,3);
		// Unlike in the rational case,	 here the square factor is stored in vi[2]
		// (instead vi[0]) 
		v1[0] = dxy[0];
		v2[0] = dxz[0];
		v3[0] = dyz[0];
		v4[0] = p0[0];
		v5[0] = p1[0];
		v6[0] = p2[0];
		v1[1] = dxy[1];
		v2[1] = dxz[1];
		v3[1] = dyz[1];
		v4[1] = p0[1];
		v5[1] = p1[1];
		v6[1] = p2[1];
		v1[2] = 1; v2[2] = 1; v3[2] = 1; v4[2] = 1; v5[2] = 1; v6[2] = 1;

		// Optimize	 dxy, dxz, dyz, p1, p2, p3
		bigint tmp, tmp2;
		if (opt_level)
			{
				if (sign(v1[0],v1[1],delta) > 0)
		 			{
			 			tmp = extract_square_factors(gcd(v1[0],v1[1]),opt_level,s)[0];
			 			if (tmp != 1)
							{
								tmp2 = tmp*tmp;
								v1[0] = v1[0]/tmp2;
								v1[1] = v1[1]/tmp2;
								v1[2] = tmp;
							}
					}
				if (sign(v2[0],v2[1],delta) > 0)
					{
			 			tmp = extract_square_factors(gcd(v2[0],v2[1]),opt_level,s)[0];
			 			if (tmp != 1)
				 			{
								tmp2 = tmp*tmp;
								v2[0] = v2[0]/tmp2;
								v2[1] = v2[1]/tmp2;
								v2[2] = tmp;
				 			}
		 			}
				if (sign(v3[0],v3[1],delta) > 0)
					{
						tmp = extract_square_factors(gcd(v3[0],v3[1]),opt_level,s)[0];
						if (tmp != 1)
							{
								tmp2 = tmp*tmp;
								v3[0] = v3[0]/tmp2;
								v3[1] = v3[1]/tmp2;
								v3[2] = tmp;
							}
					}
				if (sign(v4[0],v4[1],delta) > 0)
					{
						tmp = extract_square_factors(gcd(v4[0],v4[1]),opt_level,s)[0];
						if (tmp != 1)
						 	{
								tmp2 = tmp*tmp;
								v4[0] = v4[0]/tmp2;
								v4[1] = v4[1]/tmp2;
								v4[2] = tmp;
							}
					}
				if (sign(v5[0],v5[1],delta) > 0)
					{
						tmp = extract_square_factors(gcd(v5[0],v5[1]),opt_level,s)[0];
						if (tmp != 1)
							{
								tmp2 = tmp*tmp;
								v5[0] = v5[0]/tmp2;
								v5[1] = v5[1]/tmp2;
								v5[2] = tmp;
							}
					}
				if (sign(v6[0],v6[1],delta) > 0)
					{
						tmp = extract_square_factors(gcd(v6[0],v6[1]),opt_level,s)[0];
						if (tmp != 1)
							{
								tmp2 = tmp*tmp;
								v6[0] = v6[0]/tmp2;
								v6[1] = v6[1]/tmp2;
								v6[2] = tmp;
							}
					}
			}

		// Find the "smallest square root" among dxy, dxz, dyz, p1, p2, p3
		int id;
		math_vector <bigint> v(3,3);
		// vi >0 iff vi!=v_neg
		// If dxy < 0, then p1 > 0 and conversely (initialization)
		if (sign(v1[0],v1[1],delta) < 0)
			{
				id = 4;
				v = v4;
			}
		else
			{
				id = 1;
				v = v1;
			}
		// Find the "smallest square root" among dxy, dxz, dyz, p1, p2, p3
		// We take vj=vj[2]^2*(vj[0] + sqrt(delta)*vj[1]) instead of v
		// iff v >0 and [ (v[1]!=0 and (vj[1]=0 or ||vj||<||v||)) or (v[1]=0 and
		// vj[1]=0 and ||vj||<||v||)], we take ||v|| = abs(v[0]) + abs(v[1])
		if ((sign(v2[0], v2[1], delta) > 0) && (((v[1]!=0) && ((v2[1]==0) || 
			((abs(v2[0]) + abs(v2[1])) < (abs(v[0]) + abs(v[1]))))) || 
			((v[1]==0) && (v2[1]==0) &&	((abs(v2[0]) + abs(v2[1])) < (abs(v[0]) + abs(v[1]))))))
			{
				v = v2;
				id = 2;
			}
		if ((sign(v3[0], v3[1], delta) > 0) && (((v[1]!=0) && ((v3[1]==0) || 
			((abs(v3[0]) + abs(v3[1])) < (abs(v[0]) + abs(v[1]))))) || 
			((v[1]==0) && (v3[1]==0) &&	((abs(v3[0]) + abs(v3[1])) < (abs(v[0]) + abs(v[1]))))))
			{
				v = v3;
				id = 3;
			}
		if ((id != 4) && (sign(v4[0], v4[1], delta) > 0) && (((v[1]!=0) && ((v4[1]==0) || 
			((abs(v4[0]) + abs(v4[1])) < (abs(v[0]) + abs(v[1]))))) || ((v[1]==0) && (v4[1]==0) &&	
			((abs(v4[0]) + abs(v4[1])) < (abs(v[0]) + abs(v[1]))))))
			{
				v = v4;
				id = 4;
			}
		if ((sign(v5[0], v5[1], delta) > 0) && (((v[1]!=0) && ((v5[1]==0) || 
			((abs(v5[0]) + abs(v5[1])) < (abs(v[0]) + abs(v[1]))))) || ((v[1]==0) && (v5[1]==0) &&	
				((abs(v5[0]) + abs(v5[1])) < (abs(v[0]) + abs(v[1]))))))
			{
				v = v5;
				id = 5;
			}
		if ((sign(v6[0], v6[1], delta) > 0) && (((v[1]!=0) && ((v6[1]==0) || 
			((abs(v6[0]) + abs(v6[1])) < (abs(v[0]) + abs(v[1]))))) || ((v[1]==0) && (v6[1]==0) &&	
				((abs(v6[0]) + abs(v6[1])) < (abs(v[0]) + abs(v[1]))))))
			id = 6;

		if (id == 1)
			{
				// Switch y and z
				bigint_matrix tr(3,3);
				tr.sto(0,0,1);
				tr.sto(1,2,1);
				tr.sto(2,1,1);
			
				non_rational_conic_param1(prod(tr,prod <base_matrix <bigint> >(C1,tr)), 
																	prod(tr,prod <base_matrix<bigint> >(C2,tr)), 
																	delta, D, v1, p_trans0, 
																	p_trans1,	 p_trans2,	p_trans3, par);
			
				p_trans0 = prod(tr,p_trans0);
				p_trans1 = prod(tr,p_trans1);
				p_trans2 = prod(tr,p_trans2);
				p_trans3 = prod(tr,p_trans3);
			}
		else if (id == 2)
			non_rational_conic_param1(C1, C2, delta, D, v2, p_trans0, p_trans1,	 
																p_trans2, p_trans3, par);
		else if (id == 3)
			{
				// Switch x and y
				bigint_matrix tr(3,3);
				tr.sto(0,1,1);
				tr.sto(1,0,1);
				tr.sto(2,2,1);
			
				non_rational_conic_param1(prod(tr,prod <base_matrix<bigint> >(C1, tr)), 
																	prod(tr,prod <base_matrix<bigint> >(C2, tr)), 
																	delta, D, v3, p_trans0, 
																	p_trans1,	p_trans2,	p_trans3, par);
			
				p_trans0 = prod(tr,p_trans0);
				p_trans1 = prod(tr,p_trans1);
				p_trans2 = prod(tr,p_trans2);
				p_trans3 = prod(tr,p_trans3);
			}
		else if (id == 4)
			{
				// Switch x and z
				bigint_matrix tr(3,3);
				tr.sto(0,2,1);
				tr.sto(1,1,1);
				tr.sto(2,0,1);
			
				non_rational_conic_param2(prod(tr,prod <base_matrix<bigint> >(C1, tr)),
																	prod(tr,prod <base_matrix<bigint> >(C2, tr)), 
																	delta, D,dxz,p0,v2,v4,
																	p_trans0, p_trans1,	 p_trans2,	p_trans3, par,opt_level);
			
				p_trans0 = prod(tr,p_trans0);
				p_trans1 = prod(tr,p_trans1);
				p_trans2 = prod(tr,p_trans2);
				p_trans3 = prod(tr,p_trans3);
			}
		else if (id == 5)
			{
				// Switch y and z
				bigint_matrix tr(3,3);
				tr.sto(0,0,1);
				tr.sto(1,2,1);
				tr.sto(2,1,1);
			
				non_rational_conic_param2(prod(tr,prod <base_matrix<bigint> >(C1, tr)),
																	prod(tr,prod <base_matrix<bigint> >(C2, tr)), 
																	delta, D,dxy,p1,v1,v5,
																	p_trans0, p_trans1,	p_trans2,	p_trans3, par,opt_level);
							 
				p_trans0 = prod(tr,p_trans0);
				p_trans1 = prod(tr,p_trans1);
				p_trans2 = prod(tr,p_trans2);
				p_trans3 = prod(tr,p_trans3);
			}
		else // id = 6
			non_rational_conic_param2(C1, C2, delta, D, dxz, p2, v2, v6, p_trans0, 
																p_trans1,	p_trans2,	p_trans3, par,opt_level);
	}

	#ifndef NDEBUG
	if ((D[0].is_zero()) && (D[1].is_zero()))
		s << ">> rational conic found" << endl;
	#endif
}

///////////// Optimization of line parameterizations

// Try to reparameterize a rational line by taking as "endpoints" (ie. the points
// (u,v) = (1,0) and (0,1)) points on the planes x = 0, y = 0, z = 0, w = 0. 
// In other words, one of the coordinates is 0 when u=0 and another is zero when
// v=0!!!! This assumes that the current endpoints of line_par have already been
// optimized. Otherwise p0 and end1 (say) might represent the same point but with
// different scales !!!
// Returns the value of the new parameters corresponding to the first endpoint of
// the old parameterization
math_vector <bigint> improved_line_param(curve_param <bigint> &line_par)
{
	math_vector <bigint> p[4] = {math_vector <bigint>(4), math_vector <bigint>(4), math_vector <bigint>(4), math_vector <bigint>(4)};
	bigint c[4];

	for (unsigned int i = 0; i < 4; i++)
		if (!line_par[i].is_zero())
			{
				p[i] = line_par.eval(line_par[i][0],-line_par[i][1]);
				c[i] = content(p[i]);
			
				if (c[i] != 1)
					divide(p[i],p[i],c[i]);
			}
			
	// Initialize with endpoints of line_par
	math_vector <bigint> end1 = line_par.eval(1,0);
	math_vector <bigint> end2 = line_par.eval(0,1);
	bigint sum1 = sum_of_squares(end1);
	bigint sum2 = sum_of_squares(end2);
	bigint u1 = 1, u2 = 0, v1 = 0, v2 = 1, co1 = 1, co2 = 1;

	if (sum2 < sum1)
		{
			swap(end1,end2);
			swap(sum1,sum2);
			swap(u1,u2);
			swap(v1,v2);
		}

	bigint sum_tmp;

	for (unsigned int i = 0; i < 4; i++)
		{
			sum_tmp = sum_of_squares(p[i]);

			if ((0 < sum_tmp) && (sum_tmp < sum1))
				{
					sum2 = sum1;
					end2 = end1;
					sum1 = sum_tmp;
					end1 = p[i];
					u2 = u1;
					u1 = line_par[i][0];
					v2 = v1;
					v1 = -line_par[i][1];
					co2 = co1;
					co1 = c[i];
				}
			else if ((sum1 < sum_tmp) && (sum_tmp < sum2))
				{
					sum2 = sum_tmp;
					end2 = p[i];
					u2 = line_par[i][0];
					v2 = -line_par[i][1];
					co2 = c[i];
				}
		}

	curve_param <bigint> line_new(end1,end2);

	line_par = line_new;

	math_vector <bigint> end_old_param(4,4);
	end_old_param[0] = co1*v2;
	end_old_param[1] = -co2*v1;
	end_old_param[2] = co1*u2;
	end_old_param[3] = -co2*u1;

	return end_old_param;
}

// Try to reparameterize a line c = c1+sqrt(xi). c2 such that 
// one of the coordinates is 0 when u=0 and another is zero when v=0
// Returns the parameters corresponding to the initial endpoints
math_vector <bigint> improved_line_param(curve_param <bigint> &c1, curve_param <bigint> &c2, 
																					const bigint &xi)
{
	// c[i] = c1[i] + sqrt(xi). c2[i]
	//			 = c1[i][1].u + c1[i][0].v + sqrt(xi). c2[i][1].u + sqrt(xi). c2[i][0].v
	//			 = (c1[i][1] + sqrt(xi). c2[i][1]). u + (c1[i][0] + sqrt(xi). c2[i][0]). v
	// We consider	for i from 0 to 3 the point Mi = c(u,v) where 
	// u = (c1[i][0] + sqrt(xi). c2[i][0]) and	v = -(c1[i][1] + sqrt(xi). c2[i][1])
	// The jth coordinate of this point is 
	//		-(c1[j][1] + sqrt(xi). c2[j][1]).(c1[i][0] + sqrt(xi). c2[i][0])	
	//		+ (c1[j][0] + sqrt(xi). c2[j][0]).(c1[i][1] + sqrt(xi). c2[i][1])
	// = -c1[j][1]*c1[i][0] - xi*c2[j][1]*c2[i][0] + c1[j][0]*c1[i][1] +
	// xi*c2[j][0]*c2[i][1] +sqrt(xi)*(-c1[j][1]*c2[i][0] - c2[j][1]*c1[i][0] +
	// c1[j][0]*c2[i][1] + c2[j][0]*c1[i][1]) 
	// which is zero for j=i. 
	// Let pi and qi be the points such that Mi = pi + sqrt(xi). qi. Their
	// coordinates are	pi[j] = -c1[j][1]*c1[i][0] - xi*c2[j][1]*c2[i][0] +
	// c1[j][0]*c1[i][1] + xi*c2[j][0]*c2[i][1] qi[j] = -c1[j][1]*c2[i][0] -
	// c2[j][1]*c1[i][0] + c1[j][0]*c2[i][1] + c2[j][0]*c1[i][1] 
	// pi = c1.eval(-c1[i][0], c1[i][1])	+ xi*c2.eval(-c2[i][0], c2[i][1])
	// qi = c1.eval(-c2[i][0], c2[i][1])	+			 c2.eval(-c1[i][0], c1[i][1])

	// u = (c1[i][0] + sqrt(xi). c2[i][0]) and	v = -(c1[i][1] + sqrt(xi). c2[i][1])

	bigint_matrix p(4,4), q(4,4), u(4,2), v(4,2);
	bigint gcd_tmp[4];
	for (unsigned int i = 0; i < 4; i++)
		{
			gcd_tmp[i] = 0;
			for (unsigned int j = 0; j < 4; j++) 
			if (j != i) // if i=j, p[i,j]=q[i,j]=0
				{
					// The point pi (resp. qi) is the ith line of matrix p (resp. q).
					if ((!c1[i].is_zero()) && (!c1[j].is_zero()))
						p.sto(i, j, c1[j][1]*c1[i][0] - c1[j][0]*c1[i][1]);
					if ((!c2[i].is_zero()) && (!c2[j].is_zero()))
						p.sto(i, j, p.member(i, j) + xi*c2[j][1]*c2[i][0] - xi*c2[j][0]*c2[i][1]);
					if ((!c1[j].is_zero()) && (!c2[i].is_zero()))
						q.sto(i, j, c1[j][1]*c2[i][0] - c1[j][0]*c2[i][1]);
					if ((!c1[i].is_zero()) && (!c2[j].is_zero()))
						q.sto(i, j, q.member(i, j) + c2[j][1]*c1[i][0] - c2[j][0]*c1[i][1]);
					gcd_tmp[i] = gcd(gcd_tmp[i], p.member(i, j));
					gcd_tmp[i] = gcd(gcd_tmp[i], q.member(i, j));
				}

			// gcd_tmp is the gcd of the coeffs of the points pi and qi
			if ((gcd_tmp[i] != 1) && (gcd_tmp[i] != 0))
				for (unsigned int j = 0; j < 4; j++)
					{
						p.sto(i, j, p.member(i,j)/gcd_tmp[i]);
						q.sto(i, j, q.member(i,j)/gcd_tmp[i]);
					}
			// Store the u,v parameters corresponding to the points
			if (!c1[i].is_zero())
				{
					u.sto(i,0, c1[i][0]);
					v.sto(i,0, -c1[i][1]);
				}
			if (!c2[i].is_zero())
				{
					u.sto(i,1, c2[i][0]);
					v.sto(i,1, -c2[i][1]);
				}
		}

	// If the parameterization c of the line has its ith coordinate = 0 (ie. c[i] = 0)
	// then the point Mi that we computed has all coordinates to 0 (and thus 
	// is not a projective point). 
	// Thus 2 points Mi and Mj (that have not all coordinates to 0) are equal 
	// iff Mi[j] = Mj[i]=0
	// Similarly the point (below) A (or B) on the line	 is equal to Mi iff A[i] = 0
	// that is A1[i] = A2[i] = 0
	// Hence if points A and B are used for parameterizing the line (as below)
	// we want to consider point Mi	 for reparametizing the line iff 
	// Mi != [0,0,0,0] that is sum_of_squares!=0 and if Mi is different from A and B
	// that is iff A[i] != 0 and B[i] != 0

	// Initialize with endpoints of c (Two endpoints A1 + sqrt(xi).A2 and B1 + sqrt(xi).B2)
	math_vector <bigint> A1 = c1.eval(1,0);
	math_vector <bigint> A2 = c2.eval(1,0);
	math_vector <bigint> B1 = c1.eval(0,1);
	math_vector <bigint> B2 = c2.eval(0,1);
	bigint sumA = sum_of_squares(A1) + xi*sum_of_squares(A2);
	bigint sumB = sum_of_squares(B1) + xi*sum_of_squares(B2);

	// The parameters of the endpoints: 6 values, 3 for each endpoint
	math_vector <bigint> u1(2,2), v1(2,2), u2(2,2), v2(2,2);
	u1[0] = 1;
	v2[0] = 1;
	bigint co1 = 1, co2 = 1;

	if (sumB < sumA)
		{
			swap(A1,B1);
			swap(A2,B2);
			swap(sumA,sumB);
			swap(u1,u2);
			swap(v1,v2);
		}

	for (int i = 0; i < 4; i++)
		{ // We consider replacing B by Mi only if A!=Mi that is A[i]!=0 
			// (otherwise the line AMi is not well defined)
			if ((A1[i] != 0) || (A2[i] != 0))
				{
					bigint sum_tmp = sum_of_squares((math_vector <bigint>)p.get_row_vector(i)) +
												xi*sum_of_squares((math_vector <bigint>)q.get_row_vector(i)); 
					if ((0 < sum_tmp) && (sum_tmp < sumA))
						{ // Replace B by Mi (by replacing B by A and A by Mi)
							sumB = sumA;
							sumA = sum_tmp;
							B1 = A1;
							B2 = A2;
							A1 = p.get_row_vector(i);
							A2 = q.get_row_vector(i);
							co2 = co1;
							co1 = gcd_tmp[i];
							u2 = u1;
							v2 = v1;
							u1 = u.get_row_vector(i);
							v1 = v.get_row_vector(i);
						}
					else if ((sumA < sum_tmp) && (sum_tmp < sumB))
						{ // Replace B by Mi 
							sumB = sum_tmp;
							B1 = p.get_row_vector(i);
							B2 = q.get_row_vector(i);
							co2 = gcd_tmp[i];
							u2 = u.get_row_vector(i);
							v2 = v.get_row_vector(i);				 
						}
				}
		}

	curve_param <bigint> line_new1(A1, B1);
	curve_param <bigint> line_new2(A2, B2);
	c1 = line_new1;
	c2 = line_new2;

	math_vector <bigint> end_old_param(6,6);
	if ((v1[0] == 0) && (v1[1] == 0))
		end_old_param[0] = 1;
	else
		{
			end_old_param[0] = co1*(v1[0]*v2[0]-xi*v1[1]*v2[1]);
			end_old_param[1] = co1*(v1[0]*v2[1]-v1[1]*v2[0]);
			end_old_param[2] = -co2*(v1[0]*v1[0]-xi*v1[1]*v1[1]);
		}

	if ((u1[0] == 0) && (u1[1] == 0))
		end_old_param[3] = 1;
	else
		{
			end_old_param[3] = co1*(u1[0]*u2[0]-xi*u1[1]*u2[1]);
			end_old_param[4] = co1*(u1[0]*u2[1]-u1[1]*u2[0]);
			end_old_param[5] = -co2*(u1[0]*u1[0]-xi*u1[1]*u1[1]);
		}

	return end_old_param;
}

// Try to reparameterize a line c = c1+sqrt(xi). c2 + sqrt(D). c3 + sqrt(D.xi). c4
// such that one of the coordinates is 0 when u=0 and another is zero when v=0
// (D = D1 + sqrt(xi). D2)
void improved_line_param(curve_param <bigint> &c1, curve_param <bigint> &c2, 
													curve_param <bigint> &c3, curve_param <bigint> &c4, 
													const bigint &xi, const bigint &D1, const bigint &D2)
{
	// c[i] = c1[i] + sqrt(xi). c2[i] + sqrt(D). c3[i] + sqrt(D.xi). c4[i]
	//		= (c1[i][1] + sqrt(xi). c2[i][1] + sqrt(D). c3[i][1] + sqrt(D.xi). c4[i][1]). u
	//		= (c1[i][0] + sqrt(xi). c2[i][0] + sqrt(D). c3[i][0] + sqrt(D.xi). c4[i][0]). v
	// We consider	for i from 0 to 3 the point Mi = c(u,v) such that c[i]=0, that
	// is where	 
	// u = -(c1[i][0] + sqrt(xi). c2[i][0] + sqrt(D). c3[i][0] + sqrt(D.xi). c4[i][0])	and
	// v = (c1[i][1] + sqrt(xi). c2[i][1] + sqrt(D). c3[i][1] + sqrt(D.xi). c4[i][1])
	// The jth coordinate of this point is 
	//		- (c1[j][1] + sqrt(xi). c2[j][1] + sqrt(D). c3[j][1] + sqrt(D.xi). c4[j][1]). 
	//			 (c1[i][0] + sqrt(xi). c2[i][0] + sqrt(D). c3[i][0] + sqrt(D.xi). c4[i][0])
	//		+ (c1[j][0] + sqrt(xi). c2[j][0] + sqrt(D). c3[j][0] + sqrt(D.xi). c4[j][0]).
	//			 (c1[i][1] + sqrt(xi). c2[i][1] + sqrt(D). c3[i][1] + sqrt(D.xi). c4[i][1])
	// 
	// =	(-c1[j][1]*c1[i][0] + c1[j][0]*c1[i][1] - xi*c2[j][1]*c2[i][0] +
	// xi*c2[i][1]*c2[j][0] -D1*c3[j][1]*c3[i][0] + D1*c3[i][1]*c3[j][0] 
	//			 -D1*xi*c4[j][1]*c4[i][0] + D1*xi*c4[i][1]*c4[j][0] 
	//			 -D2*xi*c3[j][1]*c4[i][0] -D2*xi*c3[i][0]*c4[j][1]
	//			 +D2*xi*c3[i][1]*c4[j][0]	 +D2*xi*c3[j][0]*c4[i][1]) 
	//				+sqrt(xi)*(-c2[j][1]*c1[i][0] -c1[j][1]*c2[i][0] + c2[i][1]*c1[j][0]
	//				+c1[i][1]*c2[j][0] -D2*c3[j][1]*c3[i][0] + D2*c3[i][1]*c3[j][0] 
	// 			 -D2*xi*c4[j][1]*c4[i][0] + D2*xi*c4[i][1]*c4[j][0] 
	// 			 -D1*c3[j][1]*c4[i][0] -D1*c3[i][0]*c4[j][1] 
	// 			 +D1*c3[i][1]*c4[j][0]	+D1*c3[j][0]*c4[i][1])
	//				+sqrt(D)*(-c3[j][1]*c1[i][0] -c1[j][1]*c3[i][0] + c3[i][1]*c1[j][0]
	//				+c1[i][1]*c3[j][0] -xi*c2[j][1]*c4[i][0] -xi*c2[i][0]*c4[j][1] +
	//				xi*c2[i][1]*c4[j][0] +xi*c2[j][0]*c4[i][1])
	//				+sqrt(xi)*sqrt(D)*(-c4[j][1]*c1[i][0] -c1[j][1]*c4[i][0] +
	//				c4[i][1]*c1[j][0] +c1[i][1]*c4[j][0] -c3[j][1]*c2[i][0]
	//				-c3[i][0]*c2[j][1] + c3[i][1]*c2[j][0] + c3[j][0]*c2[i][1]) 
	// 
	// which is zero for j=i. 
	// Let pi, qi, ri, si be the points such that Mi = pi + sqrt(xi). qi +
	// sqrt(D).ri + sqrt(D).sqrt(xi). si 
	// Their	coordinates are 
	// pi[j] = -c1[j][1]*c1[i][0] + c1[j][0]*c1[i][1] - xi*c2[j][1]*c2[i][0] +
	// xi*c2[i][1]*c2[j][0] -D1*c3[j][1]*c3[i][0] + D1*c3[i][1]*c3[j][0] 
	//			 -D1*xi*c4[j][1]*c4[i][0] + D1*xi*c4[i][1]*c4[j][0] 
	//			 -D2*xi*c3[j][1]*c4[i][0] -D2*xi*c3[i][0]*c4[j][1]
	//			 +D2*xi*c3[i][1]*c4[j][0]	 +D2*xi*c3[j][0]*c4[i][1] 
	// qi[j] = -c2[j][1]*c1[i][0] -c1[j][1]*c2[i][0] + c2[i][1]*c1[j][0] +c1[i][1]*c2[j][0]
	// 			 -D2*c3[j][1]*c3[i][0] + D2*c3[i][1]*c3[j][0] 
	// 			 -D2*xi*c4[j][1]*c4[i][0] + D2*xi*c4[i][1]*c4[j][0] 
	// 			 -D1*c3[j][1]*c4[i][0] -D1*c3[i][0]*c4[j][1] 
	// 			 +D1*c3[i][1]*c4[j][0]	+D1*c3[j][0]*c4[i][1]
	// ri[j] = -c3[j][1]*c1[i][0] -c1[j][1]*c3[i][0] + c3[i][1]*c1[j][0] +c1[i][1]*c3[j][0]
	// 		 -xi*c2[j][1]*c4[i][0] -xi*c2[i][0]*c4[j][1] +
	// 		 xi*c2[i][1]*c4[j][0] +xi*c2[j][0]*c4[i][1] 
	// si[j] = -c4[j][1]*c1[i][0] -c1[j][1]*c4[i][0] + c4[i][1]*c1[j][0] +c1[i][1]*c4[j][0]
	// 				-c3[j][1]*c2[i][0] -c3[i][0]*c2[j][1] +
	// 				c3[i][1]*c2[j][0] + c3[j][0]*c2[i][1] 

	bigint_matrix p(4,4), q(4,4), r(4,4), s(4,4);
	bigint gcd_tmp;
	for (unsigned int i = 0; i < 4; i++)
		{
			gcd_tmp = 0;
			for (unsigned int j = 0; j < 4; j++) 
			if (j != i) // if i=j, p[i,j]=q[i,j]=0
				{
					// The points pi (resp. qi) is	 the ith line of matrix p (resp. q).
					if ((!c1[i].is_zero()) && (!c1[j].is_zero()))
						p.sto(i, j, -c1[j][1]*c1[i][0] + c1[j][0]*c1[i][1]);
					if ((!c2[i].is_zero()) && (!c2[j].is_zero()))
						p.sto(i, j, p.member(i, j) - xi*c2[j][1]*c2[i][0] + xi*c2[j][0]*c2[i][1]);
					if ((!c3[i].is_zero()) && (!c3[j].is_zero()))
						p.sto(i, j, p.member(i, j) -D1*c3[j][1]*c3[i][0] + D1*c3[i][1]*c3[j][0]);
					if ((!c4[i].is_zero()) && (!c4[j].is_zero()))
						p.sto(i,j,p.member(i,j)-D1*xi*c4[j][1]*c4[i][0]+D1*xi*c4[i][1]*c4[j][0]);
					if ((!c3[j].is_zero()) && (!c4[i].is_zero()))
						p.sto(i,j,p.member(i,j)-D2*xi*c3[j][1]*c4[i][0]+D2*xi*c3[j][0]*c4[i][1]);
					if ((!c3[i].is_zero()) && (!c4[j].is_zero()))
						p.sto(i,j,p.member(i,j)-D2*xi*c3[i][0]*c4[j][1]+D2*xi*c3[i][1]*c4[j][0]);
		
					if ((!c1[j].is_zero()) && (!c2[i].is_zero()))
						q.sto(i, j, -c1[j][1]*c2[i][0] + c1[j][0]*c2[i][1]);
					if ((!c1[i].is_zero()) && (!c2[j].is_zero()))
						q.sto(i, j, q.member(i, j) - c2[j][1]*c1[i][0] + c2[j][0]*c1[i][1]);
					if ((!c3[i].is_zero()) && (!c3[j].is_zero()))
						q.sto(i, j, q.member(i, j) -D2*c3[j][1]*c3[i][0] + D2*c3[i][1]*c3[j][0]);
					if ((!c4[i].is_zero()) && (!c4[j].is_zero()))
						q.sto(i,j,q.member(i,j)-D2*xi*c4[j][1]*c4[i][0] + D2*xi*c4[i][1]*c4[j][0]);
					if ((!c3[j].is_zero()) && (!c4[i].is_zero()))
						q.sto(i, j, q.member(i, j) -D1*c3[j][1]*c4[i][0] +D1*c3[j][0]*c4[i][1]);
					if ((!c3[i].is_zero()) && (!c4[j].is_zero()))
						q.sto(i, j, q.member(i, j) -D1*c3[i][0]*c4[j][1] +D1*c3[i][1]*c4[j][0]);
		
					if ((!c1[i].is_zero()) && (!c3[j].is_zero()))
						r.sto(i, j, -c3[j][1]*c1[i][0] +c1[i][1]*c3[j][0]);
					if ((!c1[j].is_zero()) && (!c3[i].is_zero()))
						r.sto(i, j, r.member(i, j) -c1[j][1]*c3[i][0] + c3[i][1]*c1[j][0]);
					if ((!c2[i].is_zero()) && (!c4[j].is_zero()))
						r.sto(i, j, r.member(i, j) -xi*c2[i][0]*c4[j][1] + xi*c2[i][1]*c4[j][0]);
					if ((!c2[j].is_zero()) && (!c4[i].is_zero()))
						r.sto(i, j, r.member(i, j) -xi*c2[j][1]*c4[i][0] +xi*c2[j][0]*c4[i][1]);
		
					if ((!c1[i].is_zero()) && (!c4[j].is_zero()))
						s.sto(i, j, -c4[j][1]*c1[i][0] +c1[i][1]*c4[j][0]);
					if ((!c1[j].is_zero()) && (!c4[i].is_zero()))
						s.sto(i, j, s.member(i, j) -c1[j][1]*c4[i][0] + c4[i][1]*c1[j][0]);
					if ((!c2[i].is_zero()) && (!c3[j].is_zero()))
						s.sto(i, j, s.member(i, j) -c3[j][1]*c2[i][0] + c3[j][0]*c2[i][1]);
					if ((!c2[j].is_zero()) && (!c3[i].is_zero()))
						s.sto(i, j, s.member(i, j) -c3[i][0]*c2[j][1] + c3[i][1]*c2[j][0]);
		
					gcd_tmp = gcd(gcd_tmp, p.member(i, j));
					gcd_tmp = gcd(gcd_tmp, q.member(i, j));
					gcd_tmp = gcd(gcd_tmp, r.member(i, j));
					gcd_tmp = gcd(gcd_tmp, s.member(i, j));
				}
			// gcd_tmp is the gcd of the coeffs of the points pi and qi
			if (gcd_tmp != 1) 
				for (unsigned int j = 0; j < 4; j++)
		 			{
			 			p.sto(i, j, p.member(i,j)/gcd_tmp);
			 			q.sto(i, j, q.member(i,j)/gcd_tmp);
			 			r.sto(i, j, r.member(i,j)/gcd_tmp);
			 			s.sto(i, j, s.member(i,j)/gcd_tmp);
		 			}
		}

	// If the parameterization c of the line has its ith coordinate =0 (ie. c[i]=0)
	// then the point Mi that we computed has all coordinates to 0 (and thus 
	//	is not a projective point). 
	// Thus 2 points Mi and Mj (that have not all coordinates to 0) are equal 
	// iff Mi[j]=Mj[i]=0
	// Similarly the point (below) A (or B) on the line	 is equal to Mi iff A[i]=0
	// that is A1[i]=A2[i]=A3[i]=A4[i]=0
	// Hence if points A and B are used for parameterizing the line (as below)
	// we want to consider point Mi	 for reparametizing the line iff 
	// Mi!=[0,0,0,0] that is sum_of_squares!=0 and if Mi is different from A and B
	// that is iff A[i]!=0 and B[i]!=0
	
	// Initialize with "endpoints" of c (for (u,v) = (1,0) and (0,1))
	// (Two endpoints A1 + sqrt(xi).A2 + sqrt(D).A3 + sqrt(D).sqrt(xi).A4 
	// and B1 + sqrt(xi).B2 + sqrt(D).B3 + sqrt(D).sqrt(xi).B4 )
	math_vector <bigint> A1 = c1.eval(1,0);
	math_vector <bigint> A2 = c2.eval(1,0);
	math_vector <bigint> A3 = c3.eval(1,0);
	math_vector <bigint> A4 = c4.eval(1,0);
	math_vector <bigint> B1 = c1.eval(0,1);
	math_vector <bigint> B2 = c2.eval(0,1);
	math_vector <bigint> B3 = c3.eval(0,1);
	math_vector <bigint> B4 = c4.eval(0,1);
	bigint d = abs(D1)+abs(D2);
	bigint sumA = sum_of_squares(A1) + sum_of_squares(A2) + sum_of_squares(A3)
								+sum_of_squares(A4); 
	bigint sumB = sum_of_squares(B1) + sum_of_squares(B2) + sum_of_squares(B3)
								+sum_of_squares(B4); 

	if (sumB < sumA)
		{
			swap(A1,B1);
			swap(A2,B2);
			swap(A3,B3);
			swap(A4,B4);
			swap(sumA,sumB);
		}

	for (int i = 0; i < 4; i++)
		{ // We consider replacing B by Mi only if A!=Mi that is A[i]!=0 
			// (otherwise the line AMi is not well defined
			if ((A1[i] != 0) || (A2[i] != 0) || (A3[i] != 0) || (A4[i] != 0))
				{
		 			bigint sum_tmp = sum_of_squares((math_vector <bigint>)p.get_row_vector(i))
			 										+ sum_of_squares((math_vector <bigint>)q.get_row_vector(i))
			 										+ sum_of_squares((math_vector <bigint>)r.get_row_vector(i))
			 										+ sum_of_squares((math_vector <bigint>)s.get_row_vector(i)); 
					if ((0 < sum_tmp) && (sum_tmp < sumA))
						{ // Replace B by Mi (by replacing B by A and A by Mi)
							sumB = sumA;
							sumA = sum_tmp;
							B1 = A1;
							B2 = A2;
							B3 = A3;
							B4 = A4;
							A1 = p.get_row_vector(i);
							A2 = q.get_row_vector(i);
							A3 = r.get_row_vector(i);
							A4 = s.get_row_vector(i);
						}
					else if ((sumA < sum_tmp) && (sum_tmp < sumB))
						{ // Replace B by Mi 
							sumB = sum_tmp;
							B1 = p.get_row_vector(i);
							B2 = q.get_row_vector(i);
							B3 = r.get_row_vector(i);
							B4 = s.get_row_vector(i);
						}
				}
		}
	curve_param <bigint> line_new1(A1, B1);
	curve_param <bigint> line_new2(A2, B2);
	curve_param <bigint> line_new3(A3, B3);
	curve_param <bigint> line_new4(A4, B4);
	c1 = line_new1;
	c2 = line_new2;
	c3 = line_new3;
	c4 = line_new4;
}

} // end of namespace QI

