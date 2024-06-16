// Things for checking whether computed params are alright

#include <fstream>

#include <libqi/kernel/QICheck.h>
#include <libqi/kernel/QIElem.h>

using namespace rpl;

// Enter namespace QI
namespace QI {

// Print error in file ERROR
void print_error(const bigint_matrix &q1, const bigint_matrix &q2)
{
	ostream *tmp_out;
	tmp_out = new ofstream("ERROR",ios::app);

	if (!tmp_out)
		{
			cout << "Error opening file ERROR for writing" << endl;
			exit(-1);
		}

	ostream &out_stream = *tmp_out;
	out_stream << "//* FAILURE WITH THE FOLLOWING TWO QUADRICS"<< endl;

	out_stream << quadtovec(q1) << endl;
	out_stream << quadtovec(q2) << endl;

	delete tmp_out;
}

// Checks if curve_param is zero param
bool is_zero(const curve_param <bigint> &p)
{
	bool flag = 1;

	for (rpl_size_t i = 0; i < p.capacity(); i++)
		if (!p[i].is_zero())
			{
				flag = 0;
				break;
			}

	return flag;
}

// Checks if surface_param is zero param
bool is_zero(const surface_param <bigint> &p)
{
	bool flag = 1;

	for (rpl_size_t i = 0; i < p.capacity(); i++)
		if (!p[i].is_zero())
			{
				flag = 0;
				break;
			}

	return flag;
}

// Maybe the given parameters do not represent a curve...
bool wrong_param(const curve_param <bigint> &p)
{
	// Take three sample points. If they represent the same (projective) point,
	// something fishy is going on

	if ((p[0].degree() <= 0) && (p[1].degree() <= 0) && (p[2].degree() <= 0)
			&& (p[3].degree() <= 0))
		return 0;
	else
		{
			math_vector <bigint> p1 = p.eval(1,0);
			math_vector <bigint> p2 = p.eval(0,1);
			math_vector <bigint> p3 = p.eval(1,1);

			if (are_equal(p1,p2) && are_equal(p2,p3))
				return 1;
			else
				return 0;
		}
}

// Checks if the computed param is ok by replugging in initial quadrics
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
									const curve_param <bigint> &p, ostream &s)
{
	if (is_zero(p))
		{
			s << "FAILED" << endl;
			s << "			 ----> zero param" << endl;
			print_error(q1, q2);
		}
	else if (wrong_param(p))
		{
			s << "FAILED" << endl;
			s << "			 ----> wrong param" << endl;
			print_error(q1,q2);
		}
	else
		{
			hom_polynomial <bigint> hom_1 = plug_param_in_quadric(p,q1,p);
			hom_polynomial <bigint> hom_2 = plug_param_in_quadric(p,q2,p);

			if ((hom_1.is_zero()) && (hom_2.is_zero()))
				s << "ok" << endl;
			else
				{
					s << "FAILED" << endl;
					s << "				----> " << hom_1 << endl;
					s << "				----> " << hom_2 << endl;
					print_error(q1, q2);
				}
		}
}

// Checks if the computed param is ok by replugging in initial quadrics
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
									const surface_param <bigint> &p, ostream &s)
{
	s << ">> checking param: ";

	if (is_zero(p))
		{
			s << "FAILED" << endl;
			s << "			 ----> zero param" << endl;
			print_error(q1, q2);
		}
	else
		{
			hom_hom_polynomial <bigint> hom_1 = plug_param_in_quadric(p,q1,p);
			hom_hom_polynomial <bigint> hom_2 = plug_param_in_quadric(p,q2,p);

			if ((hom_1.is_zero()) && (hom_2.is_zero()))
				s << "ok" << endl;
			else
				{
					s << "FAILED" << endl;
					s << "				----> " << hom_1 << endl;
					s << "				----> " << hom_2 << endl;
					print_error(q1, q2);
				}
		}
}

// Checks if p1+sqrt(D)*p2 is ok
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
									const curve_param <bigint> &p1, const curve_param <bigint> &p2, 
									const bigint &D, ostream &s)
{
	s << ">> checking param: ";

	if ((D.is_zero()) || (is_zero(p1)) || (is_zero(p2)))
		{
			s << "FAILED" << endl;
			s << "			 ----> zero param" << endl;
			print_error(q1, q2);
		}
	else
		{
			hom_polynomial <bigint> hom2_1;
			add(hom2_1,plug_param_in_quadric(p1,q1,p1),D*plug_param_in_quadric(p2,q1,p2));
			hom_polynomial <bigint> hom3_1 = plug_param_in_quadric(p1,q1,p2);

			hom_polynomial <bigint> hom2_2;
			add(hom2_2,plug_param_in_quadric(p1,q2,p1),D*plug_param_in_quadric(p2,q2,p2));
			hom_polynomial <bigint> hom3_2 = plug_param_in_quadric(p1,q2,p2);

			if ((hom2_1.is_zero()) && (hom2_2.is_zero()) &&
		 			(hom3_1.is_zero()) && (hom3_2.is_zero()))
				s << "ok" << endl;
			else
				{
					s << "FAILED" << endl;
					s << "				----> " << hom2_1 << "	" << hom3_1 << endl;
					s << "				----> " << hom2_2 << "	" << hom3_2 << endl;
					print_error(q1, q2);
				}
		}
}	

// Checks if param c1+sqrt(a)*c2+sqrt(b)*c3+sqrt(a*b)*c4 is ok
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
									const curve_param <bigint> &c1, const curve_param <bigint> &c2, 
									const curve_param <bigint> &c3, const curve_param <bigint> &c4,
 									const bigint &a, const bigint &b, ostream &s)
{
	s << ">> checking param: ";

	if ((a.is_zero()) || (b.is_zero()) || 
			((is_zero(c1)) && (is_zero(c2)) && (is_zero(c3)) && (is_zero(c4))))
		{
			s << "FAILED" << endl;
			s << "			 ----> zero param" << endl;
			print_error(q1, q2);
		}
	else
		{
			hom_polynomial <bigint> hom1_1 = plug_param_in_quadric(c1,q1,c1);
			add(hom1_1,hom1_1,a*plug_param_in_quadric(c2,q1,c2));
			add(hom1_1,hom1_1,b*plug_param_in_quadric(c3,q1,c3));
			add(hom1_1,hom1_1,a*b*plug_param_in_quadric(c4,q1,c4));

			hom_polynomial <bigint> hom2_1;
			add(hom2_1,plug_param_in_quadric(c1,q1,c2),b*plug_param_in_quadric(c3,q1,c4));

			hom_polynomial <bigint> hom3_1;
			add(hom3_1,plug_param_in_quadric(c1,q1,c3),a*plug_param_in_quadric(c2,q1,c4));

			hom_polynomial <bigint> hom4_1;
			add(hom4_1,plug_param_in_quadric(c1,q1,c4),plug_param_in_quadric(c2,q1,c3));

			hom_polynomial <bigint> hom1_2 = plug_param_in_quadric(c1,q2,c1);
			add(hom1_2,hom1_2,a*plug_param_in_quadric(c2,q2,c2));
			add(hom1_2,hom1_2,b*plug_param_in_quadric(c3,q2,c3));
			add(hom1_2,hom1_2,a*b*plug_param_in_quadric(c4,q2,c4));

			hom_polynomial <bigint> hom2_2;
			add(hom2_2,plug_param_in_quadric(c1,q2,c2),b*plug_param_in_quadric(c3,q2,c4));

			hom_polynomial <bigint> hom3_2;
			add(hom3_2,plug_param_in_quadric(c1,q2,c3),a*plug_param_in_quadric(c2,q2,c4));

			hom_polynomial <bigint> hom4_2;
			add(hom4_2,plug_param_in_quadric(c1,q2,c4),plug_param_in_quadric(c2,q2,c3));

			if ((hom1_1.is_zero()) && (hom2_1.is_zero()) && (hom3_1.is_zero()) && 
		 			(hom4_1.is_zero()) && (hom1_2.is_zero()) && (hom2_2.is_zero()) && 
		 			(hom3_2.is_zero()) && (hom4_2.is_zero()))
				s << "ok" << endl;
			else
				{
					s << "FAILED" << endl;
					s << "				----> " << hom1_1 << "	" << hom2_1 << "	" << hom3_1 
						<< "	 " << hom4_1 << endl;
					s << "				----> " << hom1_2 << "	" << hom2_2 << "	" << hom3_2 
						<< "	 " << hom4_2 << endl;
					print_error(q1, q2);
				}
		}
}

// Checks if param c1+sqrt(a)*c2+(c3+sqrt(a)*c4)*sqrt(b+c*sqrt(a)) is ok
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
									const curve_param <bigint> &c1, const curve_param <bigint> &c2, 
									const curve_param <bigint> &c3, const curve_param <bigint> &c4,
									const bigint &a, const bigint &b, const bigint &c, ostream &s)
{
	s << ">> checking param: ";

	if ((is_zero(c1)) && (is_zero(c2)) && (is_zero(c3)) && (is_zero(c4)))
		{
			s << "FAILED" << endl;
			s << "			 ----> zero param" << endl;
			print_error(q1, q2);
		}
	else
		{
			hom_polynomial <bigint> p11 = plug_param_in_quadric(c1,q1,c1);
			hom_polynomial <bigint> p22 = plug_param_in_quadric(c2,q1,c2);
			hom_polynomial <bigint> p33 = plug_param_in_quadric(c3,q1,c3);
			hom_polynomial <bigint> p44 = plug_param_in_quadric(c4,q1,c4);
			hom_polynomial <bigint> p12 = plug_param_in_quadric(c1,q1,c2);
			hom_polynomial <bigint> p13 = plug_param_in_quadric(c1,q1,c3);
			hom_polynomial <bigint> p14 = plug_param_in_quadric(c1,q1,c4);
			hom_polynomial <bigint> p23 = plug_param_in_quadric(c2,q1,c3);
			hom_polynomial <bigint> p24 = plug_param_in_quadric(c2,q1,c4);
			hom_polynomial <bigint> p34 = plug_param_in_quadric(c3,q1,c4);

			hom_polynomial <bigint> hom1_1 = p11;
			add(hom1_1,hom1_1,a*p22);
			add(hom1_1,hom1_1,b*p33);
			add(hom1_1,hom1_1,a*b*p44);
			add(hom1_1,hom1_1,2*a*c*p34);

			hom_polynomial <bigint> hom2_1 = p12;
			add(hom2_1,hom2_1,p12);
			add(hom2_1,hom2_1,c*p33);
			add(hom2_1,hom2_1,a*c*p44);
			add(hom2_1,hom2_1,2*b*p34);
			multiply(hom2_1,hom2_1,a); // for the case where a =0

			hom_polynomial <bigint> hom3_1, hom4_1;
			if (sign(b,c,a) != 0) // sign(b+c*sqrt(a))
				{
					add(hom3_1,p13,p13);
					add(hom3_1,hom3_1,2*a*p24);
			
					add(hom4_1,p14,p14);
					add(hom4_1,hom4_1,p23);
					add(hom4_1,hom4_1,p23);
					multiply(hom4_1, hom4_1, a); // for the case where a =0
				}

			p11 = plug_param_in_quadric(c1,q2,c1);
			p22 = plug_param_in_quadric(c2,q2,c2);
			p33 = plug_param_in_quadric(c3,q2,c3);
			p44 = plug_param_in_quadric(c4,q2,c4);
			p12 = plug_param_in_quadric(c1,q2,c2);
			p13 = plug_param_in_quadric(c1,q2,c3);
			p14 = plug_param_in_quadric(c1,q2,c4);
			p23 = plug_param_in_quadric(c2,q2,c3);
			p24 = plug_param_in_quadric(c2,q2,c4);
			p34 = plug_param_in_quadric(c3,q2,c4);

			hom_polynomial <bigint> hom1_2 = p11;
			add(hom1_2,hom1_2,a*p22);
			add(hom1_2,hom1_2,b*p33);
			add(hom1_2,hom1_2,a*b*p44);
			add(hom1_2,hom1_2,2*a*c*p34);

			hom_polynomial <bigint> hom2_2 = p12;
			add(hom2_2,hom2_2,p12);
			add(hom2_2,hom2_2,c*p33);
			add(hom2_2,hom2_2,a*c*p44);
			add(hom2_2,hom2_2,2*b*p34);
			multiply(hom2_2, hom2_2, a); // for the case where a =0

			hom_polynomial <bigint> hom3_2, hom4_2;
			if (sign(b, c, a) != 0) // sign(b+c*sqrt(a))
				{
					add(hom3_2,p13,p13);
					add(hom3_2,hom3_2,2*a*p24);
			
					add(hom4_2,p14,p14);
					add(hom4_2,hom4_2,p23);
					add(hom4_2,hom4_2,p23);
					multiply(hom4_2,hom4_2,a); // for the case where a =0
				}

			if ((hom1_1.is_zero()) && (hom2_1.is_zero()) && (hom3_1.is_zero()) && 
		 			(hom4_1.is_zero()) && (hom1_2.is_zero()) && (hom2_2.is_zero()) && 
		 			(hom3_2.is_zero()) && (hom4_2.is_zero()))
				s << "ok" << endl;
			else
				{
					s << "FAILED" << endl;
					s << "				----> " << hom1_1 << ",	 " << hom2_1 << ",	" << hom3_1 << ",	 "
						<< hom4_1 << endl;
					s << "				----> " << hom1_2 << ",	 " << hom2_2 << ",	" << hom3_2 << ",	 "
						<< hom4_2 << endl;
					print_error(q1, q2);
				}
		}
}

// Used by following procedure to check if a param is ok
bool param_ok(const hom_hom_polynomial <bigint> &hom, 
							const hom_polynomial <bigint> &pol)
{
	if (hom.is_zero())
		return 1;
	else
		{
			rpl_size_t i;

			for (i = 0; i <= hom.degree(); i++)
				if (!hom[i].is_zero())
		 			break;

			rpl_size_t d = xdegree(hom[i]);

			if (pol.degree() == 4)
				{
					if ((hom[i]*pol[d]-pol*hom[i][d]).is_zero())
						return 1;
					else
						return 0;
				}
			else // degree(pol) = 3 and hom[i] = c*s*pol ou c*t*pol
				{
					hom_polynomial <bigint> hom2 = hom[i];
					hom_polynomial <bigint> div;
			
					if (hom2[4] == 0)
						div.assign_y();
					else
						div.assign_x();
			
					divide(hom2,hom2,div);
			
					rpl_size_t d = xdegree(hom2);
					
					if ((hom2*pol[d]-pol*hom2[d]).is_zero())
						return 1;
					else
						return 0;
				}
		}
}

// Checks if the computed param is ok by replugging in initial quadrics
// Second case, when instead we have a surface param and a polynomial
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
									const surface_param <bigint> &p, const hom_polynomial <bigint> &pol, 
									ostream &s)
{
	if (((p[0].is_zero()) && (p[1].is_zero()) && (p[2].is_zero()) && (p[3].is_zero())) ||
			(pol.is_zero()))
		{
			s << "FAILED" << ": zero param" << endl;
			print_error(q1, q2);
		}
	else
		{
			s << ">> checking param: ";
	
			hom_hom_polynomial <bigint> hom = plug_param_in_quadric(p,q1,p);

			if (param_ok(hom,pol))
				{
		 			hom = plug_param_in_quadric(p,q2,p);

					if (param_ok(hom,pol))
			 			s << "ok" << endl;
		 			else
			 			{
				 			s << "FAILED" << endl;
				 			s << ">>>>>> " << hom << endl;
				 			print_error(q1, q2);
			 			}
				}
			else
				{
		 			s << "FAILED" << endl;
		 			s << ">>>>>> " << hom << endl;
		 			print_error(q1, q2);
				}
		}
}

// Checks if the parameterization s1 + sqrt(det_q[1]) s2 of the quadric 22 q lies on q
void check_param(const bigint_matrix &q, const surface_param <bigint> &s1, 
									const surface_param <bigint> &s2, const bigint &delta, ostream &s)
{
	s << ">> checking that the param of the surface (2,2) lies on the quadric (2,2): ";
	hom_hom_polynomial <bigint> hhverif1, hhverif2, hh_tmp;

	// verif= s1.q.s1 + d. s2.q.s2 + sqrt(d). 2. s1.q.s2
	hhverif1= plug_param_in_quadric(s1,q,s1); 
	hh_tmp = plug_param_in_quadric(s2,q,s2);
	multiply(hh_tmp, hh_tmp, delta);
	add(hhverif1, hhverif1, hh_tmp);	
	hhverif2= plug_param_in_quadric(s1,q,s2);	 

	if ((hhverif1.is_zero()) && (hhverif2.is_zero()))
		s << "ok" << endl;
	else
		{
			s << "FAILED" << endl;
			s << hhverif1 << ",	 " << hhverif2 << endl << endl;
		}
}

// Checks if the computed quartic param c is ok by replugging in initial quadrics
//	 c = c1 + sqrt(d). c2 + eps. sqrt(Delta). (c3 + sqrt(d). c4)
// c1, c2: degree 3, c3, c4: degree 1 (generically)
// Delta = Delta1 + sqrt(d). Delta2
// Note if check_param is ok for c then it is also ok after replacing 
// sqrt(d) by -sqrt(d) and sqrt(Delta) by -sqrt(Delta)
void check_param(const bigint_matrix &q1, const bigint_matrix &q2, 
									const curve_param <bigint> &c1, const curve_param <bigint> &c2, 
									const curve_param <bigint> &c3, const curve_param <bigint> &c4, 
									const bigint &delta, const hom_polynomial <bigint> &Delta1, 
									const hom_polynomial <bigint> &Delta2, ostream &s)
{
	// Verification that the curve c is in the quadric q1 and	 q2
	// c.q.c = c1.q.c1 + d. c2.q.c2 + Delta.(c3.q.c3 + d. c4.q.c4 + 2.sqrt(d). c3.q.c4)
	// + 2. sqrt(d). c1.q.c2 + 2.eps. sqrt(Delta). (c1.q.c3 + d. c2.q.c4 +
	// sqrt(d). c1.q.c4 + sqrt(d). c2.q.c3)
	// c.q.c = 
	//	 [c1.q.c1 + d. c2.q.c2 + Delta1.(c3.q.c3 + d. c4.q.c4) + Delta2.( 2.d. c3.q.c4) ]
	//	+sqrt(d). [ 2.c1.q.c2 + Delta1.(2.c3.q.c4) + Delta2.(c3.q.c3 + d. c4.q.c4) ]
	//	+2.eps. sqrt(Delta). [ c1.q.c3 + d. c2.q.c4 ]
	//	+2.eps. sqrt(Delta). sqrt(d). [ c1.q.c4 + c2.q.c3 ] 
	hom_polynomial <bigint> verif, verif1, verif2, verif3, verif4, htmp;
	bool check = 1; 
	s << ">> checking param: ";

	for (unsigned int i = 1; i < 3; i++)
		{
			bigint_matrix qq;
			if (i == 1)
				qq = q1;
			else 
				qq = q2;

			// verif1= c1.q.c1 + d. c2.q.c2 + Delta1.(c3.q.c3 + d. c4.q.c4) + Delta2.(
			// 2.d. c3.q.c4)	
			verif= plug_param_in_quadric(c1,qq,c1);	 // verif= c1.q.c1
			htmp = plug_param_in_quadric(c2,qq,c2);
			multiply(htmp, htmp, delta);
			add(verif, verif, htmp);	// verif= c1.q.c1 + d. c2.q.c2 
			htmp = plug_param_in_quadric(c3,qq,c3);
			multiply(htmp, htmp, Delta1);
			add(verif, verif, htmp);	// verif= c1.q.c1 + d. c2.q.c2 + Delta1.(c3.q.c3)
			htmp = plug_param_in_quadric(c4,qq,c4);
			multiply(htmp, htmp, delta);
			multiply(htmp, htmp, Delta1);
			add(verif, verif, htmp); // verif= c1.q.c1+d.c2.q.c2+Delta1.(c3.q.c3+d. c4.q.c4) 
			htmp = plug_param_in_quadric(c3,qq,c4);
			multiply(htmp, htmp, delta);
			multiply(htmp, htmp, (bigint)2);
			multiply(htmp, htmp, Delta2);
			add(verif, verif, htmp); // verif= verif1
			verif1.assign(verif);

			// verif2 =	 2.c1.q.c2 + Delta1.(2.c3.q.c4) + Delta2.(c3.q.c3 + d. c4.q.c4) 
			verif= plug_param_in_quadric(c1,qq,c2);	 
			multiply(verif, verif, (bigint)2);	//verif= 2.c1.q.c2 
			htmp = plug_param_in_quadric(c3,qq,c4);
			multiply(htmp, htmp, (bigint)2);
			multiply(htmp, htmp, Delta1);
			add(verif, verif, htmp);	// verif= 2.c1.q.c2 + Delta1.(2.c3.q.c4) 
			htmp = plug_param_in_quadric(c3,qq,c3);
			multiply(htmp, htmp, Delta2);
			add(verif, verif, htmp);	// verif= 2.c1.q.c2+Delta1.(2.c3.q.c4)+Delta2.(c3.q.c3)
			htmp = plug_param_in_quadric(c4,qq,c4);
			multiply(htmp, htmp, delta);
			multiply(htmp, htmp, Delta2);
			add(verif, verif, htmp); // verif= verif2
			verif2.assign(verif);

			// verif3 = [ c1.q.c3 + d. c2.q.c4 ]
			verif= plug_param_in_quadric(c1,qq,c3);	 
			htmp = plug_param_in_quadric(c2,qq,c4);
			multiply(htmp, htmp, delta);
			add(verif, verif, htmp); // verif= verif3
			verif3.assign(verif);

			// verif4 = [ c1.q.c4 + c2.q.c3 ] 
			verif= plug_param_in_quadric(c1,qq,c4);	 
			htmp = plug_param_in_quadric(c2,qq,c3);
			add(verif, verif, htmp); // verif= verif4
			verif4.assign(verif);
			if ((!verif1.is_zero()) || (!verif2.is_zero()) || (!verif3.is_zero()) ||
		 			(!verif4.is_zero())) 
				check = 0;
		}

	if (check)
		s << "ok" << endl;
	else
		{
			s << "FAILED" << endl;
			print_error(q1, q2);
		}
}

} // end of namespace QI
