// Intersection when the determinantal equation has two double roots

#include <libqi/kernel/QITwoMult.h>
#include <libqi/kernel/QIElem.h>
#include <libqi/kernel/QINumber.h>
#include <libqi/kernel/QIInterWith22.h>

#ifndef NDEBUG
#include <libqi/kernel/QICheck.h>
#endif

using namespace rpl;

// Enter namespace QI
namespace QI {

//////////////////////////////////////////////////////////////////////////////////
// Procedure for the conic and two lines not crossing on the conic case
quad_inter <bigint> inter_conic_2lines_not_crossing(const bigint_matrix &q1, 
																										const bigint_matrix &q2,
																										const bigint_matrix &qa, 
																										const bigint_matrix &qb, 
																										const int img_cone, 
																										const rpl_size_t &s_e, 
																										const int opt_level, ostream &s)
{
	quad_inter <bigint> ic;

	// The singular point of the cone
	math_vector <bigint> cone_sing = column(QI::singular(qa), 0);

	#ifndef NDEBUG
	s << ">> optimization of coordinates of singular point" << endl;
	#endif

	optimize(cone_sing);

	if (img_cone == 3) // Imaginary cone
		{
			// Only one component here
			ic = quad_inter <bigint>(1);

			ic.set_type(13,1);

			#ifndef NDEBUG
			//print_type(ic,s);
			#endif 

			// The point is the singular point of the cone
			curve_param <bigint> cone_p(cone_sing);

			ic.cc[0].create_component(INTER_TYPE_POINT,cone_p);
			ic.cc[0].optiflag = true;
		}
	else
		{
			// Either two or three components
			if (s_e == 1)
				{
					ic = quad_inter <bigint>(3);
			
					ic.set_type(13,3);
			
					#ifndef NDEBUG
					//print_type(ic,s);
					#endif		
				}
			else
				{
					ic = quad_inter <bigint>(2);
			
					ic.set_type(13,2);
			
					#ifndef NDEBUG
					//print_type(ic,s);
					#endif		
				}

			// Preliminaries for both lines and conic

			// The singular locus of the pair of planes
			bigint_matrix plane_sing = QI::singular(qb);

			// Since the vertex of qb falls on qa outside of its singular locus, the two
			// planes of qb are rational. Let us parameterize them
			bigint D;
			bigint_matrix m1(4,3),m2(4,1);

			pair_of_planes_param(qb,plane_sing,D,m1,m2,opt_level,s);

			// The two planes are rational so we can directly replace the first column of m1
			// by m1(0) +/- m2
			bigint_matrix m1p = m1;
			column(m1p, 0) = column(m1, 0) + column(m2, 0);
			bigint_matrix m1m = m1;
			column(m1m, 0) = column(m1, 0) - column(m2, 0);

			// Now we compute the following two matrices: one has inertia [2 1]
			// and represents the conic, the other has inertia either [1 1] or [2 0] and
			// represents the (real or imaginary) pair of lines
			bigint_matrix ap = prod(base_matrix<bigint> (prod(trans(m1p), qa)), m1p);
			bigint_matrix am = prod(base_matrix<bigint> (prod(trans(m1m), qa)), m1m);

			// Swap if ap is not the conic of inertia [2 1]
			if (det(ap).is_zero())
				{
					swap(ap,am);
					swap(m1p,m1m);
				}

			// Cut parameters of the two lines - needed later on
			math_vector <bigint> cut1(2,2), cut2(2,2);

			bigint_matrix linep(4,2),linem(4,2);
			if (s_e == -1)
				{
					// First component is point
			
					// The point is the singular point of the cone
					curve_param <bigint> cone_p(cone_sing);
			
					ic.cc[0].create_component(INTER_TYPE_POINT,cone_p);
					ic.cc[0].optiflag = true;
			
					// Set D to a non-zero value so that it goes in the right loop below
					D = 1;
				}
			else // else s_e != -1
				{
					// First component is pair of lines
			
					// The (2D) singular point of the pair of lines (could be inferred from
					// the cone singularity, but...)
					math_vector <bigint> line_sing = column(QI::singular(am), 0);
			
					optimize(line_sing);
			
					// Let us parameterize the pair of lines
					bigint_matrix m1_l(3,2),m2_l(3,1);
			
					pair_of_lines_param(am,line_sing,D,m1_l,m2_l,opt_level,s);
			
					if (D.is_zero())
						{
							// The two lines are rational
							bigint_matrix l1p = m1_l;
							column(l1p, 0) = column(m1_l, 0) + column(m2_l, 0);
							bigint_matrix l1m = m1_l;
							column(l1m, 0) = column(m1_l, 0) - column(m2_l, 0);
					
							#ifndef NDEBUG
							s << ">> optimization of transformation" << endl;
							#endif
			
							optimize_trans1(l1p);
							optimize_trans1(l1m);
			
							// Stack the transformation matrices
							linep = prod(m1m, l1p);
							linem = prod(m1m, l1m);
			
							#ifndef NDEBUG
							s << ">> optimization of parameterization" << endl;
							#endif
			
							optimize_trans1(linep);
							optimize_trans1(linem);
			
							// The two lines
							curve_param <bigint> line1(column(linep, 1), column(linep, 0)),
																	 line2(column(linem, 1), column(linem, 0));
			
							// Reparameterize the lines
							#ifndef NDEBUG
							s << ">> reparameterization of lines" << endl;
							#endif
			
							cut1 = improved_line_param(line1);
							cut2 = improved_line_param(line2);
			
							optimize_by_half(cut1);
							optimize_by_half(cut2);
			
							ic.cc[0].create_component(INTER_TYPE_LINE,line1);
							ic.cc[1].create_component(INTER_TYPE_LINE,line2);
			
							#ifndef NDEBUG
							if (s_e == 1)
								s << ">> the lines cut at " << cone_sing << endl;
							#endif 
						}
		 			else
						{
							// The lines are not rational
			
							// Stack the transformations
							m1_l = prod(m1m, m1_l);
							m2_l = prod(m1m, m2_l);
							m2_l.resize(4,2);
			
							#ifndef NDEBUG
							s << ">> optimization of parameterization" << endl;
							#endif
			
							optimize_trans2(m1_l,m2_l);
			
							// The components of the lines
							curve_param <bigint> line_rat(column(m1_l, 0),column(m1_l, 1)),
																	 line_D(column(m2_l, 0),column(m2_l, 1));
			
							// Reparameterize the lines
							#ifndef NDEBUG
							s << ">> reparameterization of lines" << endl;
							#endif
			
							improved_line_param(line_rat,line_D,D);
			
							ic.create_two_components(0,INTER_TYPE_LINE,line_rat,line_D,D);
						}
			
					// Lines are optimal
					ic.cc[0].optiflag = true;
					ic.cc[1].optiflag = true;
				} // if s_e == -1 else
			
			// If the lines are rational, their intersection with the singular line of
			// the pair of planes is a rational point on the conic
				
			if (D.is_zero())
				{
					math_vector <bigint> pp = column(linear_intersect(plane_sing,linep), 0);
					math_vector <bigint> pm = column(linear_intersect(plane_sing,linem), 0);
			
					optimize(pp);
					optimize(pm);
			
					if (s_e == 1)
						{
							// Find the parameters on the lines corresponding to pp and pm
							bigint_matrix plinep(4,2),plinem(4,2);
						
							column(plinep, 0) = (ic.cc[0].c[0]).eval(1,0);
							column(plinep, 1) = (ic.cc[0].c[0]).eval(0,1);
							
							math_vector <bigint> cutp = prod(cofactors_mat(prod(trans(plinep),plinep)), 
							 																math_vector<bigint> (prod(trans(plinep), pp)));
							
							optimize(cutp);
							
							ic.cc[0].create_cut_parameters(-1,cut1,-1,cutp);
							
							column(plinem, 0) = (ic.cc[1].c[0]).eval(1,0);
							column(plinem, 1) = (ic.cc[1].c[0]).eval(0,1);
							
							math_vector <bigint> cutm = prod(cofactors_mat(prod(trans(plinem),plinem)),
							 																math_vector<bigint>(prod(trans(plinem),pm)));
							
							optimize(cutm);
							
							ic.cc[1].create_cut_parameters(-1,cut2,-1,cutm);
						}
			 		
					// The point pp is a*plane_sing(0)+b*plane_sing(1): find a and b
			 		math_vector <bigint> ab = prod(cofactors_mat(base_matrix<bigint> (prod(trans(plane_sing),plane_sing))),
			 		 												 math_vector<bigint> (prod(trans(plane_sing), pp)));
			 		
			 		optimize(ab);
			 		
			 		// (0,ab[0],ab[1]) is a rational point on the conic!
			 		math_vector <bigint> rat_point2d(3,3);
			 		rat_point2d[1] = ab[0];
			 		rat_point2d[2] = ab[1];
			 		
			 		bigint_matrix p_trans(3,3);
			 		curve_param <bigint> conic2d(3);
			 		math_vector <bigint> l(2,2);
			 		
			 		conic_param_through_ratpoint(ap,rat_point2d,p_trans,conic2d,l,s);
					
					// Stack the transformations
					p_trans = prod(m1p, p_trans);
					
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
					
					// Rescale so that l[0]*u+l[1]*v is replaced by v
					bigint_matrix scalep(3,3);
					
					scalep.sto(0,0,l[1]*l[1]);
					scalep.sto(1,0,l[0]*l[0]);
					scalep.sto(1,1,1);
					scalep.sto(1,2,-2*l[0]);
					scalep.sto(2,0,-l[0]*l[1]);
					scalep.sto(2,2,l[1]);
					
					p_trans = prod(p_trans, scalep);
					optimize_trans3(p_trans,opt_level,s);
					
					l[0] = 0;
					l[1] = -1;
					
					curve_param <bigint> conic(4);
					multiply(conic,p_trans,conic2d);
					
					if (s_e == -1)
						{
					 		ic.cc[1].create_component(INTER_TYPE_CONIC,conic);
					 		ic.cc[1].optiflag = true;
					 		ic.cc[1].create_cut_parameters(-1,-l[1],l[0],-1,0,1);
						}
					else
						{
					 		ic.cc[2].create_component(INTER_TYPE_CONIC,conic);			
					 		ic.cc[2].optiflag = true;
					 		ic.cc[2].create_cut_parameters(-1,-l[1],l[0],-1,0,1);
						}
					
					#ifndef NDEBUG
					s << ">> line 1 cuts the conic at " << pp << endl;
					s << ">> line 2 cuts the conic at " << pm << endl;
					#endif
				}
			else
				{
					// Well the lines are not rational. But we may still end up finding a
					// rational conic...
					
					// Parameterize the conic
					bigint D2;
					bigint_matrix p_trans(3,3),p_trans2(3,3);
					curve_param <bigint> par(3);
					
					conic_param(ap,D2,p_trans,p_trans2,par,opt_level,s);
					
					if (D2.is_zero())
						{
							// Lucky you: you just found a rational conic!!
						
							// Stack the transformations
							p_trans = prod(m1p, p_trans);
						
							#ifndef NDEBUG
							s << ">> optimization of parameterization" << endl;
						 	#endif
						
							optimize_trans3(p_trans,opt_level,s);
							 
							curve_param <bigint> conic(4);
							multiply(conic,p_trans,par);
						
							if (s_e == -1)
								{
									ic.cc[1].create_component(INTER_TYPE_CONIC,conic);
									ic.cc[1].optiflag = true;
								}
							else
								{
									ic.cc[2].create_component(INTER_TYPE_CONIC,conic);			 
									ic.cc[2].optiflag = true;
								}
						}
					else
						{
							// You'll have to content yourself with one square root...
							 
							// Stack the transformations
							p_trans = prod(m1p, p_trans);
							p_trans2 = prod(m1p, p_trans2);
			
							#ifndef NDEBUG
							s << ">> optimization of parameterization" << endl;
							#endif
			
							optimize_trans4(p_trans,p_trans2,opt_level,s);
						
							curve_param <bigint> conic(4),conic2(4);
							multiply(conic,p_trans,par);
							multiply(conic2,p_trans2,par);
							 
							if (s_e == -1)
								{
									ic.cc[1].create_component(INTER_TYPE_CONIC,conic,conic2,D2);
									if (opt_level)
										ic.cc[1].optiflag = true;
									else
										ic.cc[1].optiflag = false;
								}
							else
								{
									ic.cc[2].create_component(INTER_TYPE_CONIC,conic,conic2,D2);			 
									if (opt_level)
										ic.cc[2].optiflag = true;
									else
										ic.cc[2].optiflag = false;
								}
						}
	
					if (s_e != -1)
						{
							// Temp
							ic.cc[1].create_cut_parameters(-1,0,0,-1,0,0);
							ic.cc[2].create_cut_parameters(-1,0,0,-1,0,0);
						}
				}
		}

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// Procedure for the cubic and secant line case when the 2 roots are rational
quad_inter <bigint> inter_cubic_secant_two_rat(const bigint_matrix &q1, 
																								const bigint_matrix &q2,
																								const bigint_matrix &qa, 
																								const bigint_matrix &qb, 
																								const int opt_level, ostream &s)
{
	// Two components here
	quad_inter <bigint> ic(2);

	ic.set_type(12,1);

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif
	
	// Both qa and qb are rational cones. The singular point of qa is a rational
	// point of qb

	// Singular points of cones
	math_vector <bigint> singa = column(QI::singular(qa), 0);
	math_vector <bigint> singb = column(QI::singular(qb), 0);

	#ifndef NDEBUG
	s << ">> optimization of coordinates of singular and rational points" << endl;
	#endif

	optimize(singa);
	optimize(singb);

	math_vector <bigint> rat_point = singa;

	// The line param
	curve_param <bigint> line_par(singa,singb);

	// The cut parameters
	math_vector <bigint> cut(4,4);

	// Try to find a point of small height on the line by looking at the
	// parameter values such that one of the components is zero

	#ifndef NDEBUG
	s << ">> reparameterization of line" << endl;
	#endif

	cut = improved_line_param(line_par);

	optimize_by_half(cut);

	if (are_equal(singb,line_par.eval(1,0)))
		rat_point = line_par.eval(0,1);
	else
		rat_point = line_par.eval(1,0);

	#ifndef NDEBUG
	s << ">> singular point of cone: " << singb << endl;
	s << ">> rational point on cone: " << rat_point << endl;
	#endif

	// Parameterize the cone
	surface_param <bigint> par(4);
	bigint_matrix p_trans(4,4);
	math_vector <bigint> l(2,2);
	cone_param_through_ratpoint(qb,singb,rat_point,p_trans,par,l,s);
	
	#ifndef NDEBUG
	s << ">> optimization of transformation" << endl;
	#endif

	optimize_trans3(p_trans,l,opt_level,s);

	// Do that anyway, otherwise tmp might not divide res (constant factor larger)...
	optimize(l);

	// The quadric in the canonical frame of the cone
	bigint_matrix q_tmp = prod(base_matrix<bigint> (prod(trans(p_trans), qa)) , p_trans);

	// Plug in the other quadric
	hom_hom_polynomial <bigint> res = plug_param_in_quadric(par,q_tmp,par);

	// The common factor of res (see below)
	hom_polynomial <bigint> tmp;

	math_vector <bigint> m(2,2);

	if (!l[0].is_zero())
		{
			m[0] = l[0]*res[1][2];
			m[1] = l[0]*res[1][1]-l[1]*res[1][2];
		}
	else
		{
			m[0] = res[1][1];
			m[1] = res[1][0];
		}

	optimize(m);

	// res has the form:
	//			 ---> u*(l[0]*s+l[1]*t)*(c1*p1(s,t)*u + c2*p2(s,t)*v),
	// where p1(s,t) has degree 3 and p2(s,t) = m[0]*s+m[1]*t has degree 1
	
	// The two solutions (-l[1],l[0]) and (-m[1],m[0]) correspond to the two points of
	// intersection of the line with the cubic. Let us rescale the cubic so that
	// these two points correspond to the parameters (1,0) and (0,1). Since these
	// points are quite ``simple'', the cubic should have a much nicer param.
	
	// Now rewrite the par accordingly, i.e. replace l[0]*s+l[1]*t by t and 
	// m[0]*s+m[1]*t by s
	par[0][1][2] = l[1]*l[1];
	par[0][1][1] = -2*m[1]*l[1];
	par[0][1][0] = m[1]*m[1];
	par[1][1][2] = l[0]*l[0];
	par[1][1][1] = -2*m[0]*l[0];
	par[1][1][0] = m[0]*m[0];
	par[2][1][2] = -l[0]*l[1];
	par[2][1][1] = l[0]*m[1]+l[1]*m[0];
	par[2][1][0] = -m[0]*m[1];
	
	// And do the plug again
	res = plug_param_in_quadric(par,q_tmp,par);

	// At this stage, res has the form
	//			 ---> u*t*(c1*p1(s,t)*t*u + c2*s*v)
	// where p1(s,t) has degree 2
	
	tmp.assign_y();
	
	l[0] = 0;
	l[1] = -1;
	
	m[0] = 1;
	m[1] = 0;

	hom_polynomial <bigint> p1,p2;
	divide(p1,res[2],tmp);
	divide(p2,res[1],tmp);
	negate(p2,p2);

	// Expand param
	multiply(par,p_trans,par);

	curve_param <bigint> cubic(4);
	cubic = par.eval(p2,p1);

	#ifndef NDEBUG
	s << ">> optimization of parameterization" << endl;
	#endif

	optimize(cubic);

	#ifndef NDEBUG
	s << ">> cubic and line intersect at " << singa << " and " << singb << endl;
	#endif

	ic.cc[0].create_component(INTER_TYPE_CUBIC,cubic);
	ic.cc[0].create_cut_parameters(-1,-l[1],l[0],-1,-m[1],m[0]);

	ic.cc[1].create_component(INTER_TYPE_LINE,line_par);
	ic.cc[1].create_cut_parameters(-1,cut[0],cut[1],-1,cut[2],cut[3]);

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// Procedure for the four skew lines case when the roots are rational
quad_inter <bigint> inter_four_skew_lines_rat(const bigint_matrix &q1, 
																							const bigint_matrix &q2,
																							const bigint_matrix &qa, 
																							const bigint_matrix &qb, 
																							const int opt_level, ostream &s)
{
	// Four components here
	quad_inter <bigint> ic(4);

	ic.set_type(14,2);

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// The four vertices of the quadrilateral are the intersections of the singular
	// line of one pair with the other pair
	bigint_matrix singa = QI::singular(qa);
	bigint_matrix singb = QI::singular(qb);

	// Now compute the following matrices
	bigint_matrix mab = prod(base_matrix<bigint> (prod(trans(singa), qb)) ,singa);
	bigint_matrix mba = prod(base_matrix<bigint> (prod(trans(singb), qa)) ,singb);

	// The discriminants
	bigint Dab = mab(0,1)*mab(0,1)-mab(0,0)*mab(1,1);
	bigint Dba = mba(0,1)*mba(0,1)-mba(0,0)*mba(1,1);	 

	// The (1D) points: the solutions are pab1 and pab2 if the discrim is a square,
	// pab1 +/- sqrt(Dab)*pab2 otherwise
	bigint_matrix pab1(2,1),pab2(2,1),pba1(2,1),pba2(2,1);

	// The square roots of discriminants
	bigint sqrt_Dab,sqrt_Dba;

	bool sq_Dab = is_square(sqrt_Dab,Dab);
	if (sq_Dab)
		{
			if (!mab(0,0).is_zero())
				{
					 pab1.sto(0,0,-mab(0,1)+sqrt_Dab);
					 pab1.sto(1,0,mab(0,0));
					 pab2.sto(0,0,-mab(0,1)-sqrt_Dab);
					 pab2.sto(1,0,mab(0,0));
				}
			else
				{
					 pab1.sto(0,0,1);
					 pab2.sto(0,0,mab(1,1));
					 pab2.sto(1,0,-2*mab(0,1));
				}
						
			pab1 = prod(singa, pab1);
			pab2 = prod(singa, pab2);

			#ifndef NDEBUG
			s << ">> optimization of first set of points" << endl;
			#endif

			optimize_trans1(pab1);
			optimize_trans1(pab2);
		}
	else
		{
			pab1.sto(0,0,-mab(0,1));
			pab1.sto(1,0,mab(0,0));
			pab2.sto(0,0,1);

			pab1 = prod(singa, pab1);
			pab2 = prod(singa, pab2);

			// Replace Dab by e^2*Dab' - and optimize
			#ifndef NDEBUG
			extract_message(opt_level,s,"optimization of square root 1");
			#endif

			math_vector <bigint> D_fact = extract_square_factors(Dab,opt_level,s);
		 
			// Dab --> Dab'
			Dab = D_fact[1];
		 
			// Multiply pab2 by e = D_fact[0]
			pab2 = D_fact[0]*pab2;
			
			#ifndef NDEBUG
			s << ">> optimization of first set of points" << endl;
			#endif

			optimize_trans2(pab1,pab2);
		}

	bool sq_Dba = is_square(sqrt_Dba,Dba);

	if (sq_Dba)
		{
			if (!mba(0,0).is_zero())
				{
					 pba1.sto(0,0,-mba(0,1)+sqrt_Dba);
					 pba1.sto(1,0,mba(0,0));
					 pba2.sto(0,0,-mba(0,1)-sqrt_Dba);
					 pba2.sto(1,0,mba(0,0));
				}	
						else
					{
					 pba1.sto(0,0,1);
					 pba2.sto(0,0,mba(1,1));
					 pba2.sto(1,0,-2*mba(0,1));
				}

			pba1 = prod(singb,pba1);
			pba2 = prod(singb,pba2);

			#ifndef NDEBUG
			s << ">> optimization of second set of points" << endl;
			#endif

			optimize_trans1(pba1);
			optimize_trans1(pba2);
		}
	else
		{
			pba1.sto(0,0,-mba(0,1));
			pba1.sto(1,0,mba(0,0));
			pba2.sto(0,0,1);

			pba1 = prod(singb, pba1);
			pba2 = prod(singb, pba2);

			// Replace Dba by e^2*Dba' - and optimize
			#ifndef NDEBUG
			extract_message(opt_level,s,"optimization of square root 2");
			#endif

			math_vector <bigint> D_fact = extract_square_factors(Dba,opt_level,s);
		 
			// Dba --> Dba'
			Dba = D_fact[1];

			// Multiply pba2 by e = D_fact[0]
			pba2 = D_fact[0]*pba2;

			#ifndef NDEBUG
			s << ">> optimization of second set of points" << endl;
			#endif

			optimize_trans2(pba1,pba2);
		}

	// Swap if needed
	if ((sq_Dba) && (!sq_Dab))
		{
			swap(pba1,pab1);
			swap(pba2,pab2);
			swap(Dab,Dba);
			bool tmp = sq_Dba;
			sq_Dba = sq_Dab;
			sq_Dab = tmp;
		}

	if ((sq_Dab) && (sq_Dba))
		{
			// The four lines are rational
			curve_param <bigint> line1(column(pab1, 0), column(pba1, 0)),
													 line2(column(pab1, 0), column(pba2, 0));
			curve_param <bigint> line3(column(pab2, 0), column(pba1, 0)),
													 line4(column(pab2, 0), column(pba2, 0));

			// Cut parameters of the lines 
			math_vector <bigint> cut1(4,4), cut2(4,4), cut3(4,4), cut4(4,4);

			// Reparameterize the lines

			#ifndef NDEBUG
			s << ">> reparameterization of lines" << endl;
			#endif

			cut1 = improved_line_param(line1);
			cut2 = improved_line_param(line2);
			cut3 = improved_line_param(line3);
			cut4 = improved_line_param(line4);

			optimize_by_half(cut1);
			optimize_by_half(cut2);
			optimize_by_half(cut3);
			optimize_by_half(cut4);

			ic.cc[0].create_component(INTER_TYPE_LINE,line1);
			ic.cc[0].create_cut_parameters(-1,-1,cut1);

			ic.cc[1].create_component(INTER_TYPE_LINE,line2);
			ic.cc[1].create_cut_parameters(-1,-1,cut2);

			ic.cc[2].create_component(INTER_TYPE_LINE,line3);
			ic.cc[2].create_cut_parameters(-1,-1,cut3);

			ic.cc[3].create_component(INTER_TYPE_LINE,line4);
			ic.cc[3].create_cut_parameters(-1,-1,cut4);
		}
	else if (sq_Dab)
		{
			// Two points are rational, the two others are not
			math_vector <bigint> zero(4,4);
			curve_param <bigint> line1_rat(column(pab1, 0),column(pba1, 0)),
													 line2_rat(column(pab2, 0),column(pba1, 0));
			curve_param <bigint> line_sq1(zero,column(pba2, 0)),line_sq2 = line_sq1;

			math_vector <bigint> cut1(6,6), cut2(6,6);

			// Reparameterize the lines
			#ifndef NDEBUG
			s << ">> reparameterization of lines" << endl;
			#endif
		 
			cut1 = improved_line_param(line1_rat,line_sq1,Dba);
			cut2 = improved_line_param(line2_rat,line_sq2,Dba);

			optimize_by_half(cut1);
			optimize_by_half(cut2);

			ic.create_two_components(0,INTER_TYPE_LINE,line1_rat,line_sq1,Dba);
			ic.cc[0].create_cut_parameters(-1,-1,cut1,Dba);
			ic.cc[1].create_cut_parameters(-1,cut1[0],-cut1[1],cut1[2],
						 -1,cut1[3],-cut1[4],cut1[5],Dba);

			ic.create_two_components(2,INTER_TYPE_LINE,line2_rat,line_sq2,Dba);
			ic.cc[2].create_cut_parameters(-1,-1,cut2,Dba);
			ic.cc[3].create_cut_parameters(-1,cut2[0],-cut2[1],cut2[2],
						 -1,cut2[3],-cut2[4],cut2[5],Dba);
		}
	else
		{
			// None of the points are rational, but the degree might still be 2 if 
			// sqrt(Dab) and sqrt(Dba) are in the same extension

			bigint g = gcd(Dab,Dba);
			bigint sq_ab,sq_ba;
			
			if (is_square(sq_ab,Dab/g) && is_square(sq_ba,Dba/g))
				{
					 // Degree 2
					 if (Dab != Dba)
						 { 
							 pba1 = sq_ab*pba1;
							 pba2 = sq_ba*pba2;
			
							 optimize_trans2(pba1,pba2);
						 }
			
					 curve_param <bigint> line_rat1(column(pab1, 0),column(pba1, 0));
					 curve_param <bigint> line_rat2 = line_rat1;
					 curve_param <bigint> line1_sq(column(pab2, 0),column(pba2, 0));
					 curve_param <bigint> line2_sq(column(pab2, 0),column(-pba2, 0));
			
					 math_vector <bigint> cut1(6,6), cut2(6,6);
			
					 // Reparameterize the lines
					 #ifndef NDEBUG
					 s << ">> reparameterization of lines" << endl;
					 #endif
			
					 cut1 = improved_line_param(line_rat1,line1_sq,Dab);
					 cut2 = improved_line_param(line_rat2,line2_sq,Dab);
			
					 optimize_by_half(cut1);
					 optimize_by_half(cut2);
			
					 ic.create_two_components(0,INTER_TYPE_LINE,line_rat1,line1_sq,Dab);
					 ic.cc[0].create_cut_parameters(-1,-1,cut1,Dba);
					 ic.cc[1].create_cut_parameters(-1,cut1[0],-cut1[1],cut1[2],
								 -1,cut1[3],-cut1[4],cut1[5],Dba);
			
					 ic.create_two_components(2,INTER_TYPE_LINE,line_rat2,line2_sq,Dab);
					 ic.cc[2].create_cut_parameters(-1,-1,cut2,Dab);
					 ic.cc[3].create_cut_parameters(-1,cut2[0],-cut2[1],cut2[2],
								 -1,cut2[3],-cut2[4],cut2[5],Dab);
				}
			else
				{
					 // Swap? (just for presentation)
					 if (Dab > Dba)
						 {
							 swap(pba1,pab1);
							 swap(pba2,pab2);
							 swap(Dab,Dba);
						 }
			
					 math_vector <bigint> zero(4,4);
					 curve_param <bigint> line_rat(column(pab1, 0),column(pba1, 0));
					 curve_param <bigint> line1_sq(column(pab2, 0),zero),line2_sq(zero,column(pba2, 0));
					 curve_param <bigint> zer(zero);
			
								// Reparameterize the lines
					 #ifndef NDEBUG
					 s << ">> reparameterization of lines" << endl;
					 #endif
			
					 improved_line_param(line_rat,line1_sq,line2_sq,zer,Dab,Dba,0);
			
					 // Warning: zer might have changed!!!
					 ic.create_four_components(INTER_TYPE_LINE,line_rat,line1_sq,line2_sq,zer,Dab,Dba);
			
					 // Temp
					 ic.cc[0].create_cut_parameters(-1,0,0,-1,0,0);		
					 ic.cc[1].create_cut_parameters(-1,0,0,-1,0,0);		
					 ic.cc[2].create_cut_parameters(-1,0,0,-1,0,0);		
					 ic.cc[3].create_cut_parameters(-1,0,0,-1,0,0);		
				}
		}

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// Procedure for the four skew lines case when the real part is 2 points
quad_inter <bigint> inter_four_skew_lines_2pts_rat(const bigint_matrix &q1, 
																									 const bigint_matrix &q2,
																									 const bigint_matrix &qa_cp, 
																									 const bigint_matrix &qb_cp, 
																									 const int in_qa, const int opt_level, 
																									 ostream &s)
{
	// Two components here
	quad_inter <bigint> ic(2);

	ic.set_type(14,3);

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// qa_cp and qb_cp are const: make copies
	bigint_matrix qa = qa_cp, qb = qb_cp;

	// If qb is the imaginary pair of planes, swap
	if (in_qa == 1)
		swap(qa,qb);

	// Now qa is the imaginary one
	bigint_matrix singa = QI::singular(qa);

	// Find the two points by intersecting with qb
	bigint_matrix mab = prod(base_matrix<bigint> (prod(trans(singa), qb)), singa);
	
	// The discriminant
	bigint Dab = mab(0,1)*mab(0,1)-mab(0,0)*mab(1,1);

	// The (1D) points: the solutions are pab1 and pab2 if the discrim is a square,
	// pab1 +/- sqrt(Dab)*pab2 otherwise
	bigint_matrix pab1(2,1),pab2(2,1);

	// The discriminants
	bigint sqrt_Dab;

	if (is_square(sqrt_Dab,Dab))
		{
			if (!mab(0,0).is_zero())
				{
					 pab1.sto(0,0,-mab(0,1)+sqrt_Dab);
					 pab1.sto(1,0,mab(0,0));
					 pab2.sto(0,0,-mab(0,1)-sqrt_Dab);
					 pab2.sto(1,0,mab(0,0));
				}
			else
				{
					 pab1.sto(0,0,1);
					 pab2.sto(0,0,mab(1,1));
					 pab2.sto(1,0,-2*mab(0,1));
				}

			pab1 = prod(singa,pab1);
			pab2 = prod(singa,pab2);

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize_trans1(pab1);
			optimize_trans1(pab2);
			
			// The two points
			curve_param <bigint> point1(column(pab1, 0)),point2(column(pab2, 0));

			ic.cc[0].create_component(INTER_TYPE_POINT,point1);
			ic.cc[1].create_component(INTER_TYPE_POINT,point2);
		}
	else
		{
			sqrt_Dab = 0;
			pab1.sto(0,0,-mab(0,1));
			pab1.sto(1,0,mab(0,0));
			pab2.sto(0,0,1);

			pab1 = prod(singa,pab1);
			pab2 = prod(singa,pab2);

			// Replace Dab by e^2*Dab' - and optimize
			#ifndef NDEBUG
			extract_message(opt_level,s,"optimization of square root");
			#endif

			math_vector <bigint> D_fact = extract_square_factors(Dab,opt_level,s);
		 
			// Dab --> Dab'
			Dab = D_fact[1];

			// Multiply pab2 by e = D_fact[0]
			pab2 = D_fact[0]*pab2;

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize_trans2(pab1,pab2);

			curve_param <bigint> point_rat(column(pab1, 0)),point_sq(column(pab2,0));

			ic.create_two_components(0,INTER_TYPE_POINT,point_rat,point_sq,Dab);
		}

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// Procedure for the four skew lines case when the real part is four lines or two
// points and the roots are non-rational
quad_inter <bigint> inter_four_skew_lines_non_rat(const bigint_matrix &q1, 
																									const bigint_matrix &q2, 
																									const bigint_matrix &sing, 
																									const bigint &d, 
																									const rpl_size_t &s_e,
																									const int opt_level, ostream &s)
{
	quad_inter <bigint> ic;

	// sing has four columns: among them, at least two columns are independent
	// (though the others might be linear combinations of the two points).
	// Also, if the k1i + sqrt(d)*k2i span the first line, then the k1i -
	// sqrt(d)*k2i span the second! 
	bigint_matrix k11(4,1),k21(4,1),k12(4,1),k22(4,1);
	bigint_matrix tmp1(8,2),tmp2(8,2);

	sing.split_h(tmp1,tmp2);
	tmp1.split_t(k11,k12,k21,k22);

	// Pick another column if the first two are in the same span
	if (are_equal(k11,k21,k12,k22,d))
		{
			 bigint_matrix tmp(4,1);
			 tmp2.split_t(k12,tmp,k22,tmp);
		}

	// The result of plugging u*(k11+eps1*sqrt(d)*k21)+v*(k12+eps1*sqrt(d)*k22) in
	// q1 is of the form (u	 v)*(p_rat + eps1*sqrt(d)*p_sq)*(u	v)^T
	bigint_matrix p_rat(2,2),p_sq(2,2);
	p_rat.sto(0,0,(prod(base_matrix<bigint> (prod(trans(k11),q1)), k11) + 
								 prod(base_matrix<bigint> (prod(d*trans(k21),q1)), k21))(0,0));
	p_rat.sto(0,1,(prod(base_matrix<bigint> (prod(trans(k11),q1)) , k12) + 
								 prod(base_matrix<bigint> (prod(d*trans(k21), q1)), k22))(0,0));
	p_rat.sto(1,0,p_rat(0,1));
	p_rat.sto(1,1,(prod(base_matrix<bigint> (prod(trans(k12), q1)), k12) + 
								 prod(base_matrix<bigint> (prod(d*trans(k22), q1)), k22))(0,0));

	p_sq.sto(0,0,bigint(2) * prod(base_matrix<bigint> (prod(trans(k21), q1)), k11)(0,0));
	p_sq.sto(0,1,(prod(base_matrix<bigint> (prod(trans(k21), q1)), k12) + 
								prod(base_matrix<bigint> (prod(trans(k11), q1)), k22))(0,0));
	p_sq.sto(1,0,p_sq(0, 1));
	p_sq.sto(1,1,bigint(2) * prod(base_matrix<bigint> (prod(trans(k12), q1)), k22)(0,0));

	// The (reduced) discriminants of the quadratic equation are e^2*(D1 +/- sqrt(d)*D2)
	bigint e = 1;
	bigint D1 = p_rat(0,1)*p_rat(0,1)+d*p_sq(0,1)*p_sq(0,1)-p_rat(0,0)*p_rat(1,1)
							-d*p_sq(0,0)*p_sq(1,1);
	bigint D2 = 2*p_rat(0,1)*p_sq(0,1)-p_sq(0,0)*p_rat(1,1)-p_sq(1,1)*p_rat(0,0);

	#ifndef NDEBUG
	extract_message(opt_level,s,"optimization of square root");
	#endif

	math_vector <bigint> D_fact1 = extract_square_factors(D1,opt_level,s);
	math_vector <bigint> D_fact2 = extract_square_factors(D2,opt_level,s);

	e = gcd(D_fact1[0],D_fact2[0]);

	bigint e2 = e*e;

	D1 = D1/e2;
	D2 = D2/e2;

	// The points in 3D are:
	//				po_rat + eps1*sqrt(d)*po_sq + eps2*(k11+eps1*sqrt(d)*k21)*e*sqrt(Delta), 
	// where Delta = D1 + eps1*sqrt(d)*D2
	bigint_matrix po_rat = -p_rat(0,1)*k11+p_rat(0,0)*k12-d*p_sq(0,1)*k21+d*p_sq(0,0)*k22;
	bigint_matrix po_sq = -p_rat(0,1)*k21-p_sq(0,1)*k11+p_sq(0,0)*k12+p_rat(0,0)*k22;

	// Warning: the case when p_rat(0,0) = p_sq(0,0) = 0 can happen. In that case,
	// the solution is not (-b+eps2*sqrt(Delta),2*a) but (2*c,-b+eps2*sqrt(Delta)).
	// We handle it later since it implies that the solution contains rational lines

	// Can we simplify that a bit?
	bigint ct = gcd(content(column(po_rat, 0)),content(column(po_sq, 0)));
	ct = gcd(ct,e);
			
	e = e/ct;
	po_rat = po_rat/ct;
	po_sq = po_sq/ct;

	if (s_e != 1) // Two points
		{
			// Two components here
			ic = quad_inter <bigint>(2);

			ic.set_type(14,3);

			#ifndef NDEBUG
			//print_type(ic,s);
			#endif

			// One of the singular quadrics is imaginary. Which one? (hint: the one with D1
			// +/- d*D2 < 0)
			bigint eps1;
			if (sign(D1,D2,d) == 1) // D1+sqrt(d)*D2 > 0
				eps1 = 1;
			else
				eps1 = -1;

			// The points correspond to the value of eps such that Delta > 0
			curve_param <bigint> c1(column(po_rat, 0)), c2(column(eps1*po_sq, 0));
			curve_param <bigint> c3(column(k11 * e, 0)), c4(column(eps1 * e * k21, 0));

			ic.create_two_components(0,INTER_TYPE_POINT,c1,c2,c3,c4,d,D1,eps1*D2);
		}
	else
		{
			// Four components here
			ic = quad_inter <bigint>(4);

			ic.set_type(14,2);

			#ifndef NDEBUG
			//print_type(ic,s);
			#endif

			// The lines link the two points with eps1 = 1 to the two points with eps1 = -1

			// Notations between here and paper:
			//										d	 <-->	 Delta
			//	 del2 = D1^2-d*D2^2	 <-->	 c
			
			// There are two cases according to the value of d*(D1^2-d*D2^2)
			bigint del2 = D1*D1-d*D2*D2;
			bigint sq;

			if (is_square(sq,d*del2))
				{
					// The four lines are defined on the same extension
					curve_param <bigint> c1p(column(po_rat, 0),column(-e*sq*k21, 0));
					curve_param <bigint> c1m = c1p;
					curve_param <bigint> c2p(column(po_sq, 0),column(e * sq * k11 / d, 0));
					curve_param <bigint> c2m = c2p;
					curve_param <bigint> c3p(column(e*k11, 0), column(po_rat, 0));
					curve_param <bigint> c4p(column(e*k21, 0), column(-po_sq, 0));
					curve_param <bigint> c3m(column(e*k11, 0), column(-po_rat, 0));
					curve_param <bigint> c4m(column(e*k21, 0), column(po_sq, 0));
			
					// Reparameterize the lines
					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif
			
					improved_line_param(c1p,c2p,c3p,c4p,d,D1,D2);
					improved_line_param(c1m,c2m,c3m,c4m,d,D1,D2);
			
					ic.create_two_components(0,INTER_TYPE_LINE,c1p,c2p,c3p,c4p,d,D1,D2);
					ic.create_two_components(2,INTER_TYPE_LINE,c1m,c2m,c3m,c4m,d,D1,D2);
				}
			else
				{
					// Two cases according to the value of D1^2-d*D2^2
			
					if (is_square(sq,del2))
						{
							// D2 can be equal to zero here. If sqrt(D1) and sqrt(d) are in the
							// same extension, only two lines are non-rational and we proceed in
							// the next loop. Otherwise, the four lines are non-rational and we
							// go in the second loop.
			
							bigint dsq,cd,cd1,g;
			
							if (D2 == 0)
								g = gcd(d,D1);
			
							if (((D2 == 0) && (is_square(cd,D1/g)) && (is_square(cd1,d/g))) ||
									((D2 != 0) && ((is_square(dsq,2*(D1+sq)) || is_square(dsq,2*(D1-sq))))))
								{
									// Two lines are rational and the two others have one square
									// root
						
									bigint_matrix r1p(4,1),r1m(4,1),s1p(4,1),s1m(4,1);
						
									// We handle here the case we have left aside earlier
									if ((p_rat(0,0) == 0) && (p_sq(0,0) == 0))
										{
											bigint b1 = p_rat(0,1);
											bigint b2 = p_sq(0,1);
											bigint c1 = p_rat(1,1);
											bigint c2 = p_sq(1,1);
						
											r1p = k11;
											s1p = k21;
											r1m = -c1*k11-d*c2*k21+2*b1*k12+2*d*b2*k22;
											s1m = -c2*k11-c1*k21+2*b1*k22+2*b2*k12;
										}
									else
										if (dsq != 0)
											{
												bigint co1 = 2*dsq;
												bigint co2 = e*dsq*dsq;
												bigint co3 = 2*e*D2;
									
												bigint_matrix m1r = co1*po_rat;
												bigint_matrix m2r = co2*k11+d*co3*k21;
									
												r1p = m1r+m2r;
												r1m = m1r-m2r;
									
												bigint_matrix m1s = co1*po_sq;
												bigint_matrix m2s = co2*k21+co3*k11;

												s1p = m1s+m2s;
												s1m = m1s-m2s;
											}
										else
											{
												bigint co = e*cd;
														
												bigint_matrix m1r = cd1*po_rat;
												bigint_matrix m2r = d*co*k21;
									
												r1p = m1r+m2r;
												r1m = m1r-m2r;
									
												bigint_matrix m1s = cd1*po_sq;
												bigint_matrix m2s = co*k11;
									
												s1p = m1s+m2s;
												s1m = m1s-m2s;
											}
									
									// First the two rational lines
									curve_param <bigint> line1(column(r1p, 0), column(s1p, 0)),
																			 line2(column(r1m, 0), column(s1m, 0));
									
									#ifndef NDEBUG
									s << ">> optimization of parameterization" << endl;
									#endif
									
									optimize(line1);
									optimize(line2);
									
									// Reparameterize the lines

									#ifndef NDEBUG
									s << ">> reparameterization of lines" << endl;
									#endif
							
									// Cut parameter of the lines ?????
									math_vector <bigint> cut1(4,4), cut2(4,4);
							
									cut1 = improved_line_param(line1);
									cut2 = improved_line_param(line2);
							
									ic.cc[0].create_component(INTER_TYPE_LINE,line1);
									ic.cc[1].create_component(INTER_TYPE_LINE,line2);			 
							
									// Now the two non-rational lines
									curve_param <bigint> c_r(column(r1p, 0),column(r1m, 0)),
																			 c_s(column(s1p, 0),column(-s1m, 0));
							
									#ifndef NDEBUG
									s << ">> optimization of parameterization" << endl;
							 		#endif
							
									optimize(c_r,c_s);
							
									// Reparameterize the lines
							 		#ifndef NDEBUG
									s << ">> reparameterization of lines" << endl;
							 		#endif
							
									improved_line_param(c_r,c_s,d);
							
									ic.create_two_components(2,INTER_TYPE_LINE,c_r,c_s,d);
								}
							else
								{
									// D2 might be zero here. In this case, D1 = sq and dm = em = 0
									// so special care has to be taken
						
									// Need one square root per line
									bigint tmp1 = 2*D1;
									bigint tmp2 = 2*sq;
					
									bigint dp = tmp1+tmp2;
									bigint ep = 1;
									bigint dm = tmp1-tmp2;
									bigint em = 1;
						
									// The discriminants of the planes are 
									// sqrt(Delta) = 1/2*(ep*sqrt(dp) +/- em*sqrt(dm))
						
									// We have: ep*em*sqrt(dp)*sqrt(dm) = 2*D2*sqrt(d) (if D2 != 0)
						
									#ifndef NDEBUG
									extract_message(opt_level,s,"optimization of second square root");
									#endif
						
									math_vector <bigint> D_fact = extract_square_factors(dp,opt_level,s);
									ep = D_fact[0];
									dp = D_fact[1];
						
									if (D2 == 0)
										dm = d*dp;
						
									D_fact = extract_square_factors(dm,opt_level,s);
									em = D_fact[0];
									dm = D_fact[1];
						
									bigint_matrix r1p,r1m,r2p,r2m,r3p,r3m,r4p,r4m;
						
									if (D2 != 0)
										{
											bigint_matrix k = 2*D2*k11;
						
											bigint_matrix tmp_v1 = e*(em*em*dm*k21+k);
											bigint_matrix tmp_v2 = e*(ep*ep*dp*k21+k);
						
											r1p = 4*D2*po_rat;
											r2p = tmp_v2;
											r3p = ep*tmp_v1;
											r4p = 2*ep*po_sq;
						
											r1m = r1p;
											r2m = tmp_v1;
											r3m = em*tmp_v2;
											r4m = 2*em*po_sq;
										}
									else
										{
											r1p = (bigint)2*po_rat;
											r2p = (bigint)2*po_sq;
											r3p = e*ep*k11;
											r4p = e*ep*k21;
						
											r1m = r1p;
											r2m = 2*d*po_sq;
											r3m = em*e*ep*k21;
											r4m = em*e*ep*k11;
										}
						
									#ifndef NDEBUG
									s << ">> optimization of parameterization" << endl;
									#endif
						
									optimize_trans2(r1p,r3p);
									optimize_trans2(r2p,r4p);
						
									optimize_trans2(r1m,r3m);
									optimize_trans2(r2m,r4m);
						
									curve_param <bigint> c_ratp(column(r1p, 0), column(r2p, 0)),
																			 c_sqp(column(r3p, 0), column(r4p, 0));
									curve_param <bigint> c_ratm(column(r1m, 0),column(r2m, 0)),
																			 c_sqm(column(r3m, 0),column(r4m, 0));
						
									// Reparameterize the lines
						 			#ifndef NDEBUG
									s << ">> reparameterization of lines" << endl;
						 			#endif
						
									improved_line_param(c_ratp,c_sqp,dp);
									improved_line_param(c_ratm,c_sqm,dm);
						
									ic.create_two_components(0,INTER_TYPE_LINE,c_ratp,c_sqp,dp);
									ic.create_two_components(2,INTER_TYPE_LINE,c_ratm,c_sqm,dm);
								}
						}
					else
						{
							bigint f = 1, nd1 = 2*D1, nd2 = 2, g = 1;
			
							#ifndef NDEBUG
							extract_message(opt_level,s,"optimization of square root");
							#endif
			
							math_vector <bigint> D_fact = extract_square_factors(del2,opt_level,s);
			
							g = D_fact[0];
							del2 = D_fact[1];
							nd2 = nd2*g;
			
							bigint gc = gcd(nd1,nd2);
			
							D_fact = extract_square_factors(gc,opt_level,s);
			
							f = D_fact[0];
						
							bigint f2 = f*f;
							nd1 = nd1/f2;
							nd2 = nd2/f2;
			
							math_vector <bigint> zero(4,4);
			
							bigint ef = e*f;
							bigint efg = ef*g;
			
							math_vector <bigint> k1p = column(2 * D2 * po_rat, 0);
							math_vector <bigint> k3p = column(ef * (D2 * k11 + D1 * k21), 0);
							math_vector <bigint> k4p = column(-efg * k21, 0);
			
							math_vector <bigint> k1m = column(2 * d * D2 * po_sq, 0);
							math_vector <bigint> k3m = column(ef * (d * D2 * k21 + D1 * k11), 0);
							math_vector <bigint> k4m = column(-efg * k11, 0);
			
							#ifndef NDEBUG
							s << ">> optimization of parameterization" << endl;
							#endif
										
							optimize(k1p,k3p,k4p);
							optimize(k1m,k3m,k4m);			 
			
							curve_param <bigint> c1p(k1p,k1m); 
							curve_param <bigint> c1m = c1p; 
							curve_param <bigint> c2p(zero,zero);							 
							curve_param <bigint> c2m = c2p; 
							curve_param <bigint> c3p(k3p,k3m); 
							curve_param <bigint> c3m = c3p;													
							curve_param <bigint> c4p(k4p,k4m);												
							curve_param <bigint> c4m(-k4p,-k4m); 
			
							// Reparameterize the lines
							#ifndef NDEBUG
							s << ">> reparameterization of lines" << endl;
							#endif
			
							improved_line_param(c1p,c2p,c3p,c4p,del2,nd1,nd2);
							improved_line_param(c1m,c2m,c3m,c4m,del2,nd1,-nd2);
			
							ic.create_two_components(0,INTER_TYPE_LINE,c1p,c2p,c3p,c4p,del2,nd1,nd2);
							ic.create_two_components(2,INTER_TYPE_LINE,c1m,c2m,c3m,c4m,del2,nd1,-nd2);
						}
				}

			// Temp
			ic.cc[0].create_cut_parameters(-1,0,0,-1,0,0);
			ic.cc[1].create_cut_parameters(-1,0,0,-1,0,0);
			ic.cc[2].create_cut_parameters(-1,0,0,-1,0,0);
			ic.cc[3].create_cut_parameters(-1,0,0,-1,0,0);
		}

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// Main procedure when the determinantal equation has two double roots
quad_inter <bigint> inter_two_mult(const bigint_matrix &q1, const bigint_matrix &q2, 
																	 const hom_polynomial <bigint> &det_p, 
																	 const hom_polynomial <bigint> &det_p_orig,
																	 const hom_polynomial <bigint> &gcd_p, 
																	 const int opt_level, ostream &s)
{
	#ifndef NDEBUG
	s << ">> entering inter_two_mult "<< endl;
	#endif
	
	quad_inter <bigint> ic;

	// The discriminant of the (degree 2) gcd_p
	bigint delta;
	delta.assign(discriminant2(gcd_p));
	
	// The square root of delta
	bigint sqrt_delta;

	// If the discriminant is positive and is a square, the two roots are real and
	// rational so no need for the big machinery
	if (is_square(sqrt_delta,delta))
		{
			// The two rational roots of the determinantal equation
			math_vector <bigint> root1(2,2), root2(2,2);

			if (gcd_p[0].is_zero()) // The bad case
				{
		 			// root1 is [ 0 1 ]
		 			root1[1].assign_one();
		 
		 			if (gcd_p[2].is_zero()) // root2 is [ 1 0 ]
					 	root2[0].assign_one();
					else		// root2 is [ -gcd_p[1], gcd_p[2] ]
						{
							root2[0].assign(gcd_p[1]);
							root2[0].negate();
							root2[1].assign(gcd_p[2]);
						}
				}
			else
				{
					// root1 = [2a,-b-sqrt_delta]
					multiply(root1[0],gcd_p[0],2);
					add(root1[1],gcd_p[1],sqrt_delta);
					root1[1].negate();
					 
					// root2 = [2a,-b+sqrt_delta]
					multiply(root2[0],gcd_p[0],2);
					subtract(root2[1],sqrt_delta,gcd_p[1]);
				}

			#ifndef NDEBUG
			s << ">> optimization of coordinates of roots" << endl;
			#endif

			optimize(root1);
			optimize(root2);

			// The associated matrices
			bigint_matrix qa = root1[0]*q1+root1[1]*q2;
			bigint_matrix qb = root2[0]*q1+root2[1]*q2;

			// The inertias
			math_vector <int> in_qa = inertia(qa);
			math_vector <int> in_qb = inertia(qb);
			
			// The ranks
			rpl_size_t rank_qa = in_qa[0]+in_qa[1];
			rpl_size_t rank_qb = in_qb[0]+in_qb[1];

			if (rank_qa < rank_qb)
				{
					swap(rank_qa,rank_qb);
					swap(in_qa,in_qb);
					swap(root1,root2);
					swap(qa,qb);
				}

			#ifndef NDEBUG
			s << ">> ranks of singular quadrics: " << rank_qa << " and " 
				<< rank_qb << endl; 
			#endif

			#ifndef NDEBUG
			s << ">> " << "two real rational double roots" << ": " << root1 << " and " 
				<< root2 << endl;
			#endif 
	
			if ((rank_qa == 3) && (rank_qb == 3))
				ic = inter_cubic_secant_two_rat(q1,q2,qa,qb,opt_level,s);
			else
				{
					// The sign of the determinantal equation is the same for any point
					// other than the roots
			
					// A point ``outside'' the roots
					math_vector <bigint> root_p(2,2);
			
					if (!gcd_p[2].is_zero())
						root_p[0].assign_one();
					else
						{
							root_p[0] = -gcd_p[0]-1;
							root_p[1] = gcd_p[1];
						}
			
					// The sign of det_p outside the roots
					rpl_size_t s_e;
					s_e = (det_p.eval(root_p[0],root_p[1])).sign();
			
					if (rank_qa == 3) 
						// Conic and two lines not crossing on the conic
						ic = inter_conic_2lines_not_crossing(q1,q2,qa,qb,in_qa[0],s_e,opt_level,s);
					else
						{
							if (s_e == -1)
								ic = inter_four_skew_lines_2pts_rat(q1,q2,qa,qb,in_qa[0],opt_level,s);
							else
								{
									// If inertia of one of the pair of planes is [2 0], empty intersection
						
									if (in_qa[0] == 2)
										{
											ic = quad_inter <bigint>(0);
						
											ic.set_type(14,1);
						
											#ifndef NDEBUG
											//print_type(ic,s);
											#endif	
										}
									else // There is a quadric of inertia [ 2 2 ] in the pencil
										ic = inter_four_skew_lines_rat(q1,q2,qa,qb,opt_level,s);
								}
						}
				}
		}
	else
		// The two roots are either complex or algebraically conjugate, so the ranks
		// of the associated matrices are the same 
		{
			// root = [a,b,c,d]
			// projective (a,b +/- c*sqrt{d}) if real
			// projective (a,b +/- c*i*sqrt{-d}) if complex
			math_vector <bigint> root(4,4);

			multiply(root[0],gcd_p[0],2);
			negate(root[1],gcd_p[1]);
			root[2].assign(1);
			root[3].assign(delta);

			#ifndef NDEBUG
			extract_message(opt_level,s,"optimization of coordinates of roots");
			#endif

			math_vector <bigint> sq_fact = extract_square_factors(root[3],opt_level,s);
			
			root[2].assign(sq_fact[0]);

			// Not very pretty but...
			root[3].assign(sq_fact[0]);

			optimize(root);

			root[3].assign(sq_fact[1]);

			// The singular quadrics are a*q1+(b +/- c*sqrt(d))*q2 (and similarly in the
			// complex case). The singular locus is spanned by vectors k1 +/-
			// sqrt(d)*k2, where the aggregate vector (k1 k2) is in the kernel of the 
			// 8x8 matrix:
			//			 (a*q1+b*q2		c*d*q2 )
			//			 (									 )
			//			 (	c*q2		a*q1+b*q2)
			// The size of the kernel is 4 if the singular quadrics have a
			// singular line and only 2 if they are imaginary cones (because if
			// k1+sqrt(d)*k2 is solution of the above, d*k2+sqrt(d)*k1 also is)
			bigint_matrix q_rat = root[0]*q1+root[1]*q2;
			bigint_matrix q_sq = root[2]*q2;
			bigint_matrix big_mat(8,8);
			big_mat.compose_t(q_rat,root[3]*q_sq,q_sq,q_rat);

			bigint_matrix sing = QI::singular(big_mat);

			int d_sing = sing.get_no_of_columns();

			if (delta > 0) // Roots are real conjugate
				{ 
					#ifndef NDEBUG
					s << ">> " << "two real non-rational double roots" << ": ";
					print_root(root,s,1);
					s << " and ";
					print_root(root,s,-1);
					s << endl;
					#endif
			
					if (d_sing == 4)
						{
							#ifndef NDEBUG
							s << ">> ranks of singular quadrics: 2 and 2" << endl; 
							#endif	
			
							// Sign of leading coeff of det_p, necessarily non zero
							rpl_size_t s_e = det_p[0].sign();
			
							if (s_e == 1)
								{
									// Pick a lambda between the roots. If inertia is [4,0] then
									// intersection is empty
									bigint_matrix q = root[0]*q1+root[1]*q2;
			
									if ((inertia_known_rank(q,4)[0] == 4) || (inertia_known_rank(q1,4)[0] == 4))
										{
											#ifndef NDEBUG
											s << ">> found quadric of inertia [ 4 0 ] in pencil" << endl;
											#endif	
			
											ic = quad_inter <bigint>(0);
			
											ic.set_type(14,1);
			
											#ifndef NDEBUG
											// print_type(ic,s);
											#endif	
										}
									else
										ic = inter_four_skew_lines_non_rat(q1,q2,sing,root[3],s_e,opt_level,s);
								}
				 			else // (s_e != 1) - Two points
								ic = inter_four_skew_lines_non_rat(q1,q2,sing,root[3],s_e,opt_level,s);
			 			}
		 			else
			 			{
							#ifndef NDEBUG
				 			s << ">> ranks of singular quadrics: 3 and 3" << endl; 
							#endif	

				 			bigint_matrix q;
				 			q = q1;
				 			ic = intersection_with_quadric22(q1,q2,det_p_orig,q,3,opt_level,s);
			 			}
				}
			else // Roots are complex conjugate
				{
					#ifndef NDEBUG
					s << ">> " << "two complex double roots" << ": ";
					print_root(root,s,1);
					s << " and ";
					print_root(root,s,-1);
					s << endl;
					#endif

		 			if (d_sing == 4) // rank 2
			 			{
				 			#ifndef NDEBUG
				 			s << ">> ranks of singular quadrics: 2 and 2" << endl; 
				 			#endif

				 			bigint_matrix q;
				 			q = q1;
				 			ic = intersection_with_quadric22(q1,q2,det_p_orig,q,5,opt_level,s);
			 			}
		 			else // rank 3
			 			{
				 			#ifndef NDEBUG
				 			s << ">> ranks of singular quadrics: 3 and 3" << endl; 
				 			#endif

				 			bigint_matrix q;
				 			q = q1;
				 			ic = intersection_with_quadric22(q1,q2,det_p_orig,q,4,opt_level,s);
			 			}
				}
		}
		
	#ifndef NDEBUG
	s << ">> exiting inter_two_mult "<< endl;
	#endif

	return ic;
}

} // end of namespace QI
