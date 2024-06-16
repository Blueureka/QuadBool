// Intersection when the determinantal equation has one multiple root

#include <libqi/kernel/QIOneMult.h>
#include <libqi/kernel/QIElem.h>
#include <libqi/kernel/QINumber.h>
#include <libqi/kernel/QIInterWith22.h>

using namespace rpl;

// Enter namespace QI
namespace QI {

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the case of two lines crossing on a conic
quad_inter <bigint> inter_conic_2lines_crossing(const bigint_matrix &q1, 
																								const bigint_matrix &q2,
																								const bigint_matrix &q, 
																								const bigint_matrix &q_other,
																								const bigint_matrix &q_sing,
																								const bool type_flag, 
																								const int opt_level, ostream &s)
{
	// Either one or three components
	quad_inter <bigint> ic;

	if (type_flag)
		{
			ic = quad_inter <bigint>(3);
			ic.set_type(8,2);
		}
	else
		{
			ic = quad_inter <bigint>(1);
			ic.set_type(8,1);
		}

	#ifndef NDEBUG
	// print_type(ic,s);
	#endif

	// First parameterize the pair of planes
	bigint D;
	bigint_matrix m1(4,3),m2(4,1);

	pair_of_planes_param(q,q_sing,D,m1,m2,opt_level,s);

	// The two planes are rational so we can directly replace the first column of m1
	// by m1(0) +/- m2(0)
	bigint_matrix m1p = m1;
	column(m1p, 0)	= column(m1, 0) + column(m2, 0);
	bigint_matrix m1m = m1;
	column(m1m, 0)	= column(m1, 0) - column(m2, 0);

	#ifndef NDEBUG
	s << ">> optimization of transformation" << endl;
	#endif

	optimize_trans1(m1p);
	optimize_trans1(m1m);

	// Now we compute the following two matrices: one has inertia [2 1]
	// and represents the conic, the other has inertia [2 0] and
	// represents the imaginary pair of lines
	bigint_matrix ap = prod(base_matrix<bigint> (prod(trans(m1p),q_other)),m1p);
	bigint_matrix am = prod(base_matrix<bigint> (prod(trans(m1m),q_other)),m1m);

	// Swap if ap is not the conic of inertia [2 1]
	if (det(ap).is_zero()) 
		{
			swap(ap,am);
			swap(m1p,m1m);
		}

	// The rational singular point of the pair of lines
	math_vector <bigint> line_sing = column(QI::singular(am), 0);
	math_vector <bigint> line_sing3d = prod(m1m, (math_vector <bigint>)line_sing);

	#ifndef NDEBUG
	s << ">> rational point on conic: " << line_sing3d << endl;
	#endif

	// Parameterization of the conic
	bigint_matrix p_trans(3,3);
	curve_param <bigint> conic_2d(3);
	math_vector <bigint> l(2,2);

	conic_param_through_ratpoint(ap,line_sing,p_trans,conic_2d,l,s);

	// Stack the transformations
	p_trans = prod(m1p, p_trans);

	#ifndef NDEBUG
	s << ">> optimization of parameterization" << endl;
	#endif

	optimize_trans3(p_trans,opt_level,s);

	// So that u=1, v=0 corresponds to the rational point
	swap(conic_2d[0],conic_2d[1]);

	curve_param <bigint> conic(4);
	multiply(conic,p_trans,conic_2d);

	ic.cc[0].create_component(INTER_TYPE_CONIC,conic);

	if (type_flag)
		{
			// The intersection also contains two lines - let us parameterize them
			bigint D;
			bigint_matrix m1(3,2),m2(3,1);

			pair_of_lines_param(am,line_sing,D,m1,m2,opt_level,s);

			if (D.is_zero())
				{
					// The two lines are rational
					bigint_matrix l1p = m1;
					column(l1p, 0) = column(m1, 0) + column(m2, 0);
					bigint_matrix l1m = m1;
					column(l1m, 0) = column(m1, 0) - column(m2, 0);			 
					 
					// Stack the transformation matrices
					bigint_matrix linep = prod(m1m, l1p);
					bigint_matrix linem = prod(m1m, l1m);
					 
					#ifndef NDEBUG
					s << ">> optimization of transformation" << endl;
					#endif
			
					optimize_trans1(linep);
					optimize_trans1(linem);
			
					// The two lines
					curve_param <bigint> line1(column(linep, 0),column(linep, 1)),
															line2(column(linem, 0),column(linem, 1)); 
			
					// Reparameterize the lines
			
					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif
			
					// Cut parameter of the lines ?????
					math_vector <bigint> cut1(4,4), cut2(4,4);
			
					cut1 = improved_line_param(line1);
					cut2 = improved_line_param(line2);
			
					ic.cc[1].create_component(INTER_TYPE_LINE,line1);
					ic.cc[2].create_component(INTER_TYPE_LINE,line2);
				}
			else
				{
					// The lines are not rational
			
					// Stack the transformations
					m1 = prod(m1m, m1);
					m2 = prod(m1m, m2);
					m2.resize(4,2);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
			
					optimize_trans2(m1,m2);
			
					// The components of the lines
					curve_param <bigint> line_rat(column(m1, 0), column(m1, 1)),
																line_D(column(m2, 0), column(m2, 1));
			
					// Reparameterize the lines
			
					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif
			
					improved_line_param(line_rat,line_D,D);
			
					ic.create_two_components(1,INTER_TYPE_LINE,line_rat,line_D,D);
				}

			#ifndef NDEBUG
			s << ">> the lines meet the conic at " << line_sing3d << endl;
			#endif
		}

	// Temp
	if (type_flag)
		{
			ic.cc[0].create_cut_parameter(-1,bigint(0),bigint(0));
			ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));
			ic.cc[2].create_cut_parameter(-1,bigint(0),bigint(0));
		}

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the two skew lines and double line case
quad_inter <bigint> inter_two_skew_lines_double_line(const bigint_matrix &q1, 
																											const bigint_matrix &q2,
																											const bigint_matrix &q, 
																											const bigint_matrix &q_other,
																											const bigint_matrix &q_sing,
																											const bool type_flag, 
																											const int opt_level, 
																											ostream &s)
{
	// One or three components
	quad_inter <bigint> ic;

	if (type_flag)
		{
			ic = quad_inter <bigint>(3);
			ic.set_type(9,2);
		}
	else
		{
			ic = quad_inter <bigint>(1);
			ic.set_type(9,1);
		}

	#ifndef NDEBUG
	// print_type(ic,s);
	#endif

	curve_param <bigint> line_par(column(q_sing, 0),column(q_sing, 1));

	#ifndef NDEBUG
	s << ">> reparameterization of line" << endl;
	#endif

	// Cut parameter of the lines ?????
	math_vector <bigint> cut1(4,4), cut2(4,4);

	cut1 = improved_line_param(line_par);			 

	ic.cc[0].create_component(INTER_TYPE_LINE,line_par);
	ic.cc[0].mult = 2;

	if (type_flag)
		{
			// There are also two skew lines: let us parameterize them

			// First parameterize the pair of planes
			bigint D;
			bigint_matrix m1(4,3),m2(4,1);

			pair_of_planes_param(q,q_sing,D,m1,m2,opt_level,s);

			// We want to compute (m1*par1)^T*q_other*(m1*par1), but we
			// don't do it explicitly. But since the singular line of the
			// pair of planes is a factor of this plug (corresponding to
			// `s'), the matrix has a very special shape, i.e. the 2 x 2 lower
			// right matrix is zero
			bigint_matrix a = prod(base_matrix<bigint> (prod(trans(m1), q_other)), m1);
			//bigint b0 = trans(m2)*q_other*m2)(0,0);
			math_vector<bigint> tmp = column(m2, 0);
			math_vector<bigint> tmp2 = prod(q_other, tmp);
			bigint b0 = inner_prod(tmp, tmp2);
			bigint_matrix c = prod(base_matrix<bigint> (prod(trans(m1), q_other)) , m2);

			// Extract those coefficients that are needed
			bigint a0 = a(0,0);
			bigint a1 = 2*a(0,1);
			bigint a2 = 2*a(0,2);
			bigint c0 = 2*c(0,0);
			bigint c1 = 2*c(1,0);
			bigint c2 = 2*c(2,0);

			if (D.is_zero())
				{
					// Each plane of the pair is rational
					// The plugs of the plane params in q_other are
					// s*(a0+b0+/-c0)+t*(a1+/-c1)+v*(a2+/-c2) = 0
					 		
					// Thus we can extract v = p(s,t) and build a new 4 x 2
					// transformation matrix...
					// If the coefficient of v is zero, special trick
					bigint_matrix vp(4,2),vm(4,2);
					bigint ab0 = a0 + b0;
			
					bigint tmp = a2+c2;
					if (tmp.is_zero())
						 {
					 	 	column(vp, 0) = column(m1, 2);
					 	 	column(vp, 1) = (a1	+ c1) * (column(m1, 0) + column(m2, 0))
					 									-(ab0 + c0) * column(m1, 1);
						 }
					else
						 {
					 	 	column(vp, 0) =	 tmp		* (column(m1, 0) + column(m2, 0))
					 									-(ab0 + c0) * column(m1, 2);
					 	 	column(vp, 1) = tmp * column(m1, 1) - (a1 + c1) * column(m1, 2);
						 }
			
					tmp = a2-c2;
					if (tmp.is_zero())
						 {
					 	 	column(vm, 0) = column(m1, 2);
					 	 	column(vm, 1) = (a1-c1) * (column(m1, 0) - column(m2, 0)) 
			 											-(ab0-c0) * column(m1, 1);
						 }
					else
						 {
					 	 	column(vm, 0) = tmp * (column(m1, 0)- column(m2, 0))
					 		 							-(ab0 - c0) * column(m1, 2);
					 	 
					 	 	column(vm, 1) = tmp *column(m1, 1) - (a1-c1) * column(m1, 2);
						 }
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif				 
			
					optimize_trans1(vp);
					optimize_trans1(vm);
			
					// The two skew lines
					curve_param <bigint> linep(column(vp, 0), column(vp, 1)),
					 											linem(column(vm, 0), column(vm, 1));
			
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
					// Individual planes are not rational
					// The plugs of the plane params in q_other are
					// s*(a0+D*b0)+t*a1+v*a2 +/- sqrt(D)*(c0*s+c1*t+c2*v) = 0
			
					// Extract dt*v = (cs_rat*s+ct_rat*t) + sqrt(D)*(cs_D*s+ct_D*t)
					bigint ab0 = a0+D*b0;
			
					bigint dt = a2*a2-D*c2*c2;
					
					bigint cs_rat = D*c2*c0-a2*ab0;
					bigint ct_rat = D*c2*c1-a2*a1;
					bigint cs_D = ab0*c2-a2*c0;
					bigint ct_D = a1*c2-a2*c1;
			
					// Build new projection matrices
					bigint_matrix m_rat(4,2),m_D(4,2);
					 		
					column(m_rat, 0) = dt * column(m1, 0) + cs_rat * column(m1, 2);
					
					column(m_rat, 1) = dt * column(m1, 1) + ct_rat * column(m1, 2);
					
					column(m_D, 0) = dt *	column(m2, 0) + cs_D * column(m1, 2);
					
					column(m_D, 1) = ct_D * column(m1, 2);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
			
					optimize_trans2(m_rat,m_D);
			
					// The components of the skew lines
					curve_param <bigint> line_rat(column(m_rat, 0), column(m_rat, 1)),
					 											line_D(column(m_D, 0), column(m_D, 1));
			
					// Reparameterize the lines
			
					#ifndef NDEBUG
					s << ">> reparameterization of lines" << endl;
					#endif
			
					improved_line_param(line_rat,line_D,D);
			
					ic.create_two_components(1,INTER_TYPE_LINE,line_rat,line_D,D);
				}
		}

	// Temp
	if (type_flag)
		{
			ic.cc[0].create_cut_parameters(-1,0,0,-1,0,0);
			ic.cc[1].create_cut_parameter(-1,bigint(0),bigint(0));
			ic.cc[2].create_cut_parameter(-1,bigint(0),bigint(0));
		}

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the two double lines case
quad_inter <bigint> inter_two_double_lines(const bigint_matrix &q1, 
																						const bigint_matrix &q2,
																						const bigint_matrix &q, 
																						const bigint_matrix &q_other,
																						const bigint_matrix &q_sing,
																						const rpl_size_t &s_e, const int opt_level, 
																						ostream &s)
{
	// One or two components
	quad_inter <bigint> ic;

	if (s_e == 1)
		{
			ic = quad_inter <bigint>(2);
			ic.set_type(10,2);
		}
	else
		{
			ic = quad_inter <bigint>(1);
			ic.set_type(10,1);
		}

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// First parameterize the double plane q
	bigint_matrix m1(4,3);

	double_plane_param(q,q_sing,m1,s);

	#ifndef NDEBUG
	s << ">> optimization of transformation" << endl;
	#endif

	optimize_trans1(m1);

	// Now we have a 3 x 3 matrix representing a (real or imaginary) pair of lines
	bigint_matrix am = prod(base_matrix<bigint> (prod(trans(m1), q_other)), m1);

	// The singular point of this pair of lines
	math_vector <bigint> sing_l_2d = column(QI::singular(am), 0);

	if (s_e == -1)
		{
			// The intersection is reduced to the singular point of the pair of lines
			math_vector <bigint> sing_l = prod(m1, sing_l_2d);

			#ifndef NDEBUG
			s << ">> optimization of point" << endl;
			#endif

			optimize(sing_l);

			curve_param <bigint> sing_par(sing_l);

			ic.cc[0].create_component(INTER_TYPE_POINT,sing_par);
		}
	else
		{
			// Two double lines

			// Parameterize the pair of lines
			bigint D;
			bigint_matrix m1_l(3,2),m2_l(3,1);

			pair_of_lines_param(am,sing_l_2d,D,m1_l,m2_l,opt_level,s);

			if (D.is_zero())
				{
					// The two lines are rational
					bigint_matrix l1p = m1_l;
					column(l1p, 0) = column(m1_l, 0) + column(m2_l, 0);
					bigint_matrix l1m = m1_l;
					column(l1m, 0) = column(m1_l, 0) - column(m2_l, 0);
			
					// Stack the transformation matrices
					bigint_matrix linep = prod(m1, l1p);
					bigint_matrix linem = prod(m1, l1m);
			
					#ifndef NDEBUG
					s << ">> optimization of transformation" << endl;
					#endif
			
					optimize_trans1(linep);
					optimize_trans1(linem);
			
					// The two lines
					curve_param <bigint> line1(column(linep, 1), column(linep, 0)),
					 										line2(column(linem, 1), column(linem, 0)); 
			
					// Cut parameter of the lines
					math_vector <bigint> cut1(2,2), cut2(2,2);
			
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
			else
				{
					// The lines are not rational
			
					// Stack the transformations
					m1_l = prod(m1, m1_l);
					m2_l = prod(m1, m2_l);
					m2_l.resize(4,2);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
			
					optimize_trans2(m1_l,m2_l);
			
					// The components of the lines
					curve_param <bigint> line_rat(column(m1_l, 1), column(m1_l, 0)),
					 										line_D(column(m2_l, 1), column(m2_l, 0));
			
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

			ic.cc[0].mult = 2;
			ic.cc[1].mult = 2;
		}

	ic.set_optiflag(true);

	return ic;
}


//////////////////////////////////////////////////////////////////////////////////
// The procedure for the cuspidal quartic case
quad_inter <bigint> inter_cuspidal_quartic(const bigint_matrix &q1, 
																						const bigint_matrix &q2,
																						const bigint_matrix &q, 
																						const bigint_matrix &q_other,
																						const bigint_matrix &q_sing,
																						const int opt_level, ostream &s)
{
	// One component in the intersection
	quad_inter <bigint> ic(1);

	ic.set_type(4,1);

	#ifndef NDEBUG
	// print_type(ic,s);
	#endif

	// Singular point of cone	 
	math_vector <bigint> sing = column(q_sing, 0);

	#ifndef NDEBUG
	s << ">> optimization of coordinates of singular point" << endl;
	#endif

	optimize(sing);

	#ifndef NDEBUG
	s << ">> singular point: " << sing << endl;
	#endif

	// We need to compute a point on the double line that is the intersection of
	// q with the tangent plane to q_other at sing
	// First compute the ``normal'' to the tangent plane of q_other at sing
	math_vector <bigint> dir;
	multiply(dir,q_other,sing);

	// Now a rational point on the line - We can't use a classical solve because we
	// want an answer in bigints ! (classical solve might introduce rationals)
	math_vector <bigint> rat_point = solve_proj(q,dir,sing);

	// Try to find a point of small height on the line by looking at the
	// parameter values such that one of the components is zero

	#ifndef NDEBUG
	s << ">> reparameterization of line" << endl;
	#endif

	// The line param
	curve_param <bigint> line_par(sing,rat_point);

	// Cut parameter of the lines ?????
	math_vector <bigint> cut1(4,4), cut2(4,4);

	cut1 = improved_line_param(line_par);

	if (are_equal(sing,line_par.eval(1,0)))
		rat_point = line_par.eval(0,1);
	else
		rat_point = line_par.eval(1,0);

	#ifndef NDEBUG
	s << ">> rational point on cone: " << rat_point << endl;
	#endif

	// Parameterize the cone
	surface_param <bigint> par(4);
	bigint_matrix p_trans(4,4);
	math_vector <bigint> l(2,2);

	cone_param_through_ratpoint(q,sing,rat_point,p_trans,par,l,s);

	#ifndef NDEBUG
	s << ">> optimization of transformation" << endl;
	#endif
		
	optimize_trans3(p_trans,l,opt_level,s);

	// The quadric in the canonical frame of the cone
	bigint_matrix q_tmp = prod(base_matrix<bigint>(prod(trans(p_trans), q_other)) , p_trans);

	// Plug in the other quadric
	hom_hom_polynomial <bigint> res = plug_param_in_quadric(par,q_tmp,par);

	// The solution (-l[1],l[0]) corresponds to the singular point of the
	// quartic. Let us rescale the thing so that it now stands at (1,0)
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

	//			optimize(m);

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

	l[0] = 0;
	l[1] = -1;

	hom_polynomial <bigint> p1 = res[2];
	hom_polynomial <bigint> p2;
	negate(p2,res[1]);

	// Expand param
	multiply(par,p_trans,par);

	curve_param <bigint> curve(4);
	curve = par.eval(p2,p1);

	#ifndef NDEBUG
	s << ">> optimization of parameterization" << endl;
	#endif

	optimize(curve);

	ic.cc[0].create_component(INTER_TYPE_CUSPIDAL_QUARTIC,curve);
	ic.cc[0].create_cut_parameter(-1,-l[1],l[0]);

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the two tangent conics case
quad_inter <bigint> inter_two_tangent_conics(const bigint_matrix &q1, 
																							const bigint_matrix &q2,
																							const bigint_matrix &q, 
																							const bigint_matrix &q_other,
																							const bigint_matrix &q_sing, 
																							const int img_plane,
																							const int opt_level, ostream &s)
{
	// One or two components here
	quad_inter <bigint> ic;
	
	if (img_plane == 1)
		{
			ic = quad_inter <bigint>(2);
			ic.set_type(5,2);
		}
	else
		{
			ic = quad_inter <bigint>(1);
			ic.set_type(5,1);
		}

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// The singular line of q intersects q_other in a rational double point
	curve_param <bigint> sing_line(column(q_sing, 0), column(q_sing, 1));
	curve_param <bigint> par_tmp(4);
	hom_polynomial <bigint> tmp;
	multiply(par_tmp,q_other,sing_line);
	multiply(tmp,sing_line,par_tmp);

	math_vector <bigint> db_point;
	if (!tmp[2].is_zero())
		db_point = sing_line.eval(-tmp[1],2*tmp[2]);
	else
		db_point = sing_line.eval(2*tmp[0],-tmp[1]);

	#ifndef NDEBUG
	s << ">> optimization of coordinates of rational double point" << endl;
	#endif

	optimize(db_point);

	if (img_plane == 2) // The intersection is reduced to the double point
		{
			curve_param <bigint> sing_par(db_point);

			ic.cc[0].create_component(INTER_TYPE_POINT,sing_par);
		}
	else
		{
			// First parameterize the cone			

			math_vector <bigint> cone_sing = column(QI::singular(q_other), 0);
			
			math_vector <bigint> rat_point = db_point;

			// Try to find a point of small height on the line joining the singular
			// point of the cone to db_point and take that as rational point to
			// parameterize the cone

			// The line param
			curve_param <bigint> line_par(cone_sing,db_point);
			
			// Cut parameter of the lines ?????
			math_vector <bigint> cut1(4,4), cut2(4,4);

			cut1 = improved_line_param(line_par);

			if (are_equal(cone_sing,line_par.eval(1,0)))
				rat_point = line_par.eval(0,1);
			else
				rat_point = line_par.eval(1,0);

			#ifndef NDEBUG
			s << ">> rational point on cone: " << rat_point << endl;
			#endif

			surface_param <bigint> par(4);
			bigint_matrix p_trans(4,4);
			math_vector <bigint> l(2,2);

			cone_param_through_ratpoint(q_other,cone_sing,rat_point,p_trans,par,l,s);

			// The solution (-l[1],l[0]) corresponds to the line joining the
			// singular point of the cone to rat_point. Let us rescale the param so
			// that this line (and thus the point of tangency of the two conics)
			// corresponds to the parameter (1,0).
			
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

			// Stack the transformations
			p_trans = prod(p_trans, par_trans);

			#ifndef NDEBUG
			s << ">> optimization of transformation" << endl;
			#endif
		
			optimize_trans3(p_trans,l,opt_level,s);

			l[0] = 0;
			l[1] = -1;

			// Parameterize the pair of planes
			bigint D;
			bigint_matrix m1(4,3),m2(4,1);

			pair_of_planes_param(q,q_sing,D,m1,m2,opt_level,s);

			curve_param <bigint> conic2d(3);
			conic2d[0].assign_x2();
			conic2d[1].assign_y2();
			conic2d[2].assign_xy();

			if (D.is_zero())
				{
					// The two planes are rational so we can directly replace the first
					// column of m1 by m1(0) +/- m2(0)
					bigint_matrix m1p = m1;
					column(m1p, 0) = column(m1, 0) + column(m2, 0);
					bigint_matrix m1m = m1;
					column(m1m, 0) = column(m1, 0) - column(m2, 0);
			
					#ifndef NDEBUG
					s << ">> optimization of transformation" << endl;
					#endif
			
					optimize_trans1(m1p);
					optimize_trans1(m1m);
			
					bigint_matrix v1p = prod(trans(kernel(trans(m1p))), p_trans);
					bigint_matrix n1p(4,3);
					n1p.sto(0,0,-v1p(0,3));
					n1p.sto(1,1,-v1p(0,3));
					n1p.sto(2,2,-v1p(0,3));
					n1p.sto(3,0,v1p(0,0));
					n1p.sto(3,1,v1p(0,1));
					n1p.sto(3,2,v1p(0,2));
			
					bigint_matrix v1m = prod(trans(kernel(trans(m1m))), p_trans);
					bigint_matrix n1m(4,3);
					n1m.sto(0,0,-v1m(0,3));
					n1m.sto(1,1,-v1m(0,3));
					n1m.sto(2,2,-v1m(0,3));
					n1m.sto(3,0,v1m(0,0));
					n1m.sto(3,1,v1m(0,1));
					n1m.sto(3,2,v1m(0,2));
			
					bigint_matrix p_trans1 = prod(p_trans,n1p);
					bigint_matrix p_trans2 = prod(p_trans,n1m);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
					
					optimize_trans3(p_trans1,opt_level,s);
					optimize_trans3(p_trans2,opt_level,s);
			
					curve_param <bigint> conic1(4),conic2(4);
					multiply(conic1,p_trans1,conic2d);
					multiply(conic2,p_trans2,conic2d);
			
					ic.cc[0].create_component(INTER_TYPE_CONIC,conic1);
					ic.cc[1].create_component(INTER_TYPE_CONIC,conic2);
				}
			else
				{
					// The two planes are not rational
			
					m2.resize(4,3);
			
					// The plane transformations are m1 +/- sqrt(D)*m2. We want to compute
					// the kernel of (the trans of) those things
					bigint_matrix big_mat(6,8);
					big_mat.compose_t(trans(m1),D*trans(m2),trans(m2),trans(m1));
			
					bigint_matrix sing = QI::singular(big_mat);
			
					// The kernel of trans(m1 +/- sqrt(D)*m2) is generated by k1 +/- sqrt(D)*k2
					// k3 and k4 are tmp vars
					bigint_matrix k1(4,1),k2(4,1),k3(4,1),k4(4,1);
			
					sing.split_t(k1,k3,k2,k4);
			
					bigint_matrix v_rat = prod(trans(k1), p_trans);
					bigint_matrix v_sq = prod(trans(k2), p_trans);
			
					bigint_matrix n_rat(4,3),n_sq(4,3);
					n_rat.sto(0,0,-v_rat(0,3));
					n_rat.sto(1,1,-v_rat(0,3));
					n_rat.sto(2,2,-v_rat(0,3));
					n_rat.sto(3,0,v_rat(0,0));
					n_rat.sto(3,1,v_rat(0,1));
					n_rat.sto(3,2,v_rat(0,2));
					n_sq.sto(0,0,-v_sq(0,3));
					n_sq.sto(1,1,-v_sq(0,3));
					n_sq.sto(2,2,-v_sq(0,3));
					n_sq.sto(3,0,v_sq(0,0));
					n_sq.sto(3,1,v_sq(0,1));
					n_sq.sto(3,2,v_sq(0,2));
			
					bigint_matrix p_trans_rat = prod(p_trans, n_rat);
					bigint_matrix p_trans_sq = prod(p_trans, n_sq);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
					
					optimize_trans4(p_trans_rat,p_trans_sq,opt_level,s);
			
					curve_param <bigint> conic_rat(4),conic_sq(4);
					multiply(conic_rat,p_trans_rat,conic2d);
					multiply(conic_sq,p_trans_sq,conic2d);
			
					ic.create_two_components(0,INTER_TYPE_CONIC,conic_rat,conic_sq,D);
				}

			ic.cc[0].create_cut_parameter(-1,-l[1],l[0]);
			ic.cc[1].create_cut_parameter(-1,-l[1],l[0]);
		}

	#ifndef NDEBUG
	if (img_plane == 1)
		s << ">> the two conics are tangent at " << db_point << endl;
	#endif

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the nodal quartic case
quad_inter <bigint> inter_nodal_quartic(const bigint_matrix &q1, 
																				const bigint_matrix &q2,
																				const bigint_matrix &q, 
																				const bigint_matrix &q_other, 
																				const bigint_matrix &q_sing,
																				const bigint &delta_e, const rpl_size_t &s_e, 
																				const int in_q, 
																				const int opt_level, ostream &s)
{
	// One or two components (when isolated singularity) 
	quad_inter <bigint> ic;

	if ((delta_e > 0) && (s_e == -1) && (in_q == 2))
		ic = quad_inter <bigint>(2);
	else
		ic = quad_inter <bigint>(1);

	// The node of the quartic = singular point of cone
	math_vector <bigint> sing = column(q_sing, 0);

	#ifndef NDEBUG
	s << ">> optimization of coordinates of singular point" << endl;
	#endif

	optimize(sing);

	if (delta_e < 0) // Other roots of determinantal equation are complex
		ic.set_type(2,4);
	else // Other roots of determinantal equation are real
		if (s_e == 1)
			ic.set_type(2,3);
		else
			if (in_q == 2) // The cone corresponding to the double real root is real
				ic.set_type(2,2);
			else // The cone is complex
				ic.set_type(2,1);

	#ifndef NDEBUG
	//print_type(ic,s);
	s << ">> node of quartic is at " << sing << endl;
	#endif

	if ((delta_e > 0) && (s_e == -1) && (in_q == 3))
		{
			curve_param <bigint> cone_p(sing);

			ic.cc[0].create_component(INTER_TYPE_POINT,cone_p);
			ic.cc[0].optiflag = true;
		}
	else
		{
			// This boolean is true if the singularity of the nodal quartic is isolated
			bool isolated_sing = (in_q == 2) && (s_e == -1) && (delta_e > 0);

			// Parameterize the projective cone q
			bigint_matrix p_trans(4,4),p_trans2(4,4);
			bigint D;
			surface_param <bigint> par(4);

			cone_param(q,column(q_sing,0),D,p_trans,p_trans2,par,opt_level,s);

			if (D.is_zero())
				{
					// Lucky you!!!! You just found a rational parameterization of your cone
					// while nothing guaranteed you would find one....
			
					#ifndef NDEBUG
					s << ">> optimization of transformation" << endl;
					#endif
								
					optimize_trans3(p_trans,opt_level,s);
			
					// The quadric in the canonical frame of the cone
					bigint_matrix q_tmp = prod( base_matrix<bigint> (prod(trans(p_trans),q_other)), p_trans);
			
					// Plug in the other quadric
					hom_hom_polynomial <bigint> res = plug_param_in_quadric(par,q_tmp,par);
			
					// The parameter values corresponding to the node, if applicable
					math_vector <bigint> l(2,2), m(2,2);
			
					bigint sq;
					is_square(sq,res[1][1]*res[1][1]-4*res[1][0]*res[1][2]);
			
					if ((!isolated_sing) && (sq > 0))
						{
							// res has the form:
							//			 --> u*v*(l[0]*s+l[1]*t)*(m[0]*s+m[1]*t) + u^2*p1(s,t),
							// where p1(s,t) has degree 4
							 
							// The two solutions (-l[1],l[0]) and (-m[1],m[0]) correspond to the
							// two parameters of the node of the quartic.
			
							hom_polynomial <bigint> p2(res[1]);
			
							if (p2[2] == 0)
								{
									l[1] = 1;
									m[0] = p2[1];
									m[1] = p2[0];
								}
							else
								{
									l[0] = 2*p2[2];
									m[0] = 2*p2[2];
									l[1] = p2[1]-sq;
									m[1] = p2[1]+sq;
								}
						
							optimize(l);
							optimize(m);				
			
							// Let us rescale the parameterization so that the node
							// corresponds to the parameters (1,0) and (0,1).
			
							// Rewrite the par accordingly, i.e. replace l[0]*s+l[1]*t by t and 
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
			
							l[0] = 0;
							l[1] = -1;
							m[0] = 1;
							m[1] = 0;
						}

					// We first extract p1 and p2
					hom_polynomial <bigint> p1 = res[2];
					hom_polynomial <bigint> p2;
					negate(p2,res[1]);
			
					// Expand param
					multiply(par,p_trans,par);
			
					// Compute the result
					curve_param <bigint> curve(4);
					curve = par.eval(p2,p1);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
			
					optimize(curve);
			
					ic.cc[0].create_component(INTER_TYPE_NODAL_QUARTIC,curve);
					ic.cc[0].optiflag = true;
					
					// Cut_params are wrong here (in view of above)
					if ((!isolated_sing) && (sq > 0))
						 ic.cc[0].create_cut_parameters(-1,-l[1],l[0],-1,-m[1],m[0]);
				}
			else
				{
					// Expand the params
					surface_param <bigint> par2(4);
					multiply(par2,p_trans2,par);
					multiply(par,p_trans,par);
			
					// Plug in the other quadric
					// res is par*q_other*par + D*par2*q_other*par2
					// res2 is 2*par2*q_other*par
					hom_hom_polynomial <bigint> res,res2,tmp;
					res = plug_param_in_quadric(par,q_other,par);
			
					multiply(tmp,plug_param_in_quadric(par2,q_other,par2),D);
					add(res,res,tmp);
					multiply(res2,plug_param_in_quadric(par2,q_other,par),(bigint)2);
			
					// Now res and res2 are of the form u^2*p1(s,t)+u*v*p2(s,t), where p1
					// has degree 4 and p2 has degree 2. Also, par and par2 are vectors of
					// components of the form p3(s,t)*u+p4(s,t)*v, where p3 has degree 2 and
					// p4 has degree 0 
			
					// We first extract p1 and p2
					hom_polynomial <bigint> p1 = res[2];
					hom_polynomial <bigint> p2;
					negate(p2,res[1]);
					hom_polynomial <bigint> p12 = res2[2];
					hom_polynomial <bigint> p22;
					negate(p22,res2[1]);
			
					// We then extract the components p3 and p4 of par and put them in
					// curve_param's
					// The result: curve + sqrt(D)*curve2
					curve_param <bigint> curve(4),curve2(4),cu(4);
					cu = par2.eval(p22,p12);
					multiply(cu,cu,D);
					add(curve,cu,par.eval(p2,p1));
					add(curve2,par.eval(p22,p12),par2.eval(p2,p1));
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
			
					optimize(curve,curve2);
			
					ic.cc[0].create_component(INTER_TYPE_NODAL_QUARTIC,curve,curve2,D);
			
					if (opt_level)
						 ic.cc[0].optiflag = true;
					else
						 ic.cc[0].optiflag = false;

					// The polynomial p2 + sqrt(D)*p22 factors in two linear factors
					// corresponding to the parameters at which the node of the quartic is
					// reached
					
					if (!isolated_sing)
						 {
					 	 	bigint d1 = p2[1]*p2[1]+D*p22[1]*p22[1]-4*p2[2]*p2[0]-4*D*p22[2]*p22[0];
					 	 	bigint d2 = 2*(p2[1]*p22[1]-2*p2[2]*p22[0]-2*p2[0]*p22[2]);
			
					 	 	bigint sq;
					 	 	is_square(sq,d1*d1-D*d2*d2);
			
					 	 	bigint dp = 2*(d1+sq), dm = 2*(d1-sq), dsq;
			
					 	 	// Things might simplify a bit
					 	 	if ((is_square(dsq,dp)) || (is_square(dsq,dm)))
								{
					 				math_vector <bigint> p1_rat(2,2),p2_rat(2,2),p1_D(2,2),p2_D(2,2);
			
									if ((p2[2] == 0) && (p22[2] == 0))
									 	{
									 		p1_rat[0] = D*p22[1]*p22[0]-p2[1]*p2[0];
									 		p1_rat[1] = p2[1]*p2[1]-D*p22[1]*p22[1];
									 		p1_D[0] = p22[1]*p2[0]-p22[0]*p2[1];
							
									 		optimize(p1_rat,p1_D);
							
									 		ic.cc[0].create_cut_parameters(-1,p1_rat[0],p1_D[0],p1_rat[1],
									 				 								-1,1,0,0,D);
									 	}
									else
									 	{
									 		bigint tmp1 = dsq*(-2*p2[1]+dsq), tmp2 = 2*(-dsq*p22[1]+d2);
									 		bigint tmp3 = dsq*(-2*p2[1]-dsq), tmp4 = 2*(-dsq*p22[1]-d2);
									 		bigint tmp5 = D*p22[2];
							
									 		p1_rat[0] = p2[2]*tmp1-tmp5*tmp2;
									 		p1_rat[1] = 4*dsq*(p2[2]*p2[2]-tmp5*p22[2]);
									 		p2_rat[0] = p2[2]*tmp3-tmp5*tmp4;
									 		p2_rat[1] = p1_rat[1];
							
									 		p1_D[0] = -p22[2]*tmp1+p2[2]*tmp2;
									 		p2_D[0] = -p22[2]*tmp3+p2[2]*tmp4;
							
									 		optimize(p1_rat,p1_D);
									 		optimize(p2_rat,p2_D);
				
											ic.cc[0].create_cut_parameters(-1,p1_rat[0],p1_D[0],p1_rat[1],
									 									 -1,p2_rat[0],p2_D[0],p2_rat[1],D);
										}
								}
				 			else
								// No simplification in factors: general situation
								{
									bigint ep = 1, em = 1, dm_old = dm, dp_old = dp;
						
									#ifndef NDEBUG
									extract_message(opt_level,s,"optimization of second square root");
									#endif
						
									math_vector <bigint> D_fact = extract_square_factors(dp,opt_level,s);
									ep = D_fact[0];
									dp = D_fact[1];
											
									D_fact = extract_square_factors(dm,opt_level,s);
									em = D_fact[0];
									dm = D_fact[1];
						
									bigint tmp1 = 4*d2*(p22[2]*p22[1]*D-p2[2]*p2[1]);
									bigint tmp2 = 8*d2*(p2[2]*p2[2]-D*p22[2]*p22[2]);
									bigint tmp3 = ep*(2*p2[2]*d2-p22[2]*dm_old);
									bigint tmp4 = em*(2*p2[2]*d2-p22[2]*dp_old);
									bigint tmp5 = 2*(p22[2]*p2[1]-p22[1]*p2[2]);
									
									if ((D > dp) && (D > dm))
										{
											// dp and dm are the smallest
											math_vector <bigint> p_rat(2,2),p_dp(2,2),p_dm(2,2),p_dpdm(2,2);
											
											p_rat[0] = tmp1;
											p_rat[1] = tmp2;
											
											p_dp[0] = tmp3; // coeff of sqrt(dp)
											p_dm[0] = tmp4; // coeff of sqrt(dm)
											p_dpdm[0] = ep*em*tmp5; // coeff of sqrt(dp*dm);
						
											optimize(p_rat,p_dp,p_dm,p_dpdm);
						
											if (dp > dm)
												{
										 			swap(dp,dm);
										 			swap(p_dp,p_dm);
												}
						
											ic.cc[0].create_cut_parameters(-1,p_rat[0],p_dp[0],p_dm[0],p_dpdm[0],p_rat[1],
								 															-1,p_rat[0],-p_dp[0],-p_dm[0],p_dpdm[0],p_rat[1],dp,dm);
										}
									else if ((dp > D) && (dp > dm))
										{
											// dm and D are the smallest
											math_vector <bigint> p_rat(2,2),p_dm(2,2),p_D(2,2),p_dmD(2,2);
											
											p_rat[0] = dm_old*tmp1;
											p_rat[1] = dm_old*tmp2;
											
											p_dm[0] = dm_old*tmp4;
											p_D[0] = 2*d2*dm_old*tmp5;
											p_dmD[0] = 2*d2*tmp3;
											
											optimize(p_rat,p_dm,p_D,p_dmD);
											
											if (dm > D)
												ic.cc[0].create_cut_parameters(-1,p_rat[0],p_D[0],p_dm[0],p_dmD[0],p_rat[1],
															 					-1,p_rat[0],p_D[0],-p_dm[0],-p_dmD[0],p_rat[1],D,dm);
											else
												ic.cc[0].create_cut_parameters(-1,p_rat[0],p_dm[0],p_D[0],p_dmD[0],p_rat[1],
															 					-1,p_rat[0],-p_dm[0],p_D[0],-p_dmD[0],p_rat[1],dm,D);
										}
									else
										{
											// dp and D are the smallest
											math_vector <bigint> p_rat(2,2),p_dp(2,2),p_D(2,2),p_dpD(2,2);
						
											p_rat[0] = dp_old*tmp1;
											p_rat[1] = dp_old*tmp2;
											
											p_dp[0] = dp_old*tmp3;
											p_D[0] = 2*d2*dp_old*tmp5;
											p_dpD[0] = 2*d2*tmp4;
											
											optimize(p_rat,p_dp,p_D,p_dpD);
											
											if (dp > D)
												ic.cc[0].create_cut_parameters(-1,p_rat[0],p_D[0],p_dp[0],p_dpD[0],p_rat[1],
															 					-1,p_rat[0],p_D[0],-p_dp[0],-p_dpD[0],p_rat[1],D,dp);
											else
												ic.cc[0].create_cut_parameters(-1,p_rat[0],p_dp[0],p_D[0],p_dpD[0],p_rat[1],
																			 -1,p_rat[0],-p_dp[0],p_D[0],-p_dpD[0],p_rat[1],dp,D);
										}
								}
						}
				}

			if ((in_q == 2) && (s_e == -1) && (delta_e > 0))
				{
		 			// Isolated singularity
		 			curve_param <bigint> cone_p(sing);

		 			ic.cc[1].create_component(INTER_TYPE_POINT,cone_p);
		 			ic.cc[1].optiflag = true;
				}
		}

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the double conic case
quad_inter <bigint> inter_double_conic(const bigint_matrix &q1, const bigint_matrix &q2,
																				const bigint_matrix &q, const bigint_matrix &q_sing,
																				const hom_polynomial <bigint> &det_e, 
																				const int opt_level, ostream &s)
{
	quad_inter <bigint> ic;

	// First compute the matrix associated to the other (simple) root of the
	// determinantal equation

	// The (lambda,mu) of the other (single) root of the pencil
	math_vector <bigint> root_cone(2,2);
		 
	negate(root_cone[0],det_e[0]);
	root_cone[1].assign(det_e[1]);

	// Simplify the coordinates of root_cone
	#ifndef NDEBUG
	s << ">> optimization of coordinates of cone root" << endl;
	#endif

	optimize(root_cone);

	#ifndef NDEBUG
	s << ">> cone root: " << root_cone << endl;
	#endif

	// The cone
	bigint_matrix q1_tmp,q2_tmp;
	multiply(q1_tmp, q1, root_cone[0]);
	multiply(q2_tmp, q2, root_cone[1]);

	bigint_matrix cone;
	add(cone, q1_tmp,q2_tmp);

	if (inertia_known_rank(cone,3)[0] == 3) // Cone is imaginary
		{
			#ifndef NDEBUG
			s << ">> cone is imaginary" << endl;
			#endif	

			ic = quad_inter <bigint>(0);

			ic.set_type(6,1);

			#ifndef NDEBUG
			//print_type(ic,s);
			#endif	
		}
	else // Cone is real: intersection is a double conic
		{
			// One component
			ic = quad_inter <bigint>(1);
			
			#ifndef NDEBUG
			s << ">> cone is real" << endl;
			#endif

			ic.set_type(6,2);
			
			#ifndef NDEBUG
			//print_type(ic,s);
			#endif	

			// First parameterize the double plane q
			bigint_matrix m1(4,3);

			double_plane_param(q,q_sing,m1,s);

			#ifndef NDEBUG
			s << ">> optimization of transformation" << endl;
			#endif

			optimize_trans1(m1);

			// Now we have a 3 x 3 matrix representing a real non-singular conic
			bigint_matrix conic_2d = prod(base_matrix<bigint> (prod(trans(m1),cone)),m1);

			// Parameterize the conic
			bigint D;
			bigint_matrix p_trans(3,3),p_trans2(3,3);
			curve_param <bigint> par(3);

			conic_param(conic_2d,D,p_trans,p_trans2,par,opt_level,s);

			if (D.is_zero())
				{
					// Lucky you: you just found a rational conic!!
			
					// Stack the transformations
					p_trans = prod(m1, p_trans);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
			
					optimize_trans3(p_trans,opt_level,s);
			
					curve_param <bigint> conic(4);
					multiply(conic,p_trans,par);
			
					ic.cc[0].create_component(INTER_TYPE_CONIC,conic);
					ic.cc[0].optiflag = true;
				}
			else
				{
					// Have to content oneself with one square root
			
					// Stack the transformations
					p_trans = prod(m1, p_trans);
					p_trans2 = prod(m1, p_trans2);
			
					#ifndef NDEBUG
					s << ">> optimization of parameterization" << endl;
					#endif
			
					optimize_trans4(p_trans,p_trans2,opt_level,s);
			
					curve_param <bigint> conic(4),conic2(4);
					multiply(conic,p_trans,par);
					multiply(conic2,p_trans2,par);
			
					ic.cc[0].create_component(INTER_TYPE_CONIC,conic,conic2,D);
			
					if (opt_level)
						 ic.cc[0].optiflag = true;
					else
						 ic.cc[0].optiflag = false;
				}

			ic.cc[0].mult = 2;
		}

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the cubic and tangent line case
quad_inter <bigint> inter_cubic_tangent_line(const bigint_matrix &q1, 
																							const bigint_matrix &q2,
																							const bigint_matrix &q, 
																							const bigint_matrix &q_other, 
																							const bigint_matrix &q_sing, 
																							const int opt_level, ostream &s)
{
	// Two components always
	quad_inter <bigint> ic(2);

	ic.set_type(7,1);

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif	

	// Singular point of cone	 
	math_vector <bigint> sing = column(q_sing, 0);

	#ifndef NDEBUG
	s << ">> optimization of coordinates of singular point" << endl;
	#endif

	optimize(sing);

	#ifndef NDEBUG
	s << ">> singular point: " << sing << endl;
	#endif

	// The line of the intersection is the intersection of the tangent plane
	// of q_other at sing with the cone
	math_vector <bigint> dir;
	multiply(dir,q_other,sing);

	// Now a rational point on the line - We can't use a classical solve because we
	// want an answer in bigints ! (classical solve might introduce rationals)
	math_vector <bigint> rat_point = solve_proj(q,dir,sing);

	// The line param
	curve_param <bigint> line_par(sing,rat_point);

	// Cut parameter of the line
	math_vector <bigint> cut(2,2);

	// Try to find a point of small height on the line by looking at the
	// parameter values such that one of the components is zero

	#ifndef NDEBUG
	s << ">> reparameterization of line" << endl;
	#endif

	cut = improved_line_param(line_par);

	optimize_by_half(cut);

	if (are_equal(sing,line_par.eval(1,0)))
		rat_point.assign(line_par.eval(0,1));
	else
		rat_point.assign(line_par.eval(1,0));

	#ifndef NDEBUG
	s << ">> rational point on cone: " << rat_point << endl;
	#endif

	// Parameterize the cone
	surface_param <bigint> par(4);
	bigint_matrix p_trans(4,4);
	math_vector <bigint> l(2,2);

	cone_param_through_ratpoint(q,sing,rat_point,p_trans,par,l,s);

	#ifndef NDEBUG
	s << ">> optimization of transformation" << endl;
	#endif
						
	optimize_trans3(p_trans,l,opt_level,s);

	// The quadric in the canonical frame of the cone
	bigint_matrix q_tmp = prod( base_matrix<bigint> (prod(trans(p_trans), q_other)), p_trans);

	// Plug in the other quadric
	hom_hom_polynomial <bigint> res = plug_param_in_quadric(par,q_tmp,par);

	// The common factor of res (see below)
	hom_polynomial <bigint> tmp;

	// res has the form:
	//			 ---> u^2*(l[0]*s+l[1]*t)*p1(s,t)+u*v*(l[0]*s+l[1]*t)^2
	// The solution (-l[1],l[0]) corresponds to the point of tangency of the
	// line with the cubic

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
	
	tmp.assign_y();
	
	l[0] = 0;
	l[1] = -1;

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
	s << ">> cubic and line are tangent at " << sing << endl;
	#endif

	ic.cc[0].create_component(INTER_TYPE_CUBIC,cubic);
	ic.cc[0].create_cut_parameter(-1,-l[1],l[0]);

	ic.cc[1].create_component(INTER_TYPE_LINE,line_par);
	ic.cc[1].create_cut_parameter(-1,cut);

	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The pair of non-rational planes of the pencil is parameterized by 
// m1*[u v s] +/- sqrt(D)*m2*[u] (D != 0, 1)
// q_22_pos is a matrix of inertia 22
quad_inter <bigint> inter_secant_conics_pl_not_rat(const bigint_matrix &q1, 
																										const bigint_matrix &q2,	 
																										const hom_polynomial <bigint> &det_p, 
																										const bigint_matrix &q_22_pos, 
																										const bigint &D, 
																										const bigint_matrix &m1, 
																										const bigint_matrix &m2, 
																										const unsigned int rtype,
																										const int opt_level, ostream &s)
{
	// Two components
	quad_inter <bigint> ic(2);

	ic.set_type(3,rtype);

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// Compute the new quadric q of inertia 2,2; 
	// point_on_q is a point on it; det_q its determinant in the form [a,b] with det
	// = a^2. b	 (duplicate here the beginning of intersection_with_quadric22())
	math_vector <bigint> point_on_q(4,4); 
	math_vector <bigint> det_q(2,2);

	// q_22_pos is const: make a copy
	bigint_matrix q_22 = q_22_pos;
	bigint det_R;
	find_quadric_22_and_point(q1,q2,det_p,q_22_pos,q_22,point_on_q,det_R,det_q,opt_level,s);
	
	#ifndef NDEBUG
	s << ">> quadric (2,2) found: ";
	print_quadric(q_22,s);
	s << endl;
	s << ">> decomposition of its determinant [a,b] (det = a^2*b): " << det_q << endl;
	s << ">> a point on the quadric: " << point_on_q << endl;
	#endif

	surface_param <bigint> s1(4),s2(4);	 //s1 and s2 are initialized to [0,0,0,0]
	find_parameterization_quadric_22(q_22,point_on_q,det_q,det_R,s1,s2,opt_level,s);
	// Parameterization of the quadric is	 s1 + sqrt(det_q[1]) s2
	// If det_q[1]=1 then	 s2=[0,0,0,0]

	// We first compute the implicit equation of the 2 non-rational planes 
	// The equation of the plane is the determinat of the 4x4 matrix whose 
	// first three collomns are 3 points on the planes and the last one is	[x y z w].
	// Three points are :
	//			 m1[1]+/- sqrt(D)*m2 (for u=1, v=s=0)
	//			 m2[1] for (v=1, u=s=0) and m3[1] (for s=1, u=v=0). 
	
	// We compute the six 2x2 sub-determinants of the 4x2 matrix consisting of the
	// collomn 2 and 3 of m1 
	bigint det01, det02, det03, det12, det13, det23;
	det01 = m1.member(0,1)*m1.member(1,2) - m1.member(0,2)*m1.member(1,1);
	det02 = m1.member(0,1)*m1.member(2,2) - m1.member(0,2)*m1.member(2,1);
	det03 = m1.member(0,1)*m1.member(3,2) - m1.member(0,2)*m1.member(3,1);
	det12 = m1.member(1,1)*m1.member(2,2) - m1.member(1,2)*m1.member(2,1);
	det13 = m1.member(1,1)*m1.member(3,2) - m1.member(1,2)*m1.member(3,1);
	det23 = m1.member(2,1)*m1.member(3,2) - m1.member(2,2)*m1.member(3,1);
	// Note the	 first column of m1 : [a0, a1, a2, a3] and m2 = [a0', a1', a2', a3']. 
	// Then the implicit equation of the planes are (in var x0, x1, x2, x3) :
	//	x0*[a1*det23 - a2*det13 + a3*det12] -x1*[a0*det23 - a2*det03 + a3*det02] 
	// +x2*[a0*det13 - a1*det03 + a3*det01] -x3*[a0*det12 - a1*det02 + a2*det01] 
	// +/-sqrt(D)*[the same thing with [a0, a1, a2, a3] rplaced by [a0', a1', a2', a3']. 
	//	 = (A0*x0 + A1*x1 + A2*x2 + A3*x3) +/- sqrt(D)*(A0p*x0 + A1p*x1 + A2p*x2 + A3p*x3)
	bigint A0, A1, A2, A3, A0p, A1p, A2p, A3p; 
	A0 = m1.member(1,0)*det23 - m1.member(2,0)*det13 + m1.member(3,0)*det12;
	A1 = -m1.member(0,0)*det23 + m1.member(2,0)*det03 - m1.member(3,0)*det02;
	A2 = m1.member(0,0)*det13 - m1.member(1,0)*det03 + m1.member(3,0)*det01;
	A3 = -m1.member(0,0)*det12 + m1.member(1,0)*det02 - m1.member(2,0)*det01;
	A0p = m2.member(1,0)*det23 - m2.member(2,0)*det13 + m2.member(3,0)*det12;
	A1p = -m2.member(0,0)*det23 + m2.member(2,0)*det03 - m2.member(3,0)*det02;
	A2p = m2.member(0,0)*det13 - m2.member(1,0)*det03 + m2.member(3,0)*det01;
	A3p = -m2.member(0,0)*det12 + m2.member(1,0)*det02 - m2.member(2,0)*det01;
		
	// We plug the parameterization in each of the 2 non-rational planes. 
	// Recall the parameterization of the quadric is	s1 + sqrt(det_q[1]) s2
	// Compute the polynomial poly = P1 + sqrt(xi) P2 +/- sqrt(D)*(P3 + sqrt(xi) P4)
	hom_hom_polynomial <bigint> P1, P2, P3, P4, hhtmp;
	// P1 = A0*s1[0] + A1*s1[1] + A2*s1[2] + A3*s1[3];
	multiply(P1, s1[0], A0);
	multiply(hhtmp, s1[1], A1);
	add(P1, P1, hhtmp);
	multiply(hhtmp, s1[2], A2);
	add(P1, P1, hhtmp);
	multiply(hhtmp, s1[3], A3);
	add(P1, P1, hhtmp);
	// P2 = A0*s2[0] + A1*s2[1] + A2*s2[2] + A3*s2[3];
	multiply(P2, s2[0], A0);
	multiply(hhtmp, s2[1], A1);
	add(P2, P2, hhtmp);
	multiply(hhtmp, s2[2], A2);
	add(P2, P2, hhtmp);
	multiply(hhtmp, s2[3], A3);
	add(P2, P2, hhtmp);
	// P3 = A0p*s1[0] + A1p*s1[1] + A2p*s1[2] + A3p*s1[3];
	multiply(P3, s1[0], A0p);
	multiply(hhtmp, s1[1], A1p);
	add(P3, P3, hhtmp);
	multiply(hhtmp, s1[2], A2p);
	add(P3, P3, hhtmp);
	multiply(hhtmp, s1[3], A3p);
	add(P3, P3, hhtmp);
	// P4 = A0p*s2[0] + A1p*s2[1] + A2p*s2[2] + A3p*s2[3];
	multiply(P4, s2[0], A0p);
	multiply(hhtmp, s2[1], A1p);
	add(P4, P4, hhtmp);
	multiply(hhtmp, s2[2], A2p);
	add(P4, P4, hhtmp);
	multiply(hhtmp, s2[3], A3p);
	add(P4, P4, hhtmp);

	// Linear equation : poly = P1 + sqrt(xi) P2 +/- sqrt(D)*(P3 + sqrt(xi) P4)
	// = u. (P1[1] + sqrt(xi)*P2[1] +/- sqrt(D)*P3[1] +/- sqrt(D*xi)*P4[1])
	//	+v. (P1[0] + sqrt(xi)*P2[0] +/- sqrt(D)*P3[0] +/- sqrt(D*xi)*P4[0])
	// Recall the parameterization of the quadric 22 is	 s1 + sqrt(xi) s2
	// The solution is	
	// s1(u = P1[0] + sqrt(xi)*P2[0] +/- sqrt(D)*P3[0] +/- sqrt(D*xi)*P4[0], 
	//		 v = -(P1[1] + sqrt(xi)*P2[1] +/- sqrt(D)*P3[1] +/- sqrt(D*xi)*P4[1]))
	// +sqrt(xi)*s2 (....)
	// 
	// = s1(u = P1[0], v=-P1[1]) + sqrt(xi)*s1(u = P2[0], v=-P2[1]) 
	//		 +/- sqrt(D)*s1(u = P3[0], v=-P3[1]) +/- sqrt(D*xi)*s1(u = P4[0], v=-P4[1]) 
	//	 + sqrt(xi)*s2(u = P1[0], v=-P1[1]) + xi*s2(u = P2[0], v=-P2[1]) 
	//		 +/- sqrt(D*xi)*s2(u = P3[0], v=-P3[1]) +/- sqrt(D)*xi*s2(u = P4[0], v=-P4[1]) 
	//
	// = s1(u = P1[0], v=-P1[1]) + xi*s2(u = P2[0], v=-P2[1]) 
	//		+ sqrt(xi)* (s1(u = P2[0], v=-P2[1]) + s2(u = P1[0], v=-P1[1]))
	//		+/- sqrt(D)*(s1(u = P3[0], v=-P3[1]) + xi*s2(u = P4[0], v=-P4[1]))
	//		+/- sqrt(D*xi)*(s1(u = P4[0], v=-P4[1]) + s2(u = P3[0], v=-P3[1]))

	// The output	 parameterized curve is thus : 
	//										 c = c1 + sqrt(xi). c2 + eps. sqrt(D). (c3 + sqrt(xi). c4)
	// c1 = s1(u = P1[0], v=-P1[1]) + xi*s2(u = P2[0], v=-P2[1]) 
	// c2 = s1(u = P2[0], v=-P2[1]) + s2(u = P1[0], v=-P1[1])
	// c3 = s1(u = P3[0], v=-P3[1]) + xi*s2(u = P4[0], v=-P4[1])
	// c4 = s1(u = P4[0], v=-P4[1]) + s2(u = P3[0], v=-P3[1])
	//	(If xi=1 then c2=c4=0)
	curve_param <bigint> c1(4), c2(4), c3(4), c4(4), c_tmp(4);
	hom_polynomial <bigint> htmp;
	bigint xi = det_q[1];
	if (!P1.is_zero()) 
		{
			negate(htmp, P1[1]);
			c1 = s1.eval(P1[0], htmp);
			c2 = s2.eval(P1[0], htmp);
		}
	if (!P2.is_zero())
		{
			negate(htmp, P2[1]);
			c_tmp = s2.eval(P2[0], htmp);	 
			multiply(c_tmp, c_tmp,	xi);	
			add(c1, c1, c_tmp);	 
			c_tmp = s1.eval(P2[0], htmp);	 
			add(c2, c2, c_tmp);
		}
	if (!P3.is_zero()) 
		{
			negate(htmp, P3[1]);
			c3 = s1.eval(P3[0], htmp);
			c4 = s2.eval(P3[0], htmp);
		}
	if (!P4.is_zero())
		{
			negate(htmp, P4[1]);
			c_tmp = s2.eval(P4[0], htmp);	 
			multiply(c_tmp, c_tmp,	xi);	
			add(c3, c3, c_tmp);	 
			c_tmp = s1.eval(P4[0], htmp);	 
			add(c4, c4, c_tmp);
		}

	// c = c1 + sqrt(xi). c2 + eps. sqrt(D). (c3 + sqrt(xi). c4)
	// Recall D != 0, 1, xi !=0

	bigint sq;
	if (xi == 1)	// c2 and c4 are 0
		{ 
			// the 2 conics are c1 + sqrt(D) c3 and c1 - sqrt(D) c3
			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize(c1, c3);

			ic.create_two_components(0,INTER_TYPE_CONIC,c1,c3,D);
			ic.set_optiflag(true);
		 }
	else if (is_square(sq, D*xi)) // recall neither D nor xi is 1 or 0
		{
			// the 2 conics are c = c1 + sqrt(xi). c2 + eps. f. sqrt(xi). (c3 + sqrt(xi). c4)
			// = (projective) xi. c2 + eps. sq. c3 + sqrt(xi). (c1 + eps. sq. c4)
			curve_param <bigint> c1p(4), c1m(4), c2p(4), c2m(4), tmp(4);
			multiply(c1p, c3, sq);
			multiply(c2p, c4, sq);
			QI::negate(c1m,c1p);
			QI::negate(c2m,c2p);
			multiply(tmp, c2, xi);
			add(c1p, c1p, tmp);
			add(c2p, c2p, c1);
			add(c1m, c1m, tmp);
			add(c2m, c2m, c1);

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize(c1p, c2p);
			optimize(c1m, c2m);

			ic.cc[0].create_component(INTER_TYPE_CONIC,c1p,c2p,xi);
			ic.cc[1].create_component(INTER_TYPE_CONIC,c1m,c2m,xi);
			ic.set_optiflag(true);
		}
	else
		{
			// Output parameterized curve : c = c1 + sqrt(xi). c2 + eps. sqrt(D). (c3 +
			// sqrt(xi). c4) 
			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize(c1, c2, c3, c4, s);

			ic.create_two_components(0,INTER_TYPE_CONIC,c1,c2,c3,c4,xi,D,0);
			ic.set_optiflag(false);
		}

	if ((rtype == 4) || (rtype == 6))
		{
			// Temp
			ic.cc[0].create_cut_parameters(-1,0,0,-1,0,0);
			ic.cc[1].create_cut_parameters(-1,0,0,-1,0,0);
		}
 
	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the two secant conics case, when the pair of planes is real
// and the conics cut in real space (s_e == 1)
// det_p = det_e * (a*x-b*y)^2, and delta_e is the discrimant of det_e.
// q_other is a matrix that does not correspond to the double root of the pencil
// (namely q1 or q2) q_other2 is matrix of inertia 22 (one exists) 
quad_inter <bigint> inter_secant_conics_sec(const bigint_matrix &q1, 
																						const bigint_matrix &q2, 
																						const hom_polynomial <bigint> &det_p, 
																						const bigint_matrix &q, 
																						const bigint_matrix &q_other, 
																						const bigint_matrix &q_sing, 
																						const bigint &delta_e, 
																						const int in_q, 
																						const int opt_level, ostream &s)
{
	quad_inter <bigint> ic;

	unsigned int ctype = 3, rtype;

	if (in_q == 2)
		rtype = 2;
	else
		if (delta_e < 0)
			rtype = 6;
		else
			rtype = 4;

	// Compute the two points of intersection of the conics, in all cases

	// Intersect the singular line of q with q_other
	bigint_matrix plug = prod( base_matrix<bigint> (prod(trans(q_sing), q_other)), q_sing);

	// Parameterize plug
	bigint D;
	bigint_matrix m1(2,1),m2(2,1);

	two_by_two_param(plug,D,m1,m2,opt_level,s);

	bigint_matrix p1(4,1),p2(4,1);
	if (D.is_zero())
		{
			// The two points are rational
			p1 = prod(q_sing, (m1+m2));
			p2 = prod(q_sing, (m1-m2));

			optimize_trans1(p1);
			optimize_trans1(p2);	
		}
	else
		{
			// The points are not rational
			 
			// Stack the transformations
			m1 = prod(q_sing, m1);
			m2 = prod(q_sing, m2);

			#ifndef NDEBUG
			s << ">> optimization of parameterization" << endl;
			#endif

			optimize_trans2(m1,m2);
		}

	if (in_q == 2) 
		{
			// Two components
			ic = quad_inter <bigint>(2);
			
			ic.set_type(ctype,rtype);

			#ifndef NDEBUG
			//print_type(ic,s);
			#endif

			if (D.is_zero())
				{
					curve_param <bigint> c1(column(p1, 0)),c2(column(p2, 0));
			
					ic.cc[0].create_component(INTER_TYPE_POINT,c1);
					ic.cc[1].create_component(INTER_TYPE_POINT,c2);
				}
			else
				{
		 			curve_param <bigint> p_rat(column(m1, 0)),p_D(column(m2, 0));

		 			ic.create_two_components(0,INTER_TYPE_POINT,p_rat,p_D,D);
				}

			ic.set_optiflag(true);
		}
	else // in_q == 1
		{	
			#ifndef NDEBUG
			if (D.is_zero())
				{
		 			curve_param <bigint> c1(column(p1, 0)),c2(column(p2, 0));

		 			s << ">> the two conics cut at " << c1 << " and " << c2 << endl;
				}
			else
				{
		 			curve_param <bigint> p_rat(column(m1, 0)),p_D(column(m2, 0)),c3(4),c4(4);
		 			hom_polynomial <bigint> tmp;

		 			s << ">> the two conics cut at ";
					//		 print_cp(p_rat,p_D,c3,c4,D,tmp,tmp,s);
		 			s << " and ";
					//		 print_cp(p_rat,-p_D,c3,c4,D,tmp,tmp,s);
		 			s << endl;
				}
			#endif

			// Parameterize the two planes
			bigint D2;
			bigint_matrix m1(4,3),m2(4,1);

			// q_sing is const: make a copy
			bigint_matrix q_sing_cp = q_sing;

			// We make this little change so that the intersection points are easy to
			// retrieve later on
			if (D.is_zero())
				{
					column(q_sing_cp, 0) = column(p1, 0);
					column(q_sing_cp, 1) = column(p2, 0);			
				}

			// Output is m1*[s*u t*u v] +/- sqrt(D2)*m2*[s*u] if D2 is non zero
			// Output is m1*[s*u t*u v] +/- m2*[s*u] otherwise
			pair_of_planes_param(q,q_sing_cp,D2,m1,m2,opt_level,s);

			if (!D2.is_zero())
				{
					#ifndef NDEBUG
		 			s << ">> planes are not rational" << endl;
					#endif
		 
					// The individual planes are not rational
		 			// We search for a quadric of inertia 22 in the pencil
		 			if (delta_e < 0) // det_p > 0 except at the double root => q_other is 22
			 			ic = inter_secant_conics_pl_not_rat(q1,q2,det_p,q_other,D2,
														m1,m2,rtype,opt_level,s);
		 			else // (delta_e > 0)
						// s_e > 0 (s_e = the value of det_e at the double root) otherwise 
						//	 inter_secant_conics_sec was called. 
						// Thus det_p > 0 in between the double root and its adjacent root of det_e
						// Moreover det_p > 0 at any value close enough to (but different
						// from) the double root, thus we can take the singular quadric (q) as
						// the quadric of inertia 22 because any quadric passing through a
						// point close enough to q will be of inertia 22. It's a bit sportive
						// but that should work 
						ic = inter_secant_conics_pl_not_rat(q1,q2,det_p,q,D2,m1,m2,rtype,opt_level,s);
				}
			else // if D2 = 0
				{
					#ifndef NDEBUG
		 			s << ">> planes are rational" << endl;
					#endif
			
					// Two components
					ic = quad_inter <bigint>(2);
					 
					ic.set_type(ctype,rtype);
			
					#ifndef NDEBUG
					//print_type(ic,s);
					#endif
			
					// Each individual plane is rational
					bigint_matrix m1p(4,3),m1m(4,3),ap(3,3),am(3,3);
					m1p = m1;
					column(m1p, 0) = column(m1, 0) + column(m2, 0);
					m1m = m1;
					column(m1m, 0) = column(m1, 0) - column(m2, 0);
			
					#ifndef NDEBUG
					s << ">> optimization of transformation" << endl;
					#endif
			
					optimize_trans1(m1p);
					optimize_trans1(m1m);
			
					// Now we compute the following two matrices: they represent our two conics
					ap = prod(base_matrix<bigint> (prod(trans(m1p),q_other)),m1p);
					am = prod(base_matrix<bigint> (prod(trans(m1m),q_other)),m1m);
			
					if (D.is_zero())
						{
							// We have a rational point on each conic, so we know how to find a
							// rational param
							
							// The two points of intersection correspond to (0 0 1) and (0 1 0) in
							// the space of the conics
							
							curve_param <bigint> parp(3),parm(3);
							bigint_matrix p_transp(3,3),p_transm(3,3);
							math_vector <bigint> lp(2,2),lm(2,2);
							
							conic_param_through_origin(ap,p_transp,parp,lp);		
							conic_param_through_origin(am,p_transm,parm,lm);		
							
							// Stack the transformations
							p_transp = prod(m1p, p_transp);
							p_transm = prod(m1m, p_transm);
							
							#ifndef NDEBUG
							s << ">> optimization of parameterization" << endl;
							#endif
						
							// Rescale so that l[0]*u+l[1]*v is replaced by v
							bigint_matrix scalep(3,3),scalem(3,3);
						
							scalep.sto(0,0,lp[1]*lp[1]);
							scalep.sto(1,0,lp[0]*lp[0]);
							scalep.sto(1,1,1);
							scalep.sto(1,2,-2*lp[0]);
							scalep.sto(2,0,-lp[0]*lp[1]);
							scalep.sto(2,2,lp[1]);
							
							scalem.sto(0,0,lm[1]*lm[1]);
							scalem.sto(1,0,lm[0]*lm[0]);
							scalem.sto(1,1,1);
							scalem.sto(1,2,-2*lm[0]);
							scalem.sto(2,0,-lm[0]*lm[1]);
							scalem.sto(2,2,lm[1]);
				 
							p_transp = prod(p_transp, scalep);
							p_transm = prod(p_transm, scalem);
							 
							optimize_trans3(p_transp,opt_level,s);
							optimize_trans3(p_transm,opt_level,s);
							 
							curve_param <bigint> conicp(4),conicm(4);
							multiply(conicp,p_transp,parp);
							multiply(conicm,p_transm,parm);		 
			
							ic.cc[0].create_component(INTER_TYPE_CONIC,conicp);
							ic.cc[1].create_component(INTER_TYPE_CONIC,conicm);
			
							ic.set_optiflag(true);
						}
					else // if D != 0
						{
							// The two points are not rational... but maybe we can still find a
							// rational param for the conics...
			
							// Parameterize the first conic
							bigint D3;
							bigint_matrix p_transp(3,3),p_trans2p(3,3);
							curve_param <bigint> par(3);
							 
							conic_param(ap,D3,p_transp,p_trans2p,par,opt_level,s);
							 
							if (D3.is_zero())
								{
									// Lucky you: you just found a rational conic!!
									
									// Stack the transformations
									p_transp = prod(m1p, p_transp);
						
									#ifndef NDEBUG
									s << ">> optimization of parameterization" << endl;
									#endif
						
									optimize_trans3(p_transp,opt_level,s);
									
									curve_param <bigint> conic(4);
									multiply(conic,p_transp,par);
						
									ic.cc[0].create_component(INTER_TYPE_CONIC,conic);
									ic.cc[0].optiflag = true;
								}
				 			else // if D3 != 0
								{
									// You'll have to content yourself with one square root...
									
									// Stack the transformations
									p_transp = prod(m1p, p_transp);
									p_trans2p = prod(m1p, p_trans2p);
						
									#ifndef NDEBUG
									s << ">> optimization of parameterization" << endl;
									#endif
						
									optimize_trans4(p_transp,p_trans2p,opt_level,s);
											
									curve_param <bigint> conic(4),conic2(4);
									multiply(conic,p_transp,par);
									multiply(conic2,p_trans2p,par);
						
									ic.cc[0].create_component(INTER_TYPE_CONIC,conic,conic2,D3);
						
									if (opt_level)
										ic.cc[0].optiflag = true;
									else
										ic.cc[0].optiflag = false;
								}

				 			// Parameterize the second conic
				 			bigint_matrix p_transm(3,3),p_trans2m(3,3);
				 			conic_param(am,D3,p_transm,p_trans2m,par,opt_level,s);

				 			if (D3.is_zero())
								{
									// Lucky you: you just found a rational conic!!
									
									// Stack the transformations
									p_transm = prod(m1m, p_transm);
						
									#ifndef NDEBUG
									s << ">> optimization of parameterization" << endl;
									#endif
						
									optimize_trans3(p_transm,opt_level,s);
									
									curve_param <bigint> conic(4);
									multiply(conic,p_transm,par);
						
									ic.cc[1].create_component(INTER_TYPE_CONIC,conic);
									ic.cc[1].optiflag = true;
								}
				 			else // if D3 != 0
								{
									// You'll have to content yourself with one square root...
									
									// Stack the transformations
									p_transm = prod(m1m, p_transm);
									p_trans2m = prod(m1m, p_trans2m);
						
									#ifndef NDEBUG
									s << ">> optimization of parameterization" << endl;
									#endif
											
									optimize_trans4(p_transm,p_trans2m,opt_level,s);
									
									curve_param <bigint> conic(4),conic2(4);
									multiply(conic,p_transm,par);
									multiply(conic2,p_trans2m,par);
						
									ic.cc[1].create_component(INTER_TYPE_CONIC,conic,conic2,D3);
						
									if (opt_level)
										ic.cc[1].optiflag = true;
									else
										ic.cc[1].optiflag = false;
								}				
			 			} // if D = 0 then else
				} // if D2 = 0 then else 

			// Temp
			ic.cc[0].create_cut_parameters(-1,0,0,-1,0,0);
			ic.cc[1].create_cut_parameters(-1,0,0,-1,0,0);	

		} // if in_q = 2 then else in_q = 1

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// Compute the intersection which is one conic lying in a non-rational plane.
// The non-rational plane of the pencil are parameterized by one of the two planes
// m1*[u v s] +/- sqrt(D)*m2*[u] (D != 0)
// D is here already optimized under -o
// q is a matrix that does not correspond to the double root of the pencil (namely
// q1 or q2) 
quad_inter <bigint> inter_one_conic_pl_not_rat(const bigint_matrix &q1, 
																								const bigint_matrix &q2, 
																								const hom_polynomial <bigint> &det_p, 
																								const bigint_matrix &q, const bigint &D, 
																								const bigint_matrix &m1, 
																								const bigint_matrix &m2_tmp, 
																								const int opt_level, ostream &s)
{
	// One component
	quad_inter <bigint> ic(1);

	ic.set_type(3,5);

	#ifndef NDEBUG
	//print_type(ic,s);
	#endif

	// Here the determinential polynomial is always negative execpt at the rational
	// root which corresponds to the pair of non-rational planes (parameterized by
	// m1, m2 and D).	 

	// Let T be the 4x3 parameterization matrix of the plane containing the conic,
	// namely, the parameterization is (m1[1] +/- sqrt(D)*m2)*u + m1[2]*v + m1[3]*s
	// and the matrix : its 1st collumn is m1[1] +/- sqrt(D)*m2, the 2nd and 3rd
	// collumns are	 m1[i=2, 3] 
	// For convenience we rename m2 to be a 4x3 matrix whose 2nd and 3rd collum are 0
	// Then T = m1 +/- sqrt(D)*m2
	bigint_matrix m2(4,3);
	column(m2, 0) = column(m2_tmp, 0);

	// The 3x3 matrix associated to the	 equation of the conic is 
	// C=T^t.Q.T where Q is any other (ie, not the pair of planes) quadric
	// of the pencil, namely here q.
	// One conic for say +sqrt(D) is only imaginary, and the other (for -+sqrt(D))
	// is the intersection curve we look for.	 

	// Let C=C1 +/- sqrt(D)*C2 
	//	 = (m1 + sqrt(D)*m2)*Q*(m1 + sqrt(D)*m2)
	//	 = m1*Q*m1 + D*m2*Q*m2 + sqrt(D)*(m1*Q*m2 + m2*Q*m1)
	bigint_matrix C1(3,3), C2(3,3), Ctmp(3,3), tmp(4,3);
	multiply(tmp,q,m1);
	multiply(C1, math_matrix<bigint> (trans(m1)), tmp);
	multiply(tmp,q,m2);
	multiply(Ctmp, math_matrix<bigint> (trans(m2)), tmp);
	multiply(Ctmp, Ctmp, D);
	add(C1, C1, Ctmp);

	multiply(tmp,q,m2);
	multiply(C2,math_matrix<bigint>(trans(m1)), tmp);
	multiply(tmp,q,m1);
	multiply(Ctmp,math_matrix<bigint>(trans(m2)), tmp);
	add(C2, C2, Ctmp);

	if (!inertia_non_rational_conic(C1, C2, D)) // inertia of C1+sqrt(D).C2 is [3,0]
		{
			C2.negate(); // inertia of C1+sqrt(D).C2 is now [2,1]
			m2.negate();
		}

	math_vector <bigint> D2(2,2);
	bigint_matrix p_trans0(3,3), p_trans1(3,3), p_trans2(3,3), p_trans3(3,3);
	curve_param <bigint> par(3);
	non_rational_conic_param(C1,C2,D,D2,p_trans0,p_trans1,p_trans2,
														p_trans3,par,opt_level,s);

	// The parameterization of the conic	in the plane is p*par where 
	// p = (p_trans0 + sqrt(D)*p_trans1 + sqrt(D2)*p_trans2 + sqrt(D)*sqrt(D2)*p_trans3) 
	// par = [u^2 v^2 uv],	D2 = D2[0] + sqrt(D)*D2[1]
	// Thus the parameterization of the conic	 in space is 
	// T*p*par where T be the 4x3 parameterization matrix of the plane containing
	// the conic. T = m1 + sqrt(D)*m2
	// T*p = (m1 + sqrt(D)*m2) * (p_trans0 + sqrt(D)*p_trans1 +
	// sqrt(D2)*p_trans2+sqrt(D)*sqrt(D2)*p_trans3)	 

	// T*p = (m1*p_trans0 + D*m2*p_trans1)	+ sqrt(D)*(m1*p_trans1+m2*p_trans0)
	// +sqrt(D2)*(m1*p_trans2 + D*m2*p_trans3)+sqrt(D)*sqrt(D2)*(m1*p_trans3+m2*p_trans2)
	// Let T*p = p0 + sqrt(D)*p1 + sqrt(D2)*p2 + sqrt(D)*sqrt(D2)*p3
	bigint_matrix p0(4,3), p1(4,3), p2(4,3), p3(4,3), p_tmp(4,3);
	multiply(p0, m1,p_trans0);
	multiply(p_tmp,m2,p_trans1);
	multiply(p_tmp,p_tmp, D);
	add(p0, p0, p_tmp);

	multiply(p1, m1,p_trans1);
	multiply(p_tmp,m2,p_trans0);
	add(p1, p1, p_tmp);

	multiply(p2, m1,p_trans2);
	multiply(p_tmp,m2,p_trans3);
	multiply(p_tmp,p_tmp, D);
	add(p2, p2, p_tmp);
	
	multiply(p3, m1,p_trans3);
	multiply(p_tmp,m2,p_trans2);
	add(p3, p3, p_tmp);

	curve_param <bigint> c0(4), c1(4), c2(4), c3(4);
	multiply(c0, p0, par);
	multiply(c1, p1, par);
	multiply(c2, p2, par);
	multiply(c3, p3, par);
	// The parameterization is = c0 + sqrt(D)*c1 + sqrt(D2)*c2 + sqrt(D)*sqrt(D2)*c3

	bigint gcd_tmp = gcd(cont(c0), cont(c1));
	gcd_tmp = gcd(gcd_tmp, cont(c2));
	gcd_tmp = gcd(gcd_tmp, cont(c3));
	divide(c0, c0, gcd_tmp);
	divide(c1, c1, gcd_tmp);
	divide(c2, c2, gcd_tmp);
	divide(c3, c3, gcd_tmp);

	ic.cc[0].create_component(INTER_TYPE_CONIC,c0,c1,c2,c3,D,D2[0],D2[1]);
	
	ic.set_optiflag(true);

	return ic;
}

//////////////////////////////////////////////////////////////////////////////////
// The procedure for the two secant conics case, when the pair of planes is real
// and the conics do not cut in real space (s_e == -1)
// det_p = det_e * (a*x-b*y)^2, and delta_e is the discrimant of det_e.
// q_other is a matrix that does not correspond to the double root of the pencil
// (namely q1 or q2) 
// q_other2 is matrix of inertia 22 if one exists in the pencil (that is when delta_e>0)
quad_inter <bigint> inter_secant_conics_no_sec(const bigint_matrix &q1, 
																								const bigint_matrix &q2, 
																								const hom_polynomial <bigint> &det_p, 
																								const bigint_matrix &q, 
																								const bigint_matrix &q_other, 
																								const bigint_matrix &q_other2, 
																								const bigint_matrix &q_sing, 
																								const bigint &delta_e,
																								const int opt_level, ostream &s)
{
	unsigned int ctype = 3, rtype;

	if (delta_e < 0) // Other roots of determinantal equation are complex
		rtype = 5;
	else
		rtype = 3;

	// The pair of planes is always real. Let us parameterize it
	bigint D;
	bigint_matrix m1(4,3),m2(4,1);

	// Compute the parameterization of a pair of planes
	// Output is m1*[s*u t*u v] +/- sqrt(D)*m2*[s*u] if D is non zero
	// Output is m1*[s*u t*u v] +/- m2*[s*u] otherwise
	pair_of_planes_param(q,q_sing,D,m1,m2,opt_level,s);

	if (!D.is_zero())
		{
			// The individual planes are not rational
			if (delta_e < 0) 
				return inter_one_conic_pl_not_rat(q1,q2,det_p,q_other,D,m1,m2,opt_level,s);
			else 
				return inter_secant_conics_pl_not_rat(q1,q2,det_p,q_other2,D,m1,m2,rtype,opt_level,s);
		}
	else // D = 0
		{
			// One or two components
			quad_inter <bigint> ic;

			if (delta_e < 0)
				ic = quad_inter <bigint>(1);
			else
				ic = quad_inter <bigint>(2);

			ic.set_type(ctype,rtype);

			#ifndef NDEBUG
			//print_type(ic,s);
			#endif

			// Each individual plane is rational
			bigint_matrix m1p = m1;
			column(m1p, 0) = column(m1, 0) + column(m2, 0);
			bigint_matrix m1m = m1;
			column(m1m, 0) = column(m1, 0) - column(m2, 0);

			#ifndef NDEBUG
			s << ">> optimization of transformation" << endl;
			#endif

			optimize_trans1(m1p);
			optimize_trans1(m1m);

			// Now we compute the following two matrices: they represent our two conics
			bigint_matrix ap = prod(base_matrix<bigint> (prod(trans(m1p), q_other)), m1p);
			bigint_matrix am = prod(base_matrix<bigint> (prod(trans(m1m), q_other)), m1m);

			// If delta_e < 0: one of the conics is imaginary
			if (delta_e < 0)
				{
					// Find the imaginary conic
					if (inertia(ap)[0] == 3)
						{
							// Swap ap and am
							swap(ap,am);
							swap(m1p,m1m);
						}
				}

			// Parameterize the first conic
			bigint_matrix p_trans(3,3),p_trans2(3,3);
			curve_param <bigint> par(3);

			conic_param(ap,D,p_trans,p_trans2,par,opt_level,s);

			if (D.is_zero())
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
			
					ic.cc[0].create_component(INTER_TYPE_CONIC,conic);
					ic.cc[0].optiflag = true;
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
					
					ic.cc[0].create_component(INTER_TYPE_CONIC,conic,conic2,D);
			
					if (opt_level)
						 ic.cc[0].optiflag = true;
					else
						 ic.cc[0].optiflag = false;
				}

			// Parameterize the second conic, if any
			if (delta_e > 0)
				{
					bigint_matrix p_trans(3,3),p_trans2(3,3);
					curve_param <bigint> par(3);
			
					conic_param(am,D,p_trans,p_trans2,par,opt_level,s);
			
					if (D.is_zero())
						{
							// Lucky you: you just found a rational conic!!
							 
							// Stack the transformations
							p_trans = prod(m1m, p_trans);
			
							#ifndef NDEBUG
							s << ">> optimization of parameterization" << endl;
							#endif
			
							optimize_trans3(p_trans,opt_level,s);
							 
							curve_param <bigint> conic(4);
							multiply(conic,p_trans,par);
							 
							ic.cc[1].create_component(INTER_TYPE_CONIC,conic);
							ic.cc[1].optiflag = true;
						}
		 			else
			 			{
							// You'll have to content yourself with one square root...
							 
							// Stack the transformations
							p_trans = prod(m1m, p_trans);
							p_trans2 = prod(m1m, p_trans2);
			
							#ifndef NDEBUG
							s << ">> optimization of parameterization" << endl;
							#endif
						
							optimize_trans4(p_trans,p_trans2,opt_level,s);
							 
							curve_param <bigint> conic(4),conic2(4);
							multiply(conic,p_trans,par);
							multiply(conic2,p_trans2,par);
							 
							ic.cc[1].create_component(INTER_TYPE_CONIC,conic,conic2,D);
			
							if (opt_level)
								ic.cc[1].optiflag = true;
							else
								ic.cc[1].optiflag = false;
			 			}
				}

			return ic;
		}
}

//////////////////////////////////////////////////////////////////////////////////
// The main procedure when only one real multiple root is found
quad_inter <bigint> inter_one_mult(const bigint_matrix &q1, const bigint_matrix &q2, 
																		const hom_polynomial <bigint> &det_p, 
																		const hom_polynomial <bigint> &det_p_orig,
																		const hom_polynomial <bigint> &gcd_p, 
																		const int opt_level, ostream &s)
{
	#ifndef NDEBUG
	s << ">> entering inter_one_mult" << endl;
	#endif
	
	quad_inter <bigint> ic;
	
	///////////////////////////////////////////////////////////////////////////////
	// Preliminaries: compute the multiple root, its multiplicity, the inertia of
	// the associated matrix and a few additional information
	///////////////////////////////////////////////////////////////////////////////

	// The degree of the gcd
	// d = 1: double, d = 2: triple, d = 3: quadruple
	rpl_size_t d = gcd_p.degree();

	// The (lambda,mu) of the multiple root
	math_vector <bigint> root(2,2);

	// A quadric different from the one used to parameterize
	bigint_matrix q_other;

	// The usual trick to avoid the case where the root is ``at infinity''
	if (!gcd_p[0].is_zero())
		{
			multiply(root[0],gcd_p[0],d);
			root[0].negate();
			root[1].assign(gcd_p[1]);
			q_other = q2;
		}
	else // Root is [0,1]
		{
			root[1].assign_one();
			q_other = q1;
		}		 

	// Simplify the coordinates of the root
	#ifndef NDEBUG
	s << ">> optimization of coordinates of root" << endl;
	#endif			

	optimize(root);

	// The matrix associated to the multiple root
	bigint_matrix q = root[0]*q1+root[1]*q2;

	// The singular locus of q (not needed in all cases, but eases readability of code)
	bigint_matrix q_sing = QI::singular(q);

	// Its inertia
	math_vector <int> in_q = inertia(q);

	// Its rank
	rpl_size_t rank_q = in_q[0]+in_q[1];
	
	// Compute the polynomial det_e such that det_p = det_e * (root[1]*x - root[0]*y)^(1+d)
	// Compute the sign s_e of the value of det_e at the multiple root
	hom_polynomial <bigint> det_e;
	rpl_size_t s_e;

	// det_e not really needed when quadruple root (s_e is easy to compute)
	if ((d == 1) || (d == 2)) // Double or triple real root
		{
			// First build the polynomial (root[1]*x - root[0]*y)^(1+d)
			hom_polynomial <bigint> lin;
			hom_polynomial <bigint> ptmp;
			ptmp.set_degree(1);
			ptmp[0] = -root[0];
			ptmp[1].assign(root[1]);

			// Polynomial det_e
			power(lin,ptmp,d+1);
			divide(det_e,det_p,lin);

			// Sign s_e
			s_e = (det_e.eval(root[0],root[1])).sign();
		}
	else // Quadruple real root - d == 3 (no need to compute det_e here, s_e is enough)
		if (root[1].is_zero())
			s_e = det_p[0].sign();
		else
			s_e = det_p[d+1].sign();

	///////////////////////////////////////////////////////////////////////////////
	// The main switch
	///////////////////////////////////////////////////////////////////////////////
	// d = 1: double real root case - The rank of the singular quadric can be
	// either 3 or 2
	///////////////////////////////////////////////////////////////////////////////
	if (d == 1) 
		{
			#ifndef NDEBUG
			s << ">> " << "double real root" << ": " << root << endl;
			s << ">> inertia: " << in_q << endl;
			#endif			

			// Discriminant delta_e of det_e
			bigint delta_e = discriminant2(det_e);

			///////////////////////////////////////////////////////////////////////////
			// Rank of singular quadric is 3
			if (rank_q == 3) // (C) Nodal quartic
				ic = inter_nodal_quartic(q1,q2,q,q_other,q_sing,delta_e,s_e,in_q[0],opt_level,s);
			///////////////////////////////////////////////////////////////////////////
			// Rank of singular quadric is 2
			else // (C) Two secant conics 
				if (delta_e < 0) { // Other roots of determinantal equation are complex
		 			if (s_e == 1)
			 			ic = inter_secant_conics_sec(q1,q2,det_p_orig,q,q_other,q_sing,
					 										delta_e,in_q[0],opt_level,s);
		 			else
			 			ic = inter_secant_conics_no_sec(q1,q2,det_p_orig,q,q_other,q_other,
							 								q_sing,delta_e,opt_level,s);
		}							
	else	 // Other roots of determinantal equation are real
		 if (in_q[0] == 1) // Pair of planes corresponding to double root is real
			 if (s_e == 1)
		ic = inter_secant_conics_sec(q1,q2,det_p_orig,q,q_other,q_sing,delta_e,
								in_q[0],opt_level,s);
			 else // s_e = -1
				 {
		// Decide in which case we are by testing for the presence of 
		// a [ 4 0 ] in the pencil

		// A test point 
		math_vector <bigint> root_e(2,2);

		if (det_e[2] < 0) // The double root is outside the simple roots
			{
				// Pick a point between the roots

				negate(root_e[0],det_e[1]);
				multiply(root_e[1],det_e[2],2);
			}
		else if (det_e[2] > 0) // The double root is between the simple roots
			{
				// Take infinity as test point

				root_e[0] = 1;
			}
		else // det_e[2] = 0
			if (det_e[1] > 0) // Double root is outside
				{
					root_e[0] = abs(det_e[0])+1;
					root_e[1] = det_e[1];
				}
			else
				{
					root_e[0] = -abs(det_e[0])-1;
					root_e[1] = -det_e[1];
				}

		// The associated quadric
		bigint_matrix q_det_positive = root_e[0]*q1+root_e[1]*q2;

		bigint rank0 = inertia_known_rank(q_det_positive,4)[0];

		if (rank0 == 2)
			ic = inter_secant_conics_no_sec(q1,q2,det_p_orig,q,q_other,
							q_det_positive,q_sing,delta_e,opt_level,s);
		else // rank0 == 4
			{
				ic = quad_inter <bigint>(0);

				ic.set_type(3,1);

										#ifndef NDEBUG
				//print_type(ic,s);
				#endif	
			}
				 }
		 else // Singular pair of planes is imaginary
			 if (s_e == 1) 
				 {
		ic = quad_inter <bigint>(0);

		ic.set_type(3,1);

								#ifndef NDEBUG
		//print_type(ic,s);
		#endif	
				 }
			 else // Two points
				 ic = inter_secant_conics_sec(q1,q2,det_p_orig,q,q_other,
							q_sing,delta_e,in_q[0],opt_level,s);
		}
	///////////////////////////////////////////////////////////////////////////////
	// d = 2: triple real root case - Rank of singular quadric can be either 3,
	// 2 or 1
	///////////////////////////////////////////////////////////////////////////////
	else if (d == 2)
		{
			#ifndef NDEBUG
			s << ">> " << "triple real root" << ": " << root << endl;
			s << ">> inertia: " << in_q << endl;
			#endif
			///////////////////////////////////////////////////////////////////////////
			// Rank of singular quadric is 3
			if (rank_q == 3) // (C) Cuspidal quartic
				ic = inter_cuspidal_quartic(q1,q2,q,q_other,q_sing,opt_level,s); 
			///////////////////////////////////////////////////////////////////////////
			// Rank of singular quadric is 2
			else if (rank_q == 2) // (C) Two tangent conics
				{
		 			// Here, q_other should be the cone of the pencil
		 
		 			q_other = -det_e[0]*q1+det_e[1]*q2;

		 			ic = inter_two_tangent_conics(q1,q2,q,q_other,q_sing,in_q[0],opt_level,s);
				}
			///////////////////////////////////////////////////////////////////////////
			// Rank of singular quadric is 1
			else // (C) Double conic
				ic = inter_double_conic(q1,q2,q,q_sing,det_e,opt_level,s);
		}
	///////////////////////////////////////////////////////////////////////////////
	// d = 3: quadruple real root case - Rank of singular quadric can be either 3,
	// 2, 1 or 0
	///////////////////////////////////////////////////////////////////////////////
	else
		{
			#ifndef NDEBUG
			s << ">> " << "quadruple real root" << ": " << root << endl;
			s << ">> inertia: " << in_q << endl;
			#endif

			///////////////////////////////////////////////////////////////////////////
			// Rank of singular quadric is 3
			if (rank_q == 3) // (C) Cubic and tangent line
				ic = inter_cubic_tangent_line(q1,q2,q,q_other,q_sing,opt_level,s);
			///////////////////////////////////////////////////////////////////////////
			// Rank of singular quadric is 2
			else if (rank_q == 2) 
				if (in_q[0] == 2) // Pair of planes is imaginary
		 			// Two skew lines and double line, only the double line is real
		 			ic = inter_two_skew_lines_double_line(q1,q2,q,q_other,q_sing,0,opt_level,s);
				else // Pair of planes is real
		 			if (s_e == 1)
			 			{
				 			// If the singular line of the pair of planes is entirely contained
				 			// in the other quadrics of the pencil, then 2 lines & double line.
				 			// Otherwise, conic & 2 lines crossing on conic
				 			bigint_matrix tmp = prod( base_matrix<bigint> (prod(trans(q_sing),q_other)), q_sing);

				 			if ((tmp(0,0).is_zero()) && (tmp(0,1).is_zero()) && (tmp(1,1).is_zero()))
								// Two skew lines and double line, everything is real
								ic = inter_two_skew_lines_double_line(q1,q2,q,q_other,q_sing,1,opt_level,s);
				 			else
								// Conic and two lines (crossing), everything is real
								ic = inter_conic_2lines_crossing(q1,q2,q,q_other,q_sing,1,opt_level,s);
			 			}
		 			else // s_e = -1
			 			// Conic and two lines (crossing), only the conic is real
			 			ic = inter_conic_2lines_crossing(q1,q2,q,q_other,q_sing,0,opt_level,s);
			///////////////////////////////////////////////////////////////////////////
			// Rank of singular quadric is 1
			else if (rank_q == 1) // (C) Two double lines
				ic = inter_two_double_lines(q1,q2,q,q_other,q_sing,s_e,opt_level,s);
			///////////////////////////////////////////////////////////////////////////
			// Rank of singular quadric is 0
			else // Any pencil quadric
				{
		 			// To avoid picking the zero matrix
		 			if (root[0].is_zero())
			 			q = q1;
		 			else
			 			q = q2;
			
					// Its inertia
					in_q = inertia_known_rank(q,4);
			
					if (in_q[0] == 4) // Inertia [4 0]
						{
							ic = quad_inter <bigint>(0);
			
							ic.set_type(11,1);
			
							#ifndef NDEBUG
							//print_type(ic,s);
							#endif	
						}
					else
						{
							ic = quad_inter <bigint>(1);
			
							if (in_q[0] == 3)
								ic.set_type(11,2);
							else
								ic.set_type(11,3);
			
							#ifndef NDEBUG
							//print_type(ic,s);
							#endif
			
							ic.cc[0].create_component(INTER_TYPE_SMOOTH_QUADRIC,q);
						}
				}
		}
		
	#ifndef NDEBUG
	s << ">> exiting inter_one_mult" << endl;
	#endif

	return ic;
}

} // end of namespace QI
