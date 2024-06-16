// Intersection when the determinantal equation has no multiple root

#include <libqi/kernel/QINoMult.h>
#include <libqi/kernel/QIUspensky.h>
#include <libqi/kernel/QISolve.h>
#include <libqi/kernel/QINumber.h>
#include <libqi/kernel/QIElem.h>
#include <libqi/kernel/QIInterWith22.h>

using namespace rpl;

// Enter namespace QI
namespace QI {

quad_inter <bigint> inter_no_mult(const bigint_matrix &q1, const bigint_matrix &q2, 
																	const hom_polynomial <bigint> &det_p, 
																	const hom_polynomial <bigint> &det_p_orig, 
																	const int opt_level, ostream &s)
{
	#ifndef NDEBUG
	s << ">> Entering inter_no_mult" << endl;
	#endif
	
	unsigned long deg = det_p.degree();
	interval *roots;

	unsigned int nbroot, nbroot_inf, case_flag = 0;
	bigint_matrix q;

	bool infinity_flag = 0;

	// If infinity is a solution in lambda
	if (det_p[4].is_zero())
		{
			infinity_flag = 1;
			deg = 3;
		}

	Uspensky(&roots,det_p,deg,nbroot);

	nbroot_inf = nbroot;

	// Number of roots + infinity if applicable
	if (infinity_flag)
		nbroot_inf++;

	#ifndef NDEBUG
	s << ">> number of real roots: " << nbroot_inf << endl;

	if (nbroot_inf != 0)
		{
			s << ">> intervals: ";
			affiche_roots(roots,nbroot,s);
			if (infinity_flag)
				s << ", infinity";
			s << endl;
		}
	#endif

	if (nbroot_inf == 0)
		{
			case_flag = 2;

			q = q1;
		}
	else // Two or four real roots
		{
			if (nbroot_inf == 2)
	 			case_flag = 1;
				 
			// Compute a quadric of positive determinant associated to a value either to
			// the left of the first root of the determinential equation or in between
			// the two first roots 
			hom_polynomial <bigint> deriv;
			deriv.assign(derivative(det_p,'x'));

			rpl_size_t offset = 0;
			// Pick first point to left of first root
			math_vector <bigint> p;
			p = pick_point_outside_roots(roots,infinity_flag,det_p,deriv,nbroot_inf,0);

			optimize(p);

			// The determinantal equation is negative at that point: pick the next one
			if ((det_p.eval(p[0],p[1])).sign() == -1)
				{
		 			offset = 1;
		 			p = pick_point_outside_roots(roots,infinity_flag,det_p,deriv,nbroot_inf,1);

		 			optimize(p);
				}

			#ifndef NDEBUG
			s << ">> picked test lambda 1 at " << p << ", sign ";
			if (det_p.eval(p[0],p[1]) > 0)
				s << "> 0";
			else
				s << "< 0";

			if (nbroot_inf == 4)
				s << " -- ";
			else
				s << " -- inertia [ 2 2 ] found" << endl;
			#endif

			// The matrix associated to p
			bigint_matrix q1_tmp;
			multiply(q1_tmp, q1, p[0]);
			bigint_matrix q2_tmp;
			multiply(q2_tmp, q2, p[1]);
			add(q, q1_tmp,q2_tmp);
			
			if (nbroot_inf == 4) // Four real roots
				{
					// Compute the inertia of q
					math_vector <int> in_q;
					in_q = inertia(q);
			
					#ifndef NDEBUG
					s << "inertia " << in_q << " found" << endl;
					#endif
			
					if (in_q[0] == 4) // [4 0] found: empty
						{
							quad_inter <bigint> ic(0);
			
							ic.set_type(1,1);
			
							#ifndef NDEBUG
							// print_type(ic,s);
							#endif
			
							for (unsigned int i = 0; i < nbroot; i++)
					 			mpz_clear(roots[i].c);
							delete[] roots;
			
							return ic;
						}
					else // Not enough information to decide: pick a second point
						{
							// Second test lambda
							math_vector <bigint> p2;
							p2 = pick_point_outside_roots(roots,infinity_flag,det_p,
																	deriv,nbroot_inf,2+offset);
		
							optimize(p2);
		
							#ifndef NDEBUG
							s << ">> picked test lambda 2 at " << p2 << ", sign ";
							if (det_p.eval(p2[0],p2[1]) > 0)
				 				s << "> 0";
							else
				 				s << "< 0";
							s << " -- ";
							#endif
		
							// The matrix associated to p2
							multiply(q1_tmp, q1, p2[0]);
							multiply(q2_tmp, q2,p2[1]);
							bigint_matrix qp2;
							add(qp2, q1_tmp,q2_tmp);
		
							// Its inertia
							in_q = inertia(qp2);
		
							#ifndef NDEBUG
							s << "inertia " << in_q << " found" << endl;
							#endif
		
							if (in_q[0] == 4) // [4 0] found: empty
								{
									quad_inter <bigint> ic(0);
					
									ic.set_type(1,1);
					
									#ifndef NDEBUG
									// print_type(ic,s);
									#endif
					
									for (unsigned int i = 0; i < nbroot; i++)
										mpz_clear(roots[i].c);
									delete[] roots;
					
									return ic;
								}
							else // not empty
								 {
								 	case_flag = 0;
					
								 	// choose the test point (and associated quadric) of smallest height
								 	if (sum_of_squares(p2) < sum_of_squares(p))
									 	{
											#ifndef NDEBUG
											s << ">> starting with lambda " << p2 << endl;
											#endif
					
											q = qp2;
										 }
									else
										{
											#ifndef NDEBUG
											s << ">> starting with lambda " << p << endl;
											#endif
										}
								}
						}	 // matches else // Not enough information to decide: pick a second point
				}	// matches	if (nbroot_inf == 4)	// Four real roots
		}		// matches else // Two or four real roots

	for (unsigned int i = 0; i < nbroot; i++)
		mpz_clear(roots[i].c);
	delete[] roots;

	#ifndef NDEBUG
	s << ">> exiting inter_no_mult" << endl;
	#endif
	
	// The determinantal equation passed is the one with unoptimized coefficients!
	return intersection_with_quadric22(q1,q2,det_p_orig,q,case_flag,opt_level,s);
}

} // end of namespace QI
