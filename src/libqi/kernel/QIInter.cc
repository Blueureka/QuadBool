// Main intersection loop
// Parsing of cases depending on number of multiple roots

// Bit of template instantiation
//#include "QITemplateInst.cc"
#include <libqi/kernel/QIInter.h>
#include <libqi/kernel/QIElem.h>
#include <libqi/kernel/QINumber.h>
#include <libqi/kernel/QIVanishDet.h>
#include <libqi/kernel/QINoMult.h>
#include <libqi/kernel/QIOneMult.h>
#include <libqi/kernel/QITwoMult.h>
#include <libqi/kernel/QIBench.h>

short QI_CPU_TIME_MS = 0;

using namespace rpl;

// Enter namespace QI
namespace QI {
	
// The main intersection procedure
quad_inter <bigint> intersection(const bigint_matrix &q1, const bigint_matrix &q2, 
																	const int opt_level, ostream &s)
{
	#ifndef NDEBUG
	s << ">> entering main intersection procedure" << endl;
	#endif
	
	quad_inter <bigint> result;
	int start, end;

	start = qi_bench_cputime();
	
	hom_polynomial <bigint> det_p = det_pencil(q1,q2);

	if (det_p.is_zero())
		{
			#ifndef NDEBUG
		 	s << ">> determinantal equation vanishes" << endl;
			#endif
		 result =	 inter_vanish_det(q1,q2,det_p,opt_level,s);
		}
	else
		{
			// Keep a copy of the original determinantal equation for the (2,2) case
			hom_polynomial <bigint> det_p_orig = det_p;

			#ifndef NDEBUG
			s << ">> optimization of coefficients of determinantal equation" << endl;
			#endif

			optimize(det_p);	

			#ifndef NDEBUG
			s << ">> determinantal equation: ";
			det_p.print_verbose(s,'l','m');
			s << endl;
			#endif

			hom_polynomial <bigint> gcd_p = gcd(derivative(det_p,'x'),derivative(det_p,'y'));

			#ifndef NDEBUG
			s << ">> optimization of coefficients of gcd" << endl;
			#endif		 
	
			optimize(gcd_p);

			#ifndef NDEBUG
			s << ">> gcd of derivatives of determinantal equation: ";
			gcd_p.print_verbose(s,'l','m');
			s << endl;
			#endif		

			if (gcd_p.degree() == 0) // No multiple root
				result = inter_no_mult(q1,q2,det_p,det_p_orig,opt_level,s);
			else if ((gcd_p.degree() == 1) || (gcd_p.degree() == 3) ||
				((gcd_p.degree() == 2) && (discriminant2(gcd_p).is_zero()))) // One mult. root
				result = inter_one_mult(q1,q2,det_p,det_p_orig,gcd_p,opt_level,s);
			else // Two double roots
				result = inter_two_mult(q1,q2,det_p,det_p_orig,gcd_p,opt_level,s);
		 }

	end = qi_bench_cputime();

	QI_CPU_TIME_MS = end - start;
	
	#ifndef NDEBUG
	s << ">> exiting main intersection procedure" << endl;
	#endif
	
	return result;
 }

} // end of namespace QI
