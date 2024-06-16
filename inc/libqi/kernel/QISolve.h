// Polynomial system solving

#ifndef _qi_solve_h_
#define _qi_solve_h_

/** rpl */
#include <libqi/rpl/math_vector.h>

/** QI */
#include "QIUspensky.h"
#include "QIHompoly.h"

using namespace std;
using namespace rpl;

// Enter namespace QI
namespace QI {

void Uspensky(interval **roots, const hom_polynomial <bigint> &det_p, const unsigned long deg, 
	      unsigned int &nbroot);

void affiche_roots(const interval *roots, const unsigned int nbroot, std::ostream &s);

math_vector <bigint> pick_point_outside_roots(interval *roots, const bool &infinity_flag, 
					      const hom_polynomial <bigint> &det_p, 
					      const hom_polynomial <bigint> &deriv, 
					      const unsigned int nbroot, const unsigned int index);

} // end of namespace QI

#endif
