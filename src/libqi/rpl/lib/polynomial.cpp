#include <libqi/rpl/rpl.h>
#include <libqi/rpl/polynomial.h>

namespace rpl {

bigint cont(const polynomial <bigint> &a) {
	if (a.deg < 0) {
		return (bigint)0;
	}
	else if (a.deg == 0) {
		return a[0];
	}
	bigint gcd_c = gcd(a.coeff[0],a.coeff[1]); 
	for(rpl_size_t i = 2; i <= a.deg ; ++i) {
		gcd_c = gcd(gcd_c, a.coeff[i]); 
	}
	return gcd_c;
}

polynomial <bigint> pp(const polynomial <bigint> &a) {
	polynomial <bigint> res;
	if (!a.is_zero()) {
	  bigint c = cont(a);
		divide(res, a, c);
  }
	return res;
}

polynomial <bigint> gcd(const polynomial <bigint> &aa,
			 									const polynomial <bigint> &bb)
{
	if (bb.is_zero())
		return aa;

	polynomial <bigint> q, r, a = pp(aa), b = pp(bb);
	bigint cd = gcd(cont(aa), cont(bb));
	do {
		div_rem(q, r, a, b);
		a.assign(b);
		b.assign(pp(r));
	} while (b.deg >= 0);

	return cd * a;
}

};
