#ifndef _RATIONAL_FACTORIZATION_H_
#define _RATIONAL_FACTORIZATION_H_

#include <libqi/rpl/bigint.h>
#include <vector>
#include <iostream>

// \file rational_factorization.h

namespace rpl {
//! rational_factorization 
/*! this class implements the decomposition
 *  in prime factors for bigints of small size
*   At the prototype stage, we only implements
*   the trialdiv algorithm (an Eratosthen sieve) 
*   with prime number knows until 10000
*   the factorisation is implemented trough
*   the factor method and could also be called
*   trough the trialdiv method
*   In mathematical terms we want to compute for a bigint \f$n\f$ 
*   the following decomposition 
*   \f[
*     n = \prod_{i}^{d} b_i^{e_i}
*   \f]
*   with \f$b\f$ is the bases vector and \f$e\f$ is the exponent vector
*   note that in case of partial factorisation (typically if we limit 
*   the maxfactor with calling trialdiv )
*   the factor $b_d$ is not prime, and it is called cofactor 
*/
class rational_factorization {
   public:
    rational_factorization(const bigint & n = 1)
    	: number(n),
	  cofactor(n), 
	  bases(0),
	  expos(0),
	  maxfactor(10000) {} 
   
    void assign(const bigint & n) { number = n; cofactor  = n; }

    size_t no_of_comp() const { return bases.size(); }

    size_t exponent(size_t i) const { return expos[i]; }
   
    bigint base(size_t i) const { return bases[i]; }

    bigint cofact() const { return cofactor; } 
    
    void trialdiv(const bigint & maxfact) { maxfactor = maxfact; factor(); }

    void factor();

    ~rational_factorization() {}	

    friend std::ostream & operator << (std::ostream & o , const rpl::rational_factorization & r); 
   
   private : 
 
    bigint number; /*!< the bigint we want to factorize*/
    bigint cofactor; 
    std::vector<bigint> bases; /*<! the bases vector b */
    std::vector<size_t> expos; /*<! the exponents vector e*/ 
    bigint maxfactor; /*!< currently 10000 or a smaller integer */

};


};
#endif
