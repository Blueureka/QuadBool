#ifndef _qi_uspensky_h_
#define _qi_uspensky_h_

/** GMP */
#include <gmp.h>

/** Structure used to store intervals by Uspensky */
typedef struct 
{
	mpz_t c; 
	long k; 
	unsigned int isexact; 
} interval; 

#endif
