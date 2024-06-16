// Polynomial system solving

/*
	Copyright 1999, 2002, 2004 
		Guillaume Hanrot (Guillaume.Hanrot@inria.fr)
		Fabrice Rouillier (Fabrice.Rouillier@inria.fr)
		Paul Zimmermann (Paul.Zimmermann@inria.fr)
		Sylvain Petitjean (Sylvain.Petitjean@inria.fr)

	This software is free software; you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published by
	the Free Software Foundation; either version 2.1 of the License, or (at your
	option) any later version.

	This software is distributed in the hope that it will be useful, but
	WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
	or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU Lesser General Public
	License for more details.

	If interested by a copy of the GNU Lesser General Public License, write to
	the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
	MA 02110-1301, USA. 
*/

/*
	 Uspensky's algorithm : isolation of real roots of a polynomial
	 based on Descartes' rule, following the procedure described by
	 Rouillier & Zimmermann.

	 Requires GMP.
*/

#include <libqi/kernel/QISolve.h>

using namespace rpl;

// Enter namespace QI
namespace QI {

#define vali(x) mpz_scan1((x), 0)
#define ilog2(a) mpz_sizeinbase(a,2)
#define TOT_POS -1 

/* Probably out of memory long since on nontrivial examples, but anyway. */
#define DEPTH_MAX 1024

static unsigned long usign;
static int b1 = 0, b2 = 0; 

int RemoveContent(mpz_t *P, const unsigned long deg) 
{
	unsigned long cont, i, z; 

	i = 0; while (mpz_sgn(P[i]) == 0) i++; 
	cont = vali(P[i]); 

	for( ; (i <= deg) && cont; i++)
		{
			if (mpz_sgn(P[i]) != 0) 
				{
					 z = vali(P[i]); 
					 if (z < cont) cont = z; 
				}
		}		 

	if (cont == 0) return 0; 

	for (i = 0; i <= deg; i++) 
		mpz_fdiv_q_2exp(P[i], P[i], cont);

	return cont;
}

/* Johnson's bound : 2*max(abs(-ai/an)^(1/(n-i)),
	 the maximum being taken over the i with sgn(a_i) != sgn(an) */
long bound_roots(mpz_t *t, const unsigned long deg)
{
	unsigned long i;
	long maxpow, currpow, currpow2, lan, tpos = 1; 

	currpow = currpow2 = 0;
	lan = ilog2(t[deg]) - 1; /* puiss de 2 < an */

	maxpow = -lan; 

	for (i = 0; i < deg; i++)
		{
			if (mpz_sgn(t[deg]) != mpz_sgn(t[i]))
				{
					tpos = 0; 
					currpow = ilog2(t[i]); 
					currpow -= lan; /* 2^currpow >= abs(-ai/an) */
			
					if (currpow > 0) 
						currpow2 = currpow / (deg - i); 
					else
						currpow2 = ((-currpow) / (deg - i)); 
			
					// Bug fix Sylvain in parentheses
					if (currpow2 * ((long) (deg - i)) != currpow) 
					 	currpow2++;
					/* 2^currpow2 >= abs(-ai/an)^(1/(n-i)) */
					 
					if (currpow2 > maxpow) maxpow = currpow2;
				}
		}		

	if (tpos == 1) return -1; 

	/* here 2^maxpow > max(abs(-ai/an)^(1/(n-i)), add one to get the bound */
	maxpow++; 

	return(maxpow);
}

/* Computes P(X*2^k) */
int Homoth(mpz_t *P, const long k, const unsigned long deg)
{
	unsigned long i; 
	long j;
	
	if (k > 0) { 
		j = k; 
		for (i = 1; i <= deg; i++, j += k)
			mpz_mul_2exp(P[i], P[i], j); 
	}
	else
		{
			j = deg * (-k); 
			for (i = 0; i < deg; i++, j += k)
				mpz_mul_2exp(P[i], P[i], j);			
		}

	/* Remove possible large power of 2 in content */
	return RemoveContent(P, deg);
}

/* Replaces P by the polynomial P(X+1) */
void X2XP1(mpz_t *P, const unsigned long deg)
{
	for (long i = 0; i <= (long)deg-1; i++)
		for (long j = deg-1 ; j >= i; j--)
			mpz_add(P[j], P[j], P[j+1]);
}

/* Number of sign changes in the coefficients of P(1/(X+1)) (Descartes' rule) */
unsigned long Descartes(mpz_t *P, const unsigned long deg, const long sigh, long &flag)
{
	unsigned long nb = 0; 
	long s, t; 
	unsigned long i, j;
	
	/* 
		 Prune the computation if all the coefficients are of the sign of P[deg] 
		 In that case any subsequent interval shall have the same property, 
		 we put flag at 1 to point this to Uspensky_rec.
	*/
	j = deg; t = mpz_sgn(P[j]); 
	while (j >= 0 && mpz_sgn(P[j]) == t) 
		{ 
			if (j == 0)
				break;

			j--; 
		}
		
	if (j < 0) 
		{	
			flag = -1; 
			return nb; 
		}		

	mpz_t *Q = new mpz_t[deg + 1];

	for (i = 0; i <= deg; i++) 
		mpz_init_set(Q[i], P[i]); 
	
	for (j = 0; j <= deg-1; j++) 
		mpz_add(Q[j+1], Q[j+1], Q[j]);

	s = mpz_sgn(Q[deg]); 

	flag = s && (s == mpz_sgn(P[0])) && (s == -sigh); 

	for (i = 1; i <= deg-1; i++)
		{
			/* 
				 Prune the computation if all further coefficients are of the sign of 
				 Q[deg-i] 
			*/
			j = deg - i; t = s; 
			while (j >= 0 && t == 0) { t = mpz_sgn(Q[j]); j--; }
			while (j >= 0 && mpz_sgn(Q[j]) == t) 
				{ 
		 			if (j == 0)
			 			break;
			
					j--; 
				}
			if (j < 0) 
				{	
					for (i = 0; i <= deg; i++) 
						mpz_clear(Q[i]); 
					delete[] Q; 
			
					return nb; 
				}	 
			
			for (j = 0; j <= deg - i - 1; j++)
				mpz_add(Q[j+1], Q[j+1], Q[j]); 
			
			if (s == 0) 
				s = mpz_sgn(Q[deg-i]); 
			else 
				if (s == -mpz_sgn(Q[deg-i])) 
					{ 
						if ((nb == 1 && !flag) || nb == 2)	 
							{ 
								for (i = 0; i <= deg; i++) 
									mpz_clear(Q[i]); 
								delete[] Q;
			
								return (nb + 1); 
							} 
						 
						nb++; s = -s; 
					}
			}

	if (s == -mpz_sgn(Q[0]))
		nb++; 

	for (i = 0; i <= deg; i++) 
		mpz_clear(Q[i]); 
	delete[] Q;

	return nb; 
}

/* Returns the sign of P(1/2) */
long evalhalf(mpz_t *P, const unsigned long deg)
{
	int ret;
	mpz_t x, y; 

	mpz_init_set(x, P[deg]); 
	mpz_init(y); 

	for (long j = deg - 1; j >= 0; j--) 
		{
			mpz_mul_2exp(y, P[j], deg - j); 
			mpz_add(x, x, y); 
		}

	mpz_clear(y);
	ret = mpz_sgn (x);
	mpz_clear(x);
	return ret;
}

void add_root(interval *roots, const mpz_t &c, const int k, const unsigned int flag, 
							const unsigned int nbroot)
{
	int b = (usign ? b1 : b2); 

	mpz_init(roots[nbroot].c); 

	if (k <= b) 
		{
			if (usign) 
				{ 
					mpz_neg(roots[nbroot].c, c); 
					// Bug fix Sylvain when root is exact
					if (!flag)
						mpz_sub_ui(roots[nbroot].c, roots[nbroot].c, 1); 
					mpz_mul_2exp(roots[nbroot].c, roots[nbroot].c, b-k); 
				}
			else 
				mpz_mul_2exp(roots[nbroot].c, c, b-k); 

			roots[nbroot].k = k - b; 
			roots[nbroot].isexact = flag; 
		}
	else
		{
			if (usign) 
				{ 
					mpz_neg(roots[nbroot].c, c);
				 	// Bug fix Sylvain when root is exact
				 	if (!flag)
				 		mpz_sub_ui(roots[nbroot].c, roots[nbroot].c, 1); 
				}		
			else 
				mpz_set(roots[nbroot].c, c); 

			roots[nbroot].k = k - b; 
			roots[nbroot].isexact = flag; 
		}
}

void affiche_root(const interval &z, ostream &s)
{
	mpz_t tmp; 
	mpz_init(tmp); 

	if (z.isexact != 1) { s << "]"; }

	if (z.k <= 0) 
		s << bigint(z.c);
	else 
		s << bigint(z.c) << "/2^" << z.k;

	if (z.isexact != 1) 
		{ 
			s << ", "; 

			if (z.k <= 0)			
				{ 
					mpz_set_ui(tmp, 1); 
					mpz_mul_2exp(tmp, tmp, -z.k); 
					mpz_add(tmp, z.c, tmp); 
					s << bigint(tmp);
				}
			else 
				{ 
					mpz_add_ui(tmp, z.c, 1); 
					s << bigint(tmp) << "/2^" << z.k;
				}
			
			s << "["; 
		}

	mpz_clear(tmp); 
}

void affiche_roots(const interval *roots, const unsigned int nbroot, ostream &s)
{
	for (unsigned int i = 0; i < nbroot; i++) 
		{
			affiche_root(roots[i],s); 
			if (i < nbroot - 1) s << ", "; 
		}
}

/* 
	 Check interval [c/2^k, (c+1)/2^k]. The value of k is returned, this
	 is necessary to know from where we come [i.e., what exactly is the 
	 current polynomial P] when several recursive calls return in a row. 
	 In practice, this is used to update the polynomial at HERE */
long Uspensky_rec(mpz_t *P, const mpz_t &c, const unsigned long k, unsigned long &Deg, 
									interval *roots, unsigned int &nbroot)
{
	unsigned long i, j, nb;
	long oldk;
	long shalf, flag; 
	mpz_t tmp; 

	if (k > DEPTH_MAX) 
		{
			cout << "Maximal depth reached. Check that your polynomial is " 
						<< "squarefree or increase DEPTH_MAX" << endl;
			exit(-1); 
		}

	mpz_init(tmp); 

	/* Check whether c/2^k is a root */
	if (mpz_cmp_ui(P[0], 0) == 0)
		{
			i = 1; while(mpz_cmp_ui(P[i], 0) == 0) { i++; }
			 
			for (j = 0; j < i; j++) 
			 	{
			 		add_root(roots, c, k, 1, nbroot); 
			 		nbroot++; 
			 	}
		 
			Deg -= i; /* Update the polynomial */		
			for (j = 0; j <= Deg; j++, i++)
 				mpz_set(P[j], P[i]); 
		}
			
	/* 
		 Compute the sign of P(1/2) ; thus if Descartes bound is 2, 
		 whereas sign(P(0)) = sign(P(1)) = -sign(P(1/2)) we have
		 found two roots. 
	*/
	shalf = evalhalf(P, Deg); 

	/* Compute Descartes' bound */
	nb = Descartes(P, Deg, shalf, flag);
	if (flag == TOT_POS) 
		{
			mpz_clear(tmp); 
			return TOT_POS; 
		}
			
	switch (nb)
		{
			case 0: /* no root */
				mpz_clear(tmp); 
				return k; 
	
			case 1: /* exactly one root */
				add_root(roots, c, k, 0, nbroot); 
				nbroot++; 
				mpz_clear(tmp); 
				return k; 
			 
			case 2: /* if flag != 0, one root in each half of the current interval */
				if (flag) 
					{
						mpz_set(tmp, c); 
				
						mpz_mul_2exp(tmp, tmp, 1); 
						add_root(roots, tmp, k+1, 0, nbroot); 
						nbroot++; 
				
						mpz_add_ui(tmp, tmp, 1); 
						add_root(roots, tmp, k+1, 0, nbroot); 
						nbroot++; 
				
						mpz_clear(tmp); 
						return k; 
					}
	
			default: /* recursive call on each half of the interval */
				mpz_set(tmp, c); 
	
				mpz_mul_2exp(tmp, tmp, 1);
				Homoth(P, -1, Deg); 	 
				oldk = Uspensky_rec(P, tmp, k+1, Deg, roots, nbroot); 
				if (oldk == TOT_POS) 
					{
						mpz_clear(tmp); 
						return TOT_POS; 
					}
			 
				mpz_add_ui(tmp, tmp, 1); 
				X2XP1(P, Deg); 
	
				if (oldk > (long)k + 1)
					Homoth(P, oldk - (k + 1), Deg);
	
				oldk = Uspensky_rec(P, tmp, k+1, Deg, roots, nbroot);
				if (oldk == TOT_POS) 
					{
						mpz_clear(tmp); 
						return TOT_POS; 
					}
				
				mpz_clear(tmp); 
				return oldk; 
		}
}

void Uspensky_gmp(interval *roots, mpz_t *Q, const unsigned long deg, unsigned int &nbroot)
{
	unsigned long deg1 = deg;
	long i, j;

	mpz_t e; 
	int nb_z;

	mpz_init_set_ui(e, 0); 
	nbroot = 0; 

	/* Remove 0 roots in all cases */
	nb_z = 0; while (mpz_sgn(Q[nb_z]) == 0) nb_z++; 
	for (j = 0; j < nb_z; j++) 
		{
			add_root(roots, e, 0, 1, nbroot); 
			nbroot++; 
		}

	deg1 = deg - nb_z; /* Update the polynomial */
		 
	mpz_t *P = new mpz_t[deg1 + 1]; 

	for (j = 0; j <= (long)deg1; j++) 
		mpz_init_set(P[j], Q[nb_z + j]); 

	/* First work on the positive roots. */
	b2 = bound_roots(P, deg1);

	if (b2 < 0) goto NEGATIVE; 

	Homoth(P, b2, deg1); 

	usign = 0; 
	Uspensky_rec(P, e, 0, deg1, roots, nbroot); 
	
	/* Change P into P(-X) to look for negative roots */
 	NEGATIVE:
	deg1 = deg - nb_z; 
	for (i = deg1; i >= 0; i--)
		{
			if (i % 2 == 1) 
				mpz_neg(P[i], Q[nb_z + i]); 
			else mpz_set(P[i], Q[nb_z + i]); 
		}

	b1 = bound_roots(P, deg1);

	if (b1 >= 0)
		{
			Homoth(P, b1, deg1); 
			mpz_set_ui(e, 0); 
			usign = 1; 
			Uspensky_rec(P, e, 0, deg1, roots, nbroot);		 
		}

	/* Free memory. */
	for (i = deg-nb_z; i >= 0; i--) 
		mpz_clear(P[i]); 
	delete[] P;

	mpz_clear(e);
}

// The title says it all
math_vector <bigint> pick_point_outside_roots(interval *roots, const bool &infinity_flag, 
																							const hom_polynomial <bigint> &det_p, 
																							const hom_polynomial <bigint> &deriv, 
																							const unsigned int nbroot, const unsigned int index)
{
	math_vector <bigint> point(2,2);

	if (index == 0)
		{
			// Only root is infinity
			if ((infinity_flag) && (nbroot == 1))
				point[0].assign_one();
			else
			{
				if (roots[0].k <= 0)
					point[1].assign_one();
				else
					power(point[1],2,roots[0].k);
		
				if (roots[0].isexact)
					// Take lower end of interval minus 1 as point
					subtract(point[0],bigint(roots[0].c),1);
				else
					// Take lower end of interval as point
					point[0].assign(bigint(roots[0].c));
			}

			return point;
		}
	else if ((index == nbroot-1) && (infinity_flag))
		{
			// Take upper end of interval plus 1 as point
			if (roots[index-1].k > 0)
				{
					add(point[0],bigint(roots[index-1].c),2);
					power(point[1],2,roots[index-1].k);
				}
			else
				{
					power(point[0],2,-roots[index-1].k);
					add(point[0],point[0],bigint(roots[index-1].c));
					add(point[0],point[0],1);
					point[1].assign_one();
				}
			return point;
		}
	else
		{
			math_vector <bigint> lower_end(2,2);
			math_vector <bigint> higher_end(2,2);
	
			if (roots[index-1].isexact)
				{
					lower_end[0].assign(bigint(roots[index-1].c));
			
					if (roots[index-1].k <= 0)
						lower_end[1].assign_one();
					else
						power(lower_end[1],2,roots[index-1].k);
				}
			else
				if (roots[index-1].k <= 0)
					{
						power(lower_end[0],2,-roots[index-1].k);
						add(lower_end[0],lower_end[0],bigint(roots[index-1].c));
						lower_end[1].assign_one();
					}
				else
					{
						add(lower_end[0],bigint(roots[index-1].c),1);
						power(lower_end[1],2,roots[index-1].k);
					}
			
			higher_end[0].assign(bigint(roots[index].c));
			if (roots[index].k <= 0)
				higher_end[1].assign_one();
			else
				power(higher_end[1],2,roots[index].k);

			// Are the two endpoints the same?
			bigint tmp1,tmp2,equal_end;
			multiply(tmp1,lower_end[0],higher_end[1]);
			multiply(tmp2,lower_end[1],higher_end[0]);			
			subtract(equal_end,tmp1,tmp2);

			if (!equal_end.is_zero())
				// They are different: pick a point in the middle
				if (lower_end[1] == higher_end[1]) // No need to find common denominator
					{
						add(point[0],lower_end[0],higher_end[0]);
						multiply(point[1],lower_end[1],2);
						 
						return point;
					}
				else
					{
						multiply(tmp1,higher_end[0],lower_end[1]);
						multiply(tmp2,higher_end[1],lower_end[0]);
						add(point[0],tmp1,tmp2);
			
						multiply(tmp1,lower_end[1],higher_end[1]);
						multiply(point[1],tmp1,2);
					}
				else if ((!roots[index-1].isexact) && (!roots[index].isexact))
					// None of the two roots are exact: pick the higher_end
					{
		 				point.assign(higher_end);

		 				return point;
					}
				else
					// One of the two roots is exact and equal to the endpoint of the other interval
					{
						if (roots[index-1].isexact) // Left is exact, refine right interval
							{
								// Exact root
								math_vector <bigint> ex_root(2,2);
								ex_root[0].assign(bigint(roots[index-1].c));
								if (roots[index-1].k <= 0)
									ex_root[1].assign(1);
								else
									power(ex_root[1],2,roots[index-1].k);
				
								// Evaluate sign of derivative at exact root
								rpl_size_t sg_der = (deriv.eval(ex_root[0],ex_root[1])).sign();
				
								// Evaluate at midpoint of interval and loop by cutting in half
								// until the right sign (sg_der) is reached.
								bigint l_1,m;
								if (roots[index].k > 0)
									{
										l_1.assign_one();
										power(m,2,roots[index].k+1);
									}
								else
									{
										power(l_1,2,-roots[index].k);
										m.assign(2);
									}
				
								bigint l_0;
								multiply(l_0,bigint(roots[index].c),2);
				
								add(point[0],l_0,l_1);
								point[1].assign(m);

				 				while ((det_p.eval(point[0],point[1])).sign() != sg_der)
									// Recurse
									{
										multiply(l_0,l_0,2);
										multiply(m,m,2);
										
										add(point[0],l_0,l_1);
										point[1].assign(m);			 
									}
							}
		 				else // Right is exact, refine left interval
			 				{
				 				// Exact root
								math_vector <bigint> ex_root(2,2);
								ex_root[0].assign(bigint(roots[index].c));
								if (roots[index].k <= 0)
									ex_root[1].assign(1);
								else
									power(ex_root[1],2,roots[index].k);
				
								// Evaluate sign of derivative at exact root
								rpl_size_t sg_der = (deriv.eval(ex_root[0],ex_root[1])).sign();
				
								// Evaluate at midpoint of interval and loop by cutting in half
								// until the right sign (sg_der) is reached.
								bigint l_0,l_1,m;
								if (roots[index-1].k > 0)
									{
										add(l_0,bigint(roots[index-1].c),1);
										multiply(l_0,l_0,2);
										l_1.assign_one();
										l_1.negate();
										power(m,2,roots[index-1].k+1);
									}
								else
									{
										power(l_0,2,-roots[index-1].k);
										l_1.assign(l_0);
										l_1.negate();
										add(l_0,l_0,bigint(roots[index-1].c));
										multiply(l_0,l_0,2);
										m.assign(2);
									}
							
							 	add(point[0],l_0,l_1);
							 	point[1].assign(m);
			
							 	while ((det_p.eval(point[0],point[1])).sign() != -sg_der)
									// Recurse
									{
										multiply(l_0,l_0,2);
										multiply(m,m,2);
										
										add(point[0],l_0,l_1);
										point[1].assign(m);			 
									}
			 				}
					}
			}

	return point;
}

// Output things in the right order
void clean_output(interval *roots, const unsigned long deg, const unsigned int nbroot)
{
	if ((nbroot != 0) && (mpz_sgn(roots[nbroot-1].c) < 0))
		{
			interval *roots2 = new interval[deg]; 

			int j = 0;

			for (int i = nbroot-1; i >= 0; i--)
				{
		 			if (mpz_sgn(roots[i].c) >= 0)
			 			break;

					mpz_init_set(roots2[j].c,roots[i].c);
					roots2[j].k = roots[i].k;		 
					roots2[j].isexact = roots[i].isexact;
					j++;
				}

			for (unsigned int i = 0; i < nbroot-j; i++)
				{
					 mpz_init_set(roots2[i+j].c,roots[i].c);
					 roots2[i+j].k = roots[i].k;		 
					 roots2[i+j].isexact = roots[i].isexact;
				}	

			// Copy roots2 in roots and delete roots2
			for (unsigned int i = 0; i < nbroot; i++)
				{
					 mpz_set(roots[i].c,roots2[i].c);
					 mpz_clear(roots2[i].c);
					 roots[i].k = roots2[i].k;		 
					 roots[i].isexact = roots2[i].isexact;
				}	

			delete[] roots2;
		}
}

// Convert between rpl and gmp format, then launch Uspensky
void Uspensky(interval **roots, const hom_polynomial <bigint> &det_p, const unsigned long deg, 
							unsigned int &nbroot)
{ 
	// Initialize the roots structure
	*roots = new interval[deg];

	mpz_t *P = new mpz_t[deg + 1];

	for (unsigned long i = 0; i <= deg; i++)
		mpz_init_set(P[i], det_p[i].val.get_mpz_t());

	// Do the actual computation of the roots
	Uspensky_gmp(*roots,P,deg,nbroot);

	/* Free memory. */
	for (unsigned long i = 0; i <= deg; i++) 
		mpz_clear(P[i]); 
	delete[] P;

	/* Now a bit of cleaning. */
	clean_output(*roots, deg, nbroot);
}

} // end of namespace QI
