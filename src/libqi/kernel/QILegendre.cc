// Decide if a conic with integer coefficiens has a rational solution and compute
// such a solution

// Based on the paper
//		"Solving quadratic equations using reduced unimodular quadratic forms"
// by Denis Simon
// Math. Comp. 74 (2005), 1531-1543.
// http://www.ams.org/mcom/2005-74-251/

// Written by Didier Gemmerle, with slight adjustments by Sylvain Petitjean

// There are two parts: one deals with general smooth conics, the other with Legendre equations

#include <libqi/rpl/gf_element.h> 
#include <libqi/rpl/factorization.h>

#include "QILegendre.h"
#include "QINumber.h"

using namespace rpl;

// Enter namespace QI
namespace QI {

// Compute the adjoint of a bigrational matrix
math_matrix <bigrational> Adjoint(const math_matrix <bigrational> &M)
{
	math_matrix <bigrational> HM(3,3);

	HM.sto(0,0,M(1,1)*M(2,2)-M(2,1)*M(1,2));
	HM.sto(0,1,M(0,2)*M(2,1)-M(0,1)*M(2,2));
	HM.sto(0,2,M(0,1)*M(1,2)-M(0,2)*M(1,1));
	HM.sto(1,0,M(1,2)*M(2,0)-M(1,0)*M(2,2));
	HM.sto(1,1,M(0,0)*M(2,2)-M(2,0)*M(0,2));
	HM.sto(1,2,M(0,2)*M(1,0)-M(0,0)*M(1,2));
	HM.sto(2,0,M(1,0)*M(2,1)-M(1,1)*M(2,0));
	HM.sto(2,1,M(0,1)*M(2,0)-M(0,0)*M(2,1));
	HM.sto(2,2,M(0,0)*M(1,1)-M(1,0)*M(0,1));

	return(HM);
}

// Convert a vector of rationals to a vector of integers
math_vector <bigint> RatToBigint(const math_vector <bigrational> &ratvec)
{
	math_vector <bigint> out(3,3);

	bigint mu = lcm(denominator(ratvec[0]),lcm(denominator(ratvec[1]),denominator(ratvec[2])));

	out[0] = numerator(ratvec[0]*mu);
	out[1] = numerator(ratvec[1]*mu);
	out[2] = numerator(ratvec[2]*mu);

	return out;
}

// Compute LLL for indefinite quadratic forms
//					(Algorithms 1.1 and 1.3 in D. Simon.)
// Either outputs a vector or a matrix of solutions
// The LLL constant is 3/4 here and by default
// Return code is the dimension of the solution found
lidia_size_t LLL(const bigint_matrix &G, math_vector <bigrational> &vecteurSolution, 
		 							bigint_matrix &matriceSolution, bigint_matrix &H,
		 							const bigrational &constanteLLL = bigrational(3,4))
{
	// Make a copy of the matrix G, on which we shall work
	math_matrix <bigrational> A(3,3);	 

	// The Mu coefficients in the Gram algo
	math_matrix <bigrational> M(3,3);		

	// Set to identity matrix
	M.diag(1,0);

	// Compute Gram-Schmidt	 
	for (lidia_size_t i = 0; i < 3; i++)
		{
			A.sto(i,i,G(i,i));
			if (A(i,i) == 0)
				{
		 			vecteurSolution = M.get_column_vector(i);

		 			return(1);
				}

			for (lidia_size_t j = 0; j < i ; j++)
				{
		 			A.sto(i,j,G(i,j));
		 			for (lidia_size_t k = 0; k < j; k++)
			 			A.sto(i,j,A(i,j)-M(j,k)*A(i,k));

		 			M.sto(i,j,A(i,j)/A(j,j));
		 			A.sto(i,i,A(i,i)-M(i,j)*A(i,j));
		 			if (A(i,i) == 0)
			 			{
				 			// Transpose of Adjoint of M 
				 			math_matrix <bigrational> IM = trans(Adjoint(M)); 
 
				 			// Solution 
				 			vecteurSolution = IM.get_column_vector(i); 

				 			return(1);
			 			}
				}
		}

	// A is now updated Gram matrix
	// M is the matrix of Mu coefficients of Gram-Schmidt

	// Identity matrix
	H.diag(1,0);
	
	// LLL loop
	
	bool swap;
	lidia_size_t k = 1;
	bigint q;

	while (k <= 2)
		{		 
			swap = 1;
			while (swap != 0) 
				{
		 			swap = 0;
		 
		 			// red(k,k-1);
		 			q = round(M(k,k-1));

		 			if (q != 0)
			 			{
				 			for (lidia_size_t i = 0; i < k-1; i++)
								M.sto(k,i,M(k,i)-q*M(k-1,i));

				 			M.sto(k,k-1,M(k,k-1)-q);

				 			for (lidia_size_t i = 0; i < 3; i++)
								{
									A.sto(k,i,A(k,i)-q*A(k-1,i));
									H.sto(i,k,H(i,k)-q*H(i,k-1));
								}
			 			}

		 			// Preparing to swap(k,k-1)
		 			bigrational di = -A(k-1,k-1)*A(k,k), racineCarree;
	
		 			// Is di a square? (di is determinant)
		 			if (square_root(racineCarree,di))
			 			{
				 			math_matrix <bigrational> IM = trans(Adjoint(M));

				 			bigint aux1, aux2;
				 			sqrt(aux1,numerator(di));
				 			sqrt(aux2,denominator(di));

				 			vecteurSolution[0] = aux1*IM(0,k-1)+aux2*A(k-1,k-1)*IM(0,k);
				 			vecteurSolution[1] = aux1*IM(1,k-1)+aux2*A(k-1,k-1)*IM(1,k);
				 			vecteurSolution[2] = aux1*IM(2,k-1)+aux2*A(k-1,k-1)*IM(2,k);

				 			math_matrix <bigrational> Hr(3,3);
				 			for (lidia_size_t i = 0; i < 3; i++)
								for (lidia_size_t j = 0; j < 3; j++)
									Hr.sto(i,j,(bigrational)H(i,j));

				 			vecteurSolution = Hr*vecteurSolution;

				 			return(1);
			 			}

		 			// Now we reduce [k,k-1]
		 			bigrational Mkk1 = M(k,k-1), bk1new = Mkk1*Mkk1*A(k-1,k-1)+A(k,k);
		 			bigrational Mkk1new;

					if (abs(bk1new) < constanteLLL*abs(A(k-1,k-1)))
						{
							swap = 1;
							Mkk1new = -Mkk1*A(k-1,k-1)/bk1new;
						}
					else
						swap = 0;
			
					// Compute the new matrices after the swap
					 
					// Need to invert the coordinates
					if (swap)
						{
							for (lidia_size_t j = 0; j < 3; j++) 
								{
									bigint aux = H(j,k-1);
									H.sto(j,k-1,H(j,k));
									H.sto(j,k,-aux);
								}
							for (lidia_size_t j = 0; j < k-1; j++) 
								{ 
									bigrational aux = M(k-1,j); 
									M.sto(k-1,j,M(k,j)); 
									M.sto(k,j,-aux);
								}
							for (lidia_size_t j = k+1; j < 3; j++)
								{
									bigrational aux = M(j,k); 
									M.sto(j,k,-M(j,k-1)+Mkk1*aux);
									M.sto(j,k-1,aux+Mkk1new*M(j,k));
								}
							for (lidia_size_t j = 0; j < 3; j++ )
								{
									if ((j != k) && (j != k-1))
										{
											bigrational aux = A(k-1,j); 
											A.sto(k-1,j,A(k,j)); 
											A.sto(k,j,-aux);
											aux = A(j,k-1); 
											A.sto(j,k-1,Mkk1*aux+A(j,k)); 
											A.sto(j,k,-aux-Mkk1new*A(j,k-1));
										}
								}
													 
							bigrational aux1 = A(k-1,k-1), aux2 = A(k,k-1);
							A.sto(k,k-1,-A(k-1,k)-Mkk1*aux1);
							A.sto(k-1,k-1,A(k,k)+Mkk1*aux2);
							A.sto(k,k,aux1-Mkk1new*A(k,k-1));
							A.sto(k-1,k,-aux2-Mkk1new*A(k-1,k-1));
							 
							M.sto(k,k-1,Mkk1new);
							 
							if (k != 1)
								k--; 
						}
				}
			
			for (lidia_size_t l = k-2; l >= 0; l--) 
				{
					// red(k,l)
					q = round(M(k,l));
					if (q != 0) 
						{
							for (lidia_size_t i = 1; i <= l; i++) 
								M.sto(k,i-1,M(k,i-1)-q*M(l,i-1));
			
							M.sto(k,l,M(k,l)-q);
							for (lidia_size_t i = 1; i <= 3; i++) 
								{
									A.sto(k,i-1,A(k,i-1)-q*A(l,i-1));
									H.sto(i-1,k,H(i-1,k)-q*H(i-1,l)); 
								}
						}
				}
			k++;
		} // end while

	// Did not find a solution

	matriceSolution = trans(H)*G*H;

	return(2);
}

// Compute the square root of bigint A modulo bigint B
// Returns 0 if it succeeded, -1 if A is not a square root modulo a prime number
// or its power in the decomposition, -2 if one of the Mi's is not invertible
lidia_size_t SquareRootModulo(const bigint &A, const bigint &B,
						 									const rational_factorization &bFactorise,
						 									bigint &squareRootAModuloB)
{
	// Temp variable
	bigint squareRootAModuloBtmp;

	// Find r in [0,..,|B|/2] such that r^2=A (mod |B|)
	// For that, we decompose B in prime numbers
	//			if B=p_1^n1 * p_2^n2 * ... * p_k^k

	// Compute |B|
	bigint absB = abs(B);

	// If A is negative we get back to A being in [0,B]
	bigint ANormalise;
	remainder(ANormalise,A,absB);
	// listeProduitMi and listeInverseMi are used to store information necessary to
	// apply the Chinese Remainder theorem
	bigint produitMi, inverseMi;
	// We assume B = m_1*m_2*...*m_n
	// We let M_i = B/m_i										listeProduitMi
	//				y_i = M_i^{-1} mod m_i				listeInverseMi
	//				x = sqrt{A}_1*M_1*y_1+sqrt{A}__2*M_2*y_2+...+sqrt{A}__n*M_n*y_n

	// racinesCarrees to store A modulo p_^n_i
	bigint racinesCarrees;

	// Loop on all factors of B to find there square root modulo p^n_i
	for (lidia_size_t i = 0; i < bFactorise.no_of_comp(); i++)
		{
			// Arithmetic of the Galois group GF(p^d)
			bigint p = bFactorise.base(i); 
			lidia_size_t d = bFactorise.exponent(i);

			// The field in which we are going to work
			galois_field field(p,d);

			bigint ppowerd = field.number_of_elements();

			// Variables for solving the equation ygf = sqrt (xgf)
			gf_element xgf(field);
			gf_element ygf(field);

			// Initialize xgf = ANormalise
			xgf.assign(ANormalise);

			// Test if A is a square modulo p^d
			if (xgf.is_square())
				{
		 			// Extract the square root
		 			ygf = sqrt(xgf);
			
					racinesCarrees = ygf.polynomial_rep().const_term();
					// M_i = B/m_i
					produitMi = absB/ppowerd;
			
					// Compute M_i^{-1} mod p^d
					gf_element mgf(field);
			
					mgf.assign(produitMi);
					if (!mgf.is_zero())
						{ 
					 	 	// Mi is invertible
					 	 	gf_element minversegf(field);
			
					 	 	minversegf = inverse(mgf);
					 	 	inverseMi = minversegf.polynomial_rep().const_term();
					 	 
					 	 	// Chinese theorem
					 	 	squareRootAModuloBtmp = squareRootAModuloBtmp+racinesCarrees*produitMi*inverseMi;
			 			}
		 			else // Mi is not invertible
			 			return(-2);
				} 
			else // xgf is not a square root modulo p^d
				return(-1);
		}

	remainder(squareRootAModuloB,squareRootAModuloBtmp,absB);

	if (absB-squareRootAModuloB < squareRootAModuloB)
		squareRootAModuloB = absB-squareRootAModuloB;
	
	return(0);
}

// This function computes the U matrix as described in Simon's paper
// Returns 0 if eveything went fine
//				 1 if Xa = sqrt(-c/b) mod a does not exist
//				 2 if Xb = sqrt(-c/a) mod b does not exist
//				 3 if Xc = sqrt(-b/a) mod c does not exist
lidia_size_t CalculDeU(const bigint &a, const bigint &b, const bigint &c,
												const rational_factorization &aFactorise, 
												const rational_factorization &bFactorise, 
												const rational_factorization &cFactorise,
												bigint_matrix &U)
{
	lidia_size_t codeRetour = 0;

	// Here, a, b and c are assumed squarefree and pairwise coprime

	// Warning: below, the modulus is assigned globally for all bigmods

	// Compute -c/b mod a and Xa = sqrt(-c/b) mod a
	bigint Xa,Xb,Xc;
	if (abs(a) != 1) 
		{		 
			bigmod::set_modulus(a);
			
			// First compute 1/b mod a
			bigmod bInverseTmpBigmod = b;
			bInverseTmpBigmod.invert(0);

			bigmod cTmpBigmod1 = c;
			bigmod resultTmpBigmod = - cTmpBigmod1*bInverseTmpBigmod;

			bigint moinsCsurBmoduloA = mantissa(resultTmpBigmod);

			if (SquareRootModulo(moinsCsurBmoduloA,a,aFactorise,Xa)) 
				codeRetour = 1;
		}

	if (!codeRetour)
		{
			// Compute -c/a mod b and Xb = sqrt(-c/a) mod b
			if (abs(b) != 1) 
				{	 
					bigmod::set_modulus(b);
			
					// First compute 1/a mod b
					bigmod aInverseTmpBigmod = a;
					aInverseTmpBigmod.invert(0);
			
					bigmod cTmpBigmod = c;
					bigmod resultTmpBigmod = -cTmpBigmod*aInverseTmpBigmod;
					 
					bigint moinsCsurAmoduloB = resultTmpBigmod.mantissa();
			
					if (SquareRootModulo(moinsCsurAmoduloB,b,bFactorise,Xb))
						codeRetour = 2;
				}
		}

	if (!codeRetour)
		{
			// Compute -b/a mod c and Xc = sqrt(-b/a) mod c
			if (abs(c) != 1) 
				{ 
					bigmod::set_modulus(c);
					 
					// First compute 1/a mod b
					bigmod aInverseTmpBigmod = a;
					aInverseTmpBigmod.invert(0);
			
					bigmod bTmpBigmod = b;
					bigmod resultTmpBigmod = -bTmpBigmod*aInverseTmpBigmod ;
			
					bigint moinsBsurAmoduloC = resultTmpBigmod.mantissa();
			
					if (SquareRootModulo(moinsBsurAmoduloC,c,cFactorise,Xc))
						codeRetour = 3;
				}
		}

	if (!codeRetour)
		{
			// Now we look for u,v such that b*u+c*v = 1
			bigint u,v;
			bigint pgcd = xgcd(u,v,b,c);

			// Now build U
			// First line
			U.sto(0,0,b*c); 
			U.sto(0,1,a*b*u*Xc); 
			U.sto(0,2,Xa*b*u*Xc+Xb*c*v); 
			// Second line
			U.sto(1,1,a); 
			U.sto(1,2,Xa);	
			// Third line
			U.sto(2,2,1);	 
		}

	return(codeRetour);
}

// Compute the Q' matrix as in Simon's paper
lidia_size_t CalculQprime(const bigint &a, const bigint &b, const bigint &c, 
													const rational_factorization &aFactorise, 
			  									const rational_factorization &bFactorise,
			  									const rational_factorization &cFactorise,
			  									bigint_matrix &Qprime,
			  									bigint_matrix &U)
{
	lidia_size_t codeRetour = CalculDeU(a,b,c,aFactorise,bFactorise,cFactorise,U);
	
	if (codeRetour) 
		return(codeRetour);
	else
		{
			bigint_matrix Fmatrice(3,3);

			Fmatrice.sto(0,0,a); 
			Fmatrice.sto(1,1,b);
			Fmatrice.sto(2,2,c);	
	
			Qprime = (trans(U)*Fmatrice*U)/(a*b*c);

			return (0);
		}
}

// Find a small rational solution of a*x^2+b*y^2+c*z^2 = 0, where a, b, c are
// square-free and pairwise coprime
// Returns 1 if everything went fine
lidia_size_t SolveLegendreEquationSFCP(const bigint &a, const bigint &b, const bigint &c,
																				const rational_factorization &aFactorise,
																				const rational_factorization &bFactorise,
																				const rational_factorization &cFactorise,
																				math_vector <bigrational> &sol)
{
	// Computation of Qprime
	bigint_matrix Qprime(3,3), U(3,3);

	lidia_size_t codeRetour = CalculQprime(a,b,c,aFactorise,bFactorise,cFactorise,Qprime,U);
			
	if (codeRetour)
		return(0);
	else
		{
			// Apply LLL to Qprime -- dimension of solution is necessarily 1

			bigint_matrix ms(3,3), H(3,3);

			LLL(Qprime,sol,ms,H);

			// Make U bigrat
			math_matrix <bigrational> Ur(3,3);
			for (lidia_size_t i = 0; i < 3; i++)
				for (lidia_size_t j = 0; j < 3; j++)
		 			Ur.sto(i,j,(bigrational)U(i,j));

			sol = Ur*sol;
			
			return(1);
		}
}

// Output is b,c such that a = b^2 c, c is factored
void ExtractSQ(const bigint &a, bigint &b, rational_factorization &c)
{
	// If a is a perfect square, no need to do the big thing 
	if (is_square(b,a*a.sign())) 
		{
			c.assign(a.sign());
			c.factor();
		}
	else // a is not a perfect square 
		{ 
			b = 1;

			c.assign(a);
			c.factor();

			bigint h,tmp_b;									
			lidia_size_t tmp_e, cmax = c.no_of_comp(), j = 0;

			for (lidia_size_t i = 0; i < cmax; i++) 
				{ 
					tmp_e = c.exponent(i-j); 
					tmp_b = c.base(i-j); 
			
					// Odd exponent 
					if (tmp_e % 2)												 
						{																										 
							if (tmp_e == 1) 
								h.assign_one(); 
							else 
								power(h,tmp_b,(tmp_e-1)/2); 
			
							c.set_exponent(i-j,1);
						} 
					else // Even exponent 
						{
							if (tmp_e == 2) 
								h.assign(tmp_b); 
							else 
								power(h,tmp_b,tmp_e/2); 
			
							c.set_exponent(i-j,0);
							j++;
			 			}

		 			multiply(b,b,h);										
				}												
		}
}

// Find common factors
void CommonFactors(rational_factorization &af, rational_factorization &bf,
			 							rational_factorization &cf, bigint &cx)
{
	lidia_size_t ia = 0, ib = 0, na = af.no_of_comp(), nb = bf.no_of_comp();
	bigint ta, tb;

	cx = 1;

	while ((ia < na) && (ib < nb))
		{
			ta = af.base(ia);
			tb = bf.base(ib);

			if (ta == 1)
				break;
			else if (tb == 1)
				break;
			else if (ta < tb)
				ia++;
			else if (tb < ta)
				ib++;
			else if (ta == tb)
				{
					cx = cx*ta;
			
					// Multiply cf by appropriate factorization
					rational_factorization tmp;
					tmp.assign(ta,1);
					multiply(cf,cf,tmp);
				
					af.set_exponent(ia,0);
					bf.set_exponent(ib,0);
					 	 
					na--;
					nb--;
				}
		}
}

// Find a small rational solution of a*x^2+b*y^2+c*z^2 = 0
// Returns 1 if everything went fine
lidia_size_t SolveLegendreEquation(const bigint &a, const bigint &b, const bigint &c,
					 													math_vector <bigrational> &sol)
{
	// First divide by the gcd of a,b,c
	bigint pgcd = gcd(a,gcd(b,c));
	bigint aa = a/pgcd, bb = b/pgcd, cc = c/pgcd;

	// Extract the square part of aa,bb,cc
	bigint al,bl,cl;
	rational_factorization af,bf,cf;

	ExtractSQ(aa,al,af);
	ExtractSQ(bb,bl,bf);
	ExtractSQ(cc,cl,cf);

	// Find common factors of the factorizations
	bigint ax,bx,cx;

	CommonFactors(af,bf,cf,cx);
	CommonFactors(af,cf,bf,bx);
	CommonFactors(bf,cf,af,ax);

	bigint aaa = evaluate_to_bigint(af);
	bigint bbb = evaluate_to_bigint(bf);
	bigint ccc = evaluate_to_bigint(cf);
	
	// Now solve
	if (SolveLegendreEquationSFCP(aaa,bbb,ccc,af,bf,cf,sol))
		{
			// Divide by the square factors
			sol[0] = sol[0]/al;
			sol[1] = sol[1]/bl;
			sol[2] = sol[2]/cl;

			// Multiply by constant extracted above
			sol[0] = sol[0]*ax;
			sol[1] = sol[1]*bx;
			sol[2] = sol[2]*cx;

			return(1);
		}
	else
		return(0);
}

// Computes an integer solution on a conic with integer coefficients
// Returns 1 if one such point has been found
lidia_size_t FindRationalPointOnConic(const bigint_matrix &q, math_vector <bigint> &sol,
																			std::ostream &s)
{
	#ifndef NDEBUG
	s << ">> looking for rational point on conic -- ";
	#endif

	// We first run LLL on the matrix. You never know...
	bigint_matrix q2(3,3), H(3,3);
	math_vector <bigrational> vec(3,3);

	lidia_size_t dimensionSolution = LLL(q,vec,q2,H), code;

	if (dimensionSolution == 1)
		{
			// A rational point was found by LLL
			#ifndef NDEBUG
			s << "LLL says yes :-)" << endl;
			#endif

			sol = RatToBigint(vec);
			optimize(sol);

			code = 1;
		}
	else
		{
			// LLL did not find a rational point, proceed
			// Don't forget to multiply by matrix H at the end
			
			// Transform to a diagonal matrix

			bigint A = q2(0,0), B = q2(0,1), C = q2(0,2), D = q2(1,1), E = q2(1,2), F = q2(2,2);

			bigint Dp = A*D-B*B, Ep = A*E-B*C, Fp = A*F-C*C;
			
			// Transformation matrix
			math_matrix <bigrational> tr(3,3);

			tr.sto(0,0,bigrational(1,A));
			tr.sto(1,1,bigrational(1,Dp));
			tr.sto(2,2,1);
			tr.sto(1,2,-Ep*tr(1,1));
			tr.sto(0,1,-B*tr(0,0)*tr(1,1));
			tr.sto(0,2,-Ep*tr(0,1)-C*tr(0,0));

			// With this transformation, the Legendre equation is
			//			 Dp x^2 + y^2 + (Dp*Fp-Ep^2) z^2 = 0

			bigint zco = Dp*Fp-Ep*Ep;

			if (SolveLegendreEquation(Dp,1,zco,vec))
				{
					#ifndef NDEBUG
					s << "Legendre says yes :-)" << endl;
					#endif
			
					// Make H bigrational
					math_matrix <bigrational> Hr(3,3);
					for (lidia_size_t i = 0; i < 3; i++)
						for (lidia_size_t j = 0; j < 3; j++)
							Hr.sto(i,j,(bigrational)H(i,j));
			
					// Final solution
					sol = RatToBigint(Hr*tr*vec);
			
					code = 1;
				}
			else
				{ 
					#ifndef NDEBUG
		 			s << "Legendre says no :-(" << endl;
					#endif

		 			code = 0;
				}	 
		}

	return(code);
}

} // end of namespace QI
