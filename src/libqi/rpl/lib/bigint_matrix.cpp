#include <libqi/rpl/bigint_matrix.h>
#include <libqi/rpl/math_vector.h>
#include <iostream>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace rpl {
void bigint_matrix::split_t(bigint_matrix & B,
														bigint_matrix & C,
														bigint_matrix & D,
														bigint_matrix & E) const {
	assert(B.size1() + D.size1() == this->size1());
	assert(B.size2() + C.size2() == this->size2());
	assert(D.size2() + E.size2() == this->size2());
	assert(E.size1() + C.size1() == this->size1());

	B = project(*this, range(0, B.size1()), range( 0, B.size2()));
	C = project(*this, range(0, B.size1()), range( B.size2(), this->size2()));
	D = project(*this, range(B.size1(), this->size1()), range(0, B.size2()));
	E = project(*this, range(B.size1(), this->size1()), range(B.size2(),this->size2())); 
}

void bigint_matrix::split_h(bigint_matrix & B,
														bigint_matrix & C) const {
	assert(B.size1() == this->size1());
	assert(C.size1() == this->size1());
	assert(B.size2() + C.size2() == this->size2());

	B = project(*this, range(0, B.size1()), range(0, B.size2()));
	C = project(*this, range(0, C.size1()), range(B.size2(), this->size2()));
}

void bigint_matrix::compose_t(const bigint_matrix & B,
															const bigint_matrix & C,
															const bigint_matrix & D,
															const bigint_matrix & E) {
	assert(B.size1() + D.size1() == this->size1());
	assert(B.size2() + C.size2() == this->size2());
	assert(D.size2() + E.size2() == this->size2());
	assert(E.size1() + C.size1() == this->size1());

	project(*this, range(0, B.size1()), range( 0, B.size2())) = B;

	project(*this, range(0, B.size1()), range( B.size2(), this->size2())) = C;
	project(*this, range(B.size1(), this->size1()), range(0, B.size2())) = D;
	project(*this, range(B.size1(), this->size1()), range(B.size2(),this->size2())) = E; 
}

void bigint_matrix::compose_h(const bigint_matrix & B,
															const bigint_matrix & C) {
	assert(B.size1() == this->size1());
	assert(C.size1() == this->size1());
	assert(B.size2() + C.size2() == this->size2());

	project(*this, range(0, B.size1()), range( 0, B.size2())) = B;
	project(*this, range(0, B.size1()), range( B.size2(), this->size2())) = C;
}

/////////// Code for LLL
/////////// Porting of Paul Zimmermann's code, to be found here:
/////////// http://www.loria.fr/~zimmerma/software/lll.c

//double Norm_vec(const math_vector <bigint> &x, const long n)
//{
//	long i;
//	double res;
//	bigint t, s;
//
//	for (i = 0; i <= n; i++) {
//		multiply(t, x[i], x[i]);
//		add(s, s, t);
//	}
//	mpz_sqrt(s.val.get_mpz_t(), s.val.get_mpz_t());
//	res = mpz_get_d(s.val.get_mpz_t());
//	return res;
//}
//
//void ident(bigint_matrix &X, const long m, const long n)
//{	 
//	long i, j;
////	X = bigint_matrix(m,n);
//	
//	for (i = 1; i <= m; i++)	 
//		for (j = 1; j <= n; j++)	
//			if (i == j)
//				X(i,j) = 1;
//			else	 
//				X(i,j) = 0;
//} 
//
//void InnerProduct(bigint &x, const math_vector <bigint> &a, const math_vector <bigint> &b, const long n, bigint &t1)
//{
//	long i;
//	
//	cout << n << endl;
//	cout << a << endl;
//	cout << b << endl;
//	
//
//	x = 0;
//	for (i = 1; i <= n; i++) {
//		multiply(t1, a[i], b[i]);
//		add(x, x, t1);
//	}
//}
//
//void IncrementalGS(const bigint_matrix &B, math_vector <long> &P, math_vector <bigint> &D, 
//									 bigint_matrix &lam, long &s, const long k)
//{
//		cout << "in IGS" << endl;
//cout << k << endl;
//// à corriger
//	long n = B.size2()-1;
//	bigint u, t1, t2;
//	long i, j, posj;
//
//	for (j = 1; j <= k-1; j++) {
//		posj = P[j];
//		if (posj == 0) 
//			continue;
//
//		InnerProduct(u, row(B,k), row(B,j), n, t1);
//		for (i = 1; i <= posj-1; i++) {
//			multiply(t1, D[i], u);
//			multiply(t2, lam(k,i), lam(j,i));
//			subtract(t1, t1, t2);
//			mpz_div(t1.val.get_mpz_t(), t1.val.get_mpz_t(), D[i-1].val.get_mpz_t());
//			u = t1;
//		}
//
//		lam(k,posj) = u;
//	}
//cout << "h" << endl;
//	InnerProduct(u, row(B,k), row(B,k), n, t1);
//cout << "h" << endl;
//	for (i = 1; i <= s; i++) {
//		multiply(t1, D[i], u);
//		multiply(t2, lam(k,i), lam(k,i));
//		subtract(t1, t1, t2);
//		mpz_div(t1.val.get_mpz_t(), t1.val.get_mpz_t(), D[i-1].val.get_mpz_t());
//		u = t1;
//	}
//cout << "h" << endl;
//
//	if (u.is_zero())
//		P[k] = 0;
//	else
//		{
//			s++;
//			P[k] = s;
//			D[s] = u;
//		}
//		cout << "done IGS" << endl;
//}
//
///* rounds a/d to nearest integer, breaking ties
//	 by rounding towards zero.	Assumes d > 0. */
//void BalDiv(bigint &q, const bigint &a, const bigint &d, bigint &r)
//{
//	long cmp;
//
//	mpz_fdiv_qr(q.val.get_mpz_t(), r.val.get_mpz_t(), a.val.get_mpz_t(), d.val.get_mpz_t());
//
//	mpz_mul_2exp(r.val.get_mpz_t(), r.val.get_mpz_t(), 1);
//
//	cmp = mpz_cmp(r.val.get_mpz_t(), d.val.get_mpz_t());
//	if (cmp > 0 || (cmp == 0 && mpz_cmp_ui(q.val.get_mpz_t(), 0) < 0))
//		mpz_add_ui(q.val.get_mpz_t(), q.val.get_mpz_t(), 1);
//}
//
///* c = c1 - x*c2 */
//void MulSub(bigint &c, const bigint &c1, const bigint &c2, const bigint &x, bigint &tmp)
//{
//	multiply(tmp, x, c2);
//	subtract(c, c1, tmp);
//}
//
///* c = c - x*c2 */
//// No reference ???
//void MulSubN(math_vector <bigint> c, math_vector <bigint> c2, const bigint &x, 
//						 const long n, bigint &tmp)
//{
//	long i;
//	signed long int x0;
//
//	x0 = mpz_get_si(x.val.get_mpz_t());
//	if (mpz_cmp_si(x.val.get_mpz_t(), x0) == 0 && 
//			x0 != ((signed long int) 1 << (mp_bits_per_limb - 1))) {
//		if (x0 > 0)
//			for (i = 1; i <= n; i++) {
//	 			mpz_mul_ui(tmp.val.get_mpz_t(), c2[i].val.get_mpz_t(), x0);
//	 			subtract(c[i], c[i], tmp);
//			 }
//		else if (x0 < 0) {
//			x0 = -x0;
//			for (i = 1; i <= n; i++)
//	 			mpz_addmul_ui(c[i].val.get_mpz_t(), c2[i].val.get_mpz_t(), x0);
//		}
//	}
//	else
//	{
//		for (i = 1; i <= n; i++)
//			{
//				multiply(tmp, c2[i], x);
//				subtract(c[i], c[i], tmp);
//			}
//	}
//}
//
///* ????? t1 and r are temporary variables */
//void reduce(const long k, const long l, bigint_matrix &B, math_vector <long> &P, 
//						math_vector <bigint> &D, bigint_matrix &lam, bigint_matrix &U, bigint &t1, bigint &r)
//{
//	long j;
//
//	if (P[l] == 0) 
//		return;
//
//	mpz_mul_2exp(t1.val.get_mpz_t(), lam(k,P[l]).val.get_mpz_t(), 1);
//	mpz_abs(t1.val.get_mpz_t(), t1.val.get_mpz_t());
//	if (mpz_cmp(t1.val.get_mpz_t(), D[P[l]].val.get_mpz_t()) <= 0)
//		return;
//
//	BalDiv(r, lam(k,P[l]), D[P[l]], t1);
//	// a corriger
//	MulSubN(row(B,k), row(B,l), r, B.size2()-1, t1);
//
////	if (U)
//// a corriger
//		MulSubN(row(U,k), row(U,l), r, B.size2()-1, t1);
//
//	for (j = 1; j <= l-1; j++)
//		if (P[j] != 0)
//			MulSub(lam(k,P[j]), lam(k,P[j]), lam(l,P[j]), r, t1);
//
//	MulSub(lam(k,P[l]), lam(k,P[l]), D[P[l]], r, t1);
//}
//
///* test if a*d1^2 > b*(d0*d2 + la^2), t1 and t2 are temporary variables */
//long SwapTest(const bigint &d0, const bigint &d1, const bigint &d2, const bigint &la,
//							const bigint &a, const bigint &b, bigint &t1, bigint &t2)
//{
//	multiply(t1, d0, d2);
//	multiply(t2, la, la);
//	add(t1, t1, t2);
//	multiply(t1, t1, b);
//
//	multiply(t2, d1, d1);
//	multiply(t2, t2, a);
//
//	return (mpz_cmp(t2.val.get_mpz_t(), t1.val.get_mpz_t()) > 0);
//}
//
///* c = (x*c1 + y*c2)/z, warning: c and z can be the same variable
//	 t1 and t2 are temporary variables */
//void MulAddDiv(bigint &c, const bigint &c1, const bigint &c2, 
//							 const bigint &x, const bigint &y, const bigint &z, bigint &t1, bigint &t2)
//{
//	multiply(t1, x, c1);
//	multiply(t2, y, c2);
//	add(t1, t1, t2);
//	divide(c, t1, z);
//}
//
///* c = (x*c1 - y*c2)/z, t1 is a temporary variable */
//void MulSubDiv(bigint &c, const bigint &c1, const bigint &c2, 
//							 const bigint &x, const bigint &y, const bigint &z, bigint &t1)
//{
//	multiply(t1, x, c1);
//	multiply(c, y, c2);
//	subtract(t1, t1, c);
//	divide(c, t1, z);
//}
//
///* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
//void RowTransform(bigint &c1, bigint &c2, const bigint &x, const bigint &y, 
//									const bigint &u, const bigint &v)
//{
//	bigint t1, t2;
//
//	multiply(t1, x, c1);
//	multiply(t2, y, c2);
//	add(t1, t1, t2);
//
//	multiply(t2, u, c1);
//	c1 = t1;
//	multiply(t1, v, c2);
//	add(c2, t1, t2);
//}
//
///* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
//void RowTransformN(math_vector <bigint> c1, math_vector <bigint> c2,
//									 const bigint &x, const bigint &y, const bigint &u, const bigint &v, const long n)
//{
//	bigint t1, t2;
//	long i;
//
//	for (i = 1; i <= n; i++)
//		{
//			multiply(t1, x, c1[i]);
//			multiply(t2, y, c2[i]);
//			add(t1, t1, t2);
//
//			multiply(t2, u, c1[i]);
//			c1[i] = t1;
//			multiply(t1, v, c2[i]);
//			add(c2[i], t1, t2);
//		}
//}
//
//#define lllswap(x, y) { long _tmp = (x); (x) = (y); (y) = _tmp; }
//#define lllmpz_swap_n(x, y) { math_vector <bigint> _tmp = (x); (x) = (y); (y) = _tmp; }
//
///* swaps vectors k-1 and k;	 assumes P(k-1) != 0 */
//void swapLLL (long k, bigint_matrix &B, math_vector <long> &P, math_vector <bigint> &D, 
//				 bigint_matrix &lam, bigint_matrix &U, long m, long verbose)
//{
//	 long i, j;
//	 bigint t1, t2, t3, e, x, y;
//
//	 if (P[k] != 0) {
//			if (verbose) fprintf(stderr, "swap case 1: %ld\n", k);
//
//			lllmpz_swap_n(row(B,k-1), row(B,k));
////			if (U)
//				lllmpz_swap_n(row(U,k-1), row(U,k-1));
//			
//			for (j = 1; j <= k-2; j++)
//				 if (P[j] != 0)
//						mpz_swap(lam(k-1,P[j]).val.get_mpz_t(), lam(k,P[j]).val.get_mpz_t());
//
//			for (i = k+1; i <= m; i++) {
//				 MulAddDiv(t1, lam(i,P[k]-1), lam(i,P[k]),
//									 lam(k,P[k]-1), D[P[k]-2], D[P[k]-1], t2, t3);
//				 MulSubDiv(lam(i,P[k]), lam(i,P[k]-1), lam(i,P[k]), 
//									 D[P[k]], lam(k,P[k]-1), D[P[k]-1], t2);
//				 mpz_set(lam(i,P[k]-1).val.get_mpz_t(), t1.val.get_mpz_t());
//			}
//
//			MulAddDiv(D[P[k]-1], D[P[k]], lam(k,P[k]-1),
//								D[P[k]-2], lam(k,P[k]-1), D[P[k]-1], t2, t3);
//	 }
//	 else if (mpz_cmp_ui(lam(k,P[k-1]).val.get_mpz_t(), 0) != 0) {
//			if (verbose) fprintf(stderr, "swap case 2: %ld\n", k);
//			mpz_gcdext(e.val.get_mpz_t(), x.val.get_mpz_t(), y.val.get_mpz_t(), lam(k,P[k-1]).val.get_mpz_t(), D[P[k-1]].val.get_mpz_t());
//
//			mpz_divexact(t1.val.get_mpz_t(), lam(k,P[k-1]).val.get_mpz_t(), e.val.get_mpz_t());
//			mpz_divexact(t2.val.get_mpz_t(), D[P[k-1]].val.get_mpz_t(), e.val.get_mpz_t());
//
//			mpz_set(t3.val.get_mpz_t(), t2.val.get_mpz_t());
//			mpz_neg(t2.val.get_mpz_t(), t2.val.get_mpz_t());
//			// a corriger
//			RowTransformN(row(B,k-1), row(B,k), t1, t2, y, x, B.size2()-1);
////			if (U)
//// a corriger
//				RowTransformN(row(U,k-1), row(U,k), t1, t2, y, x, B.size2()-1);
//			for (j = 1; j <= k-2; j++)
//				 if (P[j] != 0)
//						RowTransform(lam(k-1,P[j]), lam(k,P[j]), t1, t2, y, x);
//
//			mpz_mul(t2.val.get_mpz_t(), t2.val.get_mpz_t(), t2.val.get_mpz_t());
//			mpz_divexact(D[P[k-1]].val.get_mpz_t(), D[P[k-1]].val.get_mpz_t(), t2.val.get_mpz_t());
//
//			for (i = k+1; i <= m; i++)
//				 if (P[i] != 0) {
//						mpz_divexact(D[P[i]].val.get_mpz_t(), D[P[i]].val.get_mpz_t(), t2.val.get_mpz_t());
//						for (j = i+1; j <= m; j++) {
//							 mpz_divexact(lam(j,P[i]).val.get_mpz_t(), lam(j,P[i]).val.get_mpz_t(), t2.val.get_mpz_t());
//						}
//				 }
//
//			for (i = k+1; i <= m; i++) {
//				 mpz_divexact(lam(i,P[k-1]).val.get_mpz_t(), lam(i,P[k-1]).val.get_mpz_t(), t3.val.get_mpz_t());
//			}
//
//			lllswap(P[k-1], P[k]);
//	 }
//	 else {
//			if (verbose) fprintf(stderr, "swap case 3: %ld\n", k);
//
//			lllmpz_swap_n(row(B,k-1), row(B,k));
////			if (U)
//				lllmpz_swap_n(row(U,k-1), row(U,k));
//	 
//			for (j = 1; j <= k-2; j++)
//				 if (P[j] != 0)
//						mpz_swap(lam(k-1,P[j]).val.get_mpz_t(), lam(k,P[j]).val.get_mpz_t());
//
//			lllswap(P[k-1], P[k]);
//	 }
//}
//
//long LLL(bigint &det, bigint_matrix &B, bigint_matrix &U, 
//				 const bigint &a, const bigint &b, const long verbose)
//{
//// à corriger
//	long m = B.size1()-1, n = B.size2()-1, j, s, k, max_k;
//	bigint tmp1, tmp2;
//	
//	math_vector <long> P(m+1);
//	math_vector <bigint> D(m+1);
//	bigint_matrix lam(m+1,n+1);
//
//	for (j = 0; j <= m; j++)
//		mpz_init_set_ui(D[j].val.get_mpz_t(), j == 0);
//cout << "done" << endl;
////	if (U) 
//		ident(U, m,n);
//cout << "done" << endl;
//	s = 0;
//
//	k = 1;
//	max_k = 0;
//
//	while (k <= m) {
//	cout << "k: " << k << endl;
//		if (k > max_k)
//			{
//				IncrementalGS(B, P, D, lam, s, k);
//				max_k = k;
//			}
//
//		if (k == 1) {
//			k++;
//			continue;
//		}
//cout << "inr" << endl;
//		reduce(k, k-1, B, P, D, lam, U, tmp1, tmp2);
//cout << "doner" << endl;
//cout << P[k] << "   " << P[k]-1 << "   " << P[k]-2 << "   " << endl;
//		if (P[k-1] != 0 && 
//				(P[k] == 0 || SwapTest(D[P[k]], D[P[k]-1], D[P[k]-2], lam(k,P[k]-1), a, b, tmp1, tmp2))) {
//				cout << "swaplll" << endl;
//			swapLLL (k, B, P, D, lam, U, max_k, verbose);
//			k--;
//		}
//		else {	
//		   cout << "other" << endl;
//			for (j = k-2; j >= 1; j--) 
//				reduce(k, j, B, P, D, lam, U, tmp1, tmp2);
//			k++;
//		}
//	}
//	
//	cout << "outk" << endl;
//
//	mpz_set(det.val.get_mpz_t(), D[s].val.get_mpz_t());
//
//	return s;
//}
//




















// we use the kernel algorithm of the Cohen's Book : 
// A course in Computational Algrebraic Number
// Theory (p. 73)
// it seems there is an error on the stopping
// condition (P. Zimmermann correction)
// we replace i=l by (k=1 or i=1)
bigint_matrix kernel(const bigint_matrix &inputA) {

	bigint_matrix a = inputA;
	int m = a.size1();
	int n = a.size2();
	int j = n - 1;
	int k = n - 1;

	bigint_matrix U = math_matrix<bigint>::identity(n);
	math_vector<bigint> B1(m);
	math_vector<bigint> B2(n);
	
	// loop on rows
	for (int i = m - 1; i >= 0 && k >= 0; i--) {
		while (j != 0) { // Step 2 : check zero if j = 1 go step 4
			j--;
			
		// if a_ij =	0 go step 2
		if (!a(i,j).is_zero()) {
			// Step 3 : Euclidean Step
			bigint a_ik = a(i, k);
			bigint a_ik_c = a_ik;
			bigint a_ij = a(i, j);
			
			bigint u,v, q1, q2, d;
			d = xgcd(u, v, a_ik, a_ij);
			div_exact(q1, a_ik, d); // q1 = a_ik/d
			div_exact(q2, a_ij, d); // q2 = a_ij/d
			
			// we need to minimise u and v
			a_ik_c.abs();
			if (a_ik_c == d) {
				v.assign_zero();
				u = a_ik.sign();
			}

			B1 =	u * column(a, k) +	v * column(a, j);
			column(a, j) = q1 * column(a, j) - q2 * column(a, k);
			column(a, k) = B1;
			// same work on U
			B2 =	u * column(U,k) +	 v * column(U, j);
			column(U, j) = q1 * column(U, j) - q2 * column(U, k);	 
			column(U, k) = B2;
		}
	}	
	
	// Step 4 : j = 0	 
	bigint b = a(i, k);
	if (b < 0) {
		column(a, k) = -column(a, k);
	 	column(U, k) = -column(U, k);
		b.negate(); 
	}
		 
	// b = 0
	if (b.is_zero()) {
		k++;
	}
	// otherwise
	else {
		for (int w = k+1 ; w < n; ++w) { 
			bigint q;
			floor_divide(q, a(i, w), b);
			column(a, w) -= q * column(a, k);
			column(U, w) -= q * column(U, k);
		}
	}
	k--;
	j = k;
	}
	
	U.resize(n,k+1);	
	
	
//	
//	cout << U << endl;
//	
//	bigint_matrix UP(4,3);
//	UP(0,0) = 1; UP(1,0) = 1; UP(2,0) = -1; UP(3,0) = 3;
//	UP(0,1) = 3; UP(1,1) = 1; UP(2,1) = 2; UP(3,1) = 1;
//	UP(0,2) = 3; UP(1,2) = -2; UP(2,2) = 2; UP(3,2) = 4;
//	
//	
//	
//	bigint det;
//	bigint_matrix B(UP.size1()+1,UP.size2()+1), UU(UP.size1()+1,UP.size2()+1);
//	
//	for (unsigned int i = 0; i < UP.size1(); i++)
//		for (unsigned int j = 0; j < UP.size2(); j++)
//			B(i+1,j+1) = UP(i,j);
//		
//	B = trans(B);
//	UU = trans(UU);
//	
//	cout << B << endl;
//	
//  LLL(det, B, UU, 1, 1, 0);
//	
//	cout << B << endl;
//	cout << UU << endl;
//	
//	
	
	return U;
}

size_t rank(const bigint_matrix & a) {
	
	bigint_matrix k = kernel(a);
	// dim E = dim ker f + rank(f)
	return (a.size2() - k.size2());
}

bigint det(const bigint_matrix &a) {
	// checks that A is square
	assert(a.size1() == a.size2());
	bigint_matrix m = a;
	int n = m.size1();
	// step 1
	bigint c = 1;
	bigint s = 1;
	for(rpl_size_t k = 0; k <n ; ++k) {
		 // step 2
		 if (k == (n-1)) 
				 return s * m(k,k);
		 bigint p = m(k,k);
		 if (p.is_zero()) { // step 3
				 rpl_size_t i = k+1;
	 for(; i < n ; ++i) // search for a non-nul pivot
			 if ( !m(i,k).is_zero() ) break;
				 if ( i == n) { // all pivots are nul 
				return bigint(0); // then terminate
	 }
	 else {
					 for(rpl_size_t j = k; j < n ; ++j) {
							swap(m(i,j), m(k,j));
			}
			s.negate();
					 p.assign(m(k,k)); 
	 }
		 }
		 //step 4 : Main step
		for(rpl_size_t i = k + 1; i < n ; ++i) {
			for(rpl_size_t j = k + 1; j < n ; ++j) {
			bigint t = p * m(i, j) - m(i, k) * m(k, j);
			div_exact(m(i,j), t, c);
			}
		}
		c.assign(p);
	} //step 2
	std::cerr << "must not be here!" << std::endl;
	abort();
}

bigint_matrix adj(const bigint_matrix & a) {
	 assert(a.size1() == a.size2());
	 assert(a.size1() == 3);
	 // WE CODE THE ALGORITHM ONLY FOR n = 3, cf -> QIVanishDet.cc: 1651
	 bigint_matrix c(3,3);
	 c(0, 0) =	 a(1, 1) * a(2, 2) - a(2, 1) * a(1, 2);
	 c(0, 1) = -(a(1, 0) * a(2, 2) - a(2, 0) * a(1, 2));
	 c(0, 2) =	 a(1, 0) * a(2, 1) - a(2, 0) * a(1, 1);
	 
	 c(1, 0) = -(a(0, 1) * a(2, 2) - a(2, 1) * a(0, 2));
	 c(1, 1) =	 a(0, 0) * a(2, 2) - a(2, 0) * a(0, 2);
	 c(1, 2) = -(a(0, 0) * a(2, 1) - a(2, 0) * a(0, 1));

	 c(2, 0) =	 a(0, 1) * a(1, 2) - a(1, 1) * a(0, 2);
	 c(2, 1) = -(a(0, 0) * a(1, 2) - a(1, 0) * a(0, 2));
	 c(2, 2) =	 a(0, 0) * a(1, 1) - a(1, 0) * a(0, 1);

	 return c;
}

};
