#ifndef _POLY_H_
#define _POLY_H_

#include <boost/array.hpp>
#include <algorithm>
#include <functional>
#include <iostream>
#include <libqi/rpl/bigint.h>
#include <libqi/rpl/rpl.h>

//! \file polynomial.h

namespace rpl {

//! univariate polynomial class
/*! that structure implements the univariate polynomial class 
* all algorithms are under naive form (it is just a prototype)
* we use a array( i.e. compile-time bounded size array) for the storage
* we assume that maximal degree is 8!
*/
template <class T>
class polynomial {
	public:
		// members 
		typedef typename boost::array<T,9> storage_type;
		
		storage_type coeff; 
		
		rpl_size_t deg;			

		// constructors 
		polynomial() { this->assign_zero(); }

		polynomial(const T &x); 
		
		polynomial(const T *x, int d);

		polynomial(const polynomial <T> &f);
			 
		// destructor
		~polynomial() {}
		
		// copy assignment
		polynomial <T> & operator = (const polynomial <T> &a)
			{
				if (this != &a) // Beware of self-assignment
		 			this->assign(a);

				return *this;
			}
			
    polynomial <T> & operator = (const bigint &a)
      {
				this->assign(a);

				return *this;
      }

		// assignments
		void assign_zero();

		void assign_one();
		
		void assign_x();

		void assign_x2();

		void assign(const T &x);
	 
		void assign(const polynomial <T> &f);

		// basic methods
		void remove_leading_zeros();

		// access methods
		const rpl_size_t degree() const { return this->deg; }
		
		rpl_size_t degree() { return this->deg; }

		void set_degree(rpl_size_t n) { this->deg = n; }

		T & operator[] (rpl_size_t i) { return this->coeff[i]; }
		
		const T operator[] (rpl_size_t i) const { return this->coeff[i]; }

		T lead_coeff() const;

		T const_term() const { return this->coeff[0]; }

		// binary operations as functions
		void add(const polynomial <T> &g, const polynomial <T> &h); 

		void subtract(const polynomial <T> &g, const polynomial <T> &h); 

		void multiply(const polynomial <T> &g, const polynomial <T> &h); 
	 
		void power(const polynomial <T> &g, rpl_size_t e); 
		
		void negate(const polynomial<T> &x);

		// comparisons
		bool is_zero () const; 
		
		bool is_one () const;
};

//////////////////////////////////////// deported constructors

template <typename T>
polynomial<T>::polynomial(const T & x) {
	this->assign(x);
}

// constructor from a vector
template <typename T>
polynomial<T>::polynomial(const T *x, int d) {
	std::copy(x, x+ d+ 1, this->coeff.begin()); 
	this->deg = d;
	remove_leading_zeros();
}

template <typename T>
polynomial<T>::polynomial(const polynomial<T> &f) {
	this->assign(f);
}

//////////////////////////////////////// deported assignments

template <class T>
void polynomial<T>::assign_zero() {
	this->deg = -1;
	std::for_each(this->coeff.begin(),this->coeff.end(),std::mem_fun_ref(&T::assign_zero));
}

template <class T>
void polynomial<T>::assign_one() {
	this->assign_zero();
	this->coeff[0].assign_one();
	this->deg = 0;
}

template <class T>
void polynomial<T>::assign_x() {
	this->assign_zero();
	this->coeff[1].assign_one();
	this->deg = 1;
}
		
template <class T>
void polynomial<T>::assign_x2() {
	this->assign_zero();
	this->coeff[2].assign_one();
	this->deg = 2;
}

template <class T>
void polynomial<T>::assign(const T & x) {
	this->coeff[0].assign(x);
	if (x == 0)
 	  this->deg = -1;
 	else
 	  this->deg = 0; 
}

template <class T>
void polynomial<T>::assign(const polynomial <T> &f) {
	this->coeff = f.coeff;
	this->deg = f.deg;
}

//////////////////////////////////////// deported basic methods 

// Checks if some coefficients have become zero after an operation, and adjusts degree
template <class T>
void polynomial <T>::remove_leading_zeros() {
	// nothing to do for zero polynomials
	if (this->deg < 0) {
		// check that all coeffs are nul
		assert (find_if(this->coeff.begin(), this->coeff.end(), std::mem_fun_ref(&T::is_zero)) != this->coeff.end());
		return ;
	} 
	typename storage_type::const_iterator i = this->coeff.begin() + this->deg;
	while ( this->deg >= 0 && i->is_zero() ) {
		this->deg--; i--;
	}
}

//////////////////////////////////////// deported access methods

// First non zero coeff of polynomial
template <class T>
T polynomial <T>::lead_coeff() const
{
  T tmp(0);
  
  for (rpl_size_t i = this->degree(); i >= 0; i--)
    if (!this->coeff[i].is_zero())
      {
				tmp = this->coeff[i];
				break;
      }
  
  return tmp;
}

//////////////////////////////////////// deported operations

//////////////////////// divide

// Divide the coefficients of a by b and store in c
template <class T>
void divide(polynomial <T> &c, const polynomial <T> &a, const T &b)
{
  assert (b != 0);
	c.deg = a.deg;
  for (rpl_size_t i = 0; i <= c.deg; i++)
    div_exact(c.coeff[i], a.coeff[i], b);
}

template <class T>
inline polynomial <T> operator/ (const polynomial <T> &a, const T &b) {
	polynomial <T> res;
	divide(res, a, b);
	return res;
}

// Pseudo-division algorithm 
// finds (q,r) such that d^{m-n+1} a = b * q + r, where d 
// is the leading coefficient of b, m = deg(a), n = deg(b)
template <class T>
void div_rem(polynomial <T> &q, polynomial <T> &r, 
							const polynomial <T> &a, const polynomial <T> &b) {
	if (a.degree() < b.degree())
		{
			q.assign_zero();
			r = a;
			return;		 
		}						 

	rpl_size_t m = a.deg;
	rpl_size_t n = b.deg;
	rpl_size_t e = m - n +1;
	assert(!b.is_zero());
		 
	// b is a constant
	if (b.deg == 0)	{
		assert(!b.coeff[0].is_zero());
		bigint expo;
		power(expo, b.lead_coeff(), m);
		q.multiply(a, expo);
		r.assign_zero();
		return;
	} 
	
	r.assign(a);
	q.assign_zero();
	T d = b.lead_coeff();
	
	// non-optimized algorithm!
	while (r.deg >= b.deg) {
		// s <- l(r) * x^{deg(r)-deg(b)}
		polynomial <T> s;
		s.coeff[r.deg-b.deg].assign(r.lead_coeff());
		s.deg = r.deg-b.deg;
		// q <- d * q + s
		q.multiply(q, d);
		q.add(q,s);
		// r <- d * r - s * b
		r.multiply(r, d);
		polynomial <T> tmp;
		tmp.multiply(s, b);
		r.subtract(r, tmp);
		e--;
	}
	
	T ti_cul;
	power(ti_cul, d, e);
	q.multiply(q, ti_cul);
	r.multiply(r, ti_cul);
}

template <class T>
inline void divide(polynomial <T> &q, const polynomial <T> &a, const polynomial <T> &b) {
	polynomial <T> r;
	div_rem(q, r, a, b);
} 

//////////////////////// multiply

// Multiply the coefficients of a by b and store in c
template <class T>
void multiply(polynomial <T> &c, const polynomial <T> &a, const T &b)
{
	if (b.is_zero())
		c.assign_zero();
	else
		{
			c.deg = a.deg;
  		for (rpl_size_t i = 0; i <= c.deg; i++)
    		multiply(c.coeff[i], a.coeff[i], b);
    }
}

template <class T>
polynomial <T> operator* (const polynomial <T> &f, const T &g) {
	polynomial <T> res;
	multiply(res, f, g);
	return res;
}

template <class T>
polynomial <T> operator* (const T &f, const polynomial <T> &g) {
	polynomial <T> res;
	multiply(res, g, f);
	return res;
}

// Multiply a by b
template <class T>
void polynomial<T>::multiply(const polynomial <T> &a, const polynomial <T> &b) {
	if (a.deg < 0 || b.deg < 0) {
		this->assign_zero();
		return;
	}
	this->deg = a.deg + b.deg;

	typename polynomial<T>::storage_type tmp;
	// implements a naive convolution
	for (int i = 0 ; i <= this->deg; ++i) {
		tmp[i].assign_zero();			
		for(int j = 0; j<= i; ++j) {
			tmp[i] += a[j] * b[i-j]; 
		}				 
	}
	
	this->coeff = tmp; 
	remove_leading_zeros();
}

template <class T> 
inline void multiply(polynomial <T> &a, const polynomial <T> &b, const polynomial <T> &c) {
	a.multiply(b,c);
}

//////////////////////// negate

// Negate a polynomial
template <class T>
void polynomial <T>::negate(const polynomial <T> &h) {
	assign(h);
	std::for_each(this->coeff.begin(), this->coeff.end(), std::mem_fun_ref(&T::negate));
}

template <class T>
inline void negate(polynomial <T> &a, const polynomial <T> &b) {
	a.negate(b);
}

//////////////////////// add

// Add two polynomials
template <class T>
void polynomial <T>::add(const polynomial <T> &g, const polynomial <T> &h) {
	if (g.deg < 0) {
		this->assign(h);
	}
	else if (h.deg < 0) {
		this->assign(g);
	}
	else {
		std::transform(h.coeff.begin(), h.coeff.end(), g.coeff.begin(), this->coeff.begin(), 
											std::plus<T>());
		this->deg = std::max(h.deg, g.deg);
		remove_leading_zeros();
	}
}

template <class T> 
inline void add(polynomial <T> &a, const polynomial <T> &b, const polynomial <T> &c) {
	a.add(b,c);
}

template <class T>
polynomial <T> operator+ (const polynomial<T> & f, const polynomial<T> & g) {
	polynomial<T> res;
	res.add(f, g);
	return res;
}

//////////////////////// subtract

// Subtract two polynomials
template <class T>
void polynomial <T>::subtract(const polynomial <T> &g, const polynomial <T> &h) {
	std::transform(g.coeff.begin(), g.coeff.end(), 
									h.coeff.begin(), this->coeff.begin(), std::minus<T>());
	this->deg = std::max(h.deg, g.deg);
	remove_leading_zeros();
}		

template <class T> 
inline void subtract(polynomial <T> &a, const polynomial <T> &b, const polynomial <T> &c) {
	a.subtract(b,c);
}

template <class T>
polynomial <T> operator- (const polynomial <T> &f, const polynomial <T> &g) {
	polynomial<T> res;
	res.subtract(f, g);
	return res;
}

//////////////////////// power

template <class T>
void polynomial <T>::power(const polynomial <T> &g, rpl_size_t e) {
	switch(e) {
		case 0 : this->assign_one();
						 break;
		case 1 : this->assign(g);
						 break;
		case 2 : this->multiply(g, g);
						 break;
		case 3 : this->multiply(g, g); this->multiply(*this, g);
						 break;
	 	case 4 : { polynomial <T> tmp; tmp.multiply(g, g); this->multiply(tmp, tmp); }
						 break;
 		default : std::cerr << "problem with power : e must belong {0,...,4} "
												<< std::endl;
						 assert(0);				
						 break;		
		}
}

template <class T> 
inline void power(polynomial <T> &a, const polynomial <T> &b, rpl_size_t e) {
	a.power(b,e);
}

//////////////////////////////////////// deported comparisons

// Is a polynomial zero?
template <typename T>
bool polynomial <T>::is_zero() const { 
	return ((this->deg == -1) || (this->lead_coeff().is_zero()));
}

// Is a polynomial equal to 1?
template <typename T>
bool polynomial <T>::is_one() const { 
	return ((this->deg == 0) && (this->lead_coeff().is_one()));
}

// Are two polynomials equal?
template <class T>
bool operator == (const polynomial <T> &f, const polynomial <T> &g) {
  return (f-g).is_zero();
}  

////////////////////////////////////// pretty printing

template <class T>
std::ostream & operator << (std::ostream &o, const polynomial <T> &f) {
	for (int i=f.degree(); i > 0 ; i--) {
		if (!f[i].is_zero() ) {
		 	if (i!= f.degree() && f[i].sign() > 0) { 
				if (!f[i].is_one()) 
					o << "+" << f[i] << "*x^"<< i;
				else
					o << "+" << "x^"<< i;
			}
	 		else {
				if (!f[i].is_one()) 
					o << f[i] << "*x^"<< i;
				else
					o << "x^"<< i;
			}
		}
	}
	if ((f.degree() == -1) || (!f[0].is_zero())) {
		if (f[0].sign() > 0) 
			o << "+" ;
		o << f[0] ;
	}
	
	return o;
}

////////////////////////////////////// declarations when T = bigint
////////////////////////////////////// see polynomial.cpp for implementation
		
bigint cont(const polynomial <bigint> &a);

polynomial <bigint> pp(const polynomial <bigint> &a);

polynomial <bigint> gcd(const polynomial <bigint> &aa,
												const polynomial <bigint> &bb);

};

#endif
