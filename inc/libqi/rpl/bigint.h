#ifndef _rpl_bigint_h_
#define _rpl_bigint_h_

#include <gmpxx.h>
#include <iostream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>

//! \file bigint.h

namespace rpl {
extern boost::mt19937 gen; //!< the random number generator

//! multiprecision integer arithmetic class
/*! class implementing big integers
*  Roughly, we put a interface around the gmp integer functions \n
*  Thus we embed an mpz_class value into the struct and we just \n
*  call the corresponding gmp functions on it.
*/
struct bigint {
    mpz_class val; /*!< we aggregate an mpz in bigint*/
    
    //! default constructor
    bigint(): val(0) {} 
    
    //! destructor 
    ~bigint() {}

    //! copy constructor
    /*! we use the initialization list of mpz_class
     * to implement it 
     */
    bigint(const bigint & b): val(b.val) {}
    
    bigint(const char * c): val(c) {}
    
    //! constructor from mpz_t  
    explicit bigint(const mpz_t & z): val(mpz_class(z)) {} 
    
    //! constructor from int  
    bigint(int i); 
    
    //! constructor from long  
    bigint(long l); 
    
    //! constructor from unsigned long  
    bigint(unsigned long ul);

    //! constructor from double
    bigint(double d);

    //! assignment operator
    bigint & operator = (const bigint & a) { 
				if (this != &a) // Beware of self-assignment
		 			this->assign(a);

				return *this;
		}
 
    //Assignments 

    //!\f$a \leftarrow 0\f$ 
    /*! assign zero to the bigint  
     */
    void assign_zero() { val = 0; }
    
    //! \f$a \leftarrow 1\f$
    /*! assign one to the bigint 
     */
    void assign_one() { val = 1; }
    
    //! \f$a \leftarrow n\f$ 
    /*! assign a bigint 
     */
    void assign(const bigint & n) { val = n.val; }
    
    // objects modifiers

    //! \f$a \leftarrow -a\f$
    /*! in-place computation of the opposite 
     */
    void negate() { mpz_neg(val.get_mpz_t(),val.get_mpz_t()); }
    
    //! \f$a \leftrightarrow b\f$ 
    /*! swapping with the argument 
     */
    void swap(bigint & b) { mpz_swap(val.get_mpz_t(),b.val.get_mpz_t()); }
    
    //! in-place computation of absolute value
    /*! \f$a \leftarrow |a|\f$
     */
    void abs() { mpz_abs(val.get_mpz_t(),val.get_mpz_t()); }

    //! generate a random value  
    /*! \f{eqnarray*} a \leftarrow r \in   \{0,\cdots , b − 1\}  & \mbox{ if }  & b > 0  \\
     *                             r \in   \{b+1,\cdots, 0\}   & \mbox{ if }  & b < 0 \\
     *                                     \mbox{ throw error }  & \mbox{ if }  & b= 0 \f}
     */
    void randomize(const bigint & b); 
   
    // Arithmetical operators
    
    //! \f$ a \leftarrow a + 1 \f$
    /*! prefix increment operator
     *  \return a +1 
     */
    bigint & operator++() { val++; return *this; }
     
    //! \f$ a \leftarrow a + 1 \f$
    /*! prefix increment operator
     *  \return a +1 
     */
    // FIXME: the "int x" in input is not used, it is normal? 
   //  YES! it is a idiosyncrasy of C++ to distinguish post and pre increment operator
    bigint & operator++(int x) { bigint & tmp =  *this; val++; return tmp;} 

    // arithmetical operators
    //! \f$ a \leftarrow a -1 1 \f$
    /*! prefix decrement operator
     *  \return a +1 
     */
    bigint & operator--() { val--; return *this; }
     
    //! \f$ a \leftarrow a - 1 \f$
    /*! prefix decrement operator
     *  \return a +1 
     */
    /* FIXME: the "int x" in input is not used, it is normal? */
    bigint & operator--(int x) { bigint & tmp =  *this; val--; return tmp;} 

    //! \f$ a \leftarrow a + b \f$
    /*! addition unary operator
     */
    bigint & operator+= (const bigint & b) { val += b.val; return *this;}

    //! \f$ a \leftarrow a - b \f$
    /*! subtraction unary operator
     */
    bigint & operator-= (const bigint & b) { val -= b.val; return *this;}

    void add(const bigint & a, const bigint & b);
    void subtract(const bigint & a, const bigint & b);
    void multiply(const bigint & a, const bigint & b);
    
    //!  \f$ [a = 0] \f$ 
    bool is_zero () const { return cmp(val, 0) == 0; } 
    
    //! \f$ [ a = 1 ] \f$ 
    bool is_one () const { return cmp(val, 1) == 0; }
    
    //! \f$ [a  > 0 ] \f$ 
    bool is_gt_zero () const { return sgn(val) > 0; }

    //! \f$ [a  < 0 ] \f$ 
    bool is_lt_zero () const { return sgn(val) < 0; }
    
    //! test if the bigint fits on a unsigned int 
    bool is_uint () const { return mpz_fits_uint_p(val.get_mpz_t()); }
   
    //! test if the bigint is even 
 		bool is_even () const { return mpz_even_p(val.get_mpz_t()); }

    //! test if the bigint is a square 
  	bool is_square () const { return mpz_perfect_square_p (val.get_mpz_t()); }
    
    //! returns the sign of the bigint
    /*!
     * \return \f{eqnarray*} \frac{x}{|x|} & \mbox{ if }   x \neq 0 \\   0 & \mbox{ otherwise } \f}
     * */
    int sign() const { return sgn(val); }
};

//////////////////////////////////////////////
//  Deported members and non member functions 
//////////////////////////////////////////////


inline
bigint::bigint (int i) {
	val = mpz_class(i);
	//mpz_init_set_si(val.get_mpz_t(), i);
}

inline
bigint::bigint (long l) {
	val = mpz_class(l);
}

inline
bigint::bigint (unsigned long ul) {
	val = mpz_class(ul);
}

inline bigint::bigint (double d) {
	val = mpz_class(d);
}



inline 
void bigint::randomize(const bigint & b) {
    signed long int  b_cast = mpz_get_si(b.val.get_mpz_t());
    if (b_cast > 0) {
       boost::uniform_int<> dist(0, b_cast-1);
       val = dist(gen);
    }
    else if (b_cast < 0) {
       boost::uniform_int<> dist(b_cast+1,0);
        val = dist(gen);
    }
    else {
       std::cerr << "Error calling bigint.randomize() with 0 as argument"
                << std::endl;
       assert(b_cast != 0);		
    }
}

//! computes the remainder (modulo) of a / b
/*! \relates bigint 
 * computes r such that \f$a = b * q + r \mbox{ with } |r| < |b| \,,\, sgn(r) = sgn(a)  \f$
*
*/
template <class T>
inline 
void remainder(bigint & c, const bigint & a, const T & b) {
  mpz_tdiv_r(c.val.get_mpz_t(), a.val.get_mpz_t(), bigint(b).val.get_mpz_t());
}

template <>
inline 
void remainder(bigint & c, const bigint & a, const bigint & b) {
  mpz_tdiv_r(c.val.get_mpz_t(), a.val.get_mpz_t(), b.val.get_mpz_t());
}


//! *this \f$ \leftarrow a + b \f$ 
inline
void bigint::add(const bigint & a, const bigint & b) {
   mpz_add(this->val.get_mpz_t(),a.val.get_mpz_t(), b.val.get_mpz_t());
}

//! *this \f$ \leftarrow a - b \f$ 
inline
void bigint::subtract(const bigint & a, const bigint & b) {
   mpz_sub(this->val.get_mpz_t(),a.val.get_mpz_t(), b.val.get_mpz_t());
}

//! *this \f$ \leftarrow a * b \f$ 
inline
void bigint::multiply(const bigint & a, const bigint & b) {
   mpz_mul(this->val.get_mpz_t(),a.val.get_mpz_t(), b.val.get_mpz_t());
}

// functional versions

//!  \f$ a \leftarrow b+c \f$ 
/*! \relates bigint
 *   functional version of the addition method
 */  
inline
void add(bigint & a, const bigint & b, const bigint & c) {
   mpz_add(a.val.get_mpz_t(), b.val.get_mpz_t(), c.val.get_mpz_t());
}

//!  \f$ a  \leftarrow  b - c \f$ 
/*! \relates bigint
 *   functional version of the subtraction method
 */  
inline
void subtract(bigint & a, const bigint & b, const bigint & c) {
   mpz_sub(a.val.get_mpz_t(), b.val.get_mpz_t(), c.val.get_mpz_t());
}

//!  \f$ a \leftarrow b * c \f$ 
/*! \relates bigint
 *   functional version of the multiply method
 */  
inline
void multiply(bigint & a, const bigint & b, const bigint & c) {
   mpz_mul(a.val.get_mpz_t(), b.val.get_mpz_t(), c.val.get_mpz_t());
}

//!  \f$ a \leftarrow b / c \f$ 
/*! \relates bigint
 *   computes the quotient in the euclidean division of b by c 
 */  
inline 
void divide(bigint & a, const bigint & b, const bigint & c) {
    assert(!c.is_zero());
    mpz_tdiv_q(a.val.get_mpz_t(), b.val.get_mpz_t(), c.val.get_mpz_t());
}    

//! great common divisor
/*! \relates bigint 
* returns the greatest common divisor (gcd) of a and b as a positive bigint.
*/
inline 
bigint gcd(const bigint & a, const bigint & b) {
   bigint tmp;
   mpz_gcd(tmp.val.get_mpz_t(), a.val.get_mpz_t(), b.val.get_mpz_t());
   return tmp;
}


//! extended gcd function
/*! \relates bigint 
* \param[out] u
* \param[out] v
* \param[in] a 
* \param[in] b 
* \returns \f$c =  a \wedge  b\f$
* computes the gcd of a and  b and computes u and  v
* such that \f$ c = u \times a + v \times b \f$
*/

inline 
bigint xgcd (bigint & u, bigint & v, const bigint & a, const bigint & b) {
   bigint c;
   mpz_gcdext(c.val.get_mpz_t(), u.val.get_mpz_t(), v.val.get_mpz_t(), 
              a.val.get_mpz_t(), b.val.get_mpz_t());
   return c;
}

//! \f$q = a / b \f$ and \f$ r \cong a [b] \f$ 
/* * \relates bigint 
* computes the euclidean division of a by b 
* \param[out] q the quotient
* \param[out] r the remainder (|r| < |b|)
* \param[in]  a the dividend 
* \param[in]  b the divisor
*/
inline 
void div_rem (bigint & q, bigint & r, const bigint & a, const bigint & b) {
   mpz_tdiv_qr(q.val.get_mpz_t(), r.val.get_mpz_t(), a.val.get_mpz_t(), b.val.get_mpz_t());
}


//! \f$ q \leftarrow \lfloor a/b \rfloor \f$ 
/*! \relates bigint
  The remainder r has the same sign as the divisor b.
*/
inline 
void floor_divide (bigint & q, const bigint & a, const bigint & b) {
   mpz_fdiv_q(q.val.get_mpz_t(), a.val.get_mpz_t(), b.val.get_mpz_t());

}


//! \f$ q \leftarrow  a/b  \f$ 
/*! \relates bigint
 * computes the exact division of a by b : we assume
 * that b divides a
*/
inline 
void div_exact(bigint & q, const bigint & a, const bigint & b) {
    mpz_divexact(q.val.get_mpz_t(), a.val.get_mpz_t(), b.val.get_mpz_t());
}

//! \f$ [b | a] \f$
/*! \relates bigint 
 * returns a non zero value if the unsigned long b 
 * divides a
 */
inline
int divisible_ui_p (const bigint & a, unsigned long b) {
  return mpz_divisible_ui_p (a.val.get_mpz_t(), b);
}


//! \f$ a \leftarrow \lfloor \sqrt{b} \rfloor \f$ 
/*! \relates bigint 
 * puts in \c a the <b> truncated integer part</b> of square root of \c b 
*/
inline 
void floor_sqrt (bigint & a, const bigint & b) {
	if (mpz_sgn(b.val.get_mpz_t()) < 0) {
		std::cerr << "computing a square root of a negative integer" << std::endl;
                assert(0);
	}
	else {
		mpz_sqrt(a.val.get_mpz_t(), b.val.get_mpz_t());
	}
}

//! \f$ c \leftarrow a^e \f$ 
/*! \relates bigint 
*  exponentiation function
*/
inline
void power(bigint & c, const bigint & a, size_t e) {
    mpz_pow_ui(c.val.get_mpz_t(), a.val.get_mpz_t(), e);
}

//! \f$ return \leftarrow |a|\f$
/*! \relates bigint 
 * functional version computation of absolute value
 */
inline
bigint abs(const bigint & a) {
   bigint tmp = a; 
   tmp.abs();
   return tmp;
}

//! test if the bigint is a square:
/*! if no:  put 0 in r and return false
 * if yes: put sqrt(a) in r and return true
 */
inline
bool is_square (bigint & r, const bigint & a) { 
  if (!mpz_perfect_square_p (a.val.get_mpz_t()))
    {
      r = 0;
      return false;
    }
  else
    {
      mpz_sqrt (r.val.get_mpz_t(), a.val.get_mpz_t());
      return true;
    }
}


//! convert a bigint to a string 
/*! \relates bigint  
 * \returns a char * containing
 * the bigint written under decimal form
 */
inline char * bigint_to_string (const bigint &a)
{
   return mpz_get_str(NULL,10,a.val.get_mpz_t());
}



/////////////////////////////////////
// Arithmetical operations
//////////////////////////////////////


//! \f$ c \leftarrow -a\f$
/*! \relates bigint
 * computes the opposite of \c a
*/
inline 
bigint operator-(const bigint & a) {
   bigint tmp = a;
   tmp.negate();
   return tmp;
}

//! \f$ c \leftarrow a + b \f$
/*! \relates bigint
* \param[out] c first operand
* \param[in] a second operand
* \param[in] b third operand
*/
template <class T>
inline
void add (bigint & c, const bigint & a, const T & b) {
     c.val =  a.val + bigint(b).val;
}

template <>
inline
void add (bigint & c, const bigint & a, const bigint & b) {
     c.val =  a.val + b.val;
}

inline
void negate(bigint & c, const bigint & a) {
     c = a;
     c.negate();
}


//! computes the sum of two bigints
/*! \relates bigint
* \param a first operand
* \param b second operand
* \return the sum 
*/
inline bigint operator+(const bigint & a, const bigint & b) {
    bigint tmp(a);
    tmp.val += b.val;
    return tmp;
}

//! computes the difference of two bigints
/*! \relates bigint
* \param a first operand
* \param b second operand
* \return the difference 
*/
inline bigint operator-(const bigint & a, const bigint & b) {
    bigint tmp(a);
    tmp.val -= b.val;
    return tmp;
}

//! computes the product of two bigints
/*! \relates bigint
* \param a first operand
* \param b second operand
* \return the product
*/
inline bigint operator*(const bigint & a, const bigint & b) {
    bigint tmp(a);
    tmp.val *= b.val;
    return tmp;
}

//! \f$ q = a / b \f$ 
/*! \relates bigint
* \param a first operand
* \param b second operand
* \return the quotient
*/
inline bigint operator/(const bigint & a, const bigint & b) {
    bigint tmp;
//    mpz_tdiv_q(tmp.val.get_mpz_t(), a.val.get_mpz_t(), b.val.get_mpz_t());
    div_exact(tmp,a,b);
    return tmp;
}

//! equality operator
/*! \relates bigint
* \returns true if a=b
*/
inline bool operator==(const bigint & a, const bigint & b) {
    return (a.val == b.val);
}

//! different operator
/*! \relates bigint
* \returns true if a!=b
*/
inline bool operator!=(const bigint & a, const bigint & b) {
    return (a.val != b.val);
}

//!\f$[ a \leqslant b ]\f$ 
/*! \relates bigint
*! less or equal operator
* \returns true if a<=b
*/
inline bool operator<=(const bigint & a, const bigint & b) {
    return (a.val <= b.val);
}

//!\f$[ a < b ]\f$ 
/*! \relates bigint
* less than operator
* \returns true if a<b
*/
inline bool operator<(const bigint & a, const bigint & b) {
    return (a.val < b.val);
}

//! \f$[ a \geqslant b ]\f$
/*! \relates bigint
*! greater or equal operator
* \returns true if a>=b
*/
inline bool operator>=(const bigint & a, const bigint & b) {
    return (a.val >= b.val);
}

//! \f$[ a > b ]\f$
/*! \relates bigint
* ! greater than operator 
* \returns true if a>b
*/
inline bool operator>(const bigint & a, const bigint & b) {
    return (a.val > b.val);
}

//! ostream insertion operator 
/*! \relates bigint
*  flow insertion operator
*  \param o the flow
*  \param i bigint to display
*  \return the modified flow
*/
inline std::ostream & operator<< (std::ostream & o, const bigint & i) {
       o << i.val ;
       return o;
}       
//! istream extraction operator 
/*! \relates bigint
*  flow extraction operator
*  \param is the flow
*  \param i bigint to extract 
*  \return the modified flow
*/
inline std::istream & operator>> (std::istream & is, bigint & i) {
       mpz_t tmp;
       is >> tmp;
       i = bigint(tmp);
       return is;
}       

/*! \f$ a \leftrightarrow b \f$
 * \relates bigint
*  permutes \c a and \c b 
*/
inline void swap(bigint & a , bigint & b) {
      mpz_swap(a.val.get_mpz_t(), b.val.get_mpz_t());       
}

};
#endif

