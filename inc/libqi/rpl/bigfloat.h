#ifndef _rpl_bigfloat_h_
#define _rpl_bigfloat_h_

#include <gmpxx.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include <libqi/rpl/bigint.h>

// \file bigfloat.h

namespace rpl {
//! multiprecision floating-point  class
/*! class implementing big integers
*  Roughly, we put a interface around the gmp floating-point functions \n
*  Thus we embed an mpf_class value into the struct and we just \n
*  call the corresponding gmp functions on it.
*/
struct bigfloat {
    mpf_class val; /*!< we aggregate an mpf in bigfloat*/
   
    static long precision;

    static void set_precision(long p) { precision = static_cast<long>(ceil(p*log(10.)/log(2.)) + 3); 
                                        mpf_set_default_prec(precision); } 
    
    // Coplien's canonical form 
    //! default constructor
    bigfloat(): val() {}  //call mpf_init
    
    //! destructor 
    ~bigfloat() {}

    //! copy constructor
    /*! we use the initialization list of mpf_class
     * to implements it 
     */
    bigfloat(const bigfloat & b): val(b.val) {}
     
    bigfloat(const char * c): val(c) {}
    
    bigfloat(const bigint & i): val(i.val) {}
   
    //! template constructor 
    /*! we inherit it from mpf_class constructor 
     */
    bigfloat(int i); 
    
    bigfloat(long l); 
    
    bigfloat(unsigned long ul);

    bigfloat(double d);

    //! assignement operator
    bigfloat & operator=(const bigfloat & a) { if (this != &a) // Beware of self-assignment
        val = a.val; return *this; }

   //Assignements 

    //! affect to zero
    /*! \f$a \leftarrow 0\f$
     */
    void assign_zero() { val = 0; }
    
    //! affect to one 
    /*! \f$a \leftarrow 1\f$
     */
    void assign_one() { val = 1; }
    
    //! affect to integer
    /*! \f$a \leftarrow n\f$
     */
    void assign(const bigfloat & n) { val = n.val; }
    
    void assign(int i);

    void assign(long i);

    void assign(unsigned long l);

    void assign(double d);


    // objects modifiers

    //! in-place computation of the opposite 
    /*! \f$a \leftarrow -a\f$
     */
    void negate() { mpf_neg(val.get_mpf_t(),val.get_mpf_t()); }
    
    //! swapping with the argument 
    /*! \f$a \leftrightarrow b\f$
     */
    void swap(bigfloat & b) { mpf_swap(val.get_mpf_t(),b.val.get_mpf_t()); }
    
    //! in-place computation of absolute value
    /*! \f$a \leftarrow |a|\f$
     */
    void abs() { mpf_abs(val.get_mpf_t(),val.get_mpf_t()); }

    // arithmetical operators
    //! \f$ a \leftarrow a + 1 \f$
    /*! prefix increment operator
     *  \return a +1 
     */
    bigfloat & operator++() { val++; return *this; }
     
    //! \f$ a \leftarrow a + 1 \f$
    /*! prefix increment operator
     *  \return a +1 
     */
    bigfloat & operator++(int x) { bigfloat & tmp =  *this; val++; return tmp;} 


    //! \f$ a \leftarrow a + b \f$
    /*! addition unary operator
     */
    bigfloat & operator+= (const bigfloat & b) { val += b.val; return *this;}

    //! \f$ a \leftarrow a - b \f$
    /*! subtraction unary operator
     */
    bigfloat & operator-= (const bigfloat & b) { val -= b.val; return *this;}

    void add(const bigfloat & a, const bigfloat & b);
    void subtract(const bigfloat & a, const bigfloat & b);
    void multiply(const bigfloat & a, const bigfloat & b);

    //!  \f$ [a = 0] \f$ 
    bool is_zero () const { return cmp(val, 0) == 0; } 
    
    //! \f$ [ a = 1 ] \f$ 
    bool is_one () const { return cmp(val, 1) == 0; }
    
    //! \f$ [a  > 0 ] \f$ 
    bool is_gt_zero () const { return sgn(val) > 0; }

    //! \f$ [a  < 0 ] \f$ 
    bool is_lt_zero () const { return sgn(val) < 0; }
    
    //! test if the bigfloat fits on a unsigned int 
    bool is_uint () const { return mpf_fits_uint_p(val.get_mpf_t()); }
    
    //! returns the sign of the bigfloat
    /*!
     * \return \f{eqnarray*} \frac{x}{|x|} & \mbox{ if }   x \neq 0 \\   0 & \mbox{ otherwise } \f}
     * */
    int sign() const { return sgn(val); }
};

//////////////////////////////////////////////
//  Deported members and non member functions 
//////////////////////////////////////////////


inline
bigfloat::bigfloat (int i) {
	val = mpf_class(i);
}

inline
bigfloat::bigfloat (long l) {
	val = mpf_class(l);
}

inline
bigfloat::bigfloat (unsigned long ul) {
	val = mpf_class(ul);
}

inline bigfloat::bigfloat (double d) {
	val = mpf_class(d);
}




//! assignement method instantied for bigfloat case
inline
void bigfloat::assign(int i) { 
    val = mpf_class(i);
     
}

inline
void bigfloat::assign(long l) { 
    val = mpf_class(l);
     
}

inline
void bigfloat::assign(unsigned long ui) { 
    val = mpf_class(ui);
}

inline
void bigfloat::assign(double d) { 
    val = mpf_class(d);
     
}

//! *this \f$ \leftarrow a + b \f$ 
inline
void bigfloat::add(const bigfloat & a, const bigfloat & b) {
   mpf_add(this->val.get_mpf_t(),a.val.get_mpf_t(), b.val.get_mpf_t());
}

//! *this \f$ \leftarrow a - b \f$ 
inline
void bigfloat::subtract(const bigfloat & a, const bigfloat & b) {
   mpf_sub(this->val.get_mpf_t(),a.val.get_mpf_t(), b.val.get_mpf_t());
}

//! *this \f$ \leftarrow a * b \f$ 
inline
void bigfloat::multiply(const bigfloat & a, const bigfloat & b) {
   mpf_mul(this->val.get_mpf_t(),a.val.get_mpf_t(), b.val.get_mpf_t());
}

// functional versions

//!  \f$ a \leftarrow b+c \f$ 
/*! \relates bigfloat
 *   functional version of the addition method
 */  
inline
void add(bigfloat & a, const bigfloat & b, const bigfloat & c) {
   mpf_add(a.val.get_mpf_t(), b.val.get_mpf_t(), c.val.get_mpf_t());
}

//!  \f$ a  \leftarrow  b - c \f$ 
/*! \relates bigfloat
 *   functional version of the subtraction method
 */  
inline
void subtract(bigfloat & a, const bigfloat & b, const bigfloat & c) {
   mpf_sub(a.val.get_mpf_t(), b.val.get_mpf_t() ,c.val.get_mpf_t());
}

//!  \f$ a \leftarrow b * c \f$ 
/*! \relates bigfloat
 *   functional version of the multiply method
 */  
inline
void multiply(bigfloat & a, const bigfloat & b, const bigfloat & c) {
   mpf_mul(a.val.get_mpf_t(), b.val.get_mpf_t() ,c.val.get_mpf_t());
}

//! \f$ a \leftarrow \lfloor b \rfloor \f$ 
/*! \relates bigfloat 
 * puts in  a  the bigint just before of  b 
*
*/
inline 
void floor (bigint & a, const bigfloat & b) {
     bigfloat tmp;
     mpf_floor(tmp.val.get_mpf_t(), b.val.get_mpf_t());
     a.val = tmp.val;
}

//! \f$ a \leftarrow \lceil b \rceil \f$ 
/*! \relates bigfloat 
 * puts in  a  the integer just after  b 
*
*/
inline 
void ceil (bigint & a, const bigfloat & b) {
     bigfloat tmp;
     mpf_ceil(tmp.val.get_mpf_t(), b.val.get_mpf_t());
     a.val = tmp.val;
}
//! \f$ a \leftarrow \sqrt{b} \f$ 
/*! \relates bigfloat 
 * puts in  a  the square root of  b 
*
*/
inline 
void sqrt (bigfloat & a, const bigfloat & b) {
	if (mpf_sgn(b.val.get_mpf_t()) < 0) {
		std::cerr << "computing a square root of a negative integer" << std::endl;
                assert(0);
	}
	else {
		mpf_sqrt(a.val.get_mpf_t(), b.val.get_mpf_t());
	}
}

//! \f$ a \leftarrow \sqrt{b} \f$ 
/*! \relates bigfloat 
 * puts in  a  the square root of  b 
*
*/
inline 
bigfloat sqrt (const bigfloat & b) {
        bigfloat tmp;
	sqrt(tmp, b);
	return tmp;
}

//! \f$ c \leftarrow a ^e \f$ 
/*! \relates bigfloat 
*  exponentiation function
*/
inline
void power(bigfloat & c, const bigfloat & a, size_t e) {
    mpf_pow_ui(c.val.get_mpf_t(), a.val.get_mpf_t(), e);
}

//! \f$ return \leftarrow |a|\f$
/*! \relates bigint 
 * functional version computation of absolute value
 */
inline
bigfloat abs(const bigfloat & a) {
   bigfloat tmp = a; 
   tmp.abs();
   return tmp;
}

/////////////////////////////////////
// Arithmetical operations
//////////////////////////////////////


//! \f$ c \leftarrow -a\f$
/*! \relates bigfloat
 * computes the opposite of \c a
*/
inline 
bigfloat operator-(const bigfloat & a) {
   bigfloat tmp = a;
   tmp.negate();
   return tmp;
}

//! \f$ c \leftarrow a + b \f$
/*! \relates bigfloat
* \param[out] c first operand
* \param[in] b second operand
* \param[in] b second operand
*/
template <class T>
inline
void add (bigfloat & c, const bigfloat & a, const T & b) {
     c.val =  a.val + bigfloat(b).val;
}

template <>
inline
void add (bigfloat & c, const bigfloat & a, const bigfloat & b) {
     c.val =  a.val + b.val;
}


//! computes the sum of two bigfloats
/*! \relates bigfloat
* \param a first operand
* \param b second operand
* \return the sum 
*/
inline bigfloat operator+(const bigfloat & a, const bigfloat & b) {
    bigfloat tmp(a);
    tmp.val += b.val;
    return tmp;
}

//! computes the difference of two bigfloats
/*! \relates bigfloat
* \param a first operand
* \param b second operand
* \return the difference 
*/
inline bigfloat operator-(const bigfloat & a, const bigfloat & b) {
    bigfloat tmp(a);
    tmp.val -= b.val;
    return tmp;
}

//! computes the product of two bigfloats
/*! \relates bigfloat
* \param a first operand
* \param b second operand
* \return the product
*/
inline bigfloat operator*(const bigfloat & a, const bigfloat & b) {
    bigfloat tmp(a);
    tmp.val *= b.val;
    return tmp;
}

//! \f$ q = a / b \f$ 
/*! \relates bigfloat
* \param a first operand
* \param b second operand
* \return the quotient
*/
inline bigfloat operator/(const bigfloat & a, const bigfloat & b) {
    bigfloat tmp(a);
    tmp.val /= b.val;
    return tmp;
}

//! equality operator
/*! \relates bigfloat
* \returns true if a=b
*/
inline bool operator==(const bigfloat & a, const bigfloat & b) {
    return (a.val == b.val);
}

//! different operator
/*! \relates bigfloat
* \returns true if a!=b
*/
inline bool operator!=(const bigfloat & a, const bigfloat & b) {
    return (a.val != b.val);
}

//!\f$[ a \leqslant b ]\f$ 
/*! \relates bigfloat
*! less or equal operator
* \returns true if a<=b
*/
inline bool operator<=(const bigfloat & a, const bigfloat & b) {
    return (a.val <= b.val);
}

//!\f$[ a < b ]\f$ 
/*! \relates bigfloat
* less than operator
* \returns true if a<b
*/
inline bool operator<(const bigfloat & a, const bigfloat & b) {
    return (a.val < b.val);
}

//! \f$[ a \geqslant b ]\f$
/*! \relates bigfloat
*! greater or equal operator
* \returns true if a>=b
*/
inline bool operator>=(const bigfloat & a, const bigfloat & b) {
    return (a.val >= b.val);
}

//! \f$[ a > b ]\f$
/*! \relates bigfloat
* ! greater than operator 
* \returns true if a>b
*/
inline bool operator>(const bigfloat & a, const bigfloat & b) {
    return (a.val > b.val);
}

//! ostream insertion operator 
/*! \relates bigfloat
*  flow insertion operator
*  \param o the flow
*  \param i bigfloat to display
*  \return the modified flow
*/
inline std::ostream & operator<< (std::ostream & o, const bigfloat & i) {
       o << i.val ;
       return o;
}       

/*! \f$ a \leftrightarrow b \f$
 * \relates bigfloat
*  permutes \c a and \c b 
*/
inline void swap(bigfloat & a , bigfloat & b) {
      mpf_swap(a.val.get_mpf_t(), b.val.get_mpf_t());       
}

};
#endif

