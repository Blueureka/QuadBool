#ifndef _rpl_h_
#define _rpl_h_

/*! \file rpl.h
 *  \mainpage RPL Library 
 *
 * \section intro_sec Introduction
 *
 * RPL stands for Replacement Package for Lidia : it is an \em ad \em hoc
 * library designed to replace the LiDIA library in
 * the QI (Quadratic Intersections) Project. \n
 * It does \b not aims to be an alternative to LiDIA. It is
 * just a subset of some LiDIA classes, with the minimal
 * interface used in QI project. \n
 * To achieve this, we need to re-implement the following classes : 
 *  - base types
 *    -# \b bigint : arbitrary precision integers 
 *    -# \b bigmod : arbitrary precision \f$\mathbb{Z}/m\mathbb{Z}\f$ integers (QILegendre.cc)  
 *    -# \b bigfloat : arbitrary precision floats  (QIInterWith22.cc,QIOutputter.cc)
 *    -# \b bigrational : arbitrary precision integers (QILegendre.cc)
 *  - composite types
 *    -# \b polynomial<T> : univariate polynomials (in 13 files)  
 *    -# \b base_vector<T> : template type vectors (in 8 files)
 *    -# \b math_vector<T> : template type vectors for numeric (in 11 files)
 *    -# \b base_matrix<T> : template base matrix (QIParamStruct.h, QIQsicStruct.h)
 *    -# \b math_matrix<T> : template math matrix (QILegendre.cc, QIInterWith22.cc)
 *    -# \b bigint_matrix : bigint math matrix (> 13 files)
 *    -# \b rational_factorization : factorizations of rational numbers. (QIVanishDet.cc, QINumber.cc
 *     and QILegendre.cc )
 *    -# \b galois_field : QILegendre.cc
 *    -# \b gf_element : Galois Field elements (QILegendre.cc)
 *
 *  
 *  For rapidity of developpement, we will re-uses some classical libraries :
 *  - <a href="http://boost.org"> Boost</a> (for Matrix, Vectors, Array classes and the Testing Framework)
 *  - <a href="gmplib.org">GMP</a> (To implements  bigint, bigrational, bigfloat, and some mathematicals functions)
 *  
 *  the objective has two aspects : 
 *    - one is the independance from the unmaintened LiDIA libraries. 
 *    - the other comes from licensing. Using only the two previous quoted libraries (under LGPL or other licenses), \n
 *       one could therefore release the QI software as a commercial software.
 *    
 *  re-implementing from scratch some parts like polynomials seems to be naive, compared to the highly optimized \n 
 *  LIDIA library, but we must keep in mind that, it devoted to the QI software, which only uses polynomials of \n
 *  degree 4 at most. So using some compile-time optimizations techniques, we hope to keep an efficient re-implementation \n
 *
 *  As a first prototype, we could rewrite the code ignoring the file QLegendre.cc with is only required in the optimization case. \n
 *  In this case, we just need re-implement the following classes : \n 
 *
 *  - base types
 *    -# \c rpl::bigint : arbitrary precision integers 
 *    -# \c rpl::bigfloat : arbitrary precision floating points 
 *  - composite types
 *    -# \c rpl::polynomial : univariate polynomials 
 *    -# \c rpl::base_vector : template base type vectors 
 *    -# \c rpl::math_vector : template type vectors for numeric 
 *    -# \c rpl::base_matrix : template base type matrices 
 *    -# \c rpl::math_matrix : template math matrix for numeric 
 *    -# \c rpl::bigint_matrix : specialized math matrix for type bigint with some <em> ad hoc </em> linear algebra functions (det, kernel, rank) 
 *    -# \c rpl::rational_factorization : basical factorizations for bigints 
 *   
 *   \section bigint_sec bigint
 *
 *   the implementation of bigint is straightforward : we agregate an mpz_class member (coming from gmpxx.h, interface for
 *   c++ of the GMP library). 
 *   \n Note that, all operations all not implemented in the C++ style. In this case, we call
 *   the member .get_mpz_t() and apply it the corresping C GMP function
 *   
 */

//! the global namespace for RPL 
namespace rpl {

   typedef int rpl_size_t; //!< a typedef replacing lidia_size_t in the original code
}   


#endif
