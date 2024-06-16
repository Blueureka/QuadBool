#ifndef __BIGINT_MATRIX__
#define __BIGINT_MATRIX__


#include "bigint.h"
#include "math_matrix.h"

//! \file bigint_matrix.h

namespace rpl {

struct bigint_matrix : math_matrix<bigint> {
      bigint_matrix() 
                 : math_matrix<bigint> () {}
      bigint_matrix(size_t m, size_t n) 
                 : math_matrix<bigint> (m,n) {}
      bigint_matrix(const bigint_matrix & a)
                 : math_matrix<bigint>(a) {}
      
      bigint_matrix(const base_matrix<bigint> & a)
                 : math_matrix<bigint>(a) {}

      template <class AE>
      bigint_matrix(const matrix_expression<AE> & ae)
                 : math_matrix<bigint>(ae) {}
      
      bigint_matrix & operator=(const bigint_matrix & a);
      
      ~bigint_matrix() {}


      // compose and split functions
      void split_t(bigint_matrix & B,
                   bigint_matrix & C,
                   bigint_matrix & D,
                   bigint_matrix & E) const;

      void split_h(bigint_matrix & B,
                   bigint_matrix & C) const ;

      void compose_t(const bigint_matrix & B,
                     const bigint_matrix & C,
                     const bigint_matrix & D,
                     const bigint_matrix & E);

      void compose_h(const bigint_matrix & B,
                     const bigint_matrix & C);

      //! returns \f$ [ a(i,j) = 0 \, \forall j]  \f$
      /*! 
       *  test if the j\f$^{\mbox{th}}\f$ column is equal to zero
       * */ 
      bool is_column_zero(size_t i) const ;  

      bool is_row_zero(size_t i) const;

      bool is_zero() const;

};

inline
bigint_matrix & bigint_matrix::operator=(const bigint_matrix & a) {
   math_matrix<bigint>::operator=(a);
   return *this;
}

inline 
bool 
bigint_matrix::is_column_zero(size_t j) const {
  if (this->size2() <= 0)
      return true;
  for(size_t i = 0; i < this->size1() ; ++i)
     if (!this->operator()(i,j).is_zero())
         return false;
  return true;	 
}

inline 
bool 
bigint_matrix::is_zero() const {
  if ((this->size2() <= 0) || (this->size1() == 0))
      return true;
  for(size_t i = 0; i < this->size1() ; ++i) {
   for(size_t j = 0; j < this->size2() ; ++j) {  
     if (!this->operator()(i,j).is_zero())
         return false;
    }
  }
  return true;	 
}

inline 
bool 
bigint_matrix::is_row_zero(size_t i) const{
  if (this->size1() <= 0)
      return true;
  for(size_t j = 0; j < this->size2() ; ++j)
     if (!this->operator()(i,j).is_zero())
         return false;
  return true;	 
}

//! \f$ u \leftarrow \ker(A) \f$
/*! \relates  bigint_matrix 
 *  \param[in] A the input matrix
 *  \return the matrix containing the basis of the kernel as columns vectors
 *  computes the kernel of the matrix with coefficients in \f$\mathbb{Z}\f$
 *  by the algorithm of Cohen (cf A Course in Computational Algrebraic 
 *  Number theory p. 73)
 * */ 
bigint_matrix kernel(const bigint_matrix & A);

//! computes \f$ \dim \mathrm{im} (A) \f$
/*! \relates  bigint_matrix 
 *  \param[in] A the input matrix
 *  \return the rank of A, i.e. the dimension of the image of A 
 *  we use the kernel computation for that
 * */ 
size_t rank(const bigint_matrix & A);


//! \f$ u \leftarrow \det(A) \f$
/*! \relates  bigint_matrix 
 *  \param[in] A the input matrix
 *  \return a bigint equal to the determinant of A 
 *  the algorithm used is Gau\ss Barei\ss (see  
 *  A Course in Computational Algrebraic 
 *  Number theory p. 51)
 * */ 
bigint det(const bigint_matrix & A);

//! \f$ u \leftarrow \mathrm{com} (A) \f$
/*! \relates  bigint_matrix 
 *  \param[in] A the input matrix
 *  \return the <b> adjoint matrix </b> or the comatrix of A 
 *  the algorithm used is taken from(  
 *  A Course in Computational Algrebraic 
 *  Number theory p. 53)
 *  beware : the actual algorithm is different :
 *  it only works for n = 3 , we compute all 
 *  the 2x2 determinants <em> a la mano</em>
 * */ 
bigint_matrix adj(const bigint_matrix & A);

}; //end namespace
#endif


