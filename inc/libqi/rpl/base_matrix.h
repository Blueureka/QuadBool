#ifndef __BASE_MATRIX__
#define __BASE_MATRIX__

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

//! \file base_matrix.h

namespace rpl {
using namespace boost::numeric::ublas; 

//! base class for matrix  
/*! 
* this class is just an alias on matrix from ublas, with a static allocator
*/
template <class T>
struct base_matrix : matrix<T,row_major, bounded_array<T,64> > {
  // we need only of order 8 matrices -> 64 elements
  typedef bounded_array<T,64> storage_type;
  typedef row_major layout;
  typedef matrix<T,layout,storage_type> type;

  base_matrix() 
        : type() {}
  
  base_matrix(size_t m, size_t n)
        : type(m,n) {}
  
  base_matrix(const base_matrix<T> & a)
        : type(a) {}
  
  template <class AE> 
  base_matrix(const matrix_expression<AE> &ae)
        : type(ae) {}
  
  base_matrix & operator= (const base_matrix<T>& a) 
  	      { if (this != &a) // Beware of self-assignment 
  	        type::operator=(a); return *this;}

  ~base_matrix() {}

}; // end class
}; // end namespace
#endif
