#ifndef __MATH_MATRIX__
#define __MATH_MATRIX__


#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <libqi/rpl/base_matrix.h>
#include <libqi/rpl/math_vector.h>

namespace rpl {

template <typename T>
struct math_matrix : base_matrix<T> {
   
     typedef base_matrix<T> base_type;

     math_matrix()
              : base_type() {}
     math_matrix(size_t m, size_t n)
              : base_type(m,n) {}
     math_matrix(const math_matrix<T> & a)
              : base_type(a) {}
     math_matrix(const base_matrix<T> & a)
              : base_type(a) {}
     
     template <class AE> 
     math_matrix(const matrix_expression<AE> &ae)
              : base_type(ae) {}
     
     math_matrix & operator= (const math_matrix<T>& a) 
     	      { if (this != &a) // Beware of self-assignment
     	      base_type::operator=(a); return *this;}

     ~math_matrix() {}

     size_t get_no_of_rows() const { return this->size1(); }
     
     size_t get_no_of_columns() const { return this->size2(); }
 
     void set_no_of_rows(size_t i) { base_type::resize(i,this->size2()); }
     
     void set_no_of_columns(size_t j) { base_type::resize(this->size1(), j); }

     math_vector<T> get_row_vector(size_t i) { return row(*this,i); }

     math_vector<T> get_column_vector(size_t i) { return column(*this,i); }
   
     void swap_rows (size_t i, size_t j) { row(*this,i).swap(row(*this,j)); }

     void swap_columns (size_t i, size_t j) { column(*this,i).swap(column(*this,j)); }
     
     const T member(size_t i, size_t j) const { return base_type::operator()(i,j) ;} 

     void sto(size_t i, size_t j, const T & x) { base_type::operator()(i,j) = x ;} 

     static base_type identity(size_t m) { identity_matrix<T> tmp(m); return tmp; }
     
     static base_type zero(size_t m, size_t n) { zero_matrix<T> tmp(m,n); return tmp; }

     void diag(const T & a, const T & b);

     void negate() { *this = - (*this); };

     T trace();

};


template <typename T>
void math_matrix<T>::diag(const T & a, const T & b) {
  for(size_t i = 0; i < this->size1(); ++i) {
     for (size_t j = 0; j < this->size2() ; ++j )
          this->operator()(i,j) = b;
     this->operator()(i,i) = a;
  }
}

template <typename T>
inline
void multiply(math_vector<T> & v, const math_matrix<T> & a, const math_vector<T> & x) {
   v = prod(a, x);
}   

template <typename T>
inline
void multiply(math_matrix<T> & v, const math_matrix<T> & a, const math_matrix<T> & x) {
   v = prod(a, x);
}   

template <typename T>
inline
void multiply(math_matrix<T> & v, const math_matrix<T> & a, const T & x) {
   v = a * x;
}   

template <typename T>
inline
void add(math_matrix<T> & a, const math_matrix<T> & b, const math_matrix<T> & c) {
   a = b + c;
}   

template <typename T>
inline
T math_matrix<T>::trace() {
  assert( this->size1() == this->size2());
  size_t n = this->size1();
  matrix_vector_slice<matrix<T> > diag (*this, slice (0, 1, n), slice (0, 0, n));
  return std::accumulate(diag.begin(), diag.end(), T(0), std::plus<T>());     
}   

}; // end namespace
#endif 
