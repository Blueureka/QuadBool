#ifndef __BASE_VECTOR_
#define __BASE_VECTOR_

#include <libqi/rpl/rpl.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

//! \file base_vector.h


using namespace boost::numeric::ublas;

namespace rpl {

//! base class for math vector  
/*! 
* this class is just an alias on vectors from ublas, with a static allocator
*/
template <class T>
struct base_vector : boost::numeric::ublas::vector<T,bounded_array<T,11> >  {
  // we need only of order 8 vector -> 8 elements
     
     typedef  bounded_array<T,11> storage_type;
     typedef  boost::numeric::ublas::vector<T, storage_type > type;
     
     base_vector(rpl_size_t size = 0)
              : type(size) {}
     
     base_vector(rpl_size_t c, rpl_size_t s)
              : type(s) { if (s != c) { std::cerr << "problem in base_vector cstr size != capa" << std::endl; assert(-1); } } 
     
     base_vector(const base_vector<T> & a)
              : type(a) {}
     
     template <class VE>
     base_vector(const vector_expression<VE> & ve)
              : type(ve) {}
     
     base_vector(const T *a, size_t l) 
              : type(l) { std::copy(a, a+l, type::begin()); } 
     
     base_vector <T> & operator= (const base_vector<T>& a) 
     	      { 
     	      	if (this != &a) // Beware of self-assignment
					type::operator=(a); 
				return *this;}

     ~base_vector() {}
};

};// end namespace
#endif
