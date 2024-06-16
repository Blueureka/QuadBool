#ifndef __MATH_VECTOR_
#define __MATH_VECTOR_

#include <libqi/rpl/rpl.h>
#include <libqi/rpl/base_vector.h>
#include <numeric>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/typeof/typeof.hpp>

using namespace std;

namespace rpl {
template <typename T>
struct math_vector : base_vector<T> {

     math_vector(rpl_size_t size = 0)
              : base_vector<T>(size) {}
     math_vector(rpl_size_t c, rpl_size_t s)
              : base_vector<T>(s) { if (s != c) { std::cerr << "problem in math_vector cstr size != capa" << std::endl; assert(0); } } 
     math_vector(const math_vector<T> & a)
              : base_vector<T>(a) {}
     math_vector(const base_vector<T> & a)
              : base_vector<T>(a) {}
     template <class VE>
     math_vector(const vector_expression<VE> & ve)
              : base_vector<T>(ve) {}
     math_vector(const T *a, size_t l) 
              : base_vector<T>(l) { std::copy(a, a+l, base_vector<T>::begin()); } 
     
     math_vector & operator= (const math_vector<T>& a) 
     	      { base_vector<T>::operator=(a); return *this;}

     ~math_vector() {}

     rpl_size_t capacity() const { return base_vector<T>::size(); } 

     rpl_size_t get_capacity() const { return base_vector<T>::size(); } 
 
     void set_capacity(size_t cap)  { base_vector<T>::resize(cap); }

     void set_size(size_t size)  { base_vector<T>::resize(size); }
     
     const T operator[] (size_t i) const { return base_vector<T>::operator()(i); };

     T & operator[] (size_t i) { return base_vector<T>::operator()(i); };

     void concat(const math_vector<T> & u, const math_vector<T> & v);

     void assign(size_t at, const math_vector<T> & w, size_t from, size_t to);

     void assign(const math_vector<T> & v) { base_vector<T>::operator=(v); } 

     void add(const math_vector<T> & v1, const math_vector<T> & v2) { * this = v1 + v2; }
     
     void divide(const math_vector<T> & v1, const T & v2) { *this = v1/v2; }

};



template <typename T>
inline
void 
multiply(T & d, const math_vector<T> & v1, const math_vector<T> & v2) {
   d = inner_prod(v1, v2);
}

template <typename T>
inline
T 
sum_of_squares(const math_vector<T> & v) {
   return inner_prod(v, v);
}


template <typename T>
inline 
bool operator==(const math_vector<T> & v1, const math_vector<T> & v2) {
   return std::inner_product(v1.begin(), v1.end(), v2.begin(), true, std::logical_and<bool>(), std::equal_to<T>());     
}

template <typename T>
inline 
void divide(math_vector <T> &res, math_vector<T> & v1, const T & v2) {
   res.divide(v1,v2);
}

template <typename T>
inline 
void negate(math_vector<T> & res, const math_vector<T> & orig) {
  res.resize(orig.size());
  BOOST_AUTO (src, orig.begin());
  BOOST_AUTO (dest, res.begin());
  for(; src != orig.end(); ++src, ++dest) {
        negate(*dest,*src);
  }
}

template <typename T> 
inline
void math_vector<T>::concat(const math_vector<T> & u, const math_vector<T> & v) {
    this->resize(u.size() + v.size());
    project(*this, range(0, u.size())) = u;
    project(*this, range(u.size(), this->size())) = v;
}

template <typename T> 
inline
void math_vector<T>::assign(size_t at, const math_vector<T> & w, size_t from, size_t to) {
    assert(to >= from);
    size_t new_end = to - from + at;
    if (new_end > this->size()) {
         this->resize(new_end);
    }
    project(*this, range(at, new_end+1)) = project(w, range(from, to+1));
}

// Overloading of cout
template <class T>
ostream & operator << (ostream &s, const math_vector <T> &a)
{
  s << "[ ";
  for (unsigned int i = 0; i < a.size(); i++)
  	s << a[i] << " ";
  s << "]";
  return s;
}

}; //end namespace
#endif 
