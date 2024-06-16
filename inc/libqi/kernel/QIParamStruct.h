// Things concerning curve_param's and surface_param's

#ifndef _qi_param_struct_h_
#define _qi_param_struct_h_

/** RPL */
#include <libqi/rpl/bigint.h>
#include <libqi/rpl/math_matrix.h>
#include <libqi/rpl/bigint_matrix.h>
#include <boost/typeof/typeof.hpp>

/** QI */
#include <libqi/kernel/QIHompoly.h>

using namespace std;
using namespace rpl;
// Enter namespace QI
namespace QI {
  

// Template class for curve parameterizations, i.e. vectors of homogeneous 
// polynomials with bigint coefficients
template <class T>
class curve_param : public math_vector < hom_polynomial <T> >
  {
  public:
    // Warning: math_vector is public virtual math_vector< T > - so get everything 
    // directly from math_vector

    // Default constructor
    curve_param() : math_vector < hom_polynomial <T> > (4,4)
      { }

    // Create a curve_param of size a
    curve_param(const rpl_size_t &a) : math_vector < hom_polynomial <T> > (a,a) 
      { }

    // Create a curve_param from a vector
    curve_param(const math_vector <T> &v) : 
                    math_vector < hom_polynomial <T> > (v.capacity(),v.capacity())
      {
				for (rpl_size_t i = 0; i < v.capacity(); i++)
	  			this->operator[](i).assign((hom_polynomial <T>)v[i]);
      }

    // Parameterize the line through two points
    curve_param(const math_vector <T> &p1, const math_vector <T> &p2);

    // Copy constructor
    curve_param(const curve_param <T> &a): 
                    math_vector < hom_polynomial <T> > (a.capacity(),a.capacity())
      { 
				this->assign(a);
      }

    // Destructor
    ~curve_param()
      { }

    // Copy assignment
    curve_param <T> & operator = (const curve_param <T> &a)
      {
				if (this != &a) // Beware of self-assignment
				  this->assign(a);
			
				return *this;
      }

    // Methods

    // Pretty printing
    void print_verbose(ostream &s, const char x = 'u', const char y = 'v') const;

    // Evaluate a curve_param at a point 
    math_vector <T> eval(const T &a, const T &b) const;

    // Arithmetic operations
    // Multiply a curve_param by a polynomial
    void multiply(const curve_param <T> &a, const hom_polynomial <T> &b);

    // Divide a curve_param by a polynomial
    void divide(const curve_param <T> &a, const hom_polynomial <T> &b);

    // Divide a curve_param by a constant
    void divide(const curve_param <T> &a, const T &b);

    // Multiply a matrix by a curve_param
    void multiply(const math_matrix<T> &a, const curve_param <T> &b);

    // Partial derivative of this curve param
    void derivative(const curve_param <T> &a, const char v);

    // typedefs for iterators
    typedef typename math_vector<hom_polynomial<T> >::iterator iterator;
    typedef typename math_vector<hom_polynomial<T> >::const_iterator const_iterator; 
  };

/////////////////////// Member functions for curve_param

// Parameterize the line through two points
template <class T>
curve_param <T>::curve_param(const math_vector <T> &p1, const math_vector <T> &p2) : 
                               math_vector < hom_polynomial <T> > (p1.capacity(),p1.capacity())
{ 
  // Polynomial 'x'
  hom_polynomial <T> polx;
  polx.assign_x();
  
  curve_param <T> l1(p1);
  l1.multiply(l1,polx);

  // Polynomial 'y'
  hom_polynomial <T> poly;
  poly.assign_y();
  
  curve_param <T> l2(p2);
  l2.multiply(l2,poly);
  
  this->add(l1, l2);
}

// Pretty printing of param
template <class T>
void curve_param <T>::print_verbose(ostream &s, const char x, const char y) const
{
  s << "[";
  for (rpl_size_t i = 0; i < this->capacity(); i++)
    {
      this->operator[](i).print_verbose(s,x,y);
      if (i != this->capacity()-1)
				s << ", ";
    }
  s << "]";
}

// Evaluate a curve_param at a point 
template <class T>
math_vector <T> curve_param <T>::eval(const T &a, const T &b) const
{
  math_vector <T> tmp(this->size(),this->size());
  hom_polynomial <T> t;
  for (size_t i = 0; i < this->size(); i++)
    {
      t = this->operator[](i);
      tmp[i] = t.eval(a,b);
    }

  return tmp;
}

// Multiply a curve_param by a polynomial
template <class T>
void curve_param <T>::multiply(const curve_param <T> &a, const hom_polynomial <T> &b)
{
  BOOST_AUTO(cp, this->begin());
  BOOST_AUTO(ap, a.begin());

  for (rpl_size_t i = 0; i < a.capacity(); i++, cp++, ap++)
    QI::multiply(*cp,*ap,b); 
}

// Divide a curve_param by a constant
template <class T>
void curve_param <T>::divide(const curve_param <T> &a, const T &b)
{
  // No cast is used because division by a polynomial is Euclidean 
  BOOST_AUTO(cp, this->begin());
  BOOST_AUTO(ap, a.begin());

  for (rpl_size_t i = 0; i < a.capacity(); i++, cp++, ap++)
    QI::divide(*cp,*ap,b);
}

// Multiply a matrix by a curve_param
template <class T>
void curve_param <T>::multiply(const math_matrix<T> &a, const curve_param <T> &b)
{
  BOOST_AUTO(cp, this->begin());
  hom_polynomial <T> tmp;

  // Making a copy will make things easier
  curve_param <T> b_cp = b;

  for (size_t i = 0; i < a.get_no_of_rows(); i++, cp++)
    {
      BOOST_AUTO(bp, b_cp.begin());

      for (size_t j = 0; j < a.get_no_of_columns(); j++, bp++)
				if (j == 0)
	  			QI::multiply(*cp,*bp,a(i,0));
 				else
			 	  {
				    QI::multiply(tmp,*bp,a(i,j));
				    QI::add(*cp,*cp,tmp);
			 	  }
    }
}

// Partial derivative of this curve param
template <class T>
void curve_param <T>::derivative(const curve_param <T> &a, const char v) {

	BOOST_AUTO(cp, this->begin());
	BOOST_AUTO(ap, a.begin());

	for (rpl_size_t i = 0; i < a.capacity(); i++, cp++, ap++ )
		cp->derivative(*ap, v);
}

/////////////////////// Functions for curve_param

// Overloading of cout
template <class T>
inline ostream & operator << (ostream &s, const curve_param <T> &a)
{
  a.print_verbose(s);

  return s;
}

// Negating a curve_param
template <class T>
inline curve_param <T> operator - (const curve_param <T> &a)
{
  curve_param <T> c(a.capacity());

  for (rpl_size_t i = 0; i < a.capacity(); i++)
    c[i].negate(a[i]);

  return c;
}

// Multiply a curve_param by a polynomial
template <class T>
inline void multiply(curve_param <T> &c, const curve_param <T> &a, const hom_polynomial <T> &b)
{
  c.multiply(a,b);
}

// Multiply a curve_param by a constant
template <class T>
inline void multiply(curve_param <T> &c, const curve_param <T> &a, const T &b) 
{
  c.multiply(a,(hom_polynomial <T>)b);
}

// Multiply a vector by a polynomial to give a curve_param
template <class T>
inline void multiply(curve_param <T> &c, const math_vector <T> &a, const hom_polynomial <T> &b)
{
  c.multiply((curve_param <T>)a,b);
}

// Divide a curve_param by a polynomial
template <class T>
inline void divide(curve_param <T> &c, const curve_param <T> &a, const hom_polynomial <T> &b)
{
  c.divide(a,b);
}

// Divide a curve_param by a constant
template <class T>
inline void divide(curve_param <T> &c, const curve_param <T> &a, const T &b)
{
  c.divide(a,b);
}

// Multiply a vector by a constant to give a curve_param
template <class T>
inline void multiply(curve_param <T> &c, const math_vector <T> &a, const T &b)
{
  c.multiply((curve_param <T>)a,(hom_polynomial <T>)b);
}

// Multiply a matrix by a curve_param
template <class T>
inline void multiply(curve_param <T> &c, const math_matrix<T> &a, const curve_param <T> &b)
{
  c.multiply(a,b);
}

// Multiply two curve params
template <class T>
inline void multiply(hom_polynomial <T> &c, const curve_param <T> &a, const curve_param <T> &b)
{
  hom_polynomial <T> tmp; 

  multiply(c,a[0],b[0]);

  for (rpl_size_t i = 1; i < a.capacity(); i++)
    {
      multiply(tmp,a[i],b[i]);
      add(c,c,tmp);
    }
}

// Add two curve params
template <class T>
inline void add(curve_param <T> &c, const curve_param <T> &a, const curve_param <T> &b)
{
  for (rpl_size_t i = 0; i < a.capacity(); i++)
    QI::add(c[i],a[i],b[i]);
}

// Negate a curve param
template <class T>
inline void negate(curve_param <T> &a, const curve_param <T> &b)
{
	rpl::negate(a,b);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Template class for surface parameterizations, i.e. vectors of homogeneous
// polynomials whose coefficients are homogeneous coefficients
template <class T>
class surface_param : public math_vector < hom_hom_polynomial <T> >
  {
  public:
    // Constructors

    // Default
    surface_param(): math_vector < hom_hom_polynomial <T> > (4,4)
      { }

    // Create a surface_param of size a
    surface_param(const rpl_size_t &a): 
                                    math_vector < hom_hom_polynomial <T> > (a,a)
      { }

    // Copy constructor
    surface_param(const surface_param <T> &a):
                     math_vector < hom_hom_polynomial <T> > (a.capacity(),a.capacity())
      {
				this->assign(a);
      } 

    // Destructor
    ~surface_param()
      { }

    // Copy assignment
    surface_param <T> & operator = (const surface_param <T> &a)
      {
				if (this != &a) // Beware of self-assignment
				  this->assign(a);
			
				return *this;
      }

    // Methods

    // Pretty printing - (x,y) are the main variables, (z,w) are for the 
    // (polynomial) coefficients of the monomials
    void print_verbose(ostream &s, const char x = 'u', const char y = 'v',
		       const char z = 's', const char w = 't') const;

    // Evaluate a surface_param at a point 
    // ie, each coordinate of surface_param which is a hom_hom_polynomial 
    // is evaluated as a homogeneous polynomial at (a,b)
    curve_param <T> eval(const hom_polynomial <T> &a, const hom_polynomial <T> &b) const;

    // Arithmetic operations
    // Multiply a curve_param by a polynomial to give a surface_param
    void multiply(const curve_param <T> &a, const hom_polynomial <T> &b);

    // Multiply a matrix by a surface_param
    void multiply(const math_matrix <T> &a, const surface_param <T> &b);

    // Known functions
    friend void hom_hom_polynomial <T>::print_verbose(ostream &s, 
						      const char x, const char y, 
						      const char z, const char w) const;

    std::vector<std::string> surface_paramToString(const char x = 'u', const char y = 'v',
		       const char z = 's', const char w = 't') const;              
  };

/////////////////////// Member functions for surface_param

// Pretty printing
template <class T>
void surface_param <T>::print_verbose(ostream &s, const char x, const char y, 
				      const char z, const char w) const
{
  s << "[";
  for (rpl_size_t i = 0; i < this->capacity(); i++)
    {
      this->operator[](i).print_verbose(s,x,y,z,w);
      if (i != this->capacity()-1)
				s << ", ";
    }
  s << "]";
}

template <class T>
std::vector<std::string> surface_param <T>::surface_paramToString(const char x, const char y, 
				      const char z, const char w) const
{
  std::vector<std::string> ans;
  for (rpl_size_t i = 0; i < this->capacity(); i++)
    {
      std::ostringstream oss;
      this->operator[](i).print_verbose(oss,x,y,z,w);
      ans.emplace_back(oss.str());
    }
  return ans;
}

// Evaluate a surface_param at a point 
// ie, each coordinate of surface_param which is a hom_hom_polynomial 
// is evaluated as a homogeneous polynomial at (a,b)
template <class T>
curve_param <T> surface_param <T>::eval(const hom_polynomial <T> &a, 
					const hom_polynomial <T> &b) const
{
  curve_param <T> tmp(this->capacity());
  hom_hom_polynomial <T> t;

  for (rpl_size_t i = 0; i < this->capacity(); i++)
    {
      t = this->operator[](i);
      tmp[i] = t.eval(a,b);
    }

  return tmp;
}

// Multiply the coefficients of a curve_param by a polynomial
template <class T>
void surface_param <T>::multiply(const curve_param <T> &a, const hom_polynomial <T> &b)
{
  BOOST_AUTO(cp, this->begin());

  for (rpl_size_t i = 0; i < a.capacity(); i++, cp++)
    {
    	if (a[i].is_zero())
        	(*cp).set_degree(-1);
      	else
      		{
      			(*cp).set_degree(b.degree());
      			for (rpl_size_t j = 0; j <= b.degree(); j++)
  						QI::multiply((*cp)[j],a[i],b[j]);
	  			}
    }
}

// Multiply a matrix by a surface_param
template <class T>
void surface_param <T>::multiply(const math_matrix<T> &a, const surface_param <T> &b)
{
  BOOST_AUTO(cp, this->begin());
  hom_hom_polynomial <T> tmp;

  // Making a copy will make things easier
  surface_param <T> b_cp = b;

  for (size_t i = 0; i < a.get_no_of_rows(); i++, cp++)
    {
      BOOST_AUTO(bp,b_cp.begin());

      for (size_t j = 0; j < a.get_no_of_columns(); j++, bp++)
				if (j == 0)
				  QI::multiply(*cp,*bp,a(i,0));
			 	else
			 	  {
				    QI::multiply(tmp,*bp,a(i,j));
				    QI::add(*cp,*cp,tmp);
			 	  }
    }
}

/////////////////////// Functions for surface_param

// Overloading of cout
template <class T>
inline ostream & operator << (ostream &s, const surface_param <T> &a)
{
  a.print_verbose(s);

  return s;
}

// Multiply a curve_param by a polynomial to give a surface_param
template <class T>
inline void multiply(surface_param <T> &c, const curve_param <T> &a, const hom_polynomial <T> &b)
{
  c.multiply(a,b);
}

// Multiply a matrix by a surface_param 
template <class T> 
inline void multiply(surface_param <T> &c, const math_matrix<T> &a, const surface_param <T> &b) 
{ 
  c.multiply(a,b); 
}

// Multiply a vector by a surface_param
/*template <class T>
inline void multiply(surface_param <T> &c, const math_vector <T> &a, const surface_param <T> &b)
{
  c.multiply((math_matrix <T>)a,b);
}*/

// Exchange (u,v) and (s,t) in a surface param
template <class T>
inline surface_param <T> exchange_uv_st(const surface_param <T> &a)
{
  surface_param <T> tmp(a.capacity());

  for (rpl_size_t i = 0; i < a.capacity(); i++)
    tmp[i] = exchange_uv_st(a[i]);

  return tmp;
}

// Multiply two surface params
template <class T>
inline void multiply(hom_hom_polynomial <T> &c, const surface_param <T> &a,
		     const surface_param <T> &b)
{
  hom_hom_polynomial <T> tmp; 

  multiply(c,a[0],b[0]);
  
  for (size_t i = 1; i < 4; i++)
    {
      QI::multiply(tmp,a[i],b[i]);
      QI::add(c,c,tmp);
    }
}

// Add two surface params
template <class T>
inline void add(surface_param <T> &c, const surface_param <T> &a, const surface_param <T> &b)
{
  for (rpl_size_t i = 0; i < 4; i++)
    QI::add(c[i],a[i],b[i]);
}

} // end of namespace QI

#endif
