// Template classes for homogeneous polynomials

#ifndef _qi_hompoly_h_
#define _qi_hompoly_h_

#include <libqi/rpl/rpl.h>
#include <libqi/rpl/bigint.h>
#include <libqi/rpl/polynomial.h>
#include <boost/typeof/typeof.hpp>

using namespace std;
using namespace rpl;

// Enter namespace QI
namespace QI {

// Template class for homogeneous polynomials with bigint coefficients
template <class T>
  class hom_polynomial : public polynomial <T>
  {
  public:
    // Constructors
    hom_polynomial(): polynomial <T> ()
      { }
    
    hom_polynomial(const polynomial <T> &p): polynomial <T> (p)
      { }

    hom_polynomial(T x): polynomial <T> (x)
      { }

    // Copy constructor
    hom_polynomial(const hom_polynomial <T> &a)
      { 
				this->assign(a);
      }

    // Destructor
    ~hom_polynomial()
      { } 

    // Copy assignment
    hom_polynomial <T> & operator = (const hom_polynomial <T> &a)
      {
				if (this != &a) // Beware of self-assignment
					this->assign(a);

				return *this;
      }

    hom_polynomial <T> & operator = (const int a)
      {
				this->assign((bigint)a);

				return *this;
      }

    // Member functions

    // Pretty printing
    void print_verbose(ostream &s, const char x = 'x', const char y = 'y') const;

    // Derivatives of a with respect to x and y
    void derivative(const hom_polynomial <T> &a, const char v);

    // Divide a by b
    void divide(const hom_polynomial <T> &a, const hom_polynomial <T> &b);

    // Evaluate polynomial at (a,b)
    T eval(const T &a, const T &b) const;

    // Set polynomial to y
    void assign_y();

    // Set polynomial to x*y
    void assign_xy();

    // Set polynomial to y^2
    void assign_y2();
  };

///////////////////// Member functions for hom_polynomial

// Set polynomial to y
template <class T>
void hom_polynomial <T>::assign_y()
{
	this->assign_zero();
	this->coeff[0] = 1;
	this->deg = 1;
}

// Set polynomial to x*y
template <class T>
void hom_polynomial <T>::assign_xy()
{
	this->assign_zero();
  this->coeff[0] = 0;
  this->coeff[1] = 1;
  this->coeff[2] = 0;
	this->deg = 2;
}

// Set polynomial to y^2
template <class T>
void hom_polynomial <T>::assign_y2()
{
	this->assign_zero();
  this->coeff[0] = 1;
  this->coeff[1] = 0;
  this->coeff[2] = 0;
	this->deg = 2;
}

// Evaluate polynomial at (a,b)
template <class T>
T hom_polynomial <T>::eval(const T &a, const T &b) const
{
  rpl_size_t d = this->deg;
  T c;
  
  if (d == -1)
    {
      c = 0;
      return c;
    }
  
  if (a == 0)
    {
      BOOST_AUTO(cp, this->coeff.begin());
      rpl::power(c,b,d);
      rpl::multiply(c,c,*cp);
    }
  else if (b == 0)
    {
      BOOST_AUTO(cp, this->coeff.begin() + d); 
      rpl::power(c, a, d); 
      rpl::multiply(c, c, *cp); 
    }
  else
    {
      c = 0;
      
      BOOST_AUTO(cp,this->coeff.begin());
      T tmp,tmp2;
      for (rpl_size_t i = 0; i <= d; i++, cp++) 
				{
				  if (*cp != 0)
				    {
				      rpl::power(tmp, a, i);
				      rpl::power(tmp2, b, d - i);
				      rpl::multiply(tmp, tmp, tmp2);
				      rpl::multiply(tmp, tmp, *cp);
				      rpl::add(c, c, tmp);
				    }
				}
    }
    
  return c;
}

// Pretty printing of homogeneous polynomials
// When y = '1', print as a non-homogeneous polynomial
template <class T>
void hom_polynomial <T>::print_verbose(ostream &s, const char x, const char y) const
{
  rpl_size_t d = this->deg;
  bool flag = 0;

  for (rpl_size_t i = d; i >= 0; i--)
    if (!this->coeff[i].is_zero())
      {
				if (flag)
				  s << " ";
				if (this->coeff[i] < 0)
				  s << "- ";
				else if (flag)
				  s << "+ ";
				if (!(abs(this->coeff[i])).is_one())
				  {
				    s << abs(this->coeff[i]);
				    if ((d != 0) && ((y != '1') || (i != 0)))
				      s << "*";
				  }
				else if ((d == 0) || ((i == 0) && (y == '1')))
				  s << abs(this->coeff[i]);
				if (i == 1)
				  s << x;
				else if (i != 0)
				  s << x << "^" << i;
				if (y != '1')
				  {
				    if ((d != 0) && (i != 0) && (i != d))
				      s << "*";
				    if (i == d-1)
				      s << y;
				    else if (i != d)
				      s << y << "^" << d-i;
				  }
				
				flag = 1;
      }

  if (!flag)
    s << 0;
}

// Derivatives with respect to first (x) or second (y) variable
template <class T>
void hom_polynomial <T>::derivative(const hom_polynomial <T> &a, const char v)
{
  rpl_size_t d = a.deg;

  if (d <= 0)
    {
      this->set_degree(-1);
      return;
    }
  
  this->set_degree(d-1);
  T temp;
  if (v == 'x')
    {
      BOOST_AUTO(cp, this->coeff.begin());
      BOOST_AUTO(ap, a.coeff.begin() + 1);
      for (rpl_size_t i = 1; i <= d; i++, cp++, ap++) 
				{
				  temp = i;
				  rpl::multiply(*cp, *ap, temp);
				}
    }
  else
    {
      BOOST_AUTO(cp, this->coeff.begin() + d-1);
      BOOST_AUTO(ap, a.coeff.begin() + d-1);
      for (rpl_size_t i = 1; i <= d; i++, cp--, ap--) 
				{
				  temp = i;
				  rpl::multiply(*cp, *ap, temp);
				}
    }
}

// Divide a by b
template <class T>
void hom_polynomial <T>::divide(const hom_polynomial <T> &a, const hom_polynomial <T> &b)
{
  rpl_size_t ya = ydegree(a), yb = ydegree(b);

  polynomial <T> c_tmp;

  polynomial <T> aa = a;
  polynomial <T> bb = b;

  aa.remove_leading_zeros();
  bb.remove_leading_zeros();
  
  rpl::divide(c_tmp, aa, bb);

  hom_polynomial <T> c_tmp2 = (hom_polynomial <T>)c_tmp;

  // The result of divide may be a constant coefficient times the true result
  // This is not a bug of rpl, behavior is as expected!
  if (!c_tmp2.is_zero())
    {
      bigint tmp = c_tmp2.lead_coeff();
      c_tmp2.multiply(c_tmp2,aa.lead_coeff());
      rpl::divide(c_tmp2,c_tmp2,tmp);       
      rpl::divide(c_tmp2,c_tmp2,bb.lead_coeff());
    }

  int degree_c = c_tmp2.degree()+ya-yb;

  if (degree_c < -1)
    this->set_degree(-1);
  else
    {
      this->set_degree(degree_c);

      T *cp = this->coeff.data();
      for (rpl_size_t i = 0; i <= this->degree(); i++, cp++)
				if (i <= c_tmp2.degree())
				  *cp = c_tmp2[i];
				else
				  *cp = 0;
		}
}

///////////////////// Functions for hom_polynomial

// Overloading of cout
template <class T>
ostream & operator << (ostream &s, const hom_polynomial <T> &a)
{
  a.print_verbose(s);
  return s;
}

// Derivative
template <class T>
hom_polynomial <T> derivative(const hom_polynomial <T> &a, const char v = 'x')
{
  hom_polynomial <T> c;
  c.derivative(a, v);
  return c;
}

// Smallest y-degree of polynomial
template <class T>
rpl_size_t ydegree(const hom_polynomial <T> &a)
{
  rpl_size_t d = a.degree();
  rpl_size_t ydeg = -1;
  
  for (rpl_size_t i = d; i >= 0; i--)
    if (a[i] != 0)
      {
				ydeg = d-i;
				break;
      }

  return ydeg;
}

// Smallest x-degree of polynomial
template <class T>
rpl_size_t xdegree(const hom_polynomial <T> &a)
{
  rpl_size_t d = a.degree();
  rpl_size_t xdeg = -1;
  
  for (rpl_size_t i = 0; i <= d; i++)
    if (a[i] != 0)
      {
				xdeg = i;
				break;
      }

  return xdeg;
}

///////////// Modified operations from the polynomial class

// Multiply two homogeneous polynomials
template <class T>
void multiply(hom_polynomial <T> &c, const hom_polynomial <T> &a,
		     			const hom_polynomial <T> &b)
{
  if ((a.is_zero()) || (b.is_zero()))
    c.set_degree(-1);
  else
    {
      rpl_size_t ya = ydegree(a), yb = ydegree(b);

      rpl::multiply(c,a,b);

      c.set_degree(c.degree()+ya+yb);
    }
}

// Multiply a homogeneous polynomial by a bigint
template <class T>
void multiply(hom_polynomial <T> &c, const hom_polynomial <T> &a,
		     			const T &b)
{
	rpl::multiply(c,a,b);
}

// Negate a hom_polynomial 
template <class T>
void negate(hom_polynomial <T> &b, const hom_polynomial <T> &a)
{
	rpl::negate(b,a);
}

// Add two homogeneous polynomials
template <class T>
void add(hom_polynomial <T> &c, const hom_polynomial <T> &a, 
					const hom_polynomial <T> &b)
{
  rpl::add(c,a,b);
  c.set_degree(max(a.degree(),b.degree()));
}

// Compute the gcd (in x,y) of two homogeneous polynomials
template <class T>
hom_polynomial <T> gcd(const hom_polynomial <T> &a, const hom_polynomial <T> &b)
{
  if (a.is_zero())
    return b;
  else if (b.is_zero())
    return a;
  else
    { 
      rpl_size_t ya = ydegree(a), yb = ydegree(b);

      polynomial <bigint> aa = a;
      polynomial <bigint> bb = b;

      aa.remove_leading_zeros();
      bb.remove_leading_zeros();

      hom_polynomial <T> c = rpl::gcd(aa, bb);

      c.set_degree(c.degree()+min(ya,yb));

      return c;
    }
}

// Power of a hom_polynomial
template <class T>
void power(hom_polynomial <T> &b, const hom_polynomial <T> &a, 
		  			const rpl_size_t &i)
{
  rpl::power(b, a, i);
  
  b.set_degree(i*a.degree());
}

// Compute the polynomial c such that c = a/b
template <class T>
void divide(hom_polynomial <T> &c, const hom_polynomial <T> &a, 
		   			const hom_polynomial <T> &b)
{
  c.divide(a,b);
}

// Divide a polynomial by a constant
template <class T>
void divide(hom_polynomial <T> &c, const hom_polynomial <T> &a, 
		   			const T &b)
{
  rpl::divide(c,a,b);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Template class for homogeneous polynomials with hom_polynomial coefficients
template <class T>
  class hom_hom_polynomial : public hom_polynomial < hom_polynomial <T> >
  {
  public:
    // Constructors
    hom_hom_polynomial(): hom_polynomial < hom_polynomial <T> > ()
      { }

    // Initialize a hom_hom_polynomial with the coeffs of a polynomial
    hom_hom_polynomial(const hom_polynomial <T> &p): 
                          hom_polynomial < hom_polynomial <T> > ()
      { 
				rpl_size_t d = p.degree();

				this->set_degree(d);
				BOOST_AUTO(cp, this->coeff.begin());

				for (rpl_size_t i = 0; i <= d; i++, cp++)
	  			(*cp).assign((hom_polynomial <T>)p[i]);
      }

    // Copy constructor
    hom_hom_polynomial(const hom_hom_polynomial <T> &a) : hom_polynomial < hom_polynomial <T> > (a)
      { }
    
    // Destructor
    ~hom_hom_polynomial()
      { }

    // Copy assignment
    hom_hom_polynomial <T> & operator = (const hom_hom_polynomial <T> &a)
      {
				if (this != &a) // Beware of self-assignment
					this->assign(a);

				return *this;
      }

    // Member functions

    // Pretty printing - (x,y) are the main variables, (z,w) are for the 
    // (polynomial) coefficients of the monomials
    void print_verbose(ostream &s, const char x = 'u', const char y = 'v',
		       const char z = 's', const char w = 't') const;

    // Evaluate hom_hom_polynomial at (a,b)
    hom_polynomial <T> eval(const hom_polynomial <T> &a, const hom_polynomial<T> &b) const;
  };

/////////////////////// Member functions for hom_hom_polynomial

// Evaluate hom_hom_polynomial at (a,b)
template <class T>
hom_polynomial <T> hom_hom_polynomial <T>::eval(const hom_polynomial <T> &a, 
						const hom_polynomial<T> &b) const
{
  rpl_size_t d = this->deg;
  hom_polynomial <T> c;
  
  if (d == -1)
    {
      c = 0;
      return c;
    }
  
  if (a.is_zero())
    {
      BOOST_AUTO(cp, this->coeff.begin());
      QI::power(c,b,d);
      QI::multiply(c,c,*cp);
    }
  else if (b.is_zero())
    {
      BOOST_AUTO(cp, this->coeff.begin() + d); 
      QI::power(c,a,d); 
      QI::multiply(c,c,*cp); 
    }
  else
    {
      c = 0;
      
      BOOST_AUTO(cp, this->coeff.begin());
      hom_polynomial <T> tmp,tmp2;
      for (rpl_size_t i = 0; i <= d; i++, cp++) 
				{
				  if (!(*cp).is_zero())
				    {
				      QI::power(tmp,a,i);
				      QI::power(tmp2,b,d-i);
				      QI::multiply(tmp,tmp,tmp2);
				      QI::multiply(tmp,tmp,*cp);
				      QI::add(c,c,tmp);
				    }
				}
    }

  return c;
}

// Pretty printing - (x,y) are the main variables, (z,w) are for the 
// (polynomial) coefficients of the monomials
template <class T>
void hom_hom_polynomial <T>::print_verbose(ostream &s, const char x, const char y, 
					   const char z, const char w) const
{
  rpl_size_t d = this->deg;
  bool flag = 0;

  for ( rpl_size_t i = d; i >= 0; i--)
    {
      if (!this->coeff[i].is_zero())
	{
	  if (flag)
	    s << " ";

	  // Special treatment when coefficient is a single monomial
	  if (xdegree(this->coeff[i])+ydegree(this->coeff[i]) == this->coeff[i].degree())
	    if (this->coeff[i].degree() == 0)
	      {
		if (this->coeff[i][0] < 0)
		  s << "- ";
		else if (flag)
		  s << "+ ";
		if (!abs(this->coeff[i][0]).is_one())
		  s << abs(this->coeff[i][0]) << "*";
	      }
	    else
	      {
		if ((flag) && (this->coeff[i].lead_coeff() > 0))
		  s << "+ ";
		if ((xdegree(this->coeff[i]) != 0) || (w != '1'))
		  {
		    this->coeff[i].print_verbose(s,z,w);
		    if (d != 0)
		      s << "*";
		  }
	      }
	  else
	    {
	      if (flag)
		s << "+ ";
	      s << "(";
	      this->coeff[i].print_verbose(s,z,w);
	      s << ")*";
	    }
	  
	  if (i != 0)
	    {
	      s << x;
	      if (i != 1)
		s << "^" << i;
	    }
	  if ((i != d) && (i != 0))
	    s << "*";
	  if (i != d)
	    {
	      s << y;
	      if (i != d-1)
		s << "^" << d-i;
	    }
	  
	  flag = 1;
	}
    }

  if (!flag)
    s << 0;
}

///////////////////// Functions for hom_hom_polynomial

// Overloading of cout
template <class T>
inline ostream & operator << (ostream &s, const hom_hom_polynomial <T> &a)
{
  a.print_verbose(s);

  return s;
}

// Multiply a hom_hom_polynomial by a constant
template <class T>
inline void multiply(hom_hom_polynomial <T> &c, const hom_hom_polynomial <T> &a, 
		     const T &b)
{
  if (b == 0)
    c.set_degree(-1);
  else
    {
      rpl_size_t da = a.degree();

      c.set_degree(da);
      for (rpl_size_t i = 0; i <= da; i++)
				rpl::multiply(c[i],a[i],b);
    }
}

// Multiply two hom_hom_polynomials 
template <class T>
inline void multiply(hom_hom_polynomial <T> &c, const hom_hom_polynomial <T> &a,
		     const hom_hom_polynomial <T> &b)
{
  if ((a.is_zero()) || (b.is_zero()))
    c.set_degree(-1);
  else
    {
      hom_polynomial <T> tmp;

      rpl_size_t da = a.degree(), db = b.degree();

      c.set_degree(da+db);

      for (rpl_size_t i = 0; i <= da+db; i++)
				c[i] = 0;

      for (rpl_size_t i = 0; i <= da; i++)
				for (rpl_size_t j = 0; j <= db; j++)
				  {
				    QI::multiply(tmp,a[i],b[j]);
				    QI::add(c[i+j],c[i+j],tmp);
				  }
    }
}

// Add two hom_hom_polynomials -- assumes they both have the same degree (or = 0)
template <class T>
inline void add(hom_hom_polynomial <T> &c, const hom_hom_polynomial <T> &a,
 		const hom_hom_polynomial <T> &b)
{ 
  if (a.is_zero())
    c = b;
  else if (b.is_zero())
    c = a;
  else
    {
      rpl_size_t da = a.degree();
  
      c.set_degree(da);

      for (rpl_size_t j = 0; j <= da; j++)
				QI::add(c[j],a[j],b[j]);
    }
}  

// Exchange (u,v) and (s,t) in a homogeneous bi-variate polynomial 
template <class T>
inline hom_hom_polynomial <T> exchange_uv_st(const hom_hom_polynomial <T> &a)
{
  // For quadratic polynomials in (u,v) (works for any degree)
  // We have : a = u^2 . a[2] + uv . a[1] + v^2 . a[0]
  // If degree(a[i]) = 2 then (for quadratic polynomials in (s,t))
  // a[i] = s^2 . a[i][2] + st . a[i][1] + t^2 . a[i][0]
  // and     a = s^2 . (u^2 . a[2][2] + uv . a[1][2] + v^2 . a[0][2])
  //              + st .  (u^2 . a[2][1] + uv . a[1][1] + v^2 . a[0][1])
  //              + t^2 .(u^2 . a[2][0] + uv . a[1][0] + v^2 . a[0][0])
  // Else if degree(a[i]) = 1 then (for linear polynomials in (s,t))
  // a[i] = s . a[i][1] + t . a[i][0]
  //           a = s . (u^2 . a[2][1] + uv . a[1][1] + v^2 . a[0][1])
  //              + t  .(u^2 . a[2][0] + uv . a[1][0] + v^2 . a[0][0])

  // the polynomial a where the variables (u,v) and (s,t) are exchanged
  hom_hom_polynomial <T> tmp; 
  int deg = -1;

  for (rpl_size_t i = 0; i <= a.degree(); i++)
    if (a[i].degree() > deg)
      deg = a[i].degree();

  tmp.set_degree(deg);

  for (rpl_size_t j = 0; j <= tmp.degree(); j++)
    tmp[j].set_degree(a.degree()); 

  for (rpl_size_t i = 0; i <= a.degree(); i++)
    if (!a[i].is_zero())
      for (rpl_size_t j = 0; j <= tmp.degree(); j++)
	tmp[j][i] = a[i][j];

  return tmp; 
}

} // end of namespace QI

#endif

