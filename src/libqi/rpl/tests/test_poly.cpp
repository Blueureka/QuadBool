#define BOOST_TEST_MODULE bigint test

#include <boost/test/unit_test.hpp>
#include <libqi/rpl/rpl.h>
#include <libqi/rpl/bigint.h>
#include <libqi/rpl/polynomial.h>
#include <iostream>

using namespace rpl;

typedef polynomial<bigint> bpoly;
typedef polynomial<bpoly> bbpoly;

template <class T>
void checkPolyEqual(const polynomial<T> & a, const polynomial<T> b){
   BOOST_CHECK_EQUAL(a.deg, b.deg); 
   for(rpl_size_t i = 0; i <= a.deg; ++i){
       BOOST_CHECK_EQUAL(a.coeff[i],b.coeff[i]);
   }
}

BOOST_AUTO_TEST_CASE (constructors_and_assign) {

       bpoly a;
       bigint coefs[2] = { "0", "1"};
       bpoly b(coefs,1); // constructor from a pointer
       a.assign_x(); // assign_x
      
       BOOST_CHECK_EQUAL(b, a);

       b[4] = 0;
       b.remove_leading_zeros();
       std::cout << b.deg << std::endl;
       BOOST_CHECK_EQUAL(b.degree(), 1); // test for remove leading zero 

       b.assign_zero();
       BOOST_CHECK_EQUAL(b.degree(), -1);  // test de assign_zero
      
//       b.assign(a);
//       BOOST_CHECK( b.is_x());

       a.assign_one();
       bigint un(1);
       
       BOOST_CHECK_EQUAL( a, bpoly(&un, 0));

       bbpoly bb;
       bb[0].assign_x();

}


BOOST_AUTO_TEST_CASE (operations) {

       // test for negate and add
       bigint coefs[5] = { 1, 2, 3, 4, 5};
       bpoly p1(coefs, 4);
       bpoly p2;
       p2.negate(p1);
       bpoly p3 = p1 + p2 ;
       BOOST_CHECK (p3.is_zero());

       // test for derivative
       bigint coefs_deriv[4] = { 2, 6, 12, 20};
       bpoly p_prime(coefs_deriv,3);
       bpoly p1p;
//       derivative(p1p,p1);
//       BOOST_CHECK_EQUAL(p_prime, p1p);

       bpoly hix;
       hix.assign_x();
       bpoly hix_quatre;
       hix_quatre.power(hix,4);
       hix_quatre[0] = bigint(-4);
       BOOST_CHECK_EQUAL(hix_quatre.const_term(), -4);
       BOOST_CHECK_EQUAL(hix_quatre.lead_coeff(), 1);
      
       
       // test multiply
       bigint coefs_4[] = {2, 0, 1};
       bigint coefs_5[] = {1, 1};
       bpoly p4(coefs_4, 2);
       bpoly p5(coefs_5, 1);
       
       bigint coefs_45[] = { 2, 2, 1, 1};
       bpoly p45(coefs_45, 3);
       bpoly tmp; 
       tmp.multiply(p4, p5);

       checkPolyEqual(p45, tmp);
    
       // test self_multiply
       p4.multiply(p4, p5);
       checkPolyEqual(p4, tmp);

       //b test div
       bpoly q;
       bpoly r;
       div_rem(q, r, p45, p5);
       tmp = bpoly(coefs_4, 2);
       bpoly zero;
       checkPolyEqual(r, zero);
       checkPolyEqual(q, tmp);

       bigint coefs_6[] = {1, 0, 2};
       bpoly p6(coefs_6,2);
       div_rem(q,r, p1, p6);
       bpoly tmp6;
       tmp6.multiply(p6,q);
       tmp6.add(tmp6,r);
       p1.multiply(p1,bigint(8));
       // we check that d^{m-n+1} A  = BQ +R
       checkPolyEqual(p1,tmp6);
       
       // test division by a polynomial of degree 0 (a non const)
       bigint coefs_7[]  = {2, 0, 2};
       bigint coefs_8[] = {2};
       bpoly p7(coefs_7, 2);
       bpoly p8(coefs_8, 0);
       div_rem(q, r, p7, p8);

}
