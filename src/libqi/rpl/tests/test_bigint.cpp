#define BOOST_TEST_MODULE bigint test

#include <boost/test/unit_test.hpp>
#include <libqi/rpl/bigint.h>
#include <iostream>

using namespace rpl;

BOOST_AUTO_TEST_CASE (constructors_and_assign) {
     
     bigint a;                         // default constructor //
     BOOST_CHECK( a.is_zero() );       // method is_zero //

     bigint b;
     b.assign_one();                     // method assign_one
     BOOST_CHECK( b.is_one() );          // method is_one 

     a.assign_zero();                  // method assign_zero//
     BOOST_CHECK( a.is_zero() );       // method is_zero //
     
     bigint c = "1" ;           
     BOOST_CHECK( c.is_one() );
     a = c = "0" ;            
                                         // operator= (const bigint &) // 
     BOOST_CHECK( a.is_zero() );
     BOOST_CHECK( c.is_zero() );

     
     bigint d;                                                 
     double trois=3.7;
     d.assign(trois);                     // assign(double)  
     c = trois;                           // operator=(const double &)
     BOOST_CHECK_EQUAL(d , bigint(3)); 
     BOOST_CHECK_EQUAL(c , bigint(3)); 
     
     unsigned long three = 3;
     d.assign(three);                     // assign(unsigned long)  
     c = three;                           // operator=(const unsigned long &)
     BOOST_CHECK_EQUAL(d , bigint(3)); 
     BOOST_CHECK_EQUAL(c , bigint(3)); 

     int hiru = -12;                        
     d.assign(hiru);                     // assign(int)  
     c =  hiru;                          // operator=(const int &)
     BOOST_CHECK_EQUAL(d , bigint(-12)); 
     BOOST_CHECK_EQUAL(c , bigint("-12")); 

     long  hamaika= 11;                        
     d.assign(hamaika);                     // assign(long)  
     c =  hamaika;                          // operator=(const long &)
     BOOST_CHECK_EQUAL(d , bigint(11)); 
     BOOST_CHECK_EQUAL(c , bigint("11")); 

}

BOOST_AUTO_TEST_CASE (object_modifiers) {

     // negate -12 -> 12
     bigint a("-12");
     a.negate();
     BOOST_CHECK_EQUAL ( a, bigint(12) );

     // negate 5 -> -5 
     bigint b(5);
     b.negate();
     BOOST_CHECK_EQUAL ( b, bigint(-5) );
     
     
     // check that |-3| = 3 
     bigint mtrois = "-3";
     mtrois.abs();
     BOOST_CHECK_EQUAL ( mtrois, bigint(3) );

     // check that |1| = 1 
     bigint un = 1; 
     un.abs(); 
     BOOST_CHECK_EQUAL( un , bigint(1) ); 


     // check swap method
     un.swap(mtrois);
     BOOST_CHECK_EQUAL( un , bigint(3) ); 
     BOOST_CHECK_EQUAL( mtrois , bigint(1) ); 
   
     // check swap function
     swap(a, b);
     BOOST_CHECK_EQUAL( b , bigint(12) ); 
     BOOST_CHECK_EQUAL( a , bigint(-5) ); 
     
     // check randomize function
     a.randomize(bigint(10)); 
     BOOST_CHECK_GE( a , bigint (0) ); 
     BOOST_CHECK_LE( a , bigint (10-1) ); 
     b.randomize(bigint(-5)); 
     BOOST_CHECK_GE( b , bigint (-5+1) ); 
     BOOST_CHECK_LE( b , bigint (0) ); 


} 

BOOST_AUTO_TEST_CASE (related_functions) {

     // remainder test a >0 
     bigint a = 15;
     bigint b = 4;
     bigint r;
     remainder (r, a, b);
     BOOST_CHECK_EQUAL( r, bigint(3));
     
     // remainder test a <0 
     a = -14;
     unsigned long uib = 4;
     remainder (r, a, uib);
     BOOST_CHECK_EQUAL( r, bigint(-2));

     // test add(bigint, bigint, uint)
     unsigned long b_uint = 128;
     bigint c;
     add(c, a, b_uint);
     BOOST_CHECK_EQUAL( c, bigint(114));

     // test add(bigint, bigint, bigint)
     b = 15;
     add(c, a, b);
     BOOST_CHECK_EQUAL( c, bigint(1));

     // test add(bigint, bigint, int)
     int b_int = -12;
     add(c, a, b_int);
     BOOST_CHECK_EQUAL( c, bigint(-26));
     
     // test add(bigint, bigint, long)
     long b_long = -1000000;
     add(c, a, b_long);
     BOOST_CHECK_EQUAL( c, -bigint(1000014));
    
     // test operator++
     a++;
     BOOST_CHECK_EQUAL(a,bigint(-13));
     
     
     // check that two prime number have a gcd of one
     bigint trois = 3;
     bigint sept = 7;
     bigint pgcd = gcd(trois, sept);
     BOOST_CHECK ( pgcd.is_one());

     // check that gcd(234,-123) = 3 , given by sage
     bigint un_nombre = "234";
     bigint un_autre = -123;
     BOOST_CHECK_EQUAL( gcd(un_nombre, un_autre) , bigint(3)); 

     // test is_square
     //a="-2";
     bigint bidon;
     BOOST_CHECK( !a.is_square());
     //BOOST_CHECK( !is_square(bidon, a));

     a ="4";
     BOOST_CHECK( a.is_square());
     BOOST_CHECK( is_square(bidon, a));


     // check that square computes floor(sqrt(.))
     // sqrt(14) ~ 3.71
     // thus we must find 3!
     a="14";
     bigint sq_a;
     floor_sqrt(sq_a, a);
     BOOST_CHECK_EQUAL( sq_a, bigint("3"));

} 






// EOF
