#define BOOST_TEST_MODULE bigfloat test

#include <boost/test/unit_test.hpp>
#include <libqi/rpl/bigfloat.h>
#include <iostream>

using namespace rpl;
BOOST_AUTO_TEST_CASE (constructors_and_assign) {
     
     bigfloat a;                         // default constructor //
     BOOST_CHECK( a.is_zero() );       // method is_zero //

     bigfloat b;
     b.assign_one();                     // method assign_one
     BOOST_CHECK( b.is_one() );          // method is_one 

     a.assign_zero();                  // method assign_zero//
     BOOST_CHECK( a.is_zero() );       // method is_zero //
     
     bigfloat c = "1" ;           
     BOOST_CHECK( c.is_one() );
     a = c = "0" ;            
                                         // operator= (const bigfloat &) // 
     BOOST_CHECK( a.is_zero() );
     BOOST_CHECK( c.is_zero() );

     
     bigfloat d;                                                 
     double trois=3.5;
     d.assign(trois);                     // assign(double)  
     c = trois;                           // operator=(const double &)
     BOOST_CHECK_EQUAL(d , bigfloat("3.5")); 
     BOOST_CHECK_EQUAL(c , bigfloat("3.5")); 
     
     unsigned long three = 3;
     d.assign(three);                     // assign(unsigned long)  
     c = three;                           // operator=(const unsigned long &)
     BOOST_CHECK_EQUAL(d , bigfloat(3)); 
     BOOST_CHECK_EQUAL(c , bigfloat(3)); 

     int hiru = -12;                        
     d.assign(hiru);                     // assign(int)  
     c =  hiru;                          // operator=(const int &)
     BOOST_CHECK_EQUAL(d , bigfloat(-12)); 
     BOOST_CHECK_EQUAL(c , bigfloat("-12")); 

     long  hamaika= 11;                        
     d.assign(hamaika);                     // assign(long)  
     c =  hamaika;                          // operator=(const long &)
     BOOST_CHECK_EQUAL(d , bigfloat(11)); 
     BOOST_CHECK_EQUAL(c , bigfloat("11")); 

}

// EOF
