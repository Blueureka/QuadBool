#define BOOST_TEST_MODULE rational_factorization test

#include <boost/test/unit_test.hpp>
#include <libqi/rpl/bigint.h>
#include <libqi/rpl/rational_factorization.h>
#include <iostream>


using namespace rpl;
void checkFactorisation(const rational_factorization & r, 
                        const bigint * bases, 
			const size_t * expos) {

    for(size_t i = 0; i < r.no_of_comp() ; ++i) {
       BOOST_CHECK_EQUAL(r.base(i),bases[i]);
       BOOST_CHECK_EQUAL(r.exponent(i),expos[i]);
    }
}			


bigint numbers_test [] = {
     "3240",
     "29400",
     "2143",
     "12482",
     "4294967297",
     "18446744073709551617",
};

bigint base_1[] = { 2, 3, 5};
size_t exp_1[] = { 3, 4, 1};

bigint base_2[] = {2, 3, 5, 7};
size_t exp_2[] = {3, 1, 2, 2};

bigint base_3[] = {2143};
size_t exp_3[] = {1};

bigint base_4[] = {2, 79};
size_t exp_4[] = {1, 2};

bigint base_5[] = { 641, 6700417};
size_t exp_5[] = { 1, 1};

bigint base_6[] = { "18446744073709551617" };
size_t exp_6[] = { 1};

bigint * bases[6] = { base_1 , base_2, base_3, base_4, base_5, base_6};
size_t * expos[6] = { exp_1 , exp_2, exp_3, exp_4, exp_5, exp_6};


BOOST_AUTO_TEST_CASE ( basic_tests) {
     
      for(size_t i = 0; i < 6; ++i) {
          rational_factorization f = numbers_test[i];
	  f.factor();
	  checkFactorisation(f, bases[i], expos[i]);
      }

}

