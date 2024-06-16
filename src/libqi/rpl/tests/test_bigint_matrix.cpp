#define BOOST_TEST_MODULE bigint_matrix test

#include <boost/test/unit_test.hpp>
#include <libqi/rpl/bigint_matrix.h>
#include <libqi/rpl/math_vector.h>
#include <iostream>


using namespace std;
using namespace rpl;
// to test matrix equality
void checkMatEqual(const bigint_matrix & a, const bigint_matrix & b) {
    for(size_t i= 0; i < a.get_no_of_rows() ; ++i) {
        for(size_t j= 0; j < a.get_no_of_columns() ; ++j) {
         BOOST_CHECK_EQUAL(a(i,j), b(i,j));
      }
    }
}

BOOST_AUTO_TEST_CASE (constructors_and_assign) {

    bigint_matrix a;
    a.resize(3,3);
    for(size_t i= 0; i < a.get_no_of_rows() ; ++i){ 
      for(size_t j = 0; j < a.get_no_of_columns() ; ++j){
          a(i,j) = i; 
      }
    }

    bigint_matrix b;
    b.resize(3,3);
    for(size_t i= 0; i < b.get_no_of_rows() ; ++i){ 
      for(size_t j = 0; j < b.get_no_of_columns() ; ++j){
          b.sto(i,j,j); 
      }
    }

    // test for + operator
    bigint_matrix c;
    c.resize(3,3);
    c = a + b;
    
    // tests for member
    for(size_t i= 0; i < c.get_no_of_rows() ; ++i){ 
      for(size_t j = 0; j < c.get_no_of_columns() ; ++j){
          BOOST_CHECK_EQUAL(c.member(i,j), i+j);   
      }
    }

    bigint rampe[3] = { "1", "2", "3" }; 
    math_vector<bigint> v1(rampe,3);
    math_vector<bigint> v2(3);
   
    // test for prod matrix_vector product
    v2 = prod(a,v1); 

    math_vector<bigint> v3= a.get_row_vector(0);

}

BOOST_AUTO_TEST_CASE (expressions) {

    bigint coeff[] = { 1, 2, 3,
                              4, 5, 6};
    bigint_matrix A(2,3);
    std::copy(coeff, coeff + 6, A.data().begin());

    bigint_matrix B = trans(A);

    bigint coeff_trans[] = { 1,  4,
                             2,  5,
			     3,  6};

    bigint_matrix checkB(3,2);
    std::copy(coeff_trans, coeff_trans + 6, checkB.data().begin());
    checkMatEqual(B, checkB);

    bigint_matrix C = trans(B) -A ;
    //C = A.trans();
    checkMatEqual(C, math_matrix<bigint>::zero(2,3)); 

}

BOOST_AUTO_TEST_CASE (kernel_function) {
    
    // test with identity
    // must return a matrix with zero column
    bigint_matrix a =  math_matrix<bigint>::identity(3);
    BOOST_CHECK_EQUAL( kernel(a).get_no_of_columns(), 0); 
   
    BOOST_CHECK( kernel(a).is_column_zero(0));

    bigint zero_col[] = {0, 0, 0, 0};
    bigint_matrix test_is_col(4,1);
    std::copy(zero_col, zero_col +4, test_is_col.data().begin());
    BOOST_CHECK(test_is_col.is_column_zero(0));


    // test with the constant matrix 1
    a.resize(3,3);
    for (size_t i = 0; i < 3; ++i) 
       for (size_t j = 0; j < 3 ; ++j) 
            a(i,j) = 1;
    bigint result[6] = { 1, 0, 0, 1, -1 , -1}; 
    bigint_matrix ker_a(3,2);
    std::copy(result, result+6 ,ker_a.data().begin());
    BOOST_CHECK(!a.is_column_zero(0));

    checkMatEqual(kernel(a), ker_a);

    // test with the null matrix
    // must return the identity
    a = math_matrix<bigint>::zero(4,4);
    checkMatEqual(kernel(a), math_matrix<bigint>::identity(4));



    // test with the planting kernel --> must find zero
    bigint coefs_p[] = { 0,  0,  0, 
                         0,  1,  0, 
                         1,  0,  0, 
                         0,  0,  1};

    a.resize(4,3);
    std::copy(coefs_p, coefs_p + 12 ,a.data().begin());
    BOOST_CHECK_EQUAL( kernel(a).get_no_of_columns(), 0); 

    bigint coefs_v[] = { 2, 2, 4, 4, 4};
    bigint_matrix v1(5,1);
    std::copy(coefs_v, coefs_v + 5 ,v1.data().begin());
    bigint_matrix proj_sur_v = prod(v1,trans(v1));

    bigint_matrix k_v = kernel(proj_sur_v);
    bigint_matrix check_v= prod(proj_sur_v, k_v);
    BOOST_CHECK( check_v.is_zero());
    
    
    
    
    // test with another bugged matrix
    a.resize(4,3);
    bigint coefs2[]  = { 1, 0, 1, 
                         0, 0, 0, 
                         0, 0, 0, 
                         0, 1, 1 }; 

    std::copy(coefs2, coefs2 + 12 ,a.data().begin());
    BOOST_CHECK_EQUAL( kernel(a).get_no_of_columns(), 1); 
    bigint_matrix check_a = prod(a,kernel(a));
    BOOST_CHECK( check_a.is_zero()); 

}

BOOST_AUTO_TEST_CASE (det_function) {

  bigint_matrix a =  math_matrix<bigint>::identity(3);
  BOOST_CHECK_EQUAL( det(a), bigint(1)); 
  
  // test matrix with columns non-linearly dependents
  bigint_matrix quatre(3,3);
  for(size_t i = 0; i < quatre.size1(); ++i)
    for(size_t j = 0; j < quatre.size2(); ++j)
          quatre(i,j) = i * 3 + j;

  BOOST_CHECK_EQUAL( det(quatre), bigint(0));

  // test Vandermonde 4x4
  bigint coeffs[] = {1, 2, 3, 4};
  bigint_matrix vander(4,4);
  
  for(size_t i = 0; i < vander.size1(); ++i)
    for(size_t j = 0; j < vander.size2(); ++j)
          power(vander(i,j), coeffs[j], i);

  BOOST_CHECK_EQUAL( det(vander), bigint(12));

  bigint coeffs_q1[] ={ 2, -1, 0, 0,
                        -1, -2, 0, -1,
			 0, 0, 2, 0,
			 0, -1, 0, 2};
  bigint_matrix q1(4,4);
  std::copy(coeffs_q1, coeffs_q1 + 16 ,q1.data().begin());
  BOOST_CHECK_EQUAL( det(q1), bigint(-24));

  // exemple qui fait planter det_pencil
  bigint_matrix q2(4,4);
  bigint coeffs_q2[] ={ 0, 0, 0, 0,
                        1, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, -3};
  std::copy(coeffs_q2, coeffs_q2 + 16 ,q2.data().begin());
  BOOST_CHECK_EQUAL( det(q2), bigint(0));
  

}



BOOST_AUTO_TEST_CASE (split_and_compose) {
   
  bigint entiers[] = { 1, 2, 3, 4, 5,
                       6, 7, 8, 9, 10,
		       11, 12, 13, 14, 15,
		       16, 17, 18, 19, 20,
		       21, 22, 23, 24, 25 };

  bigint_matrix A(5, 5);
  std::copy(entiers, entiers + 25, A.data().begin());

  bigint_matrix B(3, 3), checkB(3, 3);
  std::copy(entiers, entiers + 3, checkB.data().begin());
  std::copy(entiers + 5 , entiers + 8, checkB.data().begin() + 3);
  std::copy(entiers + 10 , entiers + 13, checkB.data().begin() + 6);

  bigint_matrix C(3, 2), checkC(3, 2);
  std::copy(entiers + 3, entiers + 5, checkC.data().begin());
  std::copy(entiers + 8 , entiers + 10, checkC.data().begin() + 2);
  std::copy(entiers + 13 , entiers + 15, checkC.data().begin() + 4);

  bigint_matrix D(2, 3), checkD(2, 3);
  std::copy(entiers + 15, entiers + 18, checkD.data().begin());
  std::copy(entiers + 20 , entiers + 23, checkD.data().begin() + 3);
  
  bigint_matrix E(2, 2), checkE(2, 2);
  std::copy(entiers + 18, entiers + 20, checkE.data().begin());
  std::copy(entiers + 23 , entiers + 25, checkE.data().begin() + 2);


  // split_t test
  A.split_t(B, C, D, E);
  checkMatEqual(B, checkB);
  checkMatEqual(C, checkC);
  checkMatEqual(D, checkD);
  checkMatEqual(E, checkE);

  // split_h test 
  bigint_matrix Ash(A);
  bigint_matrix Bsh(3, 3), Csh(3, 2);
  Ash.resize(3, 5);
  Ash.split_h(Bsh, Csh);
  checkMatEqual(Bsh, checkB);
  checkMatEqual(Csh, checkC);

  // compose_t test
  bigint_matrix AP(5, 5);
  AP.compose_t(B, C, D, E);
  checkMatEqual(AP, A);

  // compose_h test
  bigint_matrix Ach(3,5);
  Ach.compose_h(Bsh, Csh);
  checkMatEqual(Ach, Ash);

}

BOOST_AUTO_TEST_CASE( diag_test ) {
  
  bigint entiers[] = { 1, 2, 3, 4, 5,
                       6, 7, 8, 9, 10,
		       11, 12, 13, 14, 15,
		       16, 17, 18, 19, 20,
		       21, 22, 23, 24, 25 };

  bigint_matrix A(5, 5);
  std::copy(entiers, entiers + 25, A.data().begin());

  // tipical case
  A.diag(1,0);
  
  checkMatEqual(A, math_matrix<bigint>::identity(5));

}

// EOF
