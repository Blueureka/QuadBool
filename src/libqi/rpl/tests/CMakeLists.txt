# setting libraries for GMP and Boost Test Framework
set(LIBS  ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${Boost_LIBRARIES})

# bigint test
add_executable(test_bigint test_bigint.cpp)
target_link_libraries(test_bigint rpl ${LIBS})

# bigfloat test
add_executable(test_bigfloat test_bigfloat.cpp)
target_link_libraries(test_bigfloat rpl ${LIBS})

# bigint_matrix test
add_executable(test_bigint_matrix test_bigint_matrix.cpp)
target_link_libraries(test_bigint_matrix rpl ${LIBS})

# polynomial test
add_executable(test_poly test_poly.cpp)
target_link_libraries(test_poly rpl ${LIBS})

# rational_factorization test
add_executable(test_rational_factorization test_rational_factorization.cpp)
target_link_libraries(test_rational_factorization rpl ${LIBS})
