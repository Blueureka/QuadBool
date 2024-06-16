file(GLOB SRC . *.cpp)
#message(STATUS "sources= " ${SRC})
add_library(rpl STATIC ${SRC})
#target_link_libraries(rpl ${GMP_LIBRARIES} ${GMPXX_LIBRARIES} ${Boost_LIBRARIES})
