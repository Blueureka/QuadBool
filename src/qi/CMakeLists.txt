# setting libraries for GMP and Boost Test Framework
set(LIBS  ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} )
# bigint test
add_executable(qi main.cc parse_args.cc settings.cc shell.cc)
target_link_libraries(qi ${LIBS} libqi rpl)
set(fichiers HELP pick.sed pick.sh quadrics.txt rangetest.sh debilos.sh
source_example.qi s2.qi)
file(COPY ${fichiers} DESTINATION .)
