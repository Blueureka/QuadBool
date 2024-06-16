# 添加find GMP
find_package(GMP REQUIRED)

# setting libraries for GMP and Boost Test Framework
set(LIBS  ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} )
# list of sources used in the kernel part of the lib
set(KERNEL_SRC kernel/QIElem.cc kernel/QINumber.cc kernel/QIInterWith22.cc
kernel/QIParam.cc kernel/QIInter.cc kernel/QICheck.cc kernel/QIOneMult.cc
kernel/QITwoMult.cc kernel/QINoMult.cc kernel/QISolve.cc kernel/QIBench.cc
kernel/QIVanishDet.cc)
# list of sources used in the io part of the lib
set(IO_SRC io/QIConsoleWriter.cc io/QIParser.cc io/QIWriter.cc io/QIComponentLabels.c 
io/QIOutputter.cc)
# bigint test
add_library(libqi STATIC ${KERNEL_SRC} ${IO_SRC})
target_include_directories(libqi PUBLIC ${PROJECT_SOURCE_DIR}/include)
# target_link_libraries(libqi ${GMP_LIBRARIES} ${GMPXX_LIBRARIES} ${Boost_LIBRARIES} rpl)
# 取消注释
target_link_libraries(libqi ${GMP_LIBRARIES} ${GMPXX_LIBRARIES} ${Boost_LIBRARIES} rpl)
