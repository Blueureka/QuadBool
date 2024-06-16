add_subdirectory(libqi)
# add this line to add library rpl
add_subdirectory(libqi/rpl/lib)
# add_subdirectory(bvh)
# MESSAGE(STATUS "ASDASDASDAS ${CMAKE_SOURCE_DIR}="${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories(${CMAKE_SOURCE_DIR}/3rdparty)
include_directories(${PS_BASE_DIR}/gm/inc)
include_directories(${PS_BASE_DIR}/gm/src)
include_directories(${PS_BASE_DIR}/3rd_party)

set(SOURCE_FILES
    CSG/csg.cpp
    libqi/qi.cpp
    # lange/AABB.cpp
    lange/Hybrid.cpp
    lange/TQSI.cpp
    lange/Visualization.cpp
    lange/boolean.cpp
    TEST/testcase.cpp
)

# 创建库，这里假设你的库名为 mylib
add_library(hybridlib ${SOURCE_FILES})

# 如果有其他的依赖库，可以在这里添加
target_link_libraries(hybridlib PUBLIC asan PSGMKernel igl::glfw igl::imgui  libqi
 ginac
 cln)