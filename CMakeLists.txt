cmake_minimum_required(VERSION 3.10)

project(QuadBool VERSION 1.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# set CMAKE_BUILD_TYPE default to Release
if(NOT CMAKE_BUILD_TYPE)
# set(CMAKE_BUILD_TYPE Debug)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE = " ${CMAKE_BUILD_TYPE})



message(STATUS "PROJECT_NAME =  ${PROJECT_NAME}")
add_executable(${PROJECT_NAME} main.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rdparty)

set (THIRD_PARTY_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty)

set(PS_BASE_DIR ${THIRD_PARTY_BASE_DIR}/ps-kernel)
message(STATUS "PS_BASE_DIR =  ${PS_BASE_DIR}")
# set(PSGM_BASE_DIR PS_BASE_DIR/gm)
add_subdirectory(${PS_BASE_DIR})
include_directories(${PS_BASE_DIR}/gm/inc)
include_directories(${PS_BASE_DIR}/gm/src)
include_directories(${PS_BASE_DIR}/3rd_party/eigen/include)
include_directories(${PS_BASE_DIR}/3rd_party)

include_directories(${THIRD_PARTY_BASE_DIR})
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_BUILD 2)
set(PROJECT_VERSION_PATCH 0)

option(Enable_model_view "enable to view model in testcase" ON)
message("Enable_model_view is ${Enable_model_view}")
add_definitions(-DPS_KERNEL_PATH="${CMAKE_SOURCE_DIR}")
if(${Enable_model_view})
  add_definitions(-DENABLE_PS_CLOUD_VIEW)
endif()

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)


############ for test: libigl ############
# list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake)
# Libigl
include(libigl)
# Enable the target igl::glfw
igl_include(glfw)
igl_include(imgui)
##########################################

add_subdirectory(src)

target_link_libraries(${PROJECT_NAME} PUBLIC asan PSGMKernel igl::glfw igl::imgui  libqi
 ginac
 cln
 hybridlib)