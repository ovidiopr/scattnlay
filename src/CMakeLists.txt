cmake_minimum_required(VERSION 3.10)
project(scattnlay-main)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set (${PROJECT_NAME}_VERSION "2.3")
set(${HEADERS} nmie.hpp nmie-impl.cc nmie-precision.hpp)

add_executable(scattnlay-dp farfield.cc nmie.cc ${HEADERS})
add_executable(fieldnlay-dp ${HEADERS} nearfield.cc nmie.cc)
