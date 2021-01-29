cmake_minimum_required(VERSION 3.10)
project(scattnlay-main)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set (${PROJECT_NAME}_VERSION "2.3")


add_subdirectory(src)
