cmake_minimum_required(VERSION 3.10)
project(scattnlay-main)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set (${PROJECT_NAME}_VERSION "2.3")

#target_link_libraries(/.local/lib/python3.8/site-packages/pybind11)
add_subdirectory(src)
