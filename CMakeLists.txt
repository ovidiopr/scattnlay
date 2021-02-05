cmake_minimum_required(VERSION 3.15)
project(scattnlay VERSION 2.3)

cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)

message("Build type is: ${CMAKE_BUILD_TYPE}")
message("Host OS System: ${CMAKE_HOST_SYSTEM}")
message("Hostname:  ${HOSTNAME}")
message("CMake version:  ${CMAKE_VERSION}")

# Select flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall)
add_compile_options(-funroll-loops -fstrict-aliasing)

# Set options
option(ENABLE_MP "Use multiple precision" OFF)
if (${ENABLE_MP})
    add_compile_options(-DMULTI_PRECISION=100)
endif ()

# compiler details
message("  C++ Compiler: ${CMAKE_CXX_COMPILER_ID} "
        "${CMAKE_CXX_COMPILER_VERSION} "
        "${CMAKE_CXX_COMPILER_WRAPPER}")

# installation details
message("  Installation prefix: ${CMAKE_INSTALL_PREFIX}")

# Find Boost
set(BOOSTROOT $ENV{BOOST_DIR})
if (USE_STATIC_LIBRARIES)
    set(Boost_USE_STATIC_LIBS ON)
endif ()
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED)
if (Boost_INCLUDE_DIRS)
    if (${Boost_VERSION} VERSION_LESS 1.60.0)
        message(FATAL_ERROR
                "Found Boost library is too old; required is version 1.60.0 or newer!")
    endif ()
    message("Found Boost include dir: ${Boost_INCLUDE_DIR}")
    message("Found Boost library dir: ${Boost_LIBRARY_DIR}")
    message("Found Boost libraries: ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIRS})
endif ()

#Find Python, NumPy and PyBind11
find_package(Python COMPONENTS Interpreter Development)

include_directories(${Python_INCLUDE_DIRS})

message("Python_EXECUTABLE: ${Python_EXECUTABLE}")
message("Python_FOUND: ${Python_FOUND}")
message("Python_VERSION: ${Python_VERSION}")
message("Python_Development_FOUND: ${Python_Development_FOUND}")
message("Python_LIBRARIES: ${Python_LIBRARIES}")
message("Python_INCLUDE_DIRS: ${Python_INCLUDE_DIRS}")

# Ensure that numpy is installed and read its include dir
exec_program(${Python_EXECUTABLE}
        ARGS "-c \"import numpy; print(numpy.get_include())\""
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        RETURN_VALUE NUMPY_NOT_FOUND
        )
if (NUMPY_NOT_FOUND)
    message(FATAL_ERROR "NumPy headers not found")
endif ()

# Ensure that pybind11 is installed and read its include dir
exec_program(${Python_EXECUTABLE}
        ARGS "-c \"import pybind11; print(pybind11.get_include())\""
        OUTPUT_VARIABLE PYBIND11_INCLUDE_DIR
        RETURN_VALUE PYBIND11_NOT_FOUND
        )
if (PYBIND11_NOT_FOUND)
    message(FATAL_ERROR "PyBind11 headers not found")
endif ()

#include_directories(src)
add_subdirectory(src)

#
# Copy all python scripts to the build directory.
#
set(Python_SCRIPTS scattnlay/__init__.py scattnlay/main.py)

foreach (_script ${Python_SCRIPTS})
    configure_file(
            ${PROJECT_SOURCE_DIR}/${_script}
            ${PROJECT_BINARY_DIR}/${_script}
            COPYONLY
    )
endforeach ()
