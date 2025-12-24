cmake_minimum_required(VERSION 3.18) # 3.18+ recommended for FetchContent features
project(scattnlay VERSION 2.4)

cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)

message("Build type is: ${CMAKE_BUILD_TYPE}")
message("Host OS System: ${CMAKE_HOST_SYSTEM}")
message("Hostname:  ${HOSTNAME}")
message("CMake version:  ${CMAKE_VERSION}")

# Select flags.
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-W -Wall -pedantic -Werror)
add_compile_options(-funroll-loops -fstrict-aliasing)


# compiler details
message("  C++ Compiler: ${CMAKE_CXX_COMPILER_ID} "
        "${CMAKE_CXX_COMPILER_VERSION} "
        "${CMAKE_CXX_COMPILER_WRAPPER}")

# installation details
message("  Installation prefix: ${CMAKE_INSTALL_PREFIX}")

# Set options
option(ENABLE_MP "Use multiple precision" OFF)
option(WITH_HWY "Enable Google Highway SIMD" OFF)

# --- 1. SETUP PERSISTENT PATHS ---
set(LOCAL_DEPS_DIR "${CMAKE_SOURCE_DIR}/build_deps" CACHE PATH "Persistent dependencies")
list(APPEND CMAKE_PREFIX_PATH "${LOCAL_DEPS_DIR}")

include(FetchContent)

# --- 2. THE BOOTSTRAP MACRO ---
macro(smart_import name repo tag)
    # Check if exists in system or build_deps
    find_package(${name} QUIET)

    if(NOT ${name}_FOUND)
        message(STATUS "📦 ${name} not found. Bootstrapping to ${LOCAL_DEPS_DIR}...")

        FetchContent_Declare(
            ${name}_fetch
            GIT_REPOSITORY ${repo}
            GIT_TAG        ${tag}
            SOURCE_DIR     "${LOCAL_DEPS_DIR}/sources/${name}"
        )

        FetchContent_GetProperties(${name}_fetch)
        if(NOT ${name}_fetch_POPULATED)
            FetchContent_Populate(${name}_fetch)
            
            # Identify the build type to match main project
            if(NOT CMAKE_BUILD_TYPE)
                set(BUILD_TYPE_PROP "Release")
            else()
                set(BUILD_TYPE_PROP ${CMAKE_BUILD_TYPE})
            endif()

            execute_process(
                COMMAND ${CMAKE_COMMAND} -S "${${name}_fetch_SOURCE_DIR}" -B "${${name}_fetch_BINARY_DIR}"
                        -DCMAKE_INSTALL_PREFIX=${LOCAL_DEPS_DIR}
                        -DCMAKE_BUILD_TYPE=${BUILD_TYPE_PROP}
                        -DBUILD_SHARED_LIBS=OFF
                        -DINSTALL_GTEST=ON
                        -DPYBIND11_FINDPACKAGE=ON
                COMMAND ${CMAKE_COMMAND} --build "${${name}_fetch_BINARY_DIR}" --target install -j 8
                COMMAND_ERROR_IS_FATAL ANY
            )
        endif()
        find_package(${name} REQUIRED)
    endif()
    message(STATUS "✅ Found ${name}: ${${name}_DIR}")
endmacro()

# --- 3. MANAGE DEPENDENCIES ---

# Google Test (for C++ tests)
if(NOT EMSCRIPTEN)
    smart_import(GTest https://github.com/google/googletest.git v1.14.0)
endif()

# Pybind11 (for Python bindings)
# Note: We still need the Python interpreter to find NumPy headers
if(NOT EMSCRIPTEN)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    smart_import(pybind11 https://github.com/pybind/pybind11.git v2.11.1)
endif()

# Google Highway (Planned SIMD)
if(WITH_HWY)
    smart_import(hwy https://github.com/google/highway.git 1.0.7)
endif()

# Boost (Multiprecision)
# Boost is heavy to FetchContent; we prioritize system Boost
if(ENABLE_MP)
    find_package(Boost 1.60.0 REQUIRED)
    if(Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIRS})
    endif()
endif()

# --- 4. TARGET CONFIGURATION ---
add_subdirectory(src)
add_subdirectory(examples)
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

if(NOT EMSCRIPTEN)
    enable_testing()
    add_subdirectory(tests)
    
    if (NOT DEFINED ENV{GITHUB_ENV})
        add_test(NAME tox COMMAND tox)
    endif()
endif()

# Helper target to clear deps using git clean logic
add_custom_target(clean-deps
    COMMAND git clean -fdX "${LOCAL_DEPS_DIR}"
    COMMENT "Wiping persistent dependencies using git clean"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
