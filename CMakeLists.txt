cmake_minimum_required(VERSION 3.24) # Required for FetchContent FIND_PACKAGE_ARGS
project(scattnlay VERSION 2.5)

cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()

message("Build type is: ${CMAKE_BUILD_TYPE}")
message("Host OS System: ${CMAKE_HOST_SYSTEM}")
message("Hostname:  ${HOSTNAME}")
message("CMake version:  ${CMAKE_VERSION}")

# Select flags.
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# compiler details
message("  C++ Compiler: ${CMAKE_CXX_COMPILER_ID} "
        "${CMAKE_CXX_COMPILER_VERSION} "
        "${CMAKE_CXX_COMPILER_WRAPPER}")

# installation details
message("  Installation prefix: ${CMAKE_INSTALL_PREFIX}")

# Set options
option(BUILD_PYTHON_EXT "Build Python extensions" OFF)
option(STRICT_BUILD "Fail if Boost or Highway are missing" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)
option(ENABLE_ASAN "Enable Address Sanitizer" OFF)

if(ENABLE_ASAN)
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
endif()

# 1. Boost (Multi-Precision)
find_package(Boost 1.60.0 QUIET COMPONENTS headers)
if(Boost_FOUND)
    set(ENABLE_MP ON)
else()
    if(STRICT_BUILD)
        message(FATAL_ERROR "STRICT_BUILD: Boost not found.")
    endif()
    set(ENABLE_MP OFF)
endif()

# 2. Highway (SIMD)
find_package(hwy QUIET)
if(hwy_FOUND)
    set(WITH_HWY ON)
else()
    if(STRICT_BUILD)
        message(FATAL_ERROR "STRICT_BUILD: Google Highway not found.")
    endif()
    set(WITH_HWY OFF)
endif()

include(FetchContent)

# Enable modern find-or-fetch behavior
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS CACHE STRING "Try find_package before fetching")

if(POLICY CMP0169)
  cmake_policy(SET CMP0169 OLD)
endif()

# --- DEPENDENCY MANAGEMENT (Modern FetchContent with FIND_PACKAGE_ARGS) ---

# Google Test (for C++ tests)
if(NOT EMSCRIPTEN)
    FetchContent_Declare(
        GTest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.15.2
        FIND_PACKAGE_ARGS NAMES GTest
    )
    FetchContent_MakeAvailable(GTest)
endif()

# Nanobind (for Python bindings)
if(NOT EMSCRIPTEN)
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
    find_package(nanobind CONFIG QUIET)
    if(NOT nanobind_FOUND)
        FetchContent_Declare(nanobind GIT_REPOSITORY https://github.com/wjakob/nanobind.git GIT_TAG v2.4.0)
        FetchContent_MakeAvailable(nanobind)
    endif()
endif()

# Google Highway (SIMD)
if(WITH_HWY)
    FetchContent_Declare(
        hwy
        GIT_REPOSITORY https://github.com/google/highway.git
        GIT_TAG        1.2.0
        FIND_PACKAGE_ARGS NAMES hwy
    )
    FetchContent_MakeAvailable(hwy)
endif()

# Boost (Multiprecision - header-only)
if(ENABLE_MP)
    find_package(Boost 1.60.0 REQUIRED)
endif()

# --- INTERFACE LIBRARIES FOR COMPILER FLAGS (Phase 3) ---
# Create interface libraries to avoid global pollution

# Warning flags
add_library(project_warnings INTERFACE)
target_compile_options(project_warnings INTERFACE
    -W -Wall -pedantic -Werror
)

# Optimization flags
add_library(project_optimization INTERFACE)
target_compile_options(project_optimization INTERFACE
    -funroll-loops -fstrict-aliasing
)

# --- TARGET CONFIGURATION ---
add_subdirectory(src)

if(NOT BUILD_PYTHON_EXT OR BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

#
# Copy all python scripts to the build directory.
# TODO: Move this to scattnlay subdirectory for better organization (Phase 4)
#
set(Python_SCRIPTS scattnlay/__init__.py scattnlay/main.py)

foreach (_script ${Python_SCRIPTS})
    configure_file(
            ${PROJECT_SOURCE_DIR}/${_script}
            ${PROJECT_BINARY_DIR}/${_script}
            COPYONLY
    )
endforeach ()

# Install Python scripts
install(FILES scattnlay/__init__.py scattnlay/main.py DESTINATION scattnlay)

if(NOT EMSCRIPTEN)
    enable_testing()
    add_subdirectory(tests)
    
    if (NOT DEFINED ENV{GITHUB_ENV})
        add_test(NAME tox COMMAND tox)
    endif()
endif()

# --- PHASE 5: PROFESSIONAL DISTRIBUTION ---
# Export targets for other CMake projects to use scattnlay

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Install export set
install(EXPORT scattnlayTargets
    FILE scattnlayTargets.cmake
    NAMESPACE scattnlay::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scattnlay
)

# Create version file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/scattnlayConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Create config file
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/scattnlayConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/scattnlayConfig.cmake"
    @ONLY
)

# Install config and version files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/scattnlayConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/scattnlayConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scattnlay
)
