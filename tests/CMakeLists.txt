cmake_minimum_required(VERSION 3.15)
project(scattnlay_tests C CXX)
# -- Dependency (Google Test)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GTEST_LIBRARIES})
set(LIBS ${LIBS} pthread)

# -- Output tests in directory
add_executable("test_bulk_sphere" test_bulk_sphere.cc)
target_link_libraries("test_bulk_sphere" ${LIBS})
add_test(NAME "test_bulk_sphere"
        COMMAND "test_bulk_sphere")

add_executable("test_bulk_sphere_multi_precision" test_bulk_sphere_multi_precision.cc)
target_compile_options("test_bulk_sphere_multi_precision"
        PRIVATE -DMULTI_PRECISION=100)
target_link_libraries("test_bulk_sphere_multi_precision" ${LIBS})
add_test(NAME "test_bulk_sphere_multi_precision"
        COMMAND "test_bulk_sphere_multi_precision")


add_executable("test_Riccati_Bessel_logarithmic_derivative"
        test_Riccati_Bessel_logarithmic_derivative.cc)
target_link_libraries("test_Riccati_Bessel_logarithmic_derivative" ${LIBS})
add_test(NAME "test_Riccati_Bessel_logarithmic_derivative"
        COMMAND "test_Riccati_Bessel_logarithmic_derivative")

