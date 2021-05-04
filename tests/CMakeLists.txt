cmake_minimum_required(VERSION 3.15)
project(scattnlay_tests C CXX)
# -- Dependency (Google Test)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GTEST_LIBRARIES})
set(LIBS ${LIBS} pthread)

add_compile_options(-D_GLIBCXX_DEBUG)
# -- Output tests in directory

add_executable("test_near_field"
        test_near_field.cc)
target_link_libraries("test_near_field" ${LIBS})
add_test(NAME "test_near_field"
        COMMAND "test_near_field")


# In included file test_spec_functions_data.hpp there are results of multiple
# precision computation that may overflow double precision at compile time.
set_source_files_properties(test_Riccati_Bessel_logarithmic_derivative.cc
        PROPERTIES COMPILE_FLAGS -Wno-overflow)
add_executable("test_Riccati_Bessel_logarithmic_derivative"
        test_Riccati_Bessel_logarithmic_derivative.cc)
target_link_libraries("test_Riccati_Bessel_logarithmic_derivative" ${LIBS})
add_test(NAME "test_Riccati_Bessel_logarithmic_derivative"
        COMMAND "test_Riccati_Bessel_logarithmic_derivative")


add_executable("test_bulk_sphere" test_bulk_sphere.cc)
target_link_libraries("test_bulk_sphere" ${LIBS})
add_test(NAME "test_bulk_sphere"
        COMMAND "test_bulk_sphere")


add_executable("test_bulk_sphere_multi_precision" test_bulk_sphere.cc)
target_compile_options("test_bulk_sphere_multi_precision"
        PRIVATE -DMULTI_PRECISION=100)
target_link_libraries("test_bulk_sphere_multi_precision" ${LIBS})
add_test(NAME "test_bulk_sphere_multi_precision"
        COMMAND "test_bulk_sphere_multi_precision")
