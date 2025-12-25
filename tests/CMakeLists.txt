cmake_minimum_required(VERSION 3.20)
project(scattnlay_tests C CXX)

# -- Dependency (Google Test)
find_package(GTest)

if(GTest_FOUND)
    include_directories(${GTEST_INCLUDE_DIRS})
    include_directories(${CMAKE_SOURCE_DIR})
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
        PROPERTIES COMPILE_FLAGS "-Wno-overflow -Wno-literal-range")
    add_executable("test_Riccati_Bessel_logarithmic_derivative"
        test_Riccati_Bessel_logarithmic_derivative.cc)
    target_link_libraries("test_Riccati_Bessel_logarithmic_derivative" ${LIBS})
    add_test(NAME "test_Riccati_Bessel_logarithmic_derivative"
        COMMAND "test_Riccati_Bessel_logarithmic_derivative")

    if(WITH_HWY)
        set_source_files_properties(test_SIMD_Riccati_Bessel.cc
            PROPERTIES COMPILE_FLAGS "-Wno-overflow -Wno-literal-range")
        add_executable("test_SIMD_Riccati_Bessel"
            test_SIMD_Riccati_Bessel.cc)
        target_link_libraries("test_SIMD_Riccati_Bessel" ${LIBS} hwy::hwy)
        target_compile_definitions("test_SIMD_Riccati_Bessel" PRIVATE 
            WITH_HWY
            "HWY_DISABLED_TARGETS=0x4000000"
        )
        add_test(NAME "test_SIMD_Riccati_Bessel"
            COMMAND "test_SIMD_Riccati_Bessel")
    endif()

    # TODO remove -Wno
    set_source_files_properties(test_bulk_sphere.cc
        PROPERTIES COMPILE_FLAGS "-Wno-overflow -Wno-unused-parameter")
    add_executable("test_bulk_sphere" test_bulk_sphere.cc)
    target_link_libraries("test_bulk_sphere" ${LIBS})
    add_test(NAME "test_bulk_sphere"
        COMMAND "test_bulk_sphere")

    if(Boost_FOUND)
        add_executable("test_bulk_sphere_multi_precision" test_bulk_sphere.cc)
        target_compile_options("test_bulk_sphere_multi_precision"
            PRIVATE -DMULTI_PRECISION=100)
        target_link_libraries("test_bulk_sphere_multi_precision" ${LIBS})
        add_test(NAME "test_bulk_sphere_multi_precision"
            COMMAND "test_bulk_sphere_multi_precision")

        add_executable("test_near_field_multi_precision"
            test_near_field.cc)
        target_compile_options("test_near_field_multi_precision"
            PRIVATE -DMULTI_PRECISION=100)
        target_link_libraries("test_near_field_multi_precision" ${LIBS})
        add_test(NAME "test_near_field_multi_precision"
            COMMAND "test_near_field_multi_precision")
    endif()
endif()

# Benchmark
add_executable(bench_special_functions benchmarks/bench_special_functions.cc)
if(WITH_HWY)
    target_compile_definitions(bench_special_functions PRIVATE WITH_HWY)
    target_link_libraries(bench_special_functions PRIVATE hwy::hwy)
endif()

if(WITH_HWY)
    add_executable(smoke_simd smoke_simd.cc)
    target_link_libraries(smoke_simd PRIVATE hwy::hwy)
    target_compile_definitions(smoke_simd PRIVATE WITH_HWY)
    add_test(NAME smoke_simd COMMAND smoke_simd)
endif()
