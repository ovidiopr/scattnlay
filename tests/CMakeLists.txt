cmake_minimum_required(VERSION 3.24)
project(scattnlay_tests C CXX)

# -- Dependency (Google Test)
find_package(GTest)

if(GTest_FOUND)
    include_directories(${CMAKE_SOURCE_DIR})

    # -- Output tests in directory
    add_executable("test_near_field"
        test_near_field.cc)
    target_link_libraries("test_near_field" PRIVATE 
        GTest::gtest GTest::gtest_main 
        project_warnings
    )
    add_test(NAME "test_near_field"
        COMMAND "test_near_field")

    # In included file test_spec_functions_data.hpp there are results of multiple
    # precision computation that may overflow double precision at compile time.
    add_executable("test_Riccati_Bessel_logarithmic_derivative"
        test_Riccati_Bessel_logarithmic_derivative.cc)
    target_compile_options("test_Riccati_Bessel_logarithmic_derivative" 
        PRIVATE -Wno-overflow -Wno-literal-range)
    target_link_libraries("test_Riccati_Bessel_logarithmic_derivative" PRIVATE 
        GTest::gtest GTest::gtest_main
        project_optimization
    )
    add_test(NAME "test_Riccati_Bessel_logarithmic_derivative"
        COMMAND "test_Riccati_Bessel_logarithmic_derivative")

    if(WITH_HWY)
        add_executable("test_SIMD_Riccati_Bessel"
            test_SIMD_Riccati_Bessel.cc)
        target_compile_options("test_SIMD_Riccati_Bessel"
            PRIVATE -Wno-overflow -Wno-literal-range)
        target_link_libraries("test_SIMD_Riccati_Bessel" PRIVATE 
            GTest::gtest GTest::gtest_main
            hwy::hwy
            nmie_obj_dp
            project_optimization
        )
        target_compile_definitions("test_SIMD_Riccati_Bessel" PRIVATE 
            WITH_HWY
            "HWY_DISABLED_TARGETS=0x4000000"
        )
        add_test(NAME "test_SIMD_Riccati_Bessel"
            COMMAND "test_SIMD_Riccati_Bessel")

        add_executable("test_bulk_sphere_simd"
            test_bulk_sphere_simd.cc)
        target_link_libraries("test_bulk_sphere_simd" PRIVATE 
            GTest::gtest GTest::gtest_main
            hwy::hwy
            nmie_obj_dp
            project_optimization
        )
        target_compile_definitions("test_bulk_sphere_simd" PRIVATE 
            WITH_HWY
        )
        add_test(NAME "test_bulk_sphere_simd"
            COMMAND "test_bulk_sphere_simd")

    add_executable(test_farfield_simd_benchmark test_farfield_simd_benchmark.cc)
    target_link_libraries(test_farfield_simd_benchmark PRIVATE 
        GTest::gtest GTest::gtest_main hwy::hwy nmie_obj_dp project_optimization)
    target_compile_definitions(test_farfield_simd_benchmark PRIVATE 
        WITH_HWY "HWY_DISABLED_TARGETS=0x4000000")
    add_test(NAME test_farfield_simd_benchmark COMMAND test_farfield_simd_benchmark)

    add_executable(test_nearfield_simd_benchmark test_nearfield_simd_benchmark.cc)
    target_link_libraries(test_nearfield_simd_benchmark PRIVATE 
        GTest::gtest GTest::gtest_main hwy::hwy nmie_obj_dp project_optimization)
    target_compile_definitions(test_nearfield_simd_benchmark PRIVATE 
        WITH_HWY "HWY_DISABLED_TARGETS=0x4000000")
    add_test(NAME test_nearfield_simd_benchmark COMMAND test_nearfield_simd_benchmark)
    endif()

    add_executable("test_bulk_sphere" test_bulk_sphere.cc)
    target_compile_options("test_bulk_sphere" 
        PRIVATE -Wno-overflow -Wno-unused-parameter)
    target_link_libraries("test_bulk_sphere" PRIVATE 
        GTest::gtest GTest::gtest_main
        project_optimization
    )
    add_test(NAME "test_bulk_sphere"
        COMMAND "test_bulk_sphere")

    if(Boost_FOUND)
        add_executable("test_bulk_sphere_multi_precision" test_bulk_sphere.cc)
        target_compile_options("test_bulk_sphere_multi_precision"
            PRIVATE -Wno-overflow -Wno-unused-parameter)
        target_compile_definitions("test_bulk_sphere_multi_precision" PRIVATE MULTI_PRECISION=100)
        target_link_libraries("test_bulk_sphere_multi_precision" PRIVATE
            GTest::gtest GTest::gtest_main
            Boost::headers
            project_optimization
        )
        add_test(NAME "test_bulk_sphere_multi_precision"
            COMMAND "test_bulk_sphere_multi_precision")

        add_executable("test_near_field_multi_precision"
            test_near_field.cc)
        target_compile_definitions("test_near_field_multi_precision" PRIVATE MULTI_PRECISION=100)
        target_link_libraries("test_near_field_multi_precision" PRIVATE
            GTest::gtest GTest::gtest_main
            Boost::headers
            project_warnings
        )
        add_test(NAME "test_near_field_multi_precision"
            COMMAND "test_near_field_multi_precision")
    endif()
endif()

# Benchmark
add_executable(bench_special_functions benchmarks/bench_special_functions.cc)
target_link_libraries(bench_special_functions PRIVATE 
    nmie_obj_dp
    project_optimization
)
if(WITH_HWY)
    target_compile_definitions(bench_special_functions PRIVATE WITH_HWY)
endif()

if(WITH_HWY)
    add_executable(smoke_simd smoke_simd.cc)
    target_link_libraries(smoke_simd PRIVATE 
        hwy::hwy
        project_warnings
        project_optimization
    )
    target_compile_definitions(smoke_simd PRIVATE WITH_HWY)
    add_test(NAME smoke_simd COMMAND smoke_simd)
endif()
