## Sources for python extension
#set(_scattnlay_python_sources
#        ${CMAKE_CURRENT_LIST_DIR}/nmie.hpp
#        ${CMAKE_CURRENT_LIST_DIR}/nmie.cc
#        ${CMAKE_CURRENT_LIST_DIR}/nmie-precision.hpp
#        ${CMAKE_CURRENT_LIST_DIR}/special-functions-impl.hpp
#        ${CMAKE_CURRENT_LIST_DIR}/nmie-basic.hpp
#        ${CMAKE_CURRENT_LIST_DIR}/nmie-nearfield.hpp
#        ${CMAKE_CURRENT_LIST_DIR}/pb11_wrapper.cc)
#
## Define python extension
#add_library(python3-scattnlay SHARED ${_scattnlay_python_sources})
#if (${ENABLE_MP} AND ${Boost_FOUND})
#    target_link_libraries(python3-scattnlay PRIVATE Boost::headers ${Python_LIBRARIES})
#else()
#    target_link_libraries(python3-scattnlay PRIVATE ${Python_LIBRARIES})
#endif()
#target_include_directories(python3-scattnlay PRIVATE ${NUMPY_INCLUDE_DIR} ${PYBIND11_INCLUDE_DIR})
#set_target_properties(python3-scattnlay PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

#set_target_properties(
#        python3-scattnlay
#        PROPERTIES
#        PREFIX ""
#        OUTPUT_NAME "scattnlay"
#        SUFFIX "${EXT_SUFFIX}"
#        LINKER_LANGUAGE C
#)

# Sources for far field calculation
set(_scattnlay_farfield_sources
        ${CMAKE_CURRENT_LIST_DIR}/farfield.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie.cc)

# Sources for near field calculation
set(_scattnlay_nearfield_sources
        ${CMAKE_CURRENT_LIST_DIR}/nearfield.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie.cc)

# Double Precision Targets
add_executable(scattnlay-dp ${_scattnlay_farfield_sources})
set_target_properties(scattnlay-dp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(WITH_HWY)
    target_link_libraries(scattnlay-dp PRIVATE hwy::hwy)
    # 0x4000000 disables HWY_NEON_BF16 (bit 26) which causes namespace errors on Clang/macOS
    target_compile_definitions(scattnlay-dp PRIVATE 
        WITH_HWY 
        "HWY_DISABLED_TARGETS=0x4000000"
    )
endif()

add_executable(fieldnlay-dp ${_scattnlay_nearfield_sources})
set_target_properties(fieldnlay-dp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Multi Precision Targets
if (Boost_FOUND)
    add_executable(scattnlay-mp ${_scattnlay_farfield_sources})
    target_compile_definitions(scattnlay-mp PRIVATE MULTI_PRECISION=100)
    target_link_libraries(scattnlay-mp PRIVATE Boost::headers)
    set_target_properties(scattnlay-mp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

    add_executable(fieldnlay-mp ${_scattnlay_nearfield_sources})
    target_compile_definitions(fieldnlay-mp PRIVATE MULTI_PRECISION=100)
    target_link_libraries(fieldnlay-mp PRIVATE Boost::headers)
    set_target_properties(fieldnlay-mp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# WebAssembly Target
if(EMSCRIPTEN)
    add_executable(nmiejs nmie-js-wrapper.cc)
    target_compile_options(nmiejs PRIVATE -O3 -msimd128)
    target_link_options(nmiejs PRIVATE 
        --bind 
        -lm 
        -sMODULARIZE=1 
        -sWASM=1 
        -sALLOW_MEMORY_GROWTH=1 
        -sEXPORT_NAME=createNmiejsModule 
        -sENVIRONMENT=web,node
        -O3
        -msimd128
    )
    set_target_properties(nmiejs PROPERTIES SUFFIX ".js")

    add_custom_command(TARGET nmiejs POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:nmiejs>
            ${CMAKE_SOURCE_DIR}/guiapp/src/nmiejs.js
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE_DIR:nmiejs>/nmiejs.wasm
            ${CMAKE_SOURCE_DIR}/guiapp/src/nmiejs.wasm
        COMMENT "Copying nmiejs artifacts to guiapp/src"
    )
endif()

# Install targets
install(TARGETS scattnlay-dp fieldnlay-dp DESTINATION bin)
if(Boost_FOUND)
    install(TARGETS scattnlay-mp fieldnlay-mp DESTINATION bin)
endif()
if(EMSCRIPTEN)
    install(TARGETS nmiejs DESTINATION bin)
endif()
