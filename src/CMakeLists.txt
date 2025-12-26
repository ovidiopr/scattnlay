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

# --- PHASE 2: OBJECT LIBRARIES ---
# Compile nmie.cc only once per precision level

# Double Precision Object Library
add_library(nmie_obj_dp OBJECT ${CMAKE_CURRENT_LIST_DIR}/nmie.cc)
target_include_directories(nmie_obj_dp PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(nmie_obj_dp PUBLIC project_warnings project_optimization)

# Multi Precision Object Library
if(Boost_FOUND)
    add_library(nmie_obj_mp OBJECT ${CMAKE_CURRENT_LIST_DIR}/nmie.cc)
    target_compile_definitions(nmie_obj_mp PRIVATE MULTI_PRECISION=100)
    target_link_libraries(nmie_obj_mp PUBLIC Boost::headers project_warnings project_optimization)
    target_include_directories(nmie_obj_mp PUBLIC ${CMAKE_CURRENT_LIST_DIR})
endif()

# Double Precision Targets
add_executable(scattnlay-dp ${CMAKE_CURRENT_LIST_DIR}/farfield.cc)
target_link_libraries(scattnlay-dp PRIVATE nmie_obj_dp project_warnings project_optimization)
set_target_properties(scattnlay-dp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(WITH_HWY)
    target_link_libraries(nmie_obj_dp PUBLIC hwy::hwy)
    # 0x4000000 disables HWY_NEON_BF16 (bit 26) which causes namespace errors on Clang/macOS
    target_compile_definitions(nmie_obj_dp PUBLIC 
        WITH_HWY 
        "HWY_DISABLED_TARGETS=0x4000000"
    )
endif()

add_executable(fieldnlay-dp ${CMAKE_CURRENT_LIST_DIR}/nearfield.cc)
target_link_libraries(fieldnlay-dp PRIVATE nmie_obj_dp project_warnings project_optimization)
set_target_properties(fieldnlay-dp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Multi Precision Targets
if (Boost_FOUND)
    add_executable(scattnlay-mp ${CMAKE_CURRENT_LIST_DIR}/farfield.cc)
    target_link_libraries(scattnlay-mp PRIVATE nmie_obj_mp project_warnings project_optimization)
    set_target_properties(scattnlay-mp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

    add_executable(fieldnlay-mp ${CMAKE_CURRENT_LIST_DIR}/nearfield.cc)
    target_link_libraries(fieldnlay-mp PRIVATE nmie_obj_mp project_warnings project_optimization)
    set_target_properties(fieldnlay-mp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# WebAssembly Target
if(EMSCRIPTEN)
    # WASM-specific interface library for clean encapsulation
    add_library(wasm_options INTERFACE)
    target_compile_options(wasm_options INTERFACE -O3 -msimd128)
    target_link_options(wasm_options INTERFACE 
        --bind 
        -sMODULARIZE=1 
        -sEXPORT_ES6=1                # Added: Ensure standard ES import
        -sWASM=1 
        -sALLOW_MEMORY_GROWTH=1 
        -sEXPORT_NAME=createNmiejsModule 
        -sENVIRONMENT=web,node        # Added: node support for Vitest
        -O3
        -msimd128
    )

    add_executable(nmiejs nmie-js-wrapper.cc)
    target_link_libraries(nmiejs PRIVATE wasm_options)
    set_target_properties(nmiejs PROPERTIES SUFFIX ".js")

    # Update copy command to put files in public/wasm for reliable asset fetching
    add_custom_command(TARGET nmiejs POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/guiapp/public/wasm
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:nmiejs> ${CMAKE_SOURCE_DIR}/guiapp/public/wasm/nmiejs.js
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE_DIR:nmiejs>/nmiejs.wasm ${CMAKE_SOURCE_DIR}/guiapp/public/wasm/nmiejs.wasm
    )
endif()

# Install targets
install(TARGETS scattnlay-dp fieldnlay-dp 
    EXPORT scattnlayTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

if(Boost_FOUND)
    install(TARGETS scattnlay-mp fieldnlay-mp 
        EXPORT scattnlayTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
endif()

if(EMSCRIPTEN)
    install(TARGETS nmiejs 
        EXPORT scattnlayTargets
        RUNTIME DESTINATION bin
    )
endif()

# Install header files for C++ library usage
install(FILES
    nmie.hpp
    nmie-precision.hpp
    special-functions-impl.hpp
    nmie-basic.hpp
    nmie-nearfield.hpp
    DESTINATION include/scattnlay
)
