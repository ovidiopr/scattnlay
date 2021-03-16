# Sources for python extension
set(_scattnlay_python_sources
        ${CMAKE_CURRENT_LIST_DIR}/nmie.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie-pybind11.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie-pybind11.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie-precision.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie-impl.cc
        ${CMAKE_CURRENT_LIST_DIR}/pb11_wrapper.cc)

# Define python extension
add_library(python3-scattnlay SHARED ${_scattnlay_python_sources})
target_link_libraries(python3-scattnlay PRIVATE Boost::headers ${Python_LIBRARIES})
target_include_directories(python3-scattnlay PRIVATE ${NUMPY_INCLUDE_DIR} ${PYBIND11_INCLUDE_DIR})
set_target_properties(python3-scattnlay PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set_target_properties(
        python3-scattnlay
        PROPERTIES
        PREFIX ""
        OUTPUT_NAME "scattnlay"
        SUFFIX "${EXT_SUFFIX}"
        LINKER_LANGUAGE C
)

# Sources for far field calculation
set(_scattnlay_farfield_sources
        ${CMAKE_CURRENT_LIST_DIR}/farfield.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie.cc)

# Define exe for far field calculation
add_executable(farfield ${_scattnlay_farfield_sources})
target_link_libraries(farfield PRIVATE Boost::headers)
set_target_properties(farfield PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Sources for near field calculation
set(_scattnlay_nearfield_sources
        ${CMAKE_CURRENT_LIST_DIR}/nearfield.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie.cc)

# Define exe for near field calculation
add_executable(nearfield ${_scattnlay_nearfield_sources})
target_link_libraries(nearfield PRIVATE Boost::headers)
set_target_properties(nearfield PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Rename files to match precision
if (${ENABLE_MP})
    set_property(TARGET python3-scattnlay APPEND_STRING PROPERTY OUTPUT_NAME "_mp")
    set_property(TARGET farfield APPEND_STRING PROPERTY OUTPUT_NAME "farfield-mp")
    set_property(TARGET nearfield APPEND_STRING PROPERTY OUTPUT_NAME "nearfield-mp")
else ()
    set_property(TARGET python3-scattnlay APPEND_STRING PROPERTY OUTPUT_NAME "_dp")
    set_property(TARGET farfield APPEND_STRING PROPERTY OUTPUT_NAME "farfield-dp")
    set_property(TARGET nearfield APPEND_STRING PROPERTY OUTPUT_NAME "nearfield-dp")
endif ()
