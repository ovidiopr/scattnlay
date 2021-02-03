set(_scattnlay_python_sources
        ${CMAKE_CURRENT_LIST_DIR}/nmie.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie-pybind11.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie-pybind11.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie-precision.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie-impl.cc
        ${CMAKE_CURRENT_LIST_DIR}/pb11_wrapper.cc)

add_library(python3-scattnlay SHARED ${_scattnlay_python_sources})
target_link_libraries(python3-scattnlay ${PYTHON_LIBRARIES})

set_target_properties(
        python3-scattnlay
        PROPERTIES
        PREFIX ""
        OUTPUT_NAME "python3-scattnlay"
        LINKER_LANGUAGE C
)

add_executable(farfield
        ${CMAKE_CURRENT_LIST_DIR}/farfield.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie.cc)

add_executable(nearfield
        ${CMAKE_CURRENT_LIST_DIR}/nearfield.cc
        ${CMAKE_CURRENT_LIST_DIR}/nmie.hpp
        ${CMAKE_CURRENT_LIST_DIR}/nmie.cc)

if (${ENABLE_MP})
    set_property(TARGET python3-scattnlay APPEND_STRING PROPERTY OUTPUT_NAME "_mp")
    set_property(TARGET farfield APPEND_STRING PROPERTY OUTPUT_NAME "farfield-mp")
    set_property(TARGET nearfield APPEND_STRING PROPERTY OUTPUT_NAME "nearfield-mp")
else ()
    set_property(TARGET python3-scattnlay APPEND_STRING PROPERTY OUTPUT_NAME "_sp")
    set_property(TARGET farfield APPEND_STRING PROPERTY OUTPUT_NAME "farfield-sp")
    set_property(TARGET nearfield APPEND_STRING PROPERTY OUTPUT_NAME "nearfield-sp")
endif ()
