set(CHAZELLE_MODULE_SOURCES
    polygon/polygon.cpp
    submap/submap.cpp
    submap/boundary_geometry.cpp
    submap/boundary_queries.cpp
    submap/chord_inventory.cpp
    submap/granularity_queries.cpp
    submap/normalization.cpp
    submap/tree_decomposition.cpp
    submap/validation.cpp
    separator/planar_separator.cpp
    triangulation/trapezoids.cpp
    triangulation/unimonotone.cpp
    triangulation/triangulation.cpp
    merge/conformality.cpp
    merge/fusion.cpp
    merge/fusion_rebuild.cpp
    merge/granularity.cpp
    merge/merge.cpp
    merge/ray_shooting.cpp
    visibility/chain.cpp
    visibility/naive_visibility.cpp
    visibility/bounded_regions.cpp
    visibility/region_boundary.cpp
    visibility/region_oracles.cpp
    visibility/down_phase.cpp
    visibility/up_phase.cpp)

list(TRANSFORM CHAZELLE_MODULE_SOURCES PREPEND "algorithm/" OUTPUT_VARIABLE CHAZELLE_SOURCES)

function(add_chazelle_library target)
    add_library(${target} STATIC ${ARGN})
    target_include_directories(${target} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
    target_compile_features(${target} PUBLIC cxx_std_20)
    target_link_libraries(
        ${target}
        PUBLIC GMP::GMP
        PRIVATE chazelle-warnings)
endfunction()

add_chazelle_library(chazelle-core ${CHAZELLE_SOURCES})
add_library(Chazelle::Core ALIAS chazelle-core)

list(TRANSFORM CHAZELLE_MODULE_SOURCES PREPEND "animation/" OUTPUT_VARIABLE CHAZELLE_ANIMATION_SOURCES)
list(APPEND CHAZELLE_ANIMATION_SOURCES animation/trace.cpp animation/triangulate.cpp)
add_chazelle_library(chazelle-animation ${CHAZELLE_ANIMATION_SOURCES})

add_library(chazelle-visualizer STATIC algorithm/visualizer/triangulation_svg.cpp)
target_link_libraries(
    chazelle-visualizer
    PUBLIC chazelle-core
    PRIVATE chazelle-warnings)

set(CHAZELLE_ANIMATION_PYTHON "${PROJECT_SOURCE_DIR}/.venv/bin/python" CACHE FILEPATH
    "Python interpreter with src/animation/requirements.txt installed")
configure_file(animation/config.h.in animation/config.h @ONLY)
add_executable(chazelle main.cpp animation/render.cpp)
target_include_directories(chazelle PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(chazelle PRIVATE chazelle-visualizer chazelle-animation chazelle-warnings)
set_target_properties(chazelle PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

if(BUILD_TESTING)
    add_chazelle_library(chazelle-core-checked ${CHAZELLE_SOURCES})
    target_compile_options(chazelle-core-checked PUBLIC -UNDEBUG)
    target_compile_definitions(chazelle-core-checked PUBLIC CHAZELLE_EXPENSIVE_ASSERTS)

    add_chazelle_library(chazelle-animation-checked ${CHAZELLE_ANIMATION_SOURCES})
    target_compile_options(chazelle-animation-checked PUBLIC -UNDEBUG)
    target_compile_definitions(chazelle-animation-checked PUBLIC CHAZELLE_EXPENSIVE_ASSERTS)
endif()
