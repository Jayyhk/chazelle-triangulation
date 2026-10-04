set(CHAZELLE_SOURCES
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

function(add_chazelle_library target)
    add_library(${target} STATIC ${CHAZELLE_SOURCES})
    target_include_directories(${target} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
    target_compile_features(${target} PUBLIC cxx_std_20)
    target_link_libraries(
        ${target}
        PUBLIC GMP::GMP
        PRIVATE chazelle-warnings)
endfunction()

add_chazelle_library(chazelle-core)
add_library(Chazelle::Core ALIAS chazelle-core)

add_library(chazelle-visualizer STATIC visualizer/triangulation_svg.cpp)
target_link_libraries(
    chazelle-visualizer
    PUBLIC chazelle-core
    PRIVATE chazelle-warnings)

add_executable(chazelle main.cpp)
target_link_libraries(chazelle PRIVATE chazelle-visualizer chazelle-warnings)
set_target_properties(chazelle PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

if(BUILD_TESTING)
    add_chazelle_library(chazelle-core-checked)
    target_compile_options(chazelle-core-checked PUBLIC -UNDEBUG)
    target_compile_definitions(chazelle-core-checked PUBLIC CHAZELLE_EXPENSIVE_ASSERTS)
endif()
