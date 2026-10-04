cmake_minimum_required(VERSION 3.20)
project(chazelle-triangulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_path(GMP_INCLUDE_DIR gmp.h REQUIRED)
find_library(GMP_LIBRARY NAMES gmp libgmp.so.10 REQUIRED)
add_library(GMP::GMP UNKNOWN IMPORTED)
set_target_properties(GMP::GMP PROPERTIES IMPORTED_LOCATION "${GMP_LIBRARY}"
                                          INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)

add_library(chazelle-warnings INTERFACE)
target_compile_options(chazelle-warnings INTERFACE -Wall -Wextra -Wpedantic -Wconversion
                                                   -Wsign-conversion)

include(CTest)
add_subdirectory(src)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
