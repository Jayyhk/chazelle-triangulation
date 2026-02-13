cmake_minimum_required(VERSION 3.20)
project(chazelle-triangulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Gather all source files
file(GLOB_RECURSE SOURCES src/*.cpp)

add_executable(chazelle ${SOURCES})
target_include_directories(chazelle PRIVATE src)

# Warnings
target_compile_options(chazelle PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
)

# --- Tests ---
# Collect library sources (everything except main.cpp)
file(GLOB_RECURSE LIB_SOURCES src/*.cpp)
list(FILTER LIB_SOURCES EXCLUDE REGEX "src/main\\.cpp$")

add_executable(test_fm tests/test_fm.cpp ${LIB_SOURCES})
target_include_directories(test_fm PRIVATE src)
target_compile_options(test_fm PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
)

add_executable(test_geometry tests/test_geometry.cpp ${LIB_SOURCES})
target_include_directories(test_geometry PRIVATE src)
target_compile_options(test_geometry PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
)

add_executable(test_triangulate tests/test_triangulate.cpp ${LIB_SOURCES})
target_include_directories(test_triangulate PRIVATE src)
target_compile_options(test_triangulate PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
)

add_executable(test_separator tests/test_separator.cpp ${LIB_SOURCES})
target_include_directories(test_separator PRIVATE src)
target_compile_options(test_separator PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
)

add_executable(test_integration tests/test_integration.cpp ${LIB_SOURCES})
target_include_directories(test_integration PRIVATE src)
target_compile_options(test_integration PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
)
