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
enable_testing()

# Collect library sources (everything except main.cpp)
file(GLOB_RECURSE LIB_SOURCES src/*.cpp)
list(FILTER LIB_SOURCES EXCLUDE REGEX "src/main\\.cpp$")

function(add_chazelle_test NAME)
    add_executable(${NAME} tests/${NAME}.cpp ${LIB_SOURCES})
    target_include_directories(${NAME} PRIVATE src)
    target_compile_options(${NAME} PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
        -UNDEBUG
    )
    add_test(NAME ${NAME} COMMAND ${NAME})
endfunction()

add_chazelle_test(s2.0-tests)
add_chazelle_test(s2.1-tests)
add_chazelle_test(s2.2-tests)
add_chazelle_test(s2.3-tests)
add_chazelle_test(s2.4-tests)
add_chazelle_test(s2.5-tests)
add_chazelle_test(s2-e2e-tests)
add_chazelle_test(s3.0-tests)
target_compile_options(s2-e2e-tests PRIVATE -O2)
