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
