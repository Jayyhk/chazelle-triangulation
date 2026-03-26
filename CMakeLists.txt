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

add_executable(2.0-tests tests/2.0-tests.cpp ${LIB_SOURCES})
target_include_directories(2.0-tests PRIVATE src)
target_compile_options(2.0-tests PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
    -UNDEBUG
)
add_test(NAME 2.0-tests COMMAND 2.0-tests)

add_executable(2.1-tests tests/2.1-tests.cpp ${LIB_SOURCES})
target_include_directories(2.1-tests PRIVATE src)
target_compile_options(2.1-tests PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
    -UNDEBUG
)
add_test(NAME 2.1-tests COMMAND 2.1-tests)

add_executable(2.2-tests tests/2.2-tests.cpp ${LIB_SOURCES})
target_include_directories(2.2-tests PRIVATE src)
target_compile_options(2.2-tests PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
    -UNDEBUG
)
add_test(NAME 2.2-tests COMMAND 2.2-tests)

add_executable(2.3-tests tests/2.3-tests.cpp ${LIB_SOURCES})
target_include_directories(2.3-tests PRIVATE src)
target_compile_options(2.3-tests PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
    -UNDEBUG
)
add_test(NAME 2.3-tests COMMAND 2.3-tests)

add_executable(2.4-tests tests/2.4-tests.cpp ${LIB_SOURCES})
target_include_directories(2.4-tests PRIVATE src)
target_compile_options(2.4-tests PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
    -UNDEBUG
)
add_test(NAME 2.4-tests COMMAND 2.4-tests)
