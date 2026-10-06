set(CHAZELLE_TESTS
    s2.0-tests
    s2.1-tests
    s2.2-tests
    s2.3-tests
    s2.4-tests
    s2.5-tests
    s2-e2e-tests
    s3.0-tests
    s3.1-tests
    s3.2-tests
    s3.3-tests
    s3.4-tests
    s3-e2e-tests
    s4.0-tests
    s4.1-tests
    s4.2-tests
    s4-e2e-tests
    down-phase-tests
    lt79-tests
    exact-tests
    trapezoid-tests
    unimonotone-tests
    triangulation-tests
    equal-y-tests)

foreach(name IN LISTS CHAZELLE_TESTS)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE chazelle-core-checked chazelle-warnings)
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
    add_test(NAME ${name} COMMAND ${name})
endforeach()

target_compile_options(s2-e2e-tests PRIVATE -O2)

add_executable(down-phase-production-tests down-phase-production-tests.cpp)
target_link_libraries(down-phase-production-tests PRIVATE chazelle-core chazelle-warnings)
set_property(TARGET down-phase-production-tests PROPERTY RUNTIME_OUTPUT_DIRECTORY
                                                         "${PROJECT_BINARY_DIR}")
add_test(NAME down-phase-production-tests COMMAND down-phase-production-tests)

add_executable(trapezoid-production-tests trapezoid-production-tests.cpp)
target_link_libraries(trapezoid-production-tests PRIVATE chazelle-core chazelle-warnings)
set_property(TARGET trapezoid-production-tests PROPERTY RUNTIME_OUTPUT_DIRECTORY
                                                        "${PROJECT_BINARY_DIR}")
add_test(NAME trapezoid-production-tests COMMAND trapezoid-production-tests)

add_executable(unimonotone-production-tests unimonotone-production-tests.cpp)
target_link_libraries(unimonotone-production-tests PRIVATE chazelle-core chazelle-warnings)
set_property(TARGET unimonotone-production-tests PROPERTY RUNTIME_OUTPUT_DIRECTORY
                                                          "${PROJECT_BINARY_DIR}")
add_test(NAME unimonotone-production-tests COMMAND unimonotone-production-tests)

add_executable(triangulation-production-tests triangulation-production-tests.cpp)
target_link_libraries(triangulation-production-tests PRIVATE chazelle-core chazelle-warnings)
set_property(TARGET triangulation-production-tests PROPERTY RUNTIME_OUTPUT_DIRECTORY
                                                            "${PROJECT_BINARY_DIR}")
add_test(NAME triangulation-production-tests COMMAND triangulation-production-tests)

add_executable(e2e-tests e2e-tests.cpp)
target_link_libraries(e2e-tests PRIVATE chazelle-core chazelle-warnings)
set_property(TARGET e2e-tests PROPERTY RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
add_test(NAME e2e-tests COMMAND e2e-tests)

add_executable(visualizer-tests visualizer-tests.cpp)
target_link_libraries(visualizer-tests PRIVATE chazelle-visualizer chazelle-warnings)
set_property(TARGET visualizer-tests PROPERTY RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
add_test(NAME visualizer-tests COMMAND visualizer-tests)

add_test(
    NAME cli-tests
    COMMAND
        "${CMAKE_COMMAND}" "-DCHAZELLE_EXECUTABLE=$<TARGET_FILE:chazelle>"
        "-DCLI_TEST_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/cli-tests" -P
        "${CMAKE_CURRENT_SOURCE_DIR}/cli-tests.cmake")

add_executable(animation-trace-tests animation-trace-tests.cpp)
target_link_libraries(animation-trace-tests PRIVATE chazelle-core-checked chazelle-animation-checked
                                                   chazelle-warnings)
set_property(TARGET animation-trace-tests PROPERTY RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
add_test(NAME animation-trace-tests
         COMMAND animation-trace-tests "${CMAKE_CURRENT_BINARY_DIR}/animation-traces")
set_tests_properties(animation-trace-tests PROPERTIES FIXTURES_SETUP animation-traces)

find_package(Python3 3.11 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME animation-data-tests
             COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/animation-data-tests.py"
                     "${CMAKE_CURRENT_BINARY_DIR}/animation-traces")
    set_tests_properties(animation-data-tests PROPERTIES FIXTURES_REQUIRED animation-traces)
endif()
