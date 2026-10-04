file(MAKE_DIRECTORY "${CLI_TEST_DIRECTORY}")

function(run_cli name input)
    set(input_file "${CLI_TEST_DIRECTORY}/${name}.txt")
    file(WRITE "${input_file}" "${input}")
    execute_process(
        COMMAND "${CHAZELLE_EXECUTABLE}" ${ARGN}
        INPUT_FILE "${input_file}"
        WORKING_DIRECTORY "${CLI_TEST_DIRECTORY}"
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error
        RESULT_VARIABLE result
        TIMEOUT 30)
    set(cli_output
        "${output}"
        PARENT_SCOPE)
    set(cli_error
        "${error}"
        PARENT_SCOPE)
    set(cli_result
        "${result}"
        PARENT_SCOPE)
endfunction()

function(check_polygon name input vertex_count)
    run_cli("${name}" "${input}" ${ARGN})
    if(NOT cli_result STREQUAL "0" OR NOT cli_error STREQUAL "")
        message(FATAL_ERROR "${name}: expected success, got ${cli_result}: ${cli_error}")
    endif()
    math(EXPR triangle_count "${vertex_count} - 2")
    math(EXPR line_count "${triangle_count} + 1")
    string(REGEX MATCHALL "[^\r\n]+" lines "${cli_output}")
    list(LENGTH lines actual_line_count)
    list(GET lines 0 actual_triangle_count)
    if(NOT actual_line_count EQUAL line_count OR NOT actual_triangle_count STREQUAL
                                                 "${triangle_count}")
        message(FATAL_ERROR "${name}: unexpected triangle count or output lines: ${cli_output}")
    endif()
    list(REMOVE_AT lines 0)
    set(used_vertices)
    foreach(line IN LISTS lines)
        if(NOT line MATCHES "^([0-9]+) ([0-9]+) ([0-9]+)$")
            message(FATAL_ERROR "${name}: unexpected triangle format: ${line}")
        endif()
        set(first "${CMAKE_MATCH_1}")
        set(second "${CMAKE_MATCH_2}")
        set(third "${CMAKE_MATCH_3}")
        if(first GREATER_EQUAL vertex_count
           OR second GREATER_EQUAL vertex_count
           OR third GREATER_EQUAL vertex_count
           OR first EQUAL second
           OR second EQUAL third
           OR third EQUAL first)
            message(FATAL_ERROR "${name}: invalid original triangle indices: ${line}")
        endif()
        list(APPEND used_vertices "${first}" "${second}" "${third}")
    endforeach()
    list(REMOVE_DUPLICATES used_vertices)
    list(LENGTH used_vertices used_count)
    if(NOT used_count EQUAL vertex_count)
        message(FATAL_ERROR "${name}: output omits an original vertex")
    endif()
endfunction()

function(check_invalid name input)
    run_cli("${name}" "${input}" ${ARGN})
    if(NOT cli_result STREQUAL "1"
       OR NOT cli_output STREQUAL ""
       OR cli_error STREQUAL "")
        message(FATAL_ERROR "${name}: expected an input error, got ${cli_result}: ${cli_output}")
    endif()
endfunction()

check_polygon(triangle "3\n0 0\n4 0\n2 3\n" 3)
check_polygon(counterclockwise "4\n0 0\n4 0\n4 4\n0 4\n" 4)
check_polygon(clockwise "4\n0 4\n4 4\n4 0\n0 0\n" 4)
check_polygon(concave "8\n0 0\n6 0\n6 6\n4 6\n4 2\n2 2\n2 6\n0 6\n" 8)
check_polygon(collinear "6\n0 0\n2 0\n4 0\n4 2\n4 4\n0 4\n" 6)
check_polygon(decimal "3\n+.0 -0.\n4.0e+0 -.00E-4\n2.00E0 .3e1\n" 3)
check_polygon(large_integer "3\n9007199254740992 0\n9007199254740993 0\n9007199254740992 1\n" 3)
check_polygon(small_exponent "3\n1e-400 0\n2e-400 0\n1e-400 1e-400\n" 3)
check_polygon(large_exponent "3\n1e400 0\n2e400 0\n1e400 1e400\n" 3)

check_invalid(empty "")
check_invalid(negative_count "-3\n")
check_invalid(small_count "2\n")
check_invalid(noninteger_count "3.0\n")
check_invalid(overflow_count "18446744073709551616\n")
check_invalid(missing_coordinates "3\n0 0\n4 0\n")
check_invalid(extra_coordinates "3\n0 0\n4 0\n2 3\n1 1\n")
foreach(
    token IN
    ITEMS nan
          inf
          1x
          1..0
          .
          1e
          1e-
          1e+
          -
          1e18446744073709551616
          1.0e-18446744073709551615)
    check_invalid("invalid_${token}" "3\n0 0\n4 ${token}\n2 3\n")
endforeach()

file(REMOVE_RECURSE "${CLI_TEST_DIRECTORY}/images")
string(TIMESTAMP IMAGE_DATE "%Y-%m-%d")
check_polygon(default_image "4\n0 0\n4 0\n4 4\n0 4\n" 4 --visualize)
file(READ "${CLI_TEST_DIRECTORY}/images/triangulation-${IMAGE_DATE}.svg" default_image)
if(NOT default_image MATCHES "<svg xmlns=\"http://www.w3.org/2000/svg\""
   OR NOT default_image MATCHES "id=\"triangle-1\""
   OR NOT default_image MATCHES "</svg>")
    message(FATAL_ERROR "Default visualization does not contain the complete SVG triangulation")
endif()

string(TIMESTAMP IMAGE_DATE "%Y-%m-%d")
set(TRIANGLE_IMAGE_PATH "${CLI_TEST_DIRECTORY}/images/triangle image & $-${IMAGE_DATE}.svg")
file(REMOVE "${TRIANGLE_IMAGE_PATH}")
check_polygon(named_image "3\n0 0\n4 0\n2 3\n" 3 --visualize "triangle image & $.svg")
file(READ "${TRIANGLE_IMAGE_PATH}" named_image)
if(NOT named_image MATCHES "id=\"triangle-0\"" OR named_image MATCHES "id=\"triangle-1\"")
    message(FATAL_ERROR "Named visualization has an incorrect triangle count")
endif()

run_cli(help "" --help)
if(NOT cli_result STREQUAL "0"
   OR NOT cli_error STREQUAL ""
   OR NOT cli_output MATCHES "Usage: chazelle.*--visualize")
    message(FATAL_ERROR "Help does not describe the visualizer flag")
endif()

check_invalid(unknown_option "" --unknown)
check_invalid(wrong_image_extension "" --visualize triangle.png)
check_invalid(extra_arguments "" --visualize triangle.svg extra)
check_invalid(image_path_instead_of_filename "" --visualize elsewhere/image.svg)
file(REMOVE_RECURSE "${CLI_TEST_DIRECTORY}/images")
file(WRITE "${CLI_TEST_DIRECTORY}/images" "not a directory")
check_invalid(image_directory_is_file "3\n0 0\n4 0\n2 3\n" --visualize)
file(REMOVE "${CLI_TEST_DIRECTORY}/images")
