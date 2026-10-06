set(CHAZELLE_DEV_BIN "${PROJECT_SOURCE_DIR}/.venv/bin")
file(GLOB_RECURSE CHAZELLE_CPP_FORMAT_FILES CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/src/*.cpp"
     "${PROJECT_SOURCE_DIR}/src/*.h"
     "${PROJECT_SOURCE_DIR}/tests/*.cpp"
     "${PROJECT_SOURCE_DIR}/tests/*.h")

add_custom_target(format
    COMMAND "${CHAZELLE_DEV_BIN}/clang-format" -i ${CHAZELLE_CPP_FORMAT_FILES}
    COMMAND "${CHAZELLE_DEV_BIN}/ruff" check --fix src/animation tests
    COMMAND "${CHAZELLE_DEV_BIN}/ruff" format src/animation tests
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    VERBATIM)

add_custom_target(lint
    COMMAND "${CHAZELLE_DEV_BIN}/clang-format" --dry-run --Werror ${CHAZELLE_CPP_FORMAT_FILES}
    COMMAND "${CHAZELLE_DEV_BIN}/ruff" check src/animation tests
    COMMAND "${CHAZELLE_DEV_BIN}/ruff" format --check src/animation tests
    COMMAND "${CHAZELLE_DEV_BIN}/run-clang-tidy.py"
            -p "${PROJECT_BINARY_DIR}"
            -clang-tidy-binary "${CHAZELLE_DEV_BIN}/clang-tidy"
            -j 4 -hide-progress -quiet
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    VERBATIM)
