# Unit tests (doctest). doctest's own CMake project declares a pre 3.5 minimum
# that CMake 4 rejects, so only its sources are fetched and the single header
# is exposed through an interface target.

include(FetchContent)
FetchContent_Declare(doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR do-not-add)
FetchContent_MakeAvailable(doctest)

add_library(hiccpp_doctest INTERFACE)
target_include_directories(hiccpp_doctest SYSTEM INTERFACE "${doctest_SOURCE_DIR}")

add_executable(hiccpp_tests
    test_main.cpp
    test_reader.cpp
    test_writer.cpp)
target_link_libraries(hiccpp_tests PRIVATE hiccpp::hiccpp hiccpp_doctest)
# Tests read the tracked data files and write into the build tree's scratch
# directory; both paths are passed in, never derived from the working directory.
target_compile_definitions(hiccpp_tests PRIVATE
    HICCPP_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
    HICCPP_TEST_SCRATCH="${CMAKE_CURRENT_BINARY_DIR}/scratch")

add_test(NAME hiccpp_unit_tests COMMAND hiccpp_tests)

# A downstream project consumes hiccpp through find_package on an installed
# tree and through FetchContent on the source tree.
add_test(NAME hiccpp_consumer
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/consumer
        -DPREFIX_PATH=${CMAKE_PREFIX_PATH}
        -DGENERATOR=${CMAKE_GENERATOR}
        -DEXPECTED_VERSION=${hiccpp_VERSION_MAJOR}.${hiccpp_VERSION_MINOR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/consumer/run_consumer.cmake)
set_tests_properties(hiccpp_consumer PROPERTIES TIMEOUT 900)
