# Builds the downstream project twice: once against an installed hicfilecpp
# (find_package) and once pulling hicfilecpp in with FetchContent. Both builds
# run the resulting program, which reads a .hic file through the library.

cmake_minimum_required(VERSION 3.21)

foreach(var SOURCE_DIR WORK_DIR GENERATOR EXPECTED_VERSION)
    if(NOT DEFINED ${var})
        message(FATAL_ERROR "${var} is required")
    endif()
endforeach()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

function(run)
    execute_process(COMMAND ${ARGV} RESULT_VARIABLE status)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "command failed (${status}): ${ARGV}")
    endif()
endfunction()

set(data "${SOURCE_DIR}/tests/data/SRR1791297_30.juicer_tools_1.22.01.v8.hic")

# Prefix paths travel through the environment: a ;-list on a -D argument
# would be split into separate arguments by execute_process.
string(REPLACE ";" ":" prefix_env "${PREFIX_PATH}")

# 1. Install hicfilecpp into a prefix and consume it with find_package.
set(ENV{CMAKE_PREFIX_PATH} "${prefix_env}")
run(${CMAKE_COMMAND} -S "${SOURCE_DIR}" -B "${WORK_DIR}/lib-build" -G "${GENERATOR}"
    -DCMAKE_BUILD_TYPE=Release
    -DHICFILECPP_BUILD_TESTS=OFF -DHICFILECPP_BUILD_HARNESS=OFF
    "-DCMAKE_INSTALL_PREFIX=${WORK_DIR}/prefix")
run(${CMAKE_COMMAND} --build "${WORK_DIR}/lib-build" --parallel 8)
run(${CMAKE_COMMAND} --install "${WORK_DIR}/lib-build")

set(ENV{CMAKE_PREFIX_PATH} "${WORK_DIR}/prefix:${prefix_env}")
run(${CMAKE_COMMAND} -S "${SOURCE_DIR}/tests/consumer" -B "${WORK_DIR}/find-build"
    -G "${GENERATOR}" -DCMAKE_BUILD_TYPE=Release -DCONSUME=find_package
    "-DEXPECTED_VERSION=${EXPECTED_VERSION}")
run(${CMAKE_COMMAND} --build "${WORK_DIR}/find-build" --parallel 8)
run("${WORK_DIR}/find-build/consumer" "${data}")

# 2. FetchContent from the source tree.
set(ENV{CMAKE_PREFIX_PATH} "${prefix_env}")
run(${CMAKE_COMMAND} -S "${SOURCE_DIR}/tests/consumer" -B "${WORK_DIR}/fetch-build"
    -G "${GENERATOR}" -DCMAKE_BUILD_TYPE=Release -DCONSUME=FetchContent
    "-DHICFILECPP_SOURCE=${SOURCE_DIR}")
run(${CMAKE_COMMAND} --build "${WORK_DIR}/fetch-build" --parallel 8)
run("${WORK_DIR}/fetch-build/consumer" "${data}")
