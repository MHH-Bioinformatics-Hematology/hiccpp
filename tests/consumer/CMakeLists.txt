cmake_minimum_required(VERSION 3.21)
project(hicfilecpp_consumer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_RPATH_USE_LINK_PATH TRUE)

if(CONSUME STREQUAL "find_package")
    find_package(hicfilecpp 0.3 REQUIRED)
elseif(CONSUME STREQUAL "FetchContent")
    include(FetchContent)
    FetchContent_Declare(hicfilecpp SOURCE_DIR "${HICFILECPP_SOURCE}")
    FetchContent_MakeAvailable(hicfilecpp)
else()
    message(FATAL_ERROR "CONSUME must be find_package or FetchContent")
endif()

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE hicfilecpp::hicfilecpp)
