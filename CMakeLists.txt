cmake_minimum_required(VERSION 3.21)

project(hiccpp
    VERSION 0.5.0
    DESCRIPTION "C++ reader and writer of the Juicer .hic format, API-compatible with hicstraw"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(HICCPP_BUILD_TESTS "Build the doctest unit tests" ${PROJECT_IS_TOP_LEVEL})
option(HICCPP_BUILD_EXAMPLES "Build the example programs the documentation shows"
    ${PROJECT_IS_TOP_LEVEL})
option(HICCPP_BUILD_HARNESS "Build the equivalence harness driver" ${PROJECT_IS_TOP_LEVEL})

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# The only native dependency: zlib, for the deflate-compressed blocks.
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

add_library(hiccpp
    src/codec.cpp
    src/io.cpp
    src/norms.cpp
    src/output.cpp
    src/reader.cpp
    src/writer.cpp)
add_library(hiccpp::hiccpp ALIAS hiccpp)

target_compile_features(hiccpp PUBLIC cxx_std_20)
target_include_directories(hiccpp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
# zlib is an implementation detail: no public header includes zlib.h.
target_link_libraries(hiccpp PRIVATE ZLIB::ZLIB Threads::Threads)
set_target_properties(hiccpp PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(hiccpp PRIVATE -Wall -Wextra -Wpedantic)
endif()

configure_file(include/hiccpp/version.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/hiccpp/version.hpp @ONLY)
target_include_directories(hiccpp PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>)

# ---- install and package config ----
install(TARGETS hiccpp EXPORT hiccppTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/hiccpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/hiccpp/version.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hiccpp)
install(EXPORT hiccppTargets
    NAMESPACE hiccpp::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hiccpp)
configure_package_config_file(cmake/hiccppConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/hiccppConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hiccpp)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/hiccppConfigVersion.cmake
    COMPATIBILITY SameMinorVersion)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/hiccppConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/hiccppConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hiccpp)

if(HICCPP_BUILD_HARNESS)
    add_subdirectory(harness)
endif()

if(HICCPP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(HICCPP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
