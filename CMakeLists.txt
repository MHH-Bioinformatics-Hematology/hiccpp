cmake_minimum_required(VERSION 3.21)

project(hicfilecpp
    VERSION 0.2.0
    DESCRIPTION "C++ reader and writer of the Juicer .hic format, API-compatible with hicstraw"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(HICFILECPP_BUILD_TESTS "Build the doctest unit tests" ${PROJECT_IS_TOP_LEVEL})
option(HICFILECPP_BUILD_HARNESS "Build the equivalence harness driver" ${PROJECT_IS_TOP_LEVEL})

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

add_library(hicfilecpp
    src/codec.cpp
    src/io.cpp
    src/norms.cpp
    src/output.cpp
    src/reader.cpp
    src/writer.cpp)
add_library(hicfilecpp::hicfilecpp ALIAS hicfilecpp)

target_compile_features(hicfilecpp PUBLIC cxx_std_20)
target_include_directories(hicfilecpp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
# zlib is an implementation detail: no public header includes zlib.h.
target_link_libraries(hicfilecpp PRIVATE ZLIB::ZLIB Threads::Threads)
set_target_properties(hicfilecpp PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(hicfilecpp PRIVATE -Wall -Wextra -Wpedantic)
endif()

configure_file(include/hicfilecpp/version.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/hicfilecpp/version.hpp @ONLY)
target_include_directories(hicfilecpp PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>)

# ---- install and package config ----
install(TARGETS hicfilecpp EXPORT hicfilecppTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/hicfilecpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/hicfilecpp/version.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hicfilecpp)
install(EXPORT hicfilecppTargets
    NAMESPACE hicfilecpp::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hicfilecpp)
configure_package_config_file(cmake/hicfilecppConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/hicfilecppConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hicfilecpp)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/hicfilecppConfigVersion.cmake
    COMPATIBILITY SameMinorVersion)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/hicfilecppConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/hicfilecppConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hicfilecpp)

if(HICFILECPP_BUILD_HARNESS)
    add_subdirectory(harness)
endif()

if(HICFILECPP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
