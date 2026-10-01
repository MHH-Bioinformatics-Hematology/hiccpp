# Installing and building

## Conda

```
conda install -c bioconda hiccpp
```

## From source

The library needs a C++20 compiler, CMake 3.21 or newer and zlib. doctest is
fetched for the unit tests.

```
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/zlib/prefix
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build --prefix /opt/hiccpp
```

## Using it from another project

As an installed package:

```cmake
find_package(hiccpp 0.5 REQUIRED)
target_link_libraries(app PRIVATE hiccpp::hiccpp)
```

Or fetched at configure time:

```cmake
include(FetchContent)
FetchContent_Declare(hiccpp GIT_REPOSITORY <url> GIT_TAG <tag>)
FetchContent_MakeAvailable(hiccpp)
target_link_libraries(app PRIVATE hiccpp::hiccpp)
```

The `hiccpp_consumer` test builds such a project both ways, so both paths are
covered by the test suite.
