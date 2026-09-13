# The C++ side of the equivalence harness. The Python side (oracle.py,
# compare.py, run.py) drives it; see the docstring of run.py.

add_executable(hicfilecpp-harness driver.cpp)
target_link_libraries(hicfilecpp-harness PRIVATE hicfilecpp::hicfilecpp)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(hicfilecpp-harness PRIVATE -Wall -Wextra)
endif()

# Runs a command and reports its own CPU time and peak RSS (see measure.cpp).
add_executable(hicfilecpp-measure measure.cpp)
