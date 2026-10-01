# The programs the documentation shows. They are built with the library so that
# an example that stops compiling fails the build rather than the reader.
foreach(example read_hic write_hic)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE hiccpp::hiccpp)
endforeach()
