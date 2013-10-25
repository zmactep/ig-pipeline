project(algorithm C CXX)

include_directories(../../include)
include_directories(../../include/utils)
include_directories(../../include/algorithm)

add_library(algorithm STATIC
            aho_corasick.cpp)
