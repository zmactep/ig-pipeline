project(ssw C CXX)

include_directories(../../include/ssw)

add_library(ssw STATIC
            ssw.c
            ssw_cpp.cpp)
