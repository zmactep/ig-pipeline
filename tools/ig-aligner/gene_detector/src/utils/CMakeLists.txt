project(utils C CXX)

include_directories(../../include/utils)

add_library(utils STATIC
            tokenizer.cpp)

