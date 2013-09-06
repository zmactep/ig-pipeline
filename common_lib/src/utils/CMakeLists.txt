cmake_minimum_required(VERSION 2.6)
project(utils C CXX)

include_directories(../../include/io)
include_directories(../../include/utils)

add_library(utils STATIC
            tokenizer.cpp
            kmer_generator.cpp)

