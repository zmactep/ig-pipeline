cmake_minimum_required(VERSION 2.6)

project(io C CXX)

include_directories(../../include)
include_directories(../../include/io)

add_library(io STATIC
            fasta_reader.cpp)
