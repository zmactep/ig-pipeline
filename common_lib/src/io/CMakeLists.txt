cmake_minimum_required(VERSION 2.6)

project(io C CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(../../include)
include_directories(../../include/io)
include_directories(../../include/utils)

add_library(io STATIC
            fasta_reader.cpp
            kabat_reader.cpp
            nomenclature.cpp
            generic_nomenclature.cpp
            read.cpp)
