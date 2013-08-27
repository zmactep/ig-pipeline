project(utils C CXX)

include_directories(../../include/utils)
include_directories(../../include/misc)
include_directories(../../include/data_structures)

add_library(utils STATIC
            amino_acid_converter.cpp
            reverse_complement.cpp
            kmer_generator.cpp)

