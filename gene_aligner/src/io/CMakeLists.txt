project(io C CXX)

include_directories(../../include/io)
include_directories(../../include/data_structures)
include_directories(../../include/misc)
include_directories(../../include/ssw)
include_directories(../../include/utils)
include_directories(/opt/local/include)

add_library(io STATIC
            fasta_reader.cpp
            output_formatter.cpp)
