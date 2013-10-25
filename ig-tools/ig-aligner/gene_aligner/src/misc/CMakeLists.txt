project(misc C CXX)

include_directories(/opt/local/include)
include_directories(../../include/misc)
include_directories(../../include/data_structures)

add_library(misc STATIC
            config_reader.cpp)

