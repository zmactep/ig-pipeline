project(misc C CXX)

include_directories(/opt/local/include)
include_directories(../../../../common_lib/include)
include_directories(../../include/misc)

add_library(misc STATIC
            config_reader.cpp)

