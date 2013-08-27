project(model C CXX)

include_directories(../../include/model)
include_directories(../../include/misc)
include_directories(../../include/utils)

add_library(model STATIC
            alignment.cpp
            nomenclature.cpp)

