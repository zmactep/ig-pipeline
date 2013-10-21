project(model C CXX)

include_directories(../../include)
include_directories(../../include/utils)
include_directories(../../include/model)

add_library(model STATIC
            alignment.cpp
            nomenclature.cpp)

