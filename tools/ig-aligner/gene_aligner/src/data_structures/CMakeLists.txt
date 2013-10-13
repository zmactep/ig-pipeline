project(data_structures C CXX)

include_directories(../../include/data_structures)
include_directories(../../include/misc)
include_directories(../../include/utils)
include_directories(../../include/io)
link_directories(${CMAKE_SOURCE_DIR}/utils)
add_library(data_structures STATIC
            database.cpp
            read.cpp)
target_link_libraries(data_structures utils)
