project(tasks C CXX)

include_directories(../../include/data_structures)
include_directories(../../include/algorithm)
include_directories(../../include/io)
include_directories(../../include/ssw)
include_directories(../../include/utils)
include_directories(../../include/tasks)
include_directories(../../include/misc)

link_directories(${CMAKE_SOURCE_DIR}/algorithm)
link_directories(${CMAKE_SOURCE_DIR}/misc)
link_directories(${CMAKE_SOURCE_DIR}/ssw)

add_library(tasks STATIC
            regions_finder.cpp
            task_configurator.cpp)
