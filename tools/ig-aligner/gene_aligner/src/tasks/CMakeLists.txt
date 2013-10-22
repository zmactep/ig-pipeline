project(tasks C CXX)

include_directories(../../../../common_lib/include)
include_directories(../../../../common_lib/include/algorithm)
include_directories(../../../../common_lib/include/io)
include_directories(../../../../common_lib/include/utils)
include_directories(../../../../common_lib/include/ssw)
include_directories(../../include/data_structures)
include_directories(../../include/tasks)
include_directories(../../include/misc)

link_directories(${CMAKE_SOURCE_DIR}/../../../common_lib/lib)
link_directories(${CMAKE_SOURCE_DIR}/misc)

add_library(tasks STATIC
            regions_finder.cpp
            task_configurator.cpp)
