project(data_structures C CXX)

include_directories(../../../../common_lib/include)
include_directories(../../../../common_lib/include/io)
include_directories(../../../../common_lib/include/utils)
include_directories(../../../../common_lib/include/ssw)
include_directories(../../include/data_structures)
include_directories(../../include/misc)

link_directories(../../../../common_lib/lib)

add_library(data_structures STATIC
            database.cpp
            output_formatter.cpp)

target_link_libraries(data_structures utils)
