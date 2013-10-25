cmake_minimum_required(VERSION 2.6)

project(ig_combinator)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_CURRENT_BINARY_DIR "../bin")
set(CMAKE_INSTALL_PREFIX "./")

include_directories(${Boost_INCLUDE_DIRS})
include_directories(../../common_lib/include)
include_directories(../../common_lib/include/io)
include_directories(../../common_lib/include/utils)

link_directories(../../common_lib/lib)
link_directories(${Boost_LIBRARY_DIRS})

add_executable(ig_combinator
               main.cpp)

target_link_libraries(ig_combinator z io utils)

install(TARGETS ig_combinator DESTINATION ${CMAKE_CURRENT_BINARY_DIR})