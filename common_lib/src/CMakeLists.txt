cmake_minimum_required(VERSION 2.6)

project(common_lib)
set(LIB_DIR "../lib")
set(BUILD_DIR "../build")

add_subdirectory(io)

include_directories(../include)
include_directories(../include/io)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/io/libio.a
        DESTINATION ${CMAKE_HOME_DIRECTORY}/${LIB_DIR})