cmake_minimum_required(VERSION 2.8)

project(svm_data_generator)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CURRENT_BINARY_DIR "../bin")
set(CMAKE_INSTALL_PREFIX "./")

include_directories(../include)
include_directories(../../common_lib/include)
include_directories(../../common_lib/include/io)
include_directories(../../common_lib/include/utils)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../common_lib/lib)

add_executable(svm_data_generator
  main.cpp)
target_link_libraries(svm_data_generator io utils z)

install(TARGETS svm_data_generator DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
        
