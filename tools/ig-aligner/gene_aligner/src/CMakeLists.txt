cmake_minimum_required(VERSION 2.6)

project(gene_aligner)
#find_package(TBB REQUIRED)
#find_package(LOG4CXX REQUIRED)
#find_package(gtest REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ${OpenMP_CXX_FLAGS} -std=c++11")

set(CMAKE_CURRENT_BINARY_DIR "../bin")
set(CMAKE_INSTALL_PREFIX "./")

add_subdirectory(data_structures)
add_subdirectory(misc)
add_subdirectory(tasks)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(../../../common_lib/include)
include_directories(../../../common_lib/include/io)
include_directories(../../../common_lib/include/utils)
include_directories(../../../common_lib/include/ssw)
include_directories(../../../common_lib/include/algorithm)
include_directories(../include/misc)
include_directories(../include/data_structures)
include_directories(../include/tasks)

link_directories(../../../common_lib/lib)
link_directories(${CMAKE_SOURCE_DIR}/misc)
link_directories(${CMAKE_SOURCE_DIR}/tasks)
link_directories(${CMAKE_SOURCE_DIR}/data_structures)
link_directories(${Boost_LIBRARY_DIRS})

add_executable(gene_aligner
               main.cpp)
target_link_libraries(gene_aligner utils io tasks algorithm ssw misc io data_structures tbb z log4cxx ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})

install(TARGETS gene_aligner DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(DIRECTORY "./config"
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
        FILES_MATCHING 
        PATTERN "*.ini"
        PATTERN "*.template"
        PATTERN "*.properties")
install(DIRECTORY "./test_data"
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
        FILES_MATCHING 
        PATTERN "*.fasta")