cmake_minimum_required(VERSION 2.6)

project(gene_detector)
find_package(Boost COMPONENTS program_options REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CURRENT_BINARY_DIR "../bin")
set(CMAKE_INSTALL_PREFIX "./")

add_subdirectory(misc)


include_directories(${Boost_INCLUDE_DIRS})
include_directories(/opt/local/include)
include_directories(../../../common_lib/include)
include_directories(../../../common_lib/include/utils)
include_directories(../../../common_lib/include/model)
include_directories(../include/misc)

link_directories(${Boost_LIBRARY_DIRS})
link_directories(${CMAKE_SOURCE_DIR}/misc)
link_directories(../../../common_lib/lib)

add_executable(gene_detector
               main.cpp)

target_link_libraries(gene_detector log4cxx misc model utils ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})

install(TARGETS gene_detector DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(DIRECTORY "./config"
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
        FILES_MATCHING 
        PATTERN "*.ini")
        
install(DIRECTORY "./test_data"
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
        FILES_MATCHING 
        PATTERN "*.txt")
