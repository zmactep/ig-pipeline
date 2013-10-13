cmake_minimum_required(VERSION 2.6)

project(tests)

enable_testing()

cmake_policy(SET CMP0015 OLD)

find_package (Threads)

#gtest
include_directories(../../include/io)
include_directories(../../include/utils)
include_directories(../../include/third_party/gtest-1.6.0/include)

link_directories(../../lib)
link_directories(../../third_party/lib)

add_executable(io_test
  io.cpp)

add_executable(utils_test
  utils.cpp)
    
target_link_libraries(utils_test utils gtest gtest_main pthread)
target_link_libraries(io_test io utils gtest gtest_main pthread)

add_test(NAME io_test COMMAND io_test)
add_test(NAME utils_test COMMAND utils_test)

install(TARGETS io_test DESTINATION ${PROJECT_SOURCE_DIR}/../../bin)
install(TARGETS utils_test DESTINATION ${PROJECT_SOURCE_DIR}/../../bin)

install(DIRECTORY "./data"
        DESTINATION ${PROJECT_SOURCE_DIR}/../../bin
        FILES_MATCHING 
        PATTERN "*.txt")
