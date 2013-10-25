cmake_minimum_required(VERSION 2.6)

project(common_lib)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(LIB_DIR "../lib")
set(BUILD_DIR "../build")


add_subdirectory(io)
add_subdirectory(utils)
add_subdirectory(tests)
add_subdirectory(algorithm)
add_subdirectory(ssw)
add_subdirectory(model)

include_directories(../include)
include_directories(../include/io)
include_directories(../include/utils)
include_directories(../include/algorithm)
include_directories(../include/ssw)
include_directories(../include/model)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/io/libio.a
        DESTINATION ${PROJECT_SOURCE_DIR}/${LIB_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utils/libutils.a
        DESTINATION ${PROJECT_SOURCE_DIR}/${LIB_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/algorithm/libalgorithm.a
        DESTINATION ${PROJECT_SOURCE_DIR}/${LIB_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ssw/libssw.a
        DESTINATION ${PROJECT_SOURCE_DIR}/${LIB_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/model/libmodel.a
        DESTINATION ${PROJECT_SOURCE_DIR}/${LIB_DIR})
