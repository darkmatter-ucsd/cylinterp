cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(CylinterCPP)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/home/daqtest/Processor/boost_1_81_0")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/home/daqtest/Processor/boost_1_81_0/libs")

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

add_executable(test test.cc)
target_include_directories(test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc)