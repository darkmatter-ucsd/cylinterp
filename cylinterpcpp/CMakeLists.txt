cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(CylinterCPP)

find_package(OpenMP)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/home/daqtest/Processor/boost_1_81_0")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/home/daqtest/Processor/boost_1_81_0/libs")

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

add_executable(test test.cc ${sources})
target_include_directories(test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(s2_sim s2_simulator_sandix.cc ${sources})
target_include_directories(s2_sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc)
if(OpenMP_CXX_FOUND)
    target_link_libraries(s2_sim PUBLIC OpenMP::OpenMP_CXX)
endif()