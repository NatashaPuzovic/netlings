cmake_minimum_required(VERSION 3.5)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project (netlings)

#set (CMAKE_BUILD_TYPE Debug)
set (CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Weffc++ -std=c++17")

add_executable(netlings
				 src/main.cpp
				 src/Netling.cpp
				 src/Netling.h
				 src/Population.cpp
				 src/Population.h
				 src/Params.cpp
				 src/Params.h)