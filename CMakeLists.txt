cmake_minimum_required(VERSION 3.5)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(enoise)

set (CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_FLAGS "-Wall -Weffc++ -Wshadow -Wconversion -std=c++17")
set (source_dir "${PROJECT_SOURCE_DIR}/src/")

file (GLOB source_files "${source_dir}/*")

add_executable(enoise "${source_files}")
