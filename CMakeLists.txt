cmake_minimum_required(VERSION 3.14)
project(CHIC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

add_library(chic STATIC src/chic.cpp)
target_include_directories(chic PUBLIC include)
