cmake_minimum_required(VERSION 3.28)
project(DSA_Assignment)

set(CMAKE_CXX_STANDARD 11)

add_executable(DSA_Assignment
        src/kDTree.cpp
        src/kDTree.hpp
        src/Dataset.hpp
        src/Dataset.o
        src/main.hpp
        src/main.cpp
)

