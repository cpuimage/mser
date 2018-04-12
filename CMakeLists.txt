cmake_minimum_required(VERSION 3.9)
project(mser)

set(CMAKE_CXX_STANDARD 11)

add_executable(mser
        mser.c)
target_link_libraries(mser -lm )