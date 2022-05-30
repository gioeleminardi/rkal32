cmake_minimum_required(VERSION 3.22)
project(kal32)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

add_executable(kal32 src/main.cpp)
