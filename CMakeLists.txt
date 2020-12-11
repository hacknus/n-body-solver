cmake_minimum_required(VERSION 3.17)
project(nBody)

set(CMAKE_CXX_STANDARD 14)

add_executable(nBody main.cpp body.cpp body.h io.cpp io.h math_utils.cpp math_utils.h)