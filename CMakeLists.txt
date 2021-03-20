cmake_minimum_required(VERSION 3.17)
project(nbody)

set(CMAKE_CXX_STANDARD 14)
add_executable(nbody src/main.cpp src/body.cpp inc/body.h src/io.cpp inc/io.h src/math_utils.cpp inc/math_utils.h)

find_package(MPI)
include_directories(/usr/local/include)
target_compile_options(nbody PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(nbody ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})