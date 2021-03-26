cmake_minimum_required(VERSION 3.17)
project(nbody)
include_directories(inc)

set(CMAKE_CXX_STANDARD 17)
add_executable(nbody src/main.cpp src/body.cpp inc/body.h src/io.cpp inc/io.h src/math_utils.cpp inc/math_utils.h)

find_package(MPI)
set(MPI_CXX_COMPILE_FLAGS -O3 -ffast-math)
set(MPI_CXX_LINK_FLAGS -O3)
include_directories(/usr/local/include)
target_compile_options(nbody PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(nbody ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})