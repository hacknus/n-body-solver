cmake_minimum_required(VERSION 3.17)
project(nBody)

set(CMAKE_CXX_STANDARD 14)
add_executable(nBody main.cpp body.cpp body.h io.cpp io.h math_utils.cpp math_utils.h)

find_package(MPI)
include_directories(/usr/local/include)
target_compile_options(nBody PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(nBody ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})