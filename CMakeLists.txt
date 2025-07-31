cmake_minimum_required(VERSION 3.14)

project(polynomial-factorization)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/tools/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconstexpr-steps=100000000")

add_library(factorization INTERFACE)
target_include_directories(factorization INTERFACE .)

add_subdirectory(tests)
add_subdirectory(experiments)