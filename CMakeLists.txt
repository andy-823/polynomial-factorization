cmake_minimum_required(VERSION 3.14)

project(polynomial-factorization)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/tools/cmake")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconstexpr-ops-limit=1000000000")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconstexpr-steps=100000000")
endif()

add_library(factorization INTERFACE)
target_include_directories(factorization INTERFACE .)

add_subdirectory(tests)
add_subdirectory(experiments)