cmake_minimum_required(VERSION 2.8.5)

project(lib)
enable_language(Fortran)
set(CMAKE_BUILD_TYPE Release)
file(GLOB source "./src/*.f90")
add_library(qblas STATIC ${source})
