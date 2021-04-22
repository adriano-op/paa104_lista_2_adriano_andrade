cmake_minimum_required(VERSION 3.19)
project(semana2)


add_subdirectory (matplotplusplus)

add_executable(semana2 main.cpp)

target_link_libraries (semana2 PUBLIC matplot)

set(CMAKE_CXX_STANDARD 14)