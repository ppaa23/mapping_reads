cmake_minimum_required(VERSION 3.16)
project(bioinf CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(bioinf src/main.cpp)
target_include_directories(bioinf PRIVATE include)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(bioinf OpenMP::OpenMP_CXX)
endif()
