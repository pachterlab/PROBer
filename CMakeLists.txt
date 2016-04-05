cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(PROBer CXX)

add_compile_options(-std=c++11)

find_package(ZLIB REQUIRED)
if (NOT ZLIB_FOUND)
   	message(FATA_ERROR "zlib not found. Required for the SAMtools.")
endif()

add_subdirectory(ext)
add_subdirectory(src)

