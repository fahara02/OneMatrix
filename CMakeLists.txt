
cmake_minimum_required(VERSION 3.29.0) 
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_BINARY_DIR}/conan_toolchain.cmake")

set(CMAKE_CXX_STANDARD 20)
project(CVLM
        VERSION 0.0.1
        DESCRIPTION "The leading Hello World App"
        LANGUAGES CXX)  
set(BINARY_NAME ${PROJECT_NAME})  


list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(ClangFormat)

add_clang_format_target(format-code ${PROJECT_SOURCE_DIR})





add_executable(${BINARY_NAME} src/main.cpp   )

target_include_directories( ${BINARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src )


