
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
project(aero)
include_directories(matrix)

add_library(${PROJECT_NAME} STATIC aero.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} PUBLIC matrix)