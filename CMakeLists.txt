cmake_minimum_required(VERSION 2.8.4)
set(TARGET "cplex_lab")

project(${TARGET})

set (CMAKE_CXX_STANDARD 20)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Cplex)

set(SOURCES
    source/graph.h
    source/graph.cpp
    source/solver.h
    source/solver.cpp
    source/common.hpp
    source/heuristic.hpp
    source/model.hpp
    source/utils.hpp
    source/main.cpp
)

add_executable(${TARGET} ${SOURCES})

target_include_directories(${TARGET} PRIVATE ${CPLEX_INCLUDE_DIRS})
target_link_libraries(${TARGET} PRIVATE ${CPLEX_LIBRARIES})

add_custom_command(TARGET ${TARGET} PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/graphs/ ${CMAKE_CURRENT_BINARY_DIR}/graphs/)

add_custom_command(TARGET ${TARGET} PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/graphs/ $<TARGET_FILE_DIR:${TARGET}>/graphs/)

add_custom_command(TARGET ${TARGET} PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/tasks/ ${CMAKE_CURRENT_BINARY_DIR}/tasks/)

add_custom_command(TARGET ${TARGET} PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/tasks/ $<TARGET_FILE_DIR:${TARGET}>/tasks/)

add_custom_command(TARGET ${TARGET} PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${CPLEX_BIN_DIR}/cplex${CPLEX_WIN_VERSION}0.dll $<TARGET_FILE_DIR:${TARGET}>)
