cmake_minimum_required(VERSION 3.16)

# https://cmake.org/cmake/help/latest/command/project.html
project(grid-resampler VERSION 0.1.0 
    DESCRIPTION "Create grids of different resolution to the original"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

file(GLOB SOURCE_FILES "${SOURCE_DIR}/*.cxx")

configure_file("${INCLUDE_DIR}/appconfig.h.in" "${PROJECT_BINARY_DIR}/include/appconfig.h")

set(BOOST_DIR CACHE PATH "Boost installation directory path")
set(BOOST_PRG_OPTS_LIB CACHE FILEPATH "Boost::program_options static lib filepath")

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}/include" ${BOOST_DIR})
target_link_libraries(${PROJECT_NAME} "${BOOST_PRG_OPTS_LIB}")