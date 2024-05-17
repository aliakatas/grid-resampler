cmake_minimum_required(VERSION 3.15)

# https://cmake.org/cmake/help/latest/command/project.html
project(example-project VERSION 0.1.0 
    DESCRIPTION "Example project to use as template"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

file(GLOB SOURCE_FILES "${SOURCE_DIR}/*.cxx")

configure_file("${INCLUDE_DIR}/appconfig.h.in" "${PROJECT_BINARY_DIR}/include/appconfig.h")

set(BOOST_DIR CACHE PATH "Boost installation directory path")
set(BOOST_PRG_OPTS_LIB CACHE FILEPATH "Boost::program_options static lib filepath")

add_executable(example-project ${SOURCE_FILES})

target_include_directories(example-project PUBLIC "${PROJECT_BINARY_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}/include" ${BOOST_DIR})
target_link_libraries(example-project "${BOOST_PRG_OPTS_LIB}")