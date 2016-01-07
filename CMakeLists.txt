cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH ${Cantera_PREFIX} ${CMAKE_PREFIX_PATH})

#find dependencies
find_package(Trilinos REQUIRED)
# Make sure to use same compilers and flags as Trilinos
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

include(CheckCXXSourceRuns)
include(CheckCXXSourceCompiles)
include(CheckCXXSymbolExists)
include(CheckIncludeFile)
include(CMakeForceCompiler)

#cmake_force_c_compiler(${Trilinos_C_COMPILER} GNU)
#cmake_force_cxx_compiler(${Trilinos_CXX_COMPILER} GNU)

project(Cantrilbat C CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS -std=c++11 ${CMAKE_CXX_FLAGS})


include(options)

find_package(Cantera  REQUIRED)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Cantera_INCLUDE_DIRS})

add_subdirectory(util_src)
#add_subdirectory(ctpadds)

add_subdirectory(Electrode)
add_subdirectory(1Delectrode)


