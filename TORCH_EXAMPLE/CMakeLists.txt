# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(TorchExample
  DESCRIPTION "Construction of Fourier Interpolant in 1D"
  LANGUAGES Fortran)
set(VERSION 1.0)
enable_language(Fortran)
set(CMake_Fortran_Compiler "gfortran*")

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE RELEASE)
endif()

set(FTorch_DIR /home/nico/Desktop/PHD-COURSES/librariesFOR/FTorchLIB/lib/cmake/FTorch)
find_package(FTorch)


set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -I/usr/include")
set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -std=f2008 ")
set(CMAKE_Fortran_FLAGS_DEBUG   "-Wall -O0 -g3 -fbounds-check")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")


#option(BUILD_FFT_MODULE "Build the FFT module" OFF)
#option(BUILD_OPT1 "--Build the first option" OFF)


#Set output for apps and libs:
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/apps)

#Test FTorch library:
set(TARGET testTorch.exe)
set(SOURCE_FILES src/apps/test.f90)
add_executable(${TARGET} ${SOURCE_FILES})
target_link_libraries(${TARGET} PRIVATE FTorch::ftorch)
message(STATUS "Building with Fortran PyTorch coupling")
