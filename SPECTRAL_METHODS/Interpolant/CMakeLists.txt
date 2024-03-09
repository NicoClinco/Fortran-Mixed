# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(SpectralDerivative
  DESCRIPTION "Construction of Fourier Interpolant in 1D"
  LANGUAGES Fortran)
set(VERSION 1.0)
enable_language(Fortran)
set(CMake_Fortran_Compiler "gfortran*")

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE RELEASE)
endif()


set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -I/usr/include -I/home/nico/Desktop/PHD-COURSES/librariesFOR/FFTPack/include")
set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -std=f2003 ")
set(CMAKE_Fortran_FLAGS_DEBUG   "-Wall -O0 -g3 -fbounds-check")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")


option(BUILD_FFT_MODULE "Build the FFT module" OFF)
option(BUILD_FOURIER_MODULE "Build the FOURIER module" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/apps)

if(BUILD_FFT_MODULE)
  message(STATUS "Building the fft module")
  add_library(fft_module ${CMAKE_SOURCE_DIR}/src/modules/FFT_MODULE.f90)
endif()

if(BUILD_FOURIER_MODULE)
  message(STATUS "Building the Fourier module")
  add_library(fourier_module ${CMAKE_SOURCE_DIR}/src/modules/Interp1DClass.f90)
  target_link_directories(fourier_module PUBLIC "/home/nico/Desktop/PHD-COURSES/librariesFOR/FFTPack/lib")
  target_link_libraries(fourier_module fftpack)
endif()


set(TARGET main.exe)
set(SOURCE_FILES src/apps/main.f90)


add_executable(${TARGET} ${SOURCE_FILES})
target_link_directories(${TARGET} PUBLIC "/home/nico/Desktop/PHD-COURSES/librariesFOR/FFTPack/lib"
  PUBLIC ${CMAKE_SOURCE_DIR}/modules
  )
target_link_libraries(${TARGET} fftpack fourier_module)
