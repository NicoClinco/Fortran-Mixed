# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(SpectralDerivative
  DESCRIPTION "Fourier-Galerkin spectral scheme"
  LANGUAGES Fortran)
set(VERSION 1.0)
enable_language(Fortran)
set(CMake_Fortran_COMPILER "gfortran")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
  set(OPTS "-O3 -cpp")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Compiling in debug-mode")
  set(DB_OPTS "-g -Wall -fbounds-check -fbacktrace -cpp -DDEBUG")
  set(OPTS ${DB_OPTS})
  set(MOD_DEBUG MOD_DEBUG/true_debug.f90)
endif()

set(CMAKE_Fortran_FLAGS ${OPTS})


#DEPENDENCIES: My personal computer:
set(FFTPACK_INCLUDE /home/nico/Desktop/PHD-COURSES/librariesFOR/FFTPack/include)
set(FFTPACK_LIBS /home/nico/Desktop/PHD-COURSES/librariesFOR/FFTPack/lib)
set(FFTWRAPPER_INCLUDE /home/nico/Desktop/PHD-COURSES/librariesFOR/fft_wrapper/include)
set(FFTWRAPPER_LIBS /home/nico/Desktop/PHD-COURSES/librariesFOR/fft_wrapper/libs)
set(SYS_LIBS /usr/lib/x86_64-linux-gnu)

set(LAPACK_LIBS /usr/lib/x86_64-linux-gnu)
set(SCIFORTRAN_LIBS /home/nico/Desktop/PHD-COURSES/librariesFOR/SciFORTRAN/Scifor/build)
set(SCIFORTRAN_INCLUDE /home/nico/Desktop/PHD-COURSES/librariesFOR/SciFORTRAN/Scifor/build/include)
set(NLSOLVER_LIBS /home/nico/Desktop/PHD-COURSES/librariesFOR/nlesolver-fortran/build/gfortran_63E63B2712D76620/nlesolver-fortran)
set(NLSOLVER_INCLUDE /home/nico/Desktop/PHD-COURSES/librariesFOR/nlesolver-fortran/build/gfortran_2654F75F5833692A)

set(LIB_PATHS ${FFTPACK_LIBS} ${FFTWRAPPER_LIBS} ${LAPACK_LIBS} ${SCIFORTRAN_LIBS})
set(INCLUDE_PATHS ${FFTPACK_INCLUDE} ${FFTWRAPPER_INCLUDE} ${SCIFORTRAN_INCLUDE})


set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -I/usr/include")
set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -std=f2008 ")


#option(BUILD_FFT_MODULE "Build the FFT module" OFF)
option(BUILD_FOURIER_MODULE_1D "--Build the FOURIER module 1D--" OFF)
option(BUILD_FOURIER_MODULE_2D "--Build the FOURIER module 2D--" OFF)

#Set output for apps and libs:
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/apps)

#if(BUILD_FFT_MODULE)
#  message(STATUS "Building the fft module")
#  add_library(fft_module ${CMAKE_SOURCE_DIR}/src/modules/FFT_MODULE.f90)
#endif()
###############################################################################################
if(BUILD_FOURIER_MODULE_1D)
  message(STATUS "Building the Fourier module 1d")
  add_library(fourier_module1d SHARED ${CMAKE_SOURCE_DIR}/src/modules/Interp1DClass.f90)
  target_include_directories(fourier_module1d PUBLIC ${INCLUDE_PATHS})
  target_link_directories(fourier_module1d PUBLIC ${LIB_PATHS})
  target_link_libraries(fourier_module1d fftpack scifor)
endif()
################################################################################################
if(BUILD_FOURIER_MODULE_2D)
  message(STATUS "Building the Fourier module 2d")
  add_library(fourier_module2d
    SHARED ${CMAKE_SOURCE_DIR}/src/modules/Interp2DClass.f90)
  target_include_directories(fourier_module2d
    PUBLIC ${INCLUDE_PATHS}
    PUBLIC ${CMAKE_SOURCE_DIR}/modules)
  target_link_directories(fourier_module2d
    PUBLIC ${LIB_PATHS}
    PUBLIC ${CMAKE_SOURCE_DIR}/modules)
  target_link_libraries(fourier_module2d
    fftw3
    fftw3_wrapper
    fftpack
    scifor
    fourier_module1d)
endif()
###############################################################################################
#Test 1d fourier interpolant
#set(TARGET test1d.exe)
#set(SOURCE_FILES src/apps/test1d.f90)
#add_executable(${TARGET} ${SOURCE_FILES})
#target_include_directories(${TARGET}
#  PUBLIC ${INCLUDE_PATHS}
#  PUBLIC ${CMAKE_SOURCE_DIR}/modules)
#target_link_directories(${TARGET}
#  PUBLIC ${LIB_PATHS}
#  PUBLIC ${CMAKE_SOURCE_DIR}/modules)
#target_link_libraries(${TARGET} fftw3 fftw3_wrapper fftpack scifor fourier_module1d)

#TEST-WRAPPER:
set(TARGET testWrapper.exe)
set(SOURCE_FILES src/apps/testWrapper.f90)
add_executable(${TARGET} ${SOURCE_FILES})
target_include_directories(${TARGET}
  PUBLIC ${INCLUDE_PATHS})
target_link_directories(${TARGET}
  PUBLIC ${LIB_PATHS}
  PUBLIC ${CMAKE_SOURCE_DIR}/modules)
target_link_libraries(${TARGET}
  fftw3_wrapper
  fftpack
  scifor
  fourier_module1d
  fourier_module2d)



#General-Solver:
#set(TARGET test1Dsolver.exe)
#set(SOURCE_FILES src/apps/test_solver.f90)
#add_executable(${TARGET} ${SOURCE_FILES})
#target_include_directories(${TARGET} PUBLIC ${FFTPACK_INCLUDE} PUBLIC ${FFTWRAPPER_INCLUDE})
#target_link_directories(${TARGET}
#  PUBLIC ${FFTPACK_LIBS}
#  PUBLIC ${FFTWRAPPER_LIBS}
#  PUBLIC ${CMAKE_SOURCE_DIR}/modules)
#target_link_libraries(${TARGET} fftw3 fftw3_wrapper fftpack fourier_module2d)


#Pure spectral solver:
#set(TARGET testSpectralsolver.exe)
#set(SOURCE_FILES src/apps/testSpectralSolver.f90)
#add_executable(${TARGET} ${SOURCE_FILES})
#target_include_directories(${TARGET} PUBLIC ${INCLUDE_PATHS})
#target_link_directories(${TARGET}
#  PUBLIC ${LIB_PATHS}
#  PUBLIC ${CMAKE_SOURCE_DIR}/modules
#  )
#target_link_libraries(${TARGET} fftw3 fftw3_wrapper fftpack scifor)

#IMPLICIT-SOLVER
add_executable(ImplicitSolver
  src/modules/FourierVar1d.f90
  src/modules/FourierUtils.f90
  src/apps/testImplicitSolver.f90)
target_include_directories(ImplicitSolver
  PUBLIC ${INCLUDE_PATHS}
  PUBLIC ${NLSOLVER_INCLUDE}
  PUBLIC ${CMAKE_SOURCE_DIR}/modules)
 
target_link_directories(ImplicitSolver
  PUBLIC ${LIB_PATHS}
  PUBLIC ${SYS_LIBS}
  PUBLIC ${NLSOLVER_LIBS}
  PUBLIC ${CMAKE_SOURCE_DIR}/modules)
target_link_libraries(ImplicitSolver
  fftw3
  fftw3_wrapper
  fftpack
  -Wl,--start-group
  lapacke
  lapack
  blas
  nlesolver-fortran
  nlwrapper_module
  -Wl,--end-group
  scifor
  )
