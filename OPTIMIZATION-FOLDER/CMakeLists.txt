# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(SpectralDerivative
  DESCRIPTION "Construction of Fourier Interpolant in 1D"
  LANGUAGES Fortran)
set(VERSION 1.0)
enable_language(Fortran)
set(CMake_Fortran_Compiler "gfortran")

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE RELEASE)
endif()

#DEPENDENCIES#
set(FFTPACK_INCLUDE /home/nico/Desktop/PHD-COURSES/librariesFOR/FFTPack/include)
set(FFTPACK_LIBS /home/nico/Desktop/PHD-COURSES/librariesFOR/FFTPack/lib)
set(FFTWRAPPER_INCLUDE /home/nico/Desktop/PHD-COURSES/librariesFOR/fft_wrapper/include)
set(FFTWRAPPER_LIBS /home/nico/Desktop/PHD-COURSES/librariesFOR/fft_wrapper/libs)
set(SYS_LIBS /usr/lib/x86_64-linux-gnu)
set(SCIFORTRAN_LIBS /home/nico/Desktop/PHD-COURSES/librariesFOR/SciFORTRAN/Scifor/build)
set(SCIFORTRAN_INCLUDE /home/nico/Desktop/PHD-COURSES/librariesFOR/SciFORTRAN/Scifor/build/include)
set(NLSOLVER_LIBS /home/nico/Desktop/PHD-COURSES/librariesFOR/nlesolver-fortran/build/gfortran_63E63B2712D76620/nlesolver-fortran)
set(NLSOLVER_INCLUDE /home/nico/Desktop/PHD-COURSES/librariesFOR/nlesolver-fortran/build/gfortran_2654F75F5833692A)

set(LIB_PATHS ${FFTPACK_LIBS} ${FFTWRAPPER_LIBS} ${LAPACK_LIBS} ${SCIFORTRAN_LIBS})
set(INCLUDE_PATHS ${FFTPACK_INCLUDE} ${FFTWRAPPER_INCLUDE} ${SCIFORTRAN_INCLUDE})


set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -I/usr/include")
set(CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -std=f2008")
set(CMAKE_Fortran_FLAGS_DEBUG   "-Wall -O0 -g -fbounds-check -fcheck=all")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fno-stack-arrays")


#option(BUILD_FFT_MODULE "Build the FFT module" OFF)
option(BUILD_MODULE "--Build the generic module--" OFF)


#Set output for apps and libs:
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/apps)


if(BUILD_MODULE)
  message(STATUS "--->Building the Generic Module<---")
  add_library(nlwrapper_module ${CMAKE_SOURCE_DIR}/src/modules/NLsolverWrapper.f90)
  target_include_directories(nlwrapper_module PUBLIC ${INCLUDE_PATHS} PUBLIC ${NLSOLVER_INCLUDE})
  target_link_directories(nlwrapper_module
    PUBLIC ${LIB_PATHS}
    PUBLIC ${LIB_PATHS}
    PUBLIC ${NLSOLVER_LIBS}
    PUBLIC ${SYS_LIBS}
    )
  target_link_libraries(nlwrapper_module
    -Wl,--start-group
    gfortran lapacke lapack blas nlesolver-fortran
    -Wl,--end-group)
endif()


#Add also a non linear solver:
set(TestNLsolver testNLsolver.exe)
set(SOURCE_FILES src/apps/testNLsolver.f90)
add_executable(${TestNLsolver} ${SOURCE_FILES})
target_include_directories(${TestNLsolver} PUBLIC ${INCLUDE_PATHS} PUBLIC ${NLSOLVER_INCLUDE})
target_link_directories(${TestNLsolver}
  PUBLIC ${LIB_PATHS}
  PUBLIC ${NLSOLVER_LIBS}
  PUBLIC ${CMAKE_SOURCE_DIR}/modules
  PUBLIC ${SYS_LIBS}
  )
#target_link_libraries(${TestNLsolver} gfortran lapacke lapack blas gfortran nlesolver-fortran)
target_link_libraries(${TestNLsolver} -Wl,--start-group gfortran lapacke lapack blas nlesolver-fortran -Wl,--end-group nlwrapper_module)
