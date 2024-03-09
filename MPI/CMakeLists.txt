# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(MPIexample
  DESCRIPTION "Examples of the usage of MPI"
  LANGUAGES Fortran)
set(VERSION 1.0)
enable_language (Fortran)

find_package(MPI)
if(NOT MPI_Fortran_FOUND)
  message(FATAL_ERROR "Could not find Fortran MPI.  Please set MPI_Fortran_COMPILER to point to the mpifort wrapper.")
endif()

include_directories(${MPI_Fortran_INCLUDE_PATH})
add_compile_options(${MPI_Fortran_COMPILE_FLAGS})
  

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
endif()


set(CMAKE_SKIP_BUILD_RPATH            FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH    TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(TARGET main)
#set(SOURCE_FILES src/main.f90)
#add_executable(${TARGET} ${SOURCE_FILES})
#target_link_libraries(main ${MPI_Fortran_LIBRARIES})
#set_property(TARGET main APPEND_STRING PROPERTY LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

set(SOURCE_FILES src/dotproduct.f90)
add_executable(${TARGET} ${SOURCE_FILES})
target_link_libraries(main ${MPI_Fortran_LIBRARIES})
set_property(TARGET main APPEND_STRING PROPERTY LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

