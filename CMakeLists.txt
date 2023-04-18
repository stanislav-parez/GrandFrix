cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)

project(GranLayer Fortran)

find_package(Fortran)
include(CheckFortranSourceCompiles)

# this is where we will place the Fortran module files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xHOST")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native -std=f2008")
endif()

add_executable(
    GranFrixrm
	src/GranLayer.f90
    src/principalstress.f90
	src/defect.f90
	src/updatePressure.f90
	src/rotate.f90
	src/vel_ver.f90
	src/breaks.f90
	src/smoothNonExtended.f90
	src/init.f90
	src/pssubs.f90
	src/mycommons.f90
	src/rw_conf.f90
	src/linklist.f90
	src/force.f90
	src/findforces.f90
	src/generate.f90
	src/postScriptTool.f90
	src/inputReader.f90
	src/outputWriter.f90
    )
