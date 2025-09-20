cmake_minimum_required(VERSION 3.20)

set(CMAKE_PROJECT_NAME "ForAdjoint")
set(project_version "0.2.0")
project(${CMAKE_PROJECT_NAME}, VERSION ${project_version} LANGUAGES Fortran)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Find LAPACK library
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
else()
    message(FATAL_ERROR "LAPACK not found!")
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  add_definitions(-DUSE_GNU_EXTENSIONS)
  set(CMAKE_Fortran_FLAGS "-std=gnu -O3 -fimplicit-none -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused -finline-functions -Wno-argument-mismatch")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -ggdb -fbacktrace -fbounds-check")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "NVIDIA")
  set(CMAKE_Fortran_FLAGS "-std=f2008 -O3 -Mcpu=auto -implicitnone -Merror_limit=10 -Mdclchk -Ktrap=fp:invalid,zero,overflow -warn all -Minline")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -check all")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-std08 -O3 -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -diag-disable 6477 -implicitnone -gen-interfaces -warn all -warn nointerfaces")
  set(CMAKE_Fortran_FLAGS_DEBUG "-check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
  set(CMAKE_Fortran_FLAGS "-O3 -Onoaggress -Oipa0 -hfp2 -Ovector3 -Oscalar3 -Ocache2 -Ounroll2 -Ofusion2 -M 1193 -M 1438")
  set(CMAKE_Fortran_FLAGS_DEBUG "eC -eD -ec -en -eI -ea -g -G0")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang")
  set(CMAKE_Fortran_FLAGS "-O3 -fimplicit-none -ffpe-trap=invalid,zero,overflow -funroll-loops -Wall")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fbacktrace -fbounds-check")
endif()

# compile src files
file(GLOB utils_sources "src/utils/*.f90")
file(GLOB adjoint_sources "src/adjoint_source_types/*.f90")
file(GLOB mt_sources "src/mtspec/*.f90")
add_library(foradjoint STATIC ${utils_sources} ${adjoint_sources} ${mt_sources})

# Link LAPACK to the library
target_link_libraries(foradjoint ${LAPACK_LIBRARIES})

# add example executables
add_subdirectory(examples)