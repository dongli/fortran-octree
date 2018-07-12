cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(fortran-octree VERSION 0.0.1 LANGUAGES Fortran)

find_package(OpenMP)
if (OPENMP_FOUND)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif ()

add_library(fortran_octree src/octree_mod.F90)

add_executable(octree_test.exe src/octree_test.F90)
target_link_libraries(octree_test.exe fortran_octree)