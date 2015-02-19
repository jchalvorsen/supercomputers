project(hellompi)
cmake_minimum_required(VERSION 2.6)

enable_language(C)


find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()
include_directories(${MPI_C_INCLUDE_PATH})

# To use pow function
link_libraries(m)

add_executable(sum sum.c)

add_executable(sum_omp sum_omp.c)

add_executable(sum_mpi sum_mpi.c)
target_link_libraries(sum_mpi ${MPI_C_LIBRARIES})

