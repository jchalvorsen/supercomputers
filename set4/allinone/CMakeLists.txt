project(hellompi)
cmake_minimum_required(VERSION 2.6)

enable_language(C)


find_package(MPI REQUIRED)

include_directories(${MPI_C_INCLUDE_PATH})

# To use pow function
link_libraries(m)

add_executable(sum_mpi sum_mpi.c)
target_link_libraries(sum_mpi ${MPI_C_LIBRARIES})

