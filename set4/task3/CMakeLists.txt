project(hellompi)
cmake_minimum_required(VERSION 2.6)

enable_language(C)

find_package(MPI REQUIRED)

include_directories(${MPI_C_INCLUDE_PATH})

add_executable(hellompi hellompi.c)
target_link_libraries(hellompi ${MPI_C_LIBRARIES})

add_executable(integrate integrate.c)
target_link_libraries(integrate ${MPI_C_LIBRARIES})

