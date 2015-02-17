project(hellompi)
cmake_minimum_required(VERSION 2.6)

add_subdirectory(../common ../common)

include_directories(${COMMON_INCLUDE_PATHS})
list(APPEND CMAKE_C_FLAGS ${COMMON_C_FLAGS})
add_definitions(${COMMON_DEFINITIONS})



if(UNIX)
link_libraries(m)
endif(UNIX)

if(MPI_FOUND)
  add_executable(vectorsum_mpi vectorsum_mpi.c)
  target_link_libraries(vectorsum_mpi ${MPI_C_LIBRARIES} ${COMMON_LIBRARIES})
endif()
