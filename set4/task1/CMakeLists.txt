project(set4)
cmake_minimum_required(VERSION 2.6)

enable_language(C)

add_subdirectory(../common common)

add_definitions(${COMMON_DEFINITIONS})
list(APPEND CMAKE_C_FLAGS ${COMMON_C_FLAGS})
include_directories(${COMMON_INCLUDE_PATHS})

# This is a hack. Have no idea why this works
if(UNIX)
    link_libraries(m)
endif(UNIX)

add_executable(main main.c)
target_link_libraries(main ${COMMON_LIBRARIES})




