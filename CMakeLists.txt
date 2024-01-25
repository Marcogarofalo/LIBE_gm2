cmake_minimum_required(VERSION 3.13)

project(LIBE_gm2 CXX C)


IF(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all")
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment")
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -static-libasan -g   ")#-fsanitize=float-divide-by-zero
  set(LDFLAGS "${LDFLAGS}  -fsanitize=address")
  set(CMAKE_LDFLAGS_FLAGS  " ${CMAKE_LDFLAGS_FLAGS} ${LDFLAGS} ")
ENDIF()

#############################################
# load analysis_program
#############################################
find_package(analysis_program)
include_directories(${ANALYSIS_PROGRAM_INCLUDE_DIR})
MESSAGE("ANALYSIS_PROGRAM_INCLUDE_DIR   ${ANALYSIS_PROGRAM_INCLUDE_DIR}")
MESSAGE("ANALYSIS_PROGRAM_LIBRARIES: ${ANALYSIS_PROGRAM_LIBRARIES}")
MESSAGE("ANALYSIS_PROGRAM_LIBRARY_DIR: ${ANALYSIS_PROGRAM_LIBRARY_DIR}")
IF(NOT analysis_program_FOUND)
    MESSAGE("find_package faild to find yaml-cpp: you need to do it manually")
    MESSAGE("hint:  analysis_program_FOUND_DIR or CMAKE_PREFIX_PATH need to be global path or relative to the source dir")
    MESSAGE(FATAL_ERROR)
ENDIF()

find_package(Eigen3 REQUIRED) # we need eigen
MESSAGE("IF the analysis_program was compiled with OpenMP then you also need to load it")
MESSAGE("uncomment the lines after this message")
FIND_PACKAGE( OpenMP)
if(OPENMP_FOUND)
  message("OPENMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
include_directories(SYSTEM ${OpenMP_INCLUDE_PATH})
#############################################
# add a local library for this project
#############################################
add_library(
    functions_LIBE_gm2 STATIC
    functions_LIBE_gm2.hpp
    functions_LIBE_gm2.cpp
)
target_link_libraries(functions_LIBE_gm2 PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})
target_include_directories(functions_LIBE_gm2 PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(LIBE_gm2 LIBE_gm2.cpp)
target_link_libraries(LIBE_gm2 PUBLIC functions_LIBE_gm2)
IF(WITH_ARB)
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -DWITH_ARB   ")
  target_link_libraries(LIBE_gm2 PUBLIC flint-arb flint)
ENDIF()
target_link_libraries(LIBE_gm2 PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})
