project(gridMod)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
cmake_minimum_required(VERSION 2.8.3)

#add_subdirectory(ext ext_build)

include_directories(
  # include files
  ${CMAKE_CURRENT_SOURCE_DIR}/../shared_headers
  # Eigen linear algebra library
  $ENV{EIGEN_INC}

)

message("variable is $ENV{EIGEN_INC}")

# the list of source files
set(SRCS
  src/gridMod_from_CT.cpp
)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(gridMod_from_CT ${SRCS})
target_link_libraries (gridMod_from_CT Eigen3::Eigen)

