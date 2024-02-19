project(NonEuclideanSlicedOT)
cmake_minimum_required(VERSION 3.12)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP)

include(polyscope)
include(geometry-central)
include(eigen)
include(openmp)

add_subdirectory(deps/bvh)

file(GLOB SRC
     "src/**"
     "deps/CLI11.hpp"
)

message(STATUS ${SRC})

add_executable(spherical_mesh_sampling apps/mesh_sampling.cpp ${SRC})
target_link_libraries(spherical_mesh_sampling polyscope geometry-central  bvh OpenMP::OpenMP_CXX)

add_executable(spherical_bluenoise apps/spherical_bluenoise.cpp ${SRC})
target_link_libraries(spherical_bluenoise polyscope geometry-central  bvh OpenMP::OpenMP_CXX)

add_executable(hyperbolic_bluenoise apps/hyperbolic_bluenoise.cpp ${SRC})
target_link_libraries(hyperbolic_bluenoise polyscope geometry-central  bvh OpenMP::OpenMP_CXX)

add_executable(rotation_sampling apps/rotation_sampling.cpp ${SRC})
target_link_libraries(rotation_sampling polyscope geometry-central  bvh OpenMP::OpenMP_CXX)

add_executable(direction_sampling apps/direction_sampling.cpp ${SRC})
target_link_libraries(direction_sampling polyscope geometry-central bvh OpenMP::OpenMP_CXX)

add_executable(local_hyperbolic_mesh_sampling apps/local_hyperbolic_mapping.cpp ${SRC})
target_link_libraries(local_hyperbolic_mesh_sampling polyscope geometry-central bvh OpenMP::OpenMP_CXX)

add_executable(vizu_pc apps/vizu_pc.cpp)
target_link_libraries(vizu_pc polyscope geometry-central)

