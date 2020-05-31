cmake_minimum_required(VERSION 3.10)
project(pmsm_sim)

# default to release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS "-std=c++14")

# find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)

add_subdirectory(lib/common_cpp)

add_definitions(-DPARAM_DIR="${CMAKE_CURRENT_SOURCE_DIR}/param")

include_directories(include
    lib/common_cpp/include
    # ${EIGEN3_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
)

add_executable(main
    src/main.cpp
    src/pmsm.cpp
    )
target_link_libraries(main
    ${YAML_CPP_LIBRARIES}
    )

