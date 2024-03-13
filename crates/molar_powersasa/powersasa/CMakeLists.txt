cmake_minimum_required(VERSION 3.18)
project(powersasa_wrapper CXX)

# -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD  17)

# Set CPM cache location if not provide in the command line
if(NOT CPM_SOURCE_CACHE)
set(CPM_SOURCE_CACHE "${PROJECT_BINARY_DIR}/CPM_CACHE" CACHE STRING "Location of CPM cache" FORCE)
endif()
# Include CPM
include(CPM.cmake)

# Get Eigen
CPMAddPackage(
    NAME              Eigen3
    GIT_REPOSITORY    https://gitlab.com/libeigen/eigen.git
    GIT_TAG           3.4
    OPTIONS
        "EIGEN_BUILD_DOC OFF"
        "BUILD_TESTING OFF"
        "EIGEN_BUILD_TESTING OFF"
        "EIGEN_BUILD_PKGCONFIG OFF"
)

add_library(powersasa_wrapper STATIC
    ${PROJECT_SOURCE_DIR}/wrapper.cpp
    ${PROJECT_SOURCE_DIR}/wrapper.h
    ${PROJECT_SOURCE_DIR}/power_diagram.h
    ${PROJECT_SOURCE_DIR}/power_sasa.h    
)

target_include_directories(powersasa_wrapper PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(powersasa_wrapper PRIVATE Eigen3::Eigen)

install(TARGETS powersasa_wrapper
    LIBRARY DESTINATION .
    ARCHIVE DESTINATION .
)

#add_executable(sasatest sasatest.cpp)
#target_include_directories(sasatest PRIVATE ${PROJECT_SOURCE_DIR})
#target_link_libraries(sasatest PRIVATE Eigen3::Eigen powersasa_wrapper)
#install(TARGETS sasatest DESTINATION .)
