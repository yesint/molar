#
# This file is a part of
cmake_minimum_required(VERSION 3.18)
project(gromacs_wrapper CXX)

# Set CPM cache location if not provide in the command line
if(NOT CPM_SOURCE_CACHE)
    set(CPM_SOURCE_CACHE "~/.local/CPM_CACHE" CACHE STRING "Location of CPM cache" FORCE)
endif()

include(CPM.cmake)
include(ExternalProject)

CPMAddPackage(
    NAME            GROMACS
    GIT_REPOSITORY  https://gitlab.com/gromacs/gromacs.git
    GIT_TAG         main
    DOWNLOAD_ONLY   ON
)
# Variable GROMACS_SOURCE_DIR and GROMACS_BINARY_DIR are set by CMP at this point

# File to be set as a byproduct to trigger buld correctly
set(GROMACS_LIB_FILES
    ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}muparser${CMAKE_STATIC_LIBRARY_SUFFIX}
)

# Get actual gromacs version
# Since Gromacs 2022.x version is stored in main CmakeLists file
file(READ ${GROMACS_SOURCE_DIR}/CMakeLists.txt gmxinfofile)
string(REGEX MATCH "Gromacs VERSION +([0-9]+)" match ${gmxinfofile})
set(GROMACS_VERSION ${CMAKE_MATCH_1})

message(STATUS "Gromacs version used: ${GROMACS_VERSION}")

ExternalProject_add(Gromacs_external
    SOURCE_DIR ${GROMACS_SOURCE_DIR}
    BINARY_DIR ${GROMACS_BINARY_DIR}
    CMAKE_ARGS  -DGMX_MPI=OFF
                -DGMX_GPU=OFF
                -DGMX_SIMD=none
                -DGMX_FFT_LIBRARY=fftpack
                -DBUILD_TESTING=OFF
                -DGMXAPI=OFF
                -DGMX_IMD=OFF
                -DGMX_INSTALL_NBLIB_API=OFF
                -DGMX_OPENMP=OFF
                -DGMX_THREAD_MPI=OFF
                -DCMAKE_DISABLE_FIND_PACKAGE_MPI=on
                -DBUILD_SHARED_LIBS=OFF
                -DGMX_USE_TNG=OFF
                -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                -DCMAKE_BUILD_TYPE=Release
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                
    #BUILD_COMMAND ${CMAKE_COMMAND} --build ${GROMACS_BINARY_DIR} -j 12
    # Copy compiled libs to the place where cargo sees them
    #INSTALL_COMMAND ${CMAKE_COMMAND} -E copy ${GROMACS_LIB_FILES} ${CMAKE_BINARY_DIR}/..
    BUILD_BYPRODUCTS ${GROMACS_LIB_FILES}
)

set(GROMACS_INCLUDE_DIRECTORIS
    ${GROMACS_SOURCE_DIR}/src
    ${GROMACS_SOURCE_DIR}/src/gromacs/utility/include
    ${GROMACS_SOURCE_DIR}/src/gromacs/math/include
    ${GROMACS_SOURCE_DIR}/src/gromacs/topology/include
    ${GROMACS_SOURCE_DIR}/api/legacy/include
    ${GROMACS_SOURCE_DIR}/src/external
    ${GROMACS_BINARY_DIR}/api/legacy/include   # For generated headers
)
set(GROMACS_LIBRARIES ${GROMACS_LIB_FILES})

# Configure include file with Gromacs version
#configure_file(${PROJECT_SOURCE_DIR}/gromacs_version_info.h.in
#                ${CMAKE_BINARY_DIR}/gromacs_version_info.h @ONLY)

#--------------------------------------------------------------------------------
# Create a wrapper library
#--------------------------------------------------------------------------------
add_library(gromacs_wrapper STATIC wrapper.cpp)
add_dependencies(gromacs_wrapper Gromacs_external)
target_include_directories(gromacs_wrapper PRIVATE ${GROMACS_INCLUDE_DIRECTORIS})
target_link_libraries(gromacs_wrapper PRIVATE ${GROMACS_LIBRARIES})
install(TARGETS gromacs_wrapper
    LIBRARY DESTINATION .
    ARCHIVE DESTINATION .
)