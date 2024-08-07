cmake_minimum_required(VERSION 3.18)
project(gromacs_wrapper LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD  20)
set(CMAKE_C_STANDARD    99)

message(DEBUG "ENV: ${GROMACS_SOURCE_DIR}")
message(DEBUG "ENV: ${GROMACS_BINARY_DIR}")
message(DEBUG "ENV: ${GROMACS_LIB_DIR}")

# Check if external gromacs location is provided
# User need to provide source dir, binary dir and the exact library dir
if(GROMACS_SOURCE_DIR AND GROMACS_BINARY_DIR AND GROMACS_LIB_DIR)
    # Set locations of existing lib files
    set(GROMACS_LIB_FILES
        ${GROMACS_LIB_DIR}/${CMAKE_DYNAMIC_LIBRARY_PREFIX}gromacs${CMAKE_DYNAMIC_LIBRARY_SUFFIX}
        ${GROMACS_LIB_DIR}/${CMAKE_DUNAMIC_LIBRARY_PREFIX}muparser${CMAKE_DYNAMIC_LIBRARY_SUFFIX}
    )
    # Nothing else to do
else()
    # Compiling Gromacs in place
    set(IN_PLACE TRUE)
    # Set CPM cache location if not provide in the command line
    if(NOT CPM_SOURCE_CACHE)
        set(CPM_SOURCE_CACHE "${PROJECT_BINARY_DIR}/CPM_CACHE" CACHE STRING "Location of CPM cache" FORCE)
    endif()

    include(CPM.cmake)
    include(ExternalProject)

    CPMAddPackage(
        NAME            GROMACS
        GIT_REPOSITORY  https://gitlab.com/gromacs/gromacs.git
        GIT_TAG         main
        DOWNLOAD_ONLY   ON
        BINARY_DIR      ${CMAKE_BINARY_DIR}/../..
    )
    # Variable GROMACS_SOURCE_DIR and GROMACS_BINARY_DIR are set by CMP at this point
    message(INFO "Locations: ${GROMACS_SOURCE_DIR} ${GROMACS_BINARY_DIR}")

    # Files to be set later as a byproduct to trigger buld correctly
    set(GROMACS_LIB_FILES
        ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}muparser${CMAKE_STATIC_LIBRARY_SUFFIX}
    )

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

        # If you have not enough RAM uncomment the next line and set smaller number of parallel processes
        # BUILD_COMMAND ${CMAKE_COMMAND} --build ${GROMACS_BINARY_DIR} -j 12

        # Copy compiled libs to the place where cargo sees them
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy ${GROMACS_LIB_FILES} ${CMAKE_BINARY_DIR}/..
        BUILD_BYPRODUCTS ${GROMACS_LIB_FILES}
    )
endif()

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

# Ugly workaround to pass gromacs source and bin location
# to cargo by writing them to a text file
if(IN_PLACE)
    configure_file(cmake_to_cargo.txt.in ${CMAKE_BINARY_DIR}/../cmake_to_cargo.txt)
endif()

#--------------------------------------------------------------------------------
# Create a wrapper library
#--------------------------------------------------------------------------------
add_library(gromacs_wrapper STATIC wrapper.cpp)
if(IN_PLACE)
    add_dependencies(gromacs_wrapper Gromacs_external)
endif()
target_include_directories(gromacs_wrapper PUBLIC ${GROMACS_INCLUDE_DIRECTORIS})
target_link_libraries(gromacs_wrapper PUBLIC ${GROMACS_LIBRARIES})
install(TARGETS gromacs_wrapper
    LIBRARY DESTINATION .
    ARCHIVE DESTINATION .
)