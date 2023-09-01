#------------------------------
# Modified xdrfile library
#------------------------------
cmake_minimum_required(VERSION 3.10)
project(xdrfile C)

# -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(xdrfile STATIC
    xdrfile.c
    xdrfile_trr.c
    xdrfile_xtc.c
    xdr_seek.c
    xdr_utils.c
)

target_include_directories(xdrfile PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

#------------------------------------------
# Install to the output directory directly
#------------------------------------------
install(TARGETS xdrfile
    LIBRARY DESTINATION .
    ARCHIVE DESTINATION .
)
