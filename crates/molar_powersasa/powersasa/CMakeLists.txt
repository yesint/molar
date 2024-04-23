cmake_minimum_required(VERSION 3.18)
project(powersasa_wrapper CXX)

# -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD  17)

add_library(powersasa_wrapper STATIC
    ${PROJECT_SOURCE_DIR}/wrapper.cpp
    ${PROJECT_SOURCE_DIR}/wrapper.h
    ${PROJECT_SOURCE_DIR}/power_diagram.h
    ${PROJECT_SOURCE_DIR}/power_sasa.h    
)

target_include_directories(powersasa_wrapper PRIVATE ${PROJECT_SOURCE_DIR})

install(TARGETS powersasa_wrapper
    LIBRARY DESTINATION .
    ARCHIVE DESTINATION .
)