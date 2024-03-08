cmake_minimum_required(VERSION 3.18)
project(powersasa CXX)

# -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD  20)

add_library(powersasa_wrapper STATIC
${PROJECT_SOURCE_DIR}/bindings.cpp
${PROJECT_SOURCE_DIR}/power_diagram.h
${PROJECT_SOURCE_DIR}/power_sasa.h
)

install(TARGETS powersasa_wrapper
    LIBRARY DESTINATION .
    ARCHIVE DESTINATION .
)