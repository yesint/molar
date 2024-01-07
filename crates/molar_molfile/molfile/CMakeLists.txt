#---------------------
# VMD molfile plugins
#---------------------
cmake_minimum_required(VERSION 3.10)
project(molfile_plugins C)

# -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#-------------------------------------
# VMD molfile plugins build machinery
#-------------------------------------
SET(MOLFILE_PLUGINS pdbplugin dcdplugin xyzplugin)
SET(MOLFILE_PLUGINS_FILE_LIST "")

foreach(P ${MOLFILE_PLUGINS})    
    set_source_files_properties(${P}.c PROPERTIES COMPILE_DEFINITIONS "VMDPLUGIN=${P};VMDUSECONECTRECORDS")
    list(APPEND MOLFILE_PLUGINS_FILE_LIST ${P}.c)
endforeach()

add_library(molfile STATIC 
    ${MOLFILE_PLUGINS_FILE_LIST}
    register_plugins.c
)
target_include_directories(molfile PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

#------------------------------------------
# Install to the output directory directly
#------------------------------------------
install(TARGETS molfile 
    LIBRARY DESTINATION .
    ARCHIVE DESTINATION .
)
