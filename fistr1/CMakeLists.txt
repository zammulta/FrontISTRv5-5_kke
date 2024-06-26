###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

project(fistr CXX Fortran C)

set(fistr_SOURCES)
set(fistr_INCLUDE_DIRS)
set(fistr_DEFINITIONS)
set(fistr_LIBRARIES)
set(fistr1_SOURCES)

#
# Inherit definitions from fistr options
#
set(fistr_INCLUDE_DIRS ${FrontISTR_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/hecmw1)
set(fistr_DEFINITIONS ${FrontISTR_DEFINITIONS})
set(fistr_LIBRARIES ${FrontISTR_LIBRARIES})

#
# Build for libfistr.a
#
include(src/CMakeLists.txt)
add_library(fistr STATIC ${fistr_SOURCES})
target_link_libraries(fistr hecmw)
target_include_directories(fistr PUBLIC
  ${fistr_INCLUDE_DIRS} ${hecmw_INCLUDE_DIRS}
)
target_compile_definitions(fistr PUBLIC
  ${fistr_DEFINITIONS} ${hecmw_DEFINITIONS}
)

#
# Build for fistr1 (executable)
#
include(src/main/CMakeLists.txt)
add_executable(fistr1 ${fistr1_SOURCES})
target_link_libraries(fistr1
  fistr hecmw ${hecmw_LIBRARIES} ${fistr_LIBRARIES}
  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
)

set_target_properties(fistr1 PROPERTIES LINKER_LANGUAGE "CXX")

target_include_directories(fistr1 PUBLIC
  ${fistr_INCLUDE_DIRS}
  ${hecmw_INCLUDE_DIRS}
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_compile_definitions(fistr1 PUBLIC
  ${fistr_DEFINITIONS}
  ${hecmw_DEFINITIONS}
)

#
# Set installing files
#
install(TARGETS fistr1
  DESTINATION bin
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
)

#
# Build tools
#
if(WITH_TOOLS)
  add_subdirectory(tools)
endif()
