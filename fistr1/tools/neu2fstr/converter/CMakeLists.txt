###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND neu2fstr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/cconv_mat.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CConvMessage.cpp
  ${CMAKE_CURRENT_LIST_DIR}/conv_neu2hec.cpp
  ${CMAKE_CURRENT_LIST_DIR}/conv_neu2fstr_static.cpp
  ${CMAKE_CURRENT_LIST_DIR}/conv_neu2fstr_dynamic.cpp
  ${CMAKE_CURRENT_LIST_DIR}/conv_neu2fstr_heat.cpp
)

list(APPEND neu2fstr_INCLUDE_DIRS
  ${CMAKE_CURRENT_LIST_DIR}
)

