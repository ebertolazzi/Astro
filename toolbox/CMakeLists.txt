#-----------------------------------------------------------------------#
#  file: CMakeLists.txt                                                 #
#                                                                       #
#  version: 1.0   date 14/3/2022                                        #
#                                                                       #
#  Copyright (C) 2022                                                   #
#                                                                       #
#      Enrico Bertolazzi, Francesco Biral and Paolo Bosetti             #
#      Dipartimento di Ingegneria Industriale                           #
#      Universita` degli Studi di Trento                                #
#      Via Sommarive 9, I-38123, Trento, Italy                          #
#      email: enrico.bertolazzi@unitn.it                                #
#             francesco.biral@unitn.it                                  #
#             paolo.bosetti@unitn.it                                    #
#-----------------------------------------------------------------------#


cmake_minimum_required( VERSION 3.14 )

project( Astro )

find_package( Matlab REQUIRED )

message( STATUS "Matlab_ROOT_DIR = ${Matlab_ROOT_DIR}" )
message( STATUS "PROJECT_NAME    = ${PROJECT_NAME}" )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)
if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2")
else()
  set(CMAKE_CXX_FLAGS "-Ofast -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
endif()

set( SOURCES )
file( GLOB S ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc )
foreach(F ${S})
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach()
file( GLOB S ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cc )
foreach(F ${S})
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/include )

if ( (UNIX OR LINUX) AND NOT APPLE )
  set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
  set(CMAKE_SHARED_LINKER_FLAGS "-static-libgcc -static-libstdc++")
endif()

set(
  MEX_CMDS
  AstroMexWrapper
  Lambert
  DeltaV
  Check_EQ_for_consistency
)

foreach(CMD ${MEX_CMDS})

  matlab_add_mex(
    NAME ${CMD}
    SRC src_mex/mex_${CMD}.cc src_mex/GenericContainerMatlabInterface.cc ${SOURCES}
  )

  add_custom_command(
    TARGET ${CMD}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E copy
    "$<TARGET_FILE:${CMD}>"
    "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endforeach()
