############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      <mailto:enrico.bertolazzi@unitn.it>                                 #
#                                                                          #
############################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
SET( CMAKE_VERBOSE_MAKEFILE OFF )
SET( CMAKE_INSTALL_MESSAGE NEVER )

IF ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  SET( DEBUG_SUFFIX "_debug" )
ELSE()
  SET( DEBUG_SUFFIX )
ENDIF()

IF ( BUILD_SHARED )
  IF ( APPLE )
    SET( TARGET  Astro_osx${DEBUG_SUFFIX} )
  ELSEIF( LINUX OR UNIX )
    SET( TARGET  Astro_linux${DEBUG_SUFFIX} )
  ELSE()
    SET( TARGET  Astro_win_${BITS}${DEBUG_SUFFIX} )
    SET( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )
  ENDIF()
ENDIF()

IF ( APPLE )
  SET( TARGETS Astro_osx_static${DEBUG_SUFFIX} )
ELSEIF( LINUX OR UNIX )
  SET( TARGETS Astro_linux_static${DEBUG_SUFFIX} )
ELSE()
  SET( TARGETS Astro_win_${BITS}_static${DEBUG_SUFFIX} )
ENDIF()

PROJECT( Astro CXX C )

INCLUDE( ./CMakeLists-cflags.txt )
MESSAGE( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")

SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

SET( SOURCES )
FILE( GLOB S ./src/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND SOURCES ${RF} )
ENDFOREACH (F ${S})

SET( SRCS_LIBS )
#FILE( GLOB S ./submodules/quarticRootsFlocke/src/*.cc )
#FOREACH (F ${S})
#  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
#  LIST( APPEND SRCS_LIBS ${RF} )
#ENDFOREACH (F ${S})
FILE( GLOB S ./submodules/Utils/src/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND SRCS_LIBS ${RF} )
ENDFOREACH (F ${S})
FILE( GLOB S ./submodules/GenericContainer/src/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND SRCS_LIBS ${RF} )
ENDFOREACH (F ${S})

SET( HEADERS )
FILE( GLOB S ./src/*.h* )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})
SET( HEADERS1 )
FILE( GLOB S ./src/Astro/*.h* )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS1 ${RF} )
ENDFOREACH (F ${S})

#FILE( GLOB S ./submodules/quarticRootsFlocke/src/*.hh )
#FOREACH (F ${S})
#  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
#  LIST( APPEND HEADERS ${RF} )
#ENDFOREACH (F ${S})

INCLUDE_DIRECTORIES( src lib3rd/include lib3rd/include/Utils )

SET( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib )

IF ( BUILD_SHARED )
  ADD_LIBRARY( ${TARGETS} STATIC ${SOURCES} ${HEADERS} )
  ADD_LIBRARY( ${TARGET}  SHARED ${SOURCES} ${HEADERS} )
  INSTALL(
    TARGETS ${TARGET} ${TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
ELSE()
  ADD_LIBRARY( ${TARGETS} STATIC ${SOURCES} ${HEADERS} )
  INSTALL(
    TARGETS ${TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
ENDIF()

IF( BUILD_EXECUTABLE )
  SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
  SET( EXECUTABLE )
  FOREACH ( EXE ${EXECUTABLE} )
    ADD_EXECUTABLE( ${EXE} tests-cpp/${EXE}.cc ${SRCS_LIBS} ${HEADERS} )
    IF ( UNIX )
      TARGET_LINK_LIBRARIES( ${EXE} ${TARGETS} -ldl )
    ELSE()
      TARGET_LINK_LIBRARIES( ${EXE} ${TARGETS} )
    ENDIF()
  ENDFOREACH ( EXE ${EXECUTABLE} )
ENDIF()

SET_PROPERTY( TARGET ${TARGETS} PROPERTY POSITION_INDEPENDENT_CODE ON )

INSTALL( FILES ${HEADERS}  DESTINATION ${CMAKE_INSTALL_PREFIX}/include )
INSTALL( FILES ${HEADERS1} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Astro )

MESSAGE( STATUS "Using ${SSE_FLAGS} extensions")
MESSAGE( STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}" )
MESSAGE( STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}" )
MESSAGE( STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}" )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}" )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}" )
IF(APPLE)
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()
MESSAGE( STATUS "BUILD_EXECUTABLE            = ${BUILD_EXECUTABLE}" )