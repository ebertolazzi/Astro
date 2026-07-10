############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Università degli Studi di Trento                                    #
#      email: enrico.bertolazzi@unitn.it                                   #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.25)

# ----------------------------------------------------------------------------
# Estrai versione da git tag (se disponibile)
# ----------------------------------------------------------------------------
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(GIT_VERSION MATCHES "^v")
    string(SUBSTRING ${GIT_VERSION} 1 -1 GIT_VERSION)
  endif()
endif()
if(NOT GIT_VERSION)
  set(GIT_VERSION "0.0.0")
endif()

project(
  Astro
  VERSION ${GIT_VERSION}
  DESCRIPTION "A C++ library for astrodynamics and orbital mechanics"
  HOMEPAGE_URL "https://ebertolazzi.github.io/Astro/"
  LANGUAGES CXX C
)

# ----------------------------------------------------------------------------
# Opzioni principali
# ----------------------------------------------------------------------------
# Default a Release se non specificato (per generatori a singola configurazione)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

option(BUILD_SHARED_LIBS "Build Astro as a shared library" OFF)
option(BUILD_TESTING "Build the Astro test executables" ${PROJECT_IS_TOP_LEVEL})
option(ASTRO_INSTALL "Generate the install target" ${PROJECT_IS_TOP_LEVEL})
option(ASTRO_BUILD_BENCHMARKS "Build the micro-benchmarks" ${PROJECT_IS_TOP_LEVEL})
option(ASTRO_STRICT_WARNINGS "Compile with narrowing/conversion warnings" OFF)

# Directory di installazione predefinita: ./lib
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE)
endif()
set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE)
set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE)

include(GNUInstallDirs)
include(FetchContent)

# ----------------------------------------------------------------------------
# Helper per dipendenze: locale o FetchContent
# ----------------------------------------------------------------------------
function(astro_use_local_or_fetch NAME LOCAL_DIR)
  get_filename_component(
    _local_abs
    "${LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )

  string(TOUPPER "${NAME}" _name_upper)
  string(TOLOWER "${NAME}" _name_lower)

  if(FETCHCONTENT_BASE_DIR)
    set(_fetchcontent_base_dir "${FETCHCONTENT_BASE_DIR}")
  else()
    set(_fetchcontent_base_dir "${CMAKE_BINARY_DIR}/_deps")
  endif()

  # Astro must install only Astro artifacts. Dependency install rules can
  # otherwise create symlinks/copies in lib/lib when their install scripts are
  # pulled into the top-level install target.
  set(_dependency_install_stub "${_fetchcontent_base_dir}/${_name_lower}-build/cmake_install.cmake")
  set(_saved_skip_install_rules "${CMAKE_SKIP_INSTALL_RULES}")
  set(CMAKE_SKIP_INSTALL_RULES ON)

  if(EXISTS "${_local_abs}/CMakeLists.txt")
    message(STATUS "Using local ${NAME}: ${_local_abs}")
    set(FETCHCONTENT_SOURCE_DIR_${_name_upper} "${_local_abs}" CACHE PATH "Local ${NAME} checkout" FORCE)
    FetchContent_MakeAvailable(${NAME})
    set("${NAME}_SOURCE" "local: ${_local_abs}" PARENT_SCOPE)
    set("${NAME}_SOURCE_DIR" "${_local_abs}" PARENT_SCOPE)
  else()
    message(STATUS "Fetching ${NAME}")
    FetchContent_MakeAvailable(${NAME})
    set("${NAME}_SOURCE" "FetchContent" PARENT_SCOPE)
    if(DEFINED "${_name_lower}_SOURCE_DIR")
      set("${NAME}_SOURCE_DIR" "${${_name_lower}_SOURCE_DIR}" PARENT_SCOPE)
    endif()
  endif()

  file(MAKE_DIRECTORY "${_fetchcontent_base_dir}/${_name_lower}-build")
  if(NOT EXISTS "${_dependency_install_stub}")
    file(
      WRITE "${_dependency_install_stub}"
      "# Install disabled for dependency ${NAME} inside Astro.\n"
    )
  endif()

  set(CMAKE_SKIP_INSTALL_RULES "${_saved_skip_install_rules}")
endfunction()

# ----------------------------------------------------------------------------
# Percorsi locali per checkout dei dipendenti (override via cache)
# ----------------------------------------------------------------------------
set(ASTRO_GC_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../GC" CACHE PATH "Local GenericContainer checkout")
set(ASTRO_UTILSLITE_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../UtilsLite" CACHE PATH "Local UtilsLite checkout")
set(ASTRO_QUARTIC_ROOTS_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../quarticRootsFlocke" CACHE PATH "Local quarticRootsFlocke checkout")

# ----------------------------------------------------------------------------
# 1. GenericContainer
# ----------------------------------------------------------------------------
# The dependency projects are consumed only as build dependencies. Keep their
# install/test/example rules out of Astro so lib/lib remains Astro-only.
set(GENERIC_CONTAINER_INSTALL OFF CACHE BOOL "Disable GenericContainer install rules inside Astro" FORCE)
set(GENERIC_CONTAINER_BUILD_EXAMPLES OFF CACHE BOOL "Disable GenericContainer examples inside Astro" FORCE)
set(UTILS_INSTALL_LOCAL OFF CACHE BOOL "Disable dependency install rules inside Astro" FORCE)
set(UTILS_ENABLE_TESTS OFF CACHE BOOL "Disable dependency test suites inside Astro" FORCE)

FetchContent_Declare(
  GenericContainer
  GIT_REPOSITORY https://github.com/ebertolazzi/GenericContainer.git
  GIT_TAG        GC2.0
  GIT_SHALLOW    TRUE
)
astro_use_local_or_fetch(GenericContainer "${ASTRO_GC_LOCAL_DIR}")

# Patch per Clang: forzare -stdlib=libc++ se necessario
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  foreach(_gc_target GenericContainer GenericContainerYaml)
    if(TARGET ${_gc_target})
      get_target_property(_gc_type ${_gc_target} TYPE)
      if(_gc_type STREQUAL "INTERFACE_LIBRARY")
        target_compile_options(${_gc_target} INTERFACE -stdlib=libc++)
        target_link_options(${_gc_target} INTERFACE -stdlib=libc++)
      else()
        target_compile_options(${_gc_target} PUBLIC -stdlib=libc++)
        target_link_options(${_gc_target} PUBLIC -stdlib=libc++)
      endif()
    endif()
  endforeach()
endif()

# ----------------------------------------------------------------------------
# 2. UtilsLite
# ----------------------------------------------------------------------------
FetchContent_Declare(
  UtilsLite
  GIT_REPOSITORY https://github.com/ebertolazzi/UtilsLite.git
  GIT_TAG        Utils2.0
  GIT_SHALLOW    TRUE
)
astro_use_local_or_fetch(UtilsLite "${ASTRO_UTILSLITE_LOCAL_DIR}")

# Tratta gli header di UtilsLite come SYSTEM per evitare warning
foreach(_utils_tgt utils_UtilsLite_Static utils_UtilsLite)
  if(TARGET ${_utils_tgt})
    get_target_property(_utils_inc ${_utils_tgt} INTERFACE_INCLUDE_DIRECTORIES)
    if(_utils_inc)
      set_target_properties(
        ${_utils_tgt} PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_utils_inc}"
      )
    endif()
  endif()
endforeach()

# ----------------------------------------------------------------------------
# 3. quarticRootsFlocke
# ----------------------------------------------------------------------------
FetchContent_Declare(
  quarticRootsFlocke
  GIT_REPOSITORY https://github.com/ebertolazzi/quarticRootsFlocke.git
  GIT_TAG        1.2.0
  GIT_SHALLOW    TRUE
)
astro_use_local_or_fetch(quarticRootsFlocke "${ASTRO_QUARTIC_ROOTS_LOCAL_DIR}")

# Raccogli i target delle dipendenze
set(ROOTS utils_quarticRootsFlocke_Static)
set(UTILS utils::UtilsLite_Static)
set(GC GenericContainer::Yaml)

# ----------------------------------------------------------------------------
# Funzione per impostare le directory di output per un target
# ----------------------------------------------------------------------------
function(astro_set_target_output_dirs TARGET_NAME LIB_DIR)
  set_target_properties(${TARGET_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endfunction()

# ----------------------------------------------------------------------------
# Funzione per impostare le directory di output per le dipendenze (alias)
# ----------------------------------------------------------------------------
function(astro_set_dependency_output_dirs LIB_DIR)
  foreach(tgt IN LISTS ARGN)
    get_target_property(real_tgt ${tgt} ALIASED_TARGET)
    if(NOT real_tgt)
      set(real_tgt ${tgt})
    endif()
    if(TARGET ${real_tgt})
      get_target_property(tgt_type ${real_tgt} TYPE)
      if(NOT tgt_type STREQUAL "INTERFACE_LIBRARY")
        set_target_properties(${real_tgt} PROPERTIES
          ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
          LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
        )
      endif()
    endif()
  endforeach()
endfunction()

# ----------------------------------------------------------------------------
# Raccolta delle fonti
# ----------------------------------------------------------------------------
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Raccogli tutti i file .cc (ricorsivamente) escludendo hsl_fake.cc
file(GLOB_RECURSE SOURCES "${SOURCE_DIR}/*.cc")
list(FILTER SOURCES EXCLUDE REGEX ".*hsl_fake.cc$")

# Raccogli gli header pubblici
set(PUBLIC_HEADERS)
file(GLOB_RECURSE HEADERS "${SOURCE_DIR}/*.h" "${SOURCE_DIR}/*.hpp" "${SOURCE_DIR}/*.hxx")
list(APPEND PUBLIC_HEADERS ${HEADERS})
list(FILTER PUBLIC_HEADERS EXCLUDE REGEX ".*.tmpl$")

# ----------------------------------------------------------------------------
# Target Astro
# ----------------------------------------------------------------------------
add_library(Astro ${SOURCES})
add_library(Astro::Astro ALIAS Astro)

astro_set_target_output_dirs(Astro lib/lib)

target_compile_features(Astro PUBLIC cxx_std_20)

# Include directories: pubblico per gli utenti, privato per i sorgenti
target_include_directories(
  Astro
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  PRIVATE
    "${SOURCE_DIR}"
)

# Collegamento alle dipendenze
target_link_libraries(
  Astro
  PUBLIC
    ${ROOTS}
    ${GC}
    ${UTILS}
)

# Flag di compilazione
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
  target_compile_options(Astro PRIVATE -Wall -Wextra)
  if(ASTRO_STRICT_WARNINGS)
    target_compile_options(
      Astro PRIVATE
      -Wconversion -Wsign-conversion -Wshadow -Wdouble-promotion
    )
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(Astro PUBLIC -stdlib=libc++)
  target_link_options(Astro PUBLIC -stdlib=libc++)
endif()

if(MSVC)
  target_compile_options(Astro PUBLIC /utf-8 /bigobj)
endif()

# Imposta versioni
set_target_properties(
  Astro PROPERTIES
  VERSION     ${PROJECT_VERSION}
  SOVERSION   ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME Astro
)

# Imposta directory di output per le dipendenze (vanno in lib3rd/lib)
astro_set_dependency_output_dirs(lib3rd/lib
  GenericContainer
  ${UTILS}
  ${ROOTS}
)

# ----------------------------------------------------------------------------
# Alias con suffisso di piattaforma
# ----------------------------------------------------------------------------
if(APPLE)
  execute_process(COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE)
  set(CMAKE_OSX_ARCHITECTURES ${ARCHITECTURE})
  set(ASTRO_PLATFORM_SUFFIX osx)
elseif(WIN32)
  if(MINGW)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(ASTRO_PLATFORM_SUFFIX mingw64)
    else()
      set(ASTRO_PLATFORM_SUFFIX mingw32)
    endif()
  elseif(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(ASTRO_PLATFORM_SUFFIX win64)
  else()
    set(ASTRO_PLATFORM_SUFFIX win32)
  endif()
elseif(UNIX)
  set(ASTRO_PLATFORM_SUFFIX linux)
else()
  string(TOLOWER "${CMAKE_SYSTEM_NAME}" ASTRO_PLATFORM_SUFFIX)
endif()

function(astro_add_platform_library_alias TARGET_NAME BASE_NAME)
  if(NOT TARGET ${TARGET_NAME})
    return()
  endif()

  if(BUILD_SHARED_LIBS)
    set(REAL_OUTPUT_NAME "${BASE_NAME}_${ASTRO_PLATFORM_SUFFIX}")
    set(REAL_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}${REAL_OUTPUT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(GENERIC_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}")
  else()
    set(REAL_OUTPUT_NAME "${BASE_NAME}_${ASTRO_PLATFORM_SUFFIX}_static")
    if(WIN32)
      set(REAL_FILENAME "lib${REAL_OUTPUT_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}")
      set(GENERIC_FILENAME "lib${BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    else()
      set(REAL_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${REAL_OUTPUT_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}")
      set(GENERIC_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()
  endif()

  set_target_properties(
    ${TARGET_NAME} PROPERTIES
    OUTPUT_NAME         "${REAL_OUTPUT_NAME}"
    ARCHIVE_OUTPUT_NAME "${REAL_OUTPUT_NAME}"
    LIBRARY_OUTPUT_NAME "${REAL_OUTPUT_NAME}"
    RUNTIME_OUTPUT_NAME "${REAL_OUTPUT_NAME}"
  )

  if(WIN32)
    add_custom_command(
      TARGET ${TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              "$<TARGET_FILE:${TARGET_NAME}>"
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${GENERIC_FILENAME}"
      VERBATIM
    )
  else()
    add_custom_command(
      TARGET ${TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E rm -f
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${GENERIC_FILENAME}"
      COMMAND ${CMAKE_COMMAND} -E create_symlink
              "${REAL_FILENAME}"
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${GENERIC_FILENAME}"
      VERBATIM
    )
  endif()
endfunction()

astro_add_platform_library_alias(Astro Astro)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/AstroCleanPublicLibDir.cmake" [=[
set(_public_lib_dir "${PROJECT_DIR}/lib/lib")
foreach(_name GenericContainer Lua quarticRootsFlocke UtilsLite)
  foreach(_suffix "" "_${PLATFORM_SUFFIX}" "_${PLATFORM_SUFFIX}_static")
    file(REMOVE "${_public_lib_dir}/lib${_name}${_suffix}.a")
    file(REMOVE "${_public_lib_dir}/lib${_name}${_suffix}.dylib")
    file(REMOVE "${_public_lib_dir}/lib${_name}${_suffix}.so")
  endforeach()
endforeach()
]=])

add_custom_command(
  TARGET Astro POST_BUILD
  COMMAND ${CMAKE_COMMAND}
          -DPROJECT_DIR=${CMAKE_CURRENT_SOURCE_DIR}
          -DPLATFORM_SUFFIX=${ASTRO_PLATFORM_SUFFIX}
          -P "${CMAKE_CURRENT_BINARY_DIR}/AstroCleanPublicLibDir.cmake"
  VERBATIM
)

# ----------------------------------------------------------------------------
# Test
# ----------------------------------------------------------------------------
include(CTest)
if(BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
  # Elenco specifico dei test di Astro (come nell'originale)
  set(TEST_NAMES
    test1-Lambert
    test2-DeltaV
    test3-DeltaV
  )

  add_custom_target("${PROJECT_NAME}_all_tests" ALL)

  foreach(S ${TEST_NAMES})
    add_executable(${S} "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc")
    target_link_libraries(${S} PRIVATE Astro::Astro ${CMAKE_DL_LIBS})
    # Aggiungi la directory src/ per trovare eventuali header interni
    target_include_directories(${S} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    set_target_properties(
      ${S} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    )
    add_test(NAME "${S}" COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/bin/${S}" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_dependencies("${PROJECT_NAME}_all_tests" ${S})
  endforeach()
endif()

# ----------------------------------------------------------------------------
# Installazione
# ----------------------------------------------------------------------------
if(ASTRO_INSTALL)
  install(
    TARGETS Astro
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(CODE
    "file(REMOVE_RECURSE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/Users\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libGenericContainer${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libGenericContainer${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libGenericContainer_${ASTRO_PLATFORM_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libGenericContainer_${ASTRO_PLATFORM_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libGenericContainer_${ASTRO_PLATFORM_SUFFIX}_static${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libGenericContainer_${ASTRO_PLATFORM_SUFFIX}_static${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libLua${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libLua${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libLua_${ASTRO_PLATFORM_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libLua_${ASTRO_PLATFORM_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libLua_${ASTRO_PLATFORM_SUFFIX}_static${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libLua_${ASTRO_PLATFORM_SUFFIX}_static${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libquarticRootsFlocke${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libquarticRootsFlocke${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libquarticRootsFlocke_${ASTRO_PLATFORM_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libquarticRootsFlocke_${ASTRO_PLATFORM_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libquarticRootsFlocke_${ASTRO_PLATFORM_SUFFIX}_static${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libquarticRootsFlocke_${ASTRO_PLATFORM_SUFFIX}_static${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libUtilsLite${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libUtilsLite${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libUtilsLite_${ASTRO_PLATFORM_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libUtilsLite_${ASTRO_PLATFORM_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libUtilsLite_${ASTRO_PLATFORM_SUFFIX}_static${CMAKE_STATIC_LIBRARY_SUFFIX}\")\n"
    "file(REMOVE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/lib/libUtilsLite_${ASTRO_PLATFORM_SUFFIX}_static${CMAKE_SHARED_LIBRARY_SUFFIX}\")\n"
  )

  # CPack (opzionale)
  set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
  set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
  set(CPACK_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;/build.*/")
  include(CPack)
endif()

# ----------------------------------------------------------------------------
# Riepilogo
# ----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "Astro dependency summary:")
message(STATUS "  GenericContainer  : ${GenericContainer_SOURCE}")
message(STATUS "  UtilsLite         : ${UtilsLite_SOURCE}")
message(STATUS "  quarticRootsFlocke: ${quarticRootsFlocke_SOURCE}")
message(STATUS "")