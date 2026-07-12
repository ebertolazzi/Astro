############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#                                                                          #
#  Astro -- standalone CMake build.                                        #
#                                                                          #
#  This build script configures the Astro library, its local or fetched     #
#  dependencies, optional tests/benchmarks, platform-specific output names, #
#  local install layout, and packaging metadata.                           #
#                                                                          #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Università degli Studi di Trento                                    #
#      email: enrico.bertolazzi@unitn.it                                   #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.25)

# --------------------------------------------------------------------------
# Apple FAT/universal build setup
# --------------------------------------------------------------------------
# On Apple platforms a single library can contain slices for several CPU
# architectures, for example arm64 and x86_64.  This is controlled by
# CMAKE_OSX_ARCHITECTURES.
#
# The variable must be set before project(), otherwise CMake may initialize
# compilers and targets before the requested architectures are known.
#
# Default behavior:
#   - ON  on Apple hosts or Apple cross-compilation targets.
#   - OFF on Linux/Windows and other platforms, where there is no equivalent
#     single-file universal library format.
#
# User overrides:
#   -DUTILS_BUILD_FAT=OFF
#   -DUTILS_FAT_ARCHITECTURES=arm64
#   -DCMAKE_OSX_ARCHITECTURES=arm64\;x86_64
set(UTILS_BUILD_FAT_DEFAULT OFF)
if(CMAKE_HOST_APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin|iOS|tvOS|watchOS|visionOS")
  set(UTILS_BUILD_FAT_DEFAULT ON)
endif()
option(UTILS_BUILD_FAT "Build Apple universal/fat binaries when supported" ${UTILS_BUILD_FAT_DEFAULT})
set(UTILS_FAT_ARCHITECTURES "arm64;x86_64" CACHE STRING "Apple architectures to build when UTILS_BUILD_FAT is ON")

if(UTILS_BUILD_FAT AND NOT CMAKE_OSX_ARCHITECTURES)
  set(
    CMAKE_OSX_ARCHITECTURES
    "${UTILS_FAT_ARCHITECTURES}"
    CACHE STRING "Build architectures for Apple platforms"
    FORCE
  )
endif()


# --------------------------------------------------------------------------
# Project version discovery
# --------------------------------------------------------------------------
# Prefer the latest Git tag when the source tree is a Git checkout.  Tags may
# be written either as plain semantic versions, e.g. 1.2.3, or with the common
# leading "v" prefix, e.g. v1.2.3.  The leading "v" is stripped so that the
# value can be passed directly to project(VERSION ...).
#
# If Git is unavailable, or if the source archive does not contain tags, the
# build falls back to 0.0.0.  This keeps source distributions configurable even
# outside a repository clone.
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

# --------------------------------------------------------------------------
# Top-level build options and default build type
# --------------------------------------------------------------------------
# Single-configuration generators such as Unix Makefiles and Ninja use
# CMAKE_BUILD_TYPE.  Multi-configuration generators such as Visual Studio or
# Ninja Multi-Config ignore it and select the configuration at build time.
#
# When the user did not explicitly select a configuration, default to Release
# to match the historical project workflow and avoid unoptimized local builds.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

option(BUILD_SHARED_LIBS "Build Astro as a shared library" OFF)
option(BUILD_TESTING "Build the Astro test executables" ${PROJECT_IS_TOP_LEVEL})
option(ASTRO_INSTALL "Generate the install target" ${PROJECT_IS_TOP_LEVEL})
option(ASTRO_BUILD_BENCHMARKS "Build the micro-benchmarks" ${PROJECT_IS_TOP_LEVEL})
option(ASTRO_STRICT_WARNINGS "Compile with narrowing/conversion warnings" OFF)

# --------------------------------------------------------------------------
# Local install layout
# --------------------------------------------------------------------------
# The default install prefix is the project-local lib/ directory.  This keeps
# headers, libraries, and generated package artifacts inside the source tree,
# matching the existing Rake/local-build layout used by the related projects.
#
# The cache entries remain user-overridable, but they are forced to the local
# include/lib subdirectories when the prefix is still CMake's default.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE)
endif()
set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE)
set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE)

# GNUInstallDirs provides standard install directory variables such as
# CMAKE_INSTALL_INCLUDEDIR and CMAKE_INSTALL_LIBDIR.  FetchContent is used to
# consume sibling projects either from local checkouts or from their Git URLs.
include(GNUInstallDirs)
include(FetchContent)

# --------------------------------------------------------------------------
# Dependency helper: use a local checkout when available, otherwise fetch
# --------------------------------------------------------------------------
# astro_use_local_or_fetch(NAME LOCAL_DIR) implements the dependency policy
# used by the whole project:
#
#   1. Resolve LOCAL_DIR relative to the Astro source tree.
#   2. If it contains a CMakeLists.txt, expose it to FetchContent through
#      FETCHCONTENT_SOURCE_DIR_<NAME> and build that local checkout.
#   3. Otherwise, let FetchContent download/populate the dependency declared
#      with FetchContent_Declare().
#   4. Temporarily disable dependency install rules so the top-level install
#      target installs only Astro artifacts.
#   5. Write a small cmake_install.cmake stub for the dependency build tree.
#      This avoids missing-install-script issues in some generated install
#      target graphs.
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

# --------------------------------------------------------------------------
# Local dependency checkout paths
# --------------------------------------------------------------------------
# These cache variables point to the sibling checkout layout normally used for
# local development.  A developer can override any of them from the command
# line.  If the directory does not contain a CMake project, the corresponding
# FetchContent declaration below is used instead.
set(ASTRO_GC_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../GC" CACHE PATH "Local GenericContainer checkout")
set(ASTRO_UTILSLITE_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../UtilsLite" CACHE PATH "Local UtilsLite checkout")
set(ASTRO_QUARTIC_ROOTS_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../quarticRootsFlocke" CACHE PATH "Local quarticRootsFlocke checkout")

# --------------------------------------------------------------------------
# Dependency 1: GenericContainer
# --------------------------------------------------------------------------
# GenericContainer provides the recursive heterogeneous container and YAML
# bridge target used by Astro.  It is consumed as a dependency only: install,
# examples, and dependency-local tests are disabled to keep this build focused
# on Astro outputs.
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

# Clang compatibility patch.
# The related projects historically expect libc++ when Clang is used.  Apply
# the flag to GenericContainer targets only when those targets exist, and use
# INTERFACE vs PUBLIC according to the target type.
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

# --------------------------------------------------------------------------
# Dependency 2: UtilsLite
# --------------------------------------------------------------------------
# UtilsLite provides common utility code used by Astro and by other sibling
# libraries.  The FetchContent declaration is paired with the helper above, so
# a local checkout is preferred whenever present.
FetchContent_Declare(
  UtilsLite
  GIT_REPOSITORY https://github.com/ebertolazzi/UtilsLite.git
  GIT_TAG        Utils2.0
  GIT_SHALLOW    TRUE
)
astro_use_local_or_fetch(UtilsLite "${ASTRO_UTILSLITE_LOCAL_DIR}")

# Treat UtilsLite include directories as SYSTEM includes.
# This prevents warnings from dependency headers from being reported as Astro
# warnings, which is especially useful when enabling stricter warning options
# on the Astro target itself.
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

# --------------------------------------------------------------------------
# Dependency 3: quarticRootsFlocke
# --------------------------------------------------------------------------
# quarticRootsFlocke supplies polynomial-root routines used indirectly by the
# geometry/orbital algorithms.  It follows the same local-or-fetch workflow as
# the other sibling libraries.
FetchContent_Declare(
  quarticRootsFlocke
  GIT_REPOSITORY https://github.com/ebertolazzi/quarticRootsFlocke.git
  GIT_TAG        1.2.0
  GIT_SHALLOW    TRUE
)
astro_use_local_or_fetch(quarticRootsFlocke "${ASTRO_QUARTIC_ROOTS_LOCAL_DIR}")

# Canonical dependency target names used later by target_link_libraries().
# Keeping them in variables makes the main Astro target declaration shorter
# and mirrors the naming convention used by the sibling projects.
set(ROOTS utils_quarticRootsFlocke_Static)
set(UTILS utils::UtilsLite_Static)
set(GC GenericContainer::Yaml)

# --------------------------------------------------------------------------
# Output directory helper for Astro-owned targets
# --------------------------------------------------------------------------
# CMake's global output-directory variables are intentionally avoided.  Setting
# output directories per target is more predictable when this project is used
# either as a top-level build or as a subproject through FetchContent.
function(astro_set_target_output_dirs TARGET_NAME LIB_DIR)
  set_target_properties(${TARGET_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endfunction()

# --------------------------------------------------------------------------
# Output directory helper for dependency targets
# --------------------------------------------------------------------------
# Dependencies may be referenced through ALIAS targets.  This helper resolves
# aliases to their real targets and moves non-interface dependency artifacts to
# lib3rd/lib, keeping Astro's public lib/lib directory reserved for Astro's own
# library and its compatibility link.
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

# --------------------------------------------------------------------------
# Source and public-header discovery
# --------------------------------------------------------------------------
# Source files are collected recursively from src/.  The hsl_fake.cc exclusion
# mirrors the convention used in related numerical projects where HSL stubs may
# exist but should not be compiled into this library.
#
# Public headers are collected for IDE visibility and for any future install or
# packaging logic.  Template helper files are excluded because they are not
# standalone public headers.
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Collect all implementation files recursively, excluding the optional HSL
# placeholder implementation.
file(GLOB_RECURSE SOURCES "${SOURCE_DIR}/*.cc")
list(FILTER SOURCES EXCLUDE REGEX ".*hsl_fake.cc$")

# Collect public headers from src/.  They are not passed to add_library() in
# this file, but the list remains useful for documentation, IDE integration,
# and future install logic.
set(PUBLIC_HEADERS)
file(GLOB_RECURSE HEADERS "${SOURCE_DIR}/*.h" "${SOURCE_DIR}/*.hpp" "${SOURCE_DIR}/*.hxx")
list(APPEND PUBLIC_HEADERS ${HEADERS})
list(FILTER PUBLIC_HEADERS EXCLUDE REGEX ".*.tmpl$")

# --------------------------------------------------------------------------
# Main Astro library target
# --------------------------------------------------------------------------
# The Astro target is the single public library produced by this project.  It
# exports the namespace alias Astro::Astro for consumers and writes binary
# artifacts directly under lib/lib through astro_set_target_output_dirs().
add_library(Astro ${SOURCES})
add_library(Astro::Astro ALIAS Astro)

astro_set_target_output_dirs(Astro lib/lib)

target_compile_features(Astro PUBLIC cxx_std_20)

# Public include directories are visible to consumers that link Astro.  The
# source directory remains PRIVATE because it is needed only while compiling
# Astro's own implementation files.
target_include_directories(
  Astro
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  PRIVATE
    "${SOURCE_DIR}"
)

# Link the sibling dependency targets publicly.  Their usage requirements
# include headers, compile definitions, and transitive link flags needed by
# downstream consumers of Astro::Astro.
target_link_libraries(
  Astro
  PUBLIC
    ${ROOTS}
    ${GC}
    ${UTILS}
)

# Compiler warning policy for the Astro target.  The default enables a modest
# warning baseline on GCC/Clang-compatible compilers; ASTRO_STRICT_WARNINGS adds
# conversion and narrowing diagnostics when the maintainer wants a stricter
# cleanup pass.
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

# Version metadata for shared-library sonames and exported target names.  This
# is harmless for static builds and important for shared-library packaging.
set_target_properties(
  Astro PROPERTIES
  VERSION     ${PROJECT_VERSION}
  SOVERSION   ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME Astro
)

# Move dependency artifacts to lib3rd/lib.  Astro's lib/lib directory should
# contain only libAstro_<platform>[_static] and the generic compatibility name.
astro_set_dependency_output_dirs(lib3rd/lib
  GenericContainer
  ${UTILS}
  ${ROOTS}
)

# --------------------------------------------------------------------------
# Platform suffix and generic compatibility library name
# --------------------------------------------------------------------------
# The real artifact name includes a platform suffix, for example
# libAstro_osx_static.a or libAstro_linux_static.a.  The unsuffixed file name,
# for example libAstro.a, is kept only as a compatibility alias.
#
# Do not set CMAKE_OSX_ARCHITECTURES here.  Apple FAT/universal architecture
# selection is handled before project() by the dedicated block above; changing
# it here would be too late and would override user-specified universal builds.
if(APPLE)
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
  # Ignore missing targets gracefully.  This makes the helper safe to call from
  # conditional blocks and keeps future optional targets easy to wire in.
  if(NOT TARGET ${TARGET_NAME})
    return()
  endif()

  # Compute the real artifact name and the generic compatibility name.  Shared
  # and static libraries use different suffixes and file extensions, and Windows
  # uses the historical lib<name>.a convention expected by this project family.
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

  # Apply the platform-suffixed output name to every artifact class CMake may
  # generate for the target: archive libraries, shared libraries, and runtime
  # binaries.  This keeps single- and multi-config generators consistent.
  set_target_properties(
    ${TARGET_NAME} PROPERTIES
    OUTPUT_NAME         "${REAL_OUTPUT_NAME}"
    ARCHIVE_OUTPUT_NAME "${REAL_OUTPUT_NAME}"
    LIBRARY_OUTPUT_NAME "${REAL_OUTPUT_NAME}"
    RUNTIME_OUTPUT_NAME "${REAL_OUTPUT_NAME}"
  )

  # The compatibility name is implemented as a copy on Windows, where symlink
  # creation is not always available to normal users.  On Unix-like systems it
  # is implemented as a relative symlink to avoid duplicating the library file.
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

# Generate a small cleanup script for the post-build step below.  Some sibling
# dependencies may still produce artifacts during a local build; this script
# removes those public-lib duplicates from Astro's lib/lib directory so that the
# directory exposes only Astro's library and compatibility alias.
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

# Run the cleanup script after the Astro library has been produced.  The script
# receives the project directory and the current platform suffix as CMake -D
# variables, which keeps the generated script independent from configure-time
# absolute values except for its own location.
add_custom_command(
  TARGET Astro POST_BUILD
  COMMAND ${CMAKE_COMMAND}
          -DPROJECT_DIR=${CMAKE_CURRENT_SOURCE_DIR}
          -DPLATFORM_SUFFIX=${ASTRO_PLATFORM_SUFFIX}
          -P "${CMAKE_CURRENT_BINARY_DIR}/AstroCleanPublicLibDir.cmake"
  VERBATIM
)

# --------------------------------------------------------------------------
# Test executables
# --------------------------------------------------------------------------
# Tests are built only when BUILD_TESTING is enabled and Astro is the top-level
# project.  This prevents Astro's tests from being added to a parent project's
# test suite when Astro is consumed through FetchContent.
include(CTest)
if(BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
  # Explicit list of Astro test programs.  Keeping the list explicit avoids
# accidentally compiling scratch files or incomplete experiments from src_tests/.
  set(TEST_NAMES
    test1-Lambert
    test2-DeltaV
    test3-DeltaV
  )

  add_custom_target("${PROJECT_NAME}_all_tests" ALL)

  foreach(S ${TEST_NAMES})
    add_executable(${S} "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc")
    target_link_libraries(${S} PRIVATE Astro::Astro ${CMAKE_DL_LIBS})
    # Add src/ so tests can include internal headers that are not part of the
# installed/public include tree.
    target_include_directories(${S} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    set_target_properties(
      ${S} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    )
    add_test(NAME "${S}" COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/bin/${S}" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_dependencies("${PROJECT_NAME}_all_tests" ${S})
  endforeach()
endif()

# --------------------------------------------------------------------------
# Install rules
# --------------------------------------------------------------------------
# Install only the Astro target and public include tree.  Dependency install
# rules are disabled above, and the explicit cleanup removes stale dependency
# libraries that older local build layouts may have left in Astro's public
# lib/lib directory.
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

  # Optional source/binary package generation through CPack.
  set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
  set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
  set(CPACK_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;/build.*/")
  include(CPack)
endif()

# --------------------------------------------------------------------------
# Configuration summary
# --------------------------------------------------------------------------
# Print the resolved dependency sources at the end of configuration.  This is
# intentionally concise and focuses on the information that is most useful when
# diagnosing whether a build used local checkouts or FetchContent downloads.
message(STATUS "")
message(STATUS "Astro dependency summary:")
message(STATUS "  Apple FAT build   : ${UTILS_BUILD_FAT}")
if(APPLE)
  message(STATUS "  Apple archs       : ${CMAKE_OSX_ARCHITECTURES}")
endif()
message(STATUS "  GenericContainer  : ${GenericContainer_SOURCE}")
message(STATUS "  UtilsLite         : ${UtilsLite_SOURCE}")
message(STATUS "  quarticRootsFlocke: ${quarticRootsFlocke_SOURCE}")
message(STATUS "")