cmake_minimum_required (VERSION 2.8.3)
project(robust_hex_dominant_meshing)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

if(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ext/nanoguig/ext/glfw")
  message(FATAL_ERROR "The robust_hex_dominant_meshing dependency repositories (NanoGUI, GLFW, etc.) are missing! "
    "You probably did not clone the project with --recursive. It is possible to recover "
    "by calling \"git submodule update --init --recursive\"")
endif()

# Enable folders for projects in Visual Studio
if (CMAKE_GENERATOR MATCHES "Visual Studio")
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

# Sanitize build environment for static build with C++11
if (MSVC)
  add_definitions (/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions (/D "__TBB_NO_IMPLICIT_LINKAGE")

  # Parallel build on MSVC (all targets)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

  if (NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")

    # Disable Eigen vectorization for Windows 32 bit builds (issues with unaligned access segfaults)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DEIGEN_DONT_ALIGN")
  endif()

  # Static build
  set(CompilerFlags
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
  foreach(CompilerFlag ${CompilerFlags})
    string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
  endforeach()
elseif(APPLE)
  # Try to auto-detect a suitable SDK
  execute_process(COMMAND bash -c "xcodebuild -version -sdk | grep MacOSX | grep Path | head -n 1 | cut -f 2 -d ' '" OUTPUT_VARIABLE CMAKE_OSX_SYSROOT)
  string(REGEX REPLACE "(\r?\n)+$" "" CMAKE_OSX_SYSROOT "${CMAKE_OSX_SYSROOT}")
  string(REGEX REPLACE "^.*X([0-9.]*).sdk$" "\\1" CMAKE_OSX_DEPLOYMENT_TARGET "${CMAKE_OSX_SYSROOT}")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Enable C++11 mode on GCC / Clang
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

# Compile with compiler warnings
if(MSVC)
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-gnu-anonymous-struct -Wno-c99-extensions -Wno-nested-anon-types")
  endif()
endif()


# Header file directories of dependency libraries
include_directories(
  # Eigen
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/nanoguig/ext/eigen
  # GLFW
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/nanoguig/ext/glfw/include
  # GLEW
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/nanoguig/ext/glew/include
  # Pseudorandom number generator
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/pcg32
  # RPLY mesh I/O library
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/rply
  # pybind11
  #${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11/include
  # Resource file
  ${CMAKE_CURRENT_BINARY_DIR}
)


add_definitions(-DSINGLE_PRECISION)
#add_definitions(-DT_VTAG)

# Resource file (icons etc.)
set(EXTRA_SOURCE "")
if (APPLE)
  set(EXTRA_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/resources/im.icns")
elseif(WIN32)
  set(EXTRA_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/resources/im.rc")
endif()

add_library(robust_hex_dominant_meshing_lib
    src/aabb.h
  src/adjacency.h
    src/common.h
    src/quat.h
  src/timer.h
  src/hierarchy.cpp src/hierarchy.h
    src/orientations.h
    ext/rply/rply.c
)


add_executable(robust_hex_dominant_meshing MACOSX_BUNDLE
src/main.cpp
)

target_link_libraries(robust_hex_dominant_meshing robust_hex_dominant_meshing_lib)

if (WIN32)
# Quench some warnings on MSVC
if (MSVC)
  set_source_files_properties(ext/rply/rply.c PROPERTIES COMPILE_FLAGS "/wd4127")
endif()
elseif (APPLE)
# Build an application bundle on OSX
set_target_properties(robust_hex_dominant_meshing PROPERTIES MACOSX_BUNDLE_BUNDLE_NAME "robust_hex_dominant_meshing")
set_target_properties(robust_hex_dominant_meshing PROPERTIES MACOSX_BUNDLE_BUNDLE_GUI_IDENTIFIER "com.im.robust_hex_dominant_meshing")
set_target_properties(robust_hex_dominant_meshing PROPERTIES MACOSX_BUNDLE_ICON_FILE im.icns)
set_target_properties(robust_hex_dominant_meshing PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/resources/MacOSXBundleInfo.plist.in)
set_source_files_properties(resources/im.icns PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
else()
# Insulate from a few types of ABI changes by statically linking against libgcc and libstdc++
set_target_properties(robust_hex_dominant_meshing PROPERTIES LINK_FLAGS "-static-libgcc")
#set_target_properties(robust_hex_dominant_meshing PROPERTIES LINK_FLAGS "-static-libgcc -static-libstdc++")
endif()

if (UNIX AND NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
add_custom_command(TARGET robust_hex_dominant_meshing POST_BUILD COMMAND strip $<TARGET_FILE:robust_hex_dominant_meshing>)
endif()

#find_package(PythonInterp REQUIRED)
#set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6)
#find_package(PythonLibs REQUIRED)
#find_package(PythonInterp REQUIRED)
#include_directories(${PYTHON_INCLUDE_DIR} include)
#
#
#set(PYBIND11_DIR ext/pybind11 CACHE PATH "Path to pybind11")
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11/ pybind11)
#

#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11)
#pybind11_add_module(rosy src/main.cpp)
#target_link_libraries(rosy robust_hex_dominant_meshing_lib ${PYTHON_LIBRARIES})