cmake_minimum_required (VERSION 2.8.12)
project(robust_hex_dominant_meshing)

### Compilation flags: adapt to your needs ###
if(MSVC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj /w") ### Enable parallel compilation
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2") #### Libigl requires a modern C++ compiler that supports c++11
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "./" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations -Wno-unused-parameter -Wno-deprecated-register -Wno-return-type-c-linkage")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated-declarations -Wno-unused-parameter -Wno-deprecated-register -Wno-return-type-c-linkage")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Os")
endif()



#set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6)
#find_package(PythonLibs REQUIRED)
#find_package(PythonInterp REQUIRED)
#include_directories(${PYTHON_INCLUDE_DIR} include)

add_subdirectory(ext/pybind11)


if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  # set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    # "MinSizeRel" "RelWithDebInfo")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

if(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ext/nanoguig/ext/glfw")
  message(FATAL_ERROR "The robust_hex_dominant_meshing dependency repositories (NanoGUI, GLFW, etc.) are missing! "
    "You probably did not clone the project with --recursive. It is possible to recover "
    "by calling \"git submodule update --init --recursive\"")
endif()

if(APPLE)
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
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-gnu-anonymous-struct -Wno-c99-extensions -Wno-nested-anon-types")
  endif()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

# add_definitions(-DSINGLE_PRECISION)
# #add_definitions(-DT_VTAG)

# add_library(robust_hex_dominant_meshing_lib
#   src/aabb.h
#   src/adjacency.h
#   src/common.h
#   src/quat.h
#   src/timer.h
#   src/hierarchy.cpp src/hierarchy.h
#   src/orientations.h
# )


# target_link_libraries(robust_hex_dominant_meshing robust_hex_dominant_meshing_lib)

# # Quench some warnings on MSVC
# if (APPLE)
# # Build an application bundle on OSX
# else()
# # Insulate from a few types of ABI changes by statically linking against libgcc and libstdc++
# set_target_properties(robust_hex_dominant_meshing PROPERTIES LINK_FLAGS "-static-libgcc")
# #set_target_properties(robust_hex_dominant_meshing PROPERTIES LINK_FLAGS "-static-libgcc -static-libstdc++")
# endif()

# if (UNIX AND NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
# add_custom_command(TARGET robust_hex_dominant_meshing POST_BUILD COMMAND strip $<TARGET_FILE:robust_hex_dominant_meshing>)
# endif()



option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR      "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_COMISO           "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK             "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE           "Use Embree"         OFF)
option(LIBIGL_WITH_LIM              "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          OFF)
option(LIBIGL_WITH_NANOGUI          "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW      "Use GLFW"           OFF)
option(LIBIGL_WITH_PNG              "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON           "Use Python"         ON)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER           "Use OpenGL viewer"  OFF)
option(LIBIGL_WITH_XML              "Use XML"            OFF)
option(LIBIGL_WITH_LIM              "Use LIM"            OFF)
option(ROSY_WITH_REMESH "Add CGAL and remeshing" ON)

if(ROSY_WITH_REMESH)
    set(LIBIGL_WITH_CGAL ON)
endif()

find_package(LIBIGL REQUIRED QUIET)
# Header file directories of dependency libraries
include_directories(
  # Eigen
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/nanoguig/ext/eigen
  # Pseudorandom number generator
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/pcg32
  # libIGL
)
include_directories(${PYTHON_INCLUDE_DIR} include)
pybind11_add_module(rosy 
    src/bind.cpp 
    src/main.cpp 
    src/hierarchy.cpp)
#target_link_libraries(rosy ${PYTHON_LIBRARIES})
target_link_libraries(rosy PUBLIC igl::core)

if(ROSY_WITH_REMESH)
target_link_libraries(rosy PUBLIC igl::cgal)
target_compile_definitions(rosy PRIVATE REMESHER)
endif()

target_compile_definitions(rosy PRIVATE PYBIND)

#add_executable(cmdln
#   src/main.cpp
#   src/hierarchy.cpp
#)
