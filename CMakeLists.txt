# project settings
cmake_minimum_required(VERSION 3.18)

project(ExaSort LANGUAGES C CXX VERSION 0.1.0)

include(GNUInstallDirs)

# find dependencies
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

##---------------------------------------------------------------------------##
## Print the revision number to stdout
##---------------------------------------------------------------------------##
FIND_PACKAGE(Git)
IF(GIT_FOUND AND IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  EXECUTE_PROCESS(
    COMMAND           ${GIT_EXECUTABLE} log --pretty=format:%H -n 1
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE   ExaSort_GIT_COMMIT_HASH
    )
ELSE()
  SET(ExaSort_GIT_COMMIT_HASH "Not a git repository")
ENDIF()
MESSAGE(STATUS "ExaSort Revision = '${ExaSort_GIT_COMMIT_HASH}'")

##---------------------------------------------------------------------------##
## General build options that will be used by various build tools/scripts    ## 
##---------------------------------------------------------------------------##
# Build options that will be picked up by the LLNL BLT tools to help with building 
# and testing
set(ENABLE_MPI ON CACHE BOOL "Build MPI support")
set(ENABLE_CLANG_CUDA "Use Clang's native CUDA support" ON)

# Build options for libraries, disable extras
set(ENABLE_TESTS "Build tests" On)

if (ENABLE_CUDA)
  # Separable compilation is required by ExaSort, set before load BLT
  set(CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "")
  if (NOT DEFINED CUDA_ARCH)
    message(STATUS "CUDA compute architecture set to ExaSort default sm_60 since it was not specified")
    set(CUDA_ARCH "sm_60" CACHE STRING "Set CUDA_ARCH to ExaSort minimum supported" FORCE)
  endif()
endif()

##---------------------------------------------------------------------------##
## Bring up the BLT tools
## We use LLNL BLT to help with building, linking, and testing on various HPC 
## systems. Before bringing in BLT, we set a variety of flags whch tell it the 
## featurs we're interested in.
##---------------------------------------------------------------------------##
if (DEFINED BLT_SOURCE_DIR)
    # Support having a shared BLT outside of the repository if given a BLT_SOURCE_DIR
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
else()
    # Use internal BLT if no BLT_SOURCE_DIR is given
    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/thirdparty/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR
            "The BLT git submodule is not present. "
            "Either run the following two commands in your git repository: \n"
            "    git submodule init\n"
            "    git submodule update\n"
            "Or add -DBLT_SOURCE_DIR=/path/to/blt to your CMake command." )
    endif()
endif()

# Default to C++17 so we can use if constexpr and [[maybe_unused]]
if (NOT BLT_CXX_STD)
    set(BLT_CXX_STD "c++17" CACHE STRING "")
endif()

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

##---------------------------------------------------------------------------##
## Get the dependent packages we need, primarily Kokkos and ClangFormat      ##
## We will later want MPIAdvance here
##---------------------------------------------------------------------------##
# Get Kokkos to for local node parallelism and sorting
find_package(Kokkos 4.0 REQUIRED)

# Get Clang Format to use for making sure the resulting code is
# properly formatted
find_package(CLANG_FORMAT 10)

# library
add_subdirectory(src)

# examples
add_subdirectory(examples)

# tests
add_subdirectory(tests)

# Add a target for formatting the code using Clang
if(CLANG_FORMAT_FOUND)
  file(GLOB_RECURSE FORMAT_SOURCES src/*.cpp src/*.hpp tests/*.hpp tests/*.cpp examples/*.cpp examples/*.hpp)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXECUTABLE} -i -style=file ${FORMAT_SOURCES}
    DEPENDS ${FORMAT_SOURCES})
endif()
