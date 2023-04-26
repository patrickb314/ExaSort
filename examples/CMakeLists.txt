include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

blt_add_executable( NAME       exasortbench
                    SOURCES    benchmark.cpp 
                    DEPENDS_ON exasort)
