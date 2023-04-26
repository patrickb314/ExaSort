
#configure_file(Beatnik_Config.hpp.in Beatnik_Config.hpp)

set(HEADERS
  ExaSort.hpp
  KokkosSorter.hpp
  BitonicSorter.hpp
# MultiBitonicSorter.hpp
# SampleSorter.hpp
)

set(DEPENDS_ON
  Kokkos::Kokkos
  mpi
  )
 
set(INCLUDES
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
  )
       
blt_add_library( NAME exasort
                 HEADERS ${HEADERS}
                 SOURCES ${SOURCES}
                 INCLUDES ${INCLUDES}
                 DEPENDS_ON ${DEPENDS_ON} 
  )
