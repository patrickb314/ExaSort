# ExaSort - an exascale distributed sort benchmark

Sorting is a key component of many high-performance computing applications, and a wide range of sorting algorithms have been developed with different performance characteristics. The best sorting approach can thus vary dramatically based on the characteristics of the data being sorted, the performance of the local sort, and system communication characteristics. 

ExaSort is a distributed sorting benchmark designed for two purposes:
  1. Test the performance of different GPU-accelerated distributed sorting 
     algorithms under different data distributions and layouts. 
  1. Drive the development of new communication abstractions that will:
    1. Better support algorithms that require interleaved communication and computation
    1. Provide applications and algorithms communication perfomrance information to allow that to choose the best approach based on system and data characteristics.


## Overview 

ExaSort has the following key components implemented as template classes the provide the relevant entry points and concepts based on the C++ std::ranges::sort interface.
  * Driver - 
  * DistributedSort - 
std::random_access_iterator I, std::sentinel_for<I> S,

Development Plan:

Basic MVP development
  1. Create driver that creates the data to sort with the appropriate parameters of the data.
  2. Implement a first quicksort using the flecsph algorithm to flesh out the interfaces
  3. Add basic check that makes sure the data is sorted appropriately.
  4. Add in a multistep bitonic sort
  5. Add in the flecsph sample sort
  6. Add in the ability to customize the key and value data for sorting

Communication abstraction
  1. Design basic coupled compute/communication sort schedule abstraction.
  1. Design abstraction to provide communication schedules in a sort-specific way.

Sort Modeling
  1. Create simple model of the performance of each sort in terms of LogGOP parameters and other parameters
  1. Create communication abstractions to provide relevant performance data to each sort
  1. Create interface to sort performance for each algorithm 

In addition, ExaSort is designed to drive the development of  on  Crucially, ExaSort seeks to //adapt// the sorting algorithm used based on th different sot=  evaluate the performance of different distributed sorting approaches on GPU-accelerated systems. appro uses Kokkos for the local sort  
