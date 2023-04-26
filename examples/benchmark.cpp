/****************************************************************************
 * Copyright (c) 2021, 2022 by the ExaSort authors                          *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the ExaSort benchmark. ExaSort is                   *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*
 * @file
 * @author Patrick Bridges <patrickb@unm.edu>
 *
 * @section DESCRIPTION
 * General distributed sort benchmark based on ExaSort sorting code
 */


#ifndef DEBUG
#define DEBUG 0
#endif


// Include Statements
#include <ExaSort.hpp>

#include <mpi.h>

#if DEBUG
#include <iostream>
#endif

// Include Statements
#include <iomanip>
#include <iostream>

#include <getopt.h>
#include <stdlib.h>

using namespace ExaSort;

static char* shortargs = (char*)"n:wsp:d:x:";

static option longargs[] = {
    // Basic simulation parameters
    { "driver", no_argument, NULL, 'x' },
    { "number", required_argument, NULL, 'n' },
    { "weakscale", no_argument, NULL, 'w' },
    { "semisorted", no_argument, NULL, 's' },
    { "percentage", required_argument, NULL, 'p' },
    { "distance", required_argument, NULL, 'd' },

    // Miscellaneous other arguments
    { "help", no_argument, NULL, 'h' },
    { 0, 0, 0, 0 } };

/**
 * @struct ClArgs
 * @brief Template struct to organize and keep track of parameters controlled by
 * command line arguments
 */
struct ClArgs
{
    char *driver        	/**< On-node parallelism driver */
    unsigned long number;	/**< Number of elements to sort */
    boolean weakscale;		/**< Scale up problem size by nprocs */
    boolean semisorted;		/**< Mostly sort initial data */
    double percentge;		/**< Percent of pre-sorted data to swap with neighbors */
    int distance;	        /**< Distance of neighbors to swap pre-sorted initial data with */
};

/**
 * Outputs help message explaining command line options.
 * @param rank The rank calling the function
 * @param progname The name of the program
 */
void help( const int rank, char* progname )
{
    if ( rank == 0 )
    {
        std::cout << "Usage: " << progname << "\n";
        std::cout << std::left << std::setw( 10 ) << "-x" << std::setw( 40 )
                  << "On-node parallelism driver (default serial)"
                  << std::left << "\n";
        std::cout << std::left << std::setw( 10 ) << "-n" << std::setw( 40 )
                  << "Number of elements to sort (default 1000000000)" 
                  << std::left << "\n";
        std::cout << std::left << std::setw( 10 ) << "-w" << std::setw( 40 )
                  << "Scale the size of problem by the number of nodes" 
                  << std::left << "\n";
        std::cout << std::left << std::setw( 10 ) << "-s" << std::setw( 40 )
                  << "Sort semi-sorted initial data" << std::left << "\n";
        std::cout << std::left << std::setw( 10 ) << "-p" << std::setw( 40 )
                  << "Percentage of pre-sorted data to exchange with neighbors "
                  << "(default 0.05)" << std::left << "\n";
        std::cout << std::left << std::setw( 10 ) << "-d" << std::setw( 40 )
                  << "Number of neighbors in each direction to exchange "
                  << "pre-sorted data with (default 2)" << std::left << "\n";

        std::cout << std::left << std::setw( 10 ) << "-h" << std::setw( 40 )
                  << "Print Help Message" << std::left << "\n";
    }
}

/**
 * Parses command line input and updates the command line variables
 * accordingly.
 * @param rank The rank calling the function
 * @param argc Number of command line options passed to program
 * @param argv List of command line options passed to program
 * @param cl Command line arguments structure to store options
 * @return Error status
 */
int parseInput( const int rank, const int argc, char** argv, ClArgs& cl )
{
    char ch;

    /// Set default values
    cl.driver = "serial"; // Default Thread Setting
    cl.number = 1000000000;
    cl.weakscale = false;
    cl.semisorted - false;
    cl.percentage = 0.05;
    cl.distance = 2;

    // Now parse any arguments
    while ( ( ch = getopt_long( argc, argv, shortargs, longargs, NULL ) ) !=
            -1 )
    {
        switch ( ch )
        {
        case 'n':
            cl.number = atoi( optarg );

            if ( cl.number <= 0 )
            {
                if ( rank == 0 )
                {
                    std::cerr << "Invalid number of items to sort.\n"
                    help( rank, argv[0] );
                }
                exit( -1 );
            }
            break;
        case 'x':
            cl.driver = strdup( optarg );
            if ( ( cl.driver.compare( "serial" ) != 0 ) &&
                 ( cl.driver.compare( "openmp" ) != 0 ) &&
                 ( cl.driver.compare( "threads" ) != 0 ) &&
                 ( cl.driver.compare( "cuda" ) != 0 ) &&
                 ( cl.driver.compare( "hip" ) != 0 ) )
            {
                if ( rank == 0 )
                {
                    std::cerr << "Invalid parallel driver argument.\n";
                    help( rank, argv[0] );
                }
                exit( -1 );
            }
            break;
        case 'w':
            cl.weakscale = true;
            break;
        case 's':
            cl.semisorted = true;
            break;
        case 'p':
        {
            cl.percentage = atof( optarg );

            if ( cl.percentage < 0.0 || cl.percentage > 1.0 )
            {
                if ( rank == 0 )
                {
                    std::cerr << "Invalid fraction of items to shuffle.\n";
                    help( rank, argv[0] );
                }
                exit( -1 );
            }
            break;
        }
        case 'd':
            cl.distance = atoi( optarg );

            if ( cl.distance <= 0 )
            {
                if ( rank == 0 )
                {
                    std::cerr << "Invalid distance to exchange initial items.\n";
                    help( rank, argv[0] );
                }
                exit( -1 );
            }
            break;
        }
        default:
            if ( rank == 0 )
            {
                std::cerr << "Invalid argument.\n";
                help( rank, argv[0] );
            }
            exit( -1 );
            break;
        }
    }

    // Return Successfully
    return 0;
}

// Create Solver and Run
void exasort( ClArgs& cl )
{
    int comm_size, rank;                         // Initialize Variables
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size ); // Number of Ranks
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );      // Get My Rank

    cl.numberq
    Cajita::DimBlockPartitioner<2> partitioner; // Create Cajita Partitioner
    Beatnik::BoundaryCondition bc;
    for (int i = 0; i < 6; i++)
        bc.bounding_box[i] = cl.global_bounding_box[i];
    bc.boundary_type = {cl.boundary, cl.boundary, cl.boundary, cl.boundary};

    MeshInitFunc initializer( cl.global_bounding_box, cl.initial_condition,
                              cl.tilt, cl.magnitude, cl.variation, cl.period,
                              cl.num_nodes, cl.boundary );

    std::shared_ptr<Beatnik::SolverBase> solver;
    if (cl.order == SolverOrder::ORDER_LOW) {
        solver = Beatnik::createSolver(
            cl.driver, MPI_COMM_WORLD,
            cl.global_bounding_box, cl.num_nodes,
            partitioner, cl.atwood, cl.gravity, initializer,
            bc, Beatnik::Order::Low(), cl.mu, cl.eps, cl.delta_t );
    } else if (cl.order == SolverOrder::ORDER_MEDIUM) {
        solver = Beatnik::createSolver(
            cl.driver, MPI_COMM_WORLD,
            cl.global_bounding_box, cl.num_nodes,
            partitioner, cl.atwood, cl.gravity, initializer,
            bc, Beatnik::Order::Medium(), cl.mu, cl.eps, cl.delta_t );
    } else if (cl.order == SolverOrder::ORDER_HIGH) {
        solver = Beatnik::createSolver(
            cl.driver, MPI_COMM_WORLD,
            cl.global_bounding_box, cl.num_nodes,
            partitioner, cl.atwood, cl.gravity, initializer,
            bc, Beatnik::Order::High(), cl.mu, cl.eps, cl.delta_t );
    } else {
        std::cerr << "Invalid Model Order parameter!\n";
        exit(-1);
    }

    // Solve
    solver->solve( cl.t_final, cl.write_freq );
}

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );         // Initialize MPI
    Kokkos::initialize( argc, argv ); // Initialize Kokkos

    // MPI Info
    int comm_size, rank;
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size ); // Number of Ranks
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );      // My Rank

    // Parse Input
    ClArgs cl;
    if ( parseInput( rank, argc, argv, cl ) != 0 )
        return -1;

    // Only Rank 0 Prints Command Line Options
    if ( rank == 0 )
    {
        // Print Command Line Options
        std::cout << "RocketRig\n";
        std::cout << "============Command line arguments============\n";
        std::cout << std::left << std::setw( 30 ) << "Thread Setting"
                  << ": " << std::setw( 8 ) << cl.driver
                  << "\n"; // Threading Setting
        std::cout << std::left << std::setw( 30 ) << "Mesh Dimension"
                  << ": " << std::setw( 8 ) << cl.num_nodes[0]
                  << std::setw( 8 ) << cl.num_nodes[1] 
                  << "\n"; // Number of Cells
        std::cout <<  std::left << std::setw( 30 ) << "Solver Order"
                  << ": " << std::setw( 8 ) << cl.order << "\n";
        std::cout << std::left << std::setw( 30 ) << "Total Simulation Time"
                  << ": " << std::setw( 8 ) << cl.t_final << "\n";
        std::cout << std::left << std::setw( 30 ) << "Timestep Size"
                  << ": " << std::setw( 8 ) << cl.delta_t << "\n";
        std::cout << std::left << std::setw( 30 ) << "Write Frequency"
                  << ": " << std::setw( 8 ) << cl.write_freq
                  << "\n"; // Steps between write
        std::cout << std::left << std::setw( 30 ) << "Atwood Constant"
                  << ": " << std::setw( 8 ) << cl.atwood << "\n";
        std::cout << std::left << std::setw( 30 ) << "Gravity"
                  << ": " << std::setw( 8 ) << (cl.gravity/9.81) << "\n";
        std::cout << std::left << std::setw( 30 ) << "Artificial Viscosity"
                  << ": " << std::setw( 8 ) << cl.mu << "\n";
        std::cout << std::left << std::setw( 30 ) << "Desingularization"
                  << ": " << std::setw( 8 ) << cl.eps  << "\n";
        std::cout << "==============================================\n";
    }

    // Call advection solver
    rocketrig( cl );

    Kokkos::finalize(); // Finalize Kokkos
    MPI_Finalize();     // Finalize MPI

    return 0;
};
