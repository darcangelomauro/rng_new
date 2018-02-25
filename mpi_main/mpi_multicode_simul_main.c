#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include "time.h"
#include "global.h"
#include "macros.h"
#include "action.h"
#include "update.h"

int main()
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    


    // Initialize random number generator
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL) + rank);



    // Read G values from file
    char* name_gval = alloc_rank_filename(rank, "gval");  
    FILE* fG = fopen(name_gval, "r");
    if(fG == NULL)
    {
        printf("Error: %s not found\n", name_gval);
        exit(EXIT_FAILURE);
    }
    free(name_gval);
    
    double incr_G;
    int rep_G;
    int r1 = 0;

    r1 += fscanf(fG, "%lf", &incr_G);
    r1 += fscanf(fG, "%d", &rep_G);
    
    if(r1 != 2)
    {
        printf("Error: failed to read parameters\n");
        exit(EXIT_FAILURE);
    }

    fclose(fG);
    


    // Simulation
    multicode_wrapper(dirac42, delta42, incr_G, rep_G, rank, r);
    


    // Free random number generator
    gsl_rng_free(r);
    
    // Finalize the MPI environment.
    MPI_Finalize();
}
