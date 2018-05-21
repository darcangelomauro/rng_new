#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "global.h"
#include "update.h"
#include "codeop.h"
#include "dirac.h"
#include "action.h"

int main(int argc, char** argv)
{
    // ****************************************** HANDLE ARGUMENTS ******************************
    if(argc != 3)
    {
        printf("Error: Not enough or too many arguments. The arguments must be [rank] and [codename]\n");
        exit(EXIT_FAILURE);
    }
    int rank = strtol(argv[1], NULL, 10);


    // ****************************************** RNG STUFF ******************************
    // Initialize random number generator
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL) + rank); // SUBSTITUTE RANK WITH ARRAY INDEX



    // ****************************************** READ GVAL FILE ******************************
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
    




    // ****************************************** SIMULATION ******************************
    char* varG_code = alloc_rank_filename2(rank, argv[2]);
    multicode_wrapper(dirac42, delta42, incr_G, rep_G, rank, varG_code, r);
    


    // ****************************************** FREE MEMORY ******************************
    // Free codename
    free(varG_code);
    // Free random number generator
    gsl_rng_free(r);
}
