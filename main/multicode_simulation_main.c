#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "time.h"
#include "global.h"
#include "macros.h"
#include "action.h"
#include "update.h"

int main()
{
    // random number generator initialization
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // Read G values from file
    FILE* fG = fopen("gval.txt", "r");
    if(fG == NULL)
    {
        printf("Error: gval.txt not found\n");
        exit(EXIT_FAILURE);
    }
    
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
    
    
    // simulation
    multicode_wrapper(dirac42, delta42, incr_G, rep_G, r);
    
    // free random number generator
    gsl_rng_free(r);
}
