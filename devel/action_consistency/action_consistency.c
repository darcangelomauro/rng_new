#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_complex.h>
#include "time.h"
#include "math.h"
#include "global.h"
#include "matop.h"
#include "action.h"
#include "update.h"
#include "dirac.h"
#include "geom.h"

#define REP (10)
#define ITER (1)
#define ACTION_B actionD4D2_bruteforce
#define TYPE D4D2t20
#define MOD0
#define MOD1 _err
#define PASTER_NAME(x, y, z) x ## y ## z
#define FILENAME(x, y, z) PASTER_NAME(x, y, z)
#define STR(x) STR1(x)
#define STR1(x) #x

int main()
{
    // random number generator initialization
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // file
    FILE* fS = fopen(STR(FILENAME(TYPE, MOD0, MOD0)), "w");
    FILE* ferr = fopen(STR(FILENAME(TYPE, MOD1, MOD0)), "w");

    
    // simulation
    for(int idx=0; idx<REP; idx++)
    {
        if(!(idx%ITER))
            printf("Iteration %d\n", idx);

        // Initialize
        init_data("init.txt");
        FUNCTION(geom_check, CLIFF_P, CLIFF_Q)();
        init_cold(FUNCTION(dirac, 42, MOD0));
        
        // Simulation
        double ar = 0.;
        for(int i=0; i<Nsw; i++)
        {
            ar += sweep(FUNCTION(delta, 42, MOD0), r);
            
            if(doit(i, ADJ))
                adjust();
        }

        // Comparison
        double S_b = ACTION_B();
        fprintf(fS, "%.15lf %.15lf\n", S, S_b);
        fprintf(ferr, "%.15lf\n", fabs(S-S_b)/S_b);
        
        // Free
        simulation_free();
    }

    // free stuff 
    gsl_rng_free(r);
    fclose(fS);
    fclose(ferr);
}
