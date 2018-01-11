#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include "time.h"
#include "global.h"
#include "macros.h"
#include "matop.h"
#include "action.h"
#include "update.h"

#define REP 500
#define ACTION_B P_actionD2_b

int main()
{
    // random number generator initialization
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // file
    FILE* fS = fopen("D2t40.txt", "w");
    FILE* f2 = fopen("D2t40_time.txt", "w");

    init_data();
    GEOM_CHECK();
    init_cold(ACTION_B, P_gamma);
    fprintf(fS, "%.15lf\n", S);
    
    clock_t start1;
    clock_t cumul1 = 0;
    clock_t start2;
    clock_t cumul2 = 0;
    // simulation
    double Sdelta; 
    double dS = S; 
    for(int i=0; i<REP; i++)
    {

        // decide what matrix gets updated
        int uM = (int)(nHL*gsl_rng_uniform(r));
        while(uM == nHL)
            uM = (int)(nHL*gsl_rng_uniform(r));

        gsl_matrix_complex* dM = gsl_matrix_complex_alloc(dim, dim);
        if(uM < nH)
            generate_HL(dM, 0, dim, r);  
        else
            generate_HL(dM, 1, dim, r);  


        start1 = clock();
        Sdelta = delta2(dM, uM);
        dS += Sdelta;
        cumul1 += clock()-start1;
        printf("%lf %lf\n", dS, Sdelta);
        
        gsl_matrix_complex_add(MAT[uM], dM);
        free(dM);

        start2 = clock();
        gsl_complex cS = ACTION_B();
        cumul2 += clock()-start2;
        
        fprintf(fS, "%.15lf %.15lf\n", dS, GSL_REAL(cS));
    }

    fprintf(f2, "nH = %d, nL = %d, REP = %d, n = %d:\n", nH, nL, REP, dim);
    fprintf(f2, "decomp time %lf\n", (double)cumul1/(double)CLOCKS_PER_SEC);
    fprintf(f2, "brutef time %lf\n", (double)cumul2/(double)CLOCKS_PER_SEC);

    // free random number generator
    gsl_rng_free(r);
    fclose(fS);
    fclose(f2);
}
