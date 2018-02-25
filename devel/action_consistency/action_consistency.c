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
#include "macros.h"
#include "matop.h"
#include "action.h"
#include "update.h"

#define REP (10000)
#define ACTION_B P_actionD4D2_b
#define HERM 1
#define HERM_STEP 100 
#define SC 1
#define ITER 100

int main()
{
    // random number generator initialization
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // file
    FILE* fS;
    FILE* f2;
    FILE* f3;
    if(HERM)
    {
        fS = fopen("D4D2t13H.txt", "w");
        f2 = fopen("D4D2t13H_time.txt", "w");
        f3 = fopen("D4D2t13H_err.txt", "w");
    }
    else
    {
        fS = fopen("D4D2t13.txt", "w");
        f2 = fopen("D4D2t13_time.txt", "w");
        f2 = fopen("D4D2t13_err.txt", "w");
    }

    init_data();
    GEOM_CHECK();
    init_cold(ACTION_B, P_gamma);
    
    clock_t start1;
    clock_t cumul1 = 0;
    clock_t start2;
    clock_t cumul2 = 0;
    clock_t start3;
    clock_t cumul3 = 0;
    // simulation
    double dS1 = S; 
    //double dS2 = S; 
    for(int idx=0; idx<REP; idx++)
    {
        if(!(idx % ITER))
            printf("iteration: %d\n", idx);
        if(!(idx % HERM_STEP))
            adjust();

        // decide what matrix gets updated
        int uM = (int)(nHL*gsl_rng_uniform(r));
        while(uM == nHL)
            uM = (int)(nHL*gsl_rng_uniform(r));

        gsl_matrix_complex* dM = gsl_matrix_complex_calloc(dim, dim);
        //generate_HL_1(dM, dim, r);
        
        // decide matrix entry that gets updated
        int i = (int)(dim*gsl_rng_uniform(r));
        while(i == dim)
            i = (int)(dim*gsl_rng_uniform(r));

        int j = (int)(dim*gsl_rng_uniform(r));
        while(j == dim)
            j = (int)(dim*gsl_rng_uniform(r));


        double x=0;
        double y=0;
        // generate random entry
        if(i==j)
        {
            // generate x uniformly between -SC and SC
            x = -1 + 2*gsl_rng_uniform(r);
            x*=SC;
            
            gsl_matrix_complex_set(dM, i, j, gsl_complex_rect(2.*x, 0.));
        }
        else
        {
            // generate x uniformly between -SC and SC
            x = -1 + 2*gsl_rng_uniform(r);
            y = -1 + 2*gsl_rng_uniform(r);
            x*=SC;
            y*=SC;

            gsl_matrix_complex_set(dM, i, j, gsl_complex_rect(x, y));
            gsl_matrix_complex_set(dM, j, i, gsl_complex_rect(x, -y));
        }
        

        start1 = clock();
        dS1 += G*delta2(uM, i, j, gsl_complex_rect(x, y)) + delta4(uM, i, j, gsl_complex_rect(x, y));
        cumul1 += clock()-start1;
        
        gsl_matrix_complex_add(MAT[uM], dM);
        gsl_matrix_complex_free(dM);
        
        start2 = clock();
        double rS = G*dirac2() + dirac4();
        cumul2 += clock()-start2;
        
        if(!(idx%1))
        {
            start3 = clock();
            double cS = ACTION_B();
            cumul3 += clock()-start3;

            fprintf(fS, "%.15lf %.15lf %.15lf\n", dS1, rS, cS);
            fprintf(f3, "%.15lf %.15lf\n", fabs(dS1-cS)/cS, fabs(rS-cS)/cS);
        }
    
    }

    fprintf(f2, "nH = %d, nL = %d, REP = %d, n = %d:\n", nH, nL, REP, dim);
    fprintf(f2, "decomp time 1 %lf\n", (double)cumul1/(double)CLOCKS_PER_SEC);
    fprintf(f2, "decomp time 2 %lf\n", (double)cumul2/(double)CLOCKS_PER_SEC);
    fprintf(f2, "brutef time %lf\n", (double)cumul3/(double)CLOCKS_PER_SEC);

    // free random number generator
    gsl_rng_free(r);
    fclose(fS);
    fclose(f2);
    fclose(f3);
}
