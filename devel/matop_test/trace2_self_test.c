#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include "matop.h"
#include "time.h"
#include "math.h"

#define N 1000
#define REP 10

int main()
{
    // random number generator initialization
    const gsl_rng_type* T;
    gsl_rng* r;
    T = gsl_rng_ranlxd1;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));
    
    clock_t start1;
    clock_t cumul1 = 0;
    clock_t start2;
    clock_t cumul2 = 0;
    for(int i=0; i<REP; i++)
    {
        gsl_matrix_complex* m = gsl_matrix_complex_alloc(N, N);
        gsl_matrix_complex* mm = gsl_matrix_complex_alloc(N, N);
        generate_HL(m, 0, N, r);
        
        start1 = clock();
        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, m, m, GSL_COMPLEX_ZERO, mm);
        double tr2_bf = GSL_REAL(trace(mm));
        cumul1 += clock()-start1;

        start2 = clock();
        double tr2 = trace2_self(m);
        cumul2 += clock()-start2;
        double diff = fabs(tr2_bf-tr2);
        if(diff > 1e-12)
        {
            printf("iteration: %d\n", i);
            printf("trace brutf: %.15lf\n", tr2_bf);
            printf("trace smart: %.15lf\n", tr2);
            printf("difference: %.15lf\n", tr2_bf-tr2);
        }

        gsl_matrix_complex_free(m);
        gsl_matrix_complex_free(mm);
    }

    printf("time brutf %lf\n", (double)cumul1/(double)CLOCKS_PER_SEC);
    printf("time smart %lf\n", (double)cumul2/(double)CLOCKS_PER_SEC);
}





    
    

