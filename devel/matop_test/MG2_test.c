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

#define N 100
#define REP 1000

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
        gsl_matrix_complex* m1 = gsl_matrix_complex_alloc(N, N);
        gsl_matrix_complex* m2 = gsl_matrix_complex_alloc(N, N);
        gsl_matrix_complex* mm = gsl_matrix_complex_alloc(N, N);
        generate_HL(m1, 0, N, r);
        generate_HL(m2, 0, N, r);
        
        // decide matrix entry that gets taken
        int i = (int)(N*gsl_rng_uniform(r));
        while(i == N)
            i = (int)(N*gsl_rng_uniform(r));

        int j = (int)(N*gsl_rng_uniform(r));
        while(j == N)
            j = (int)(N*gsl_rng_uniform(r));
        
        start1 = clock();
        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, m1, m2, GSL_COMPLEX_ZERO, mm);
        gsl_complex el2_bf = gsl_matrix_complex_get(mm, i, j);
        cumul1 += clock()-start1;

        start2 = clock();
        gsl_complex el2 = MG2(m1, m2, i, j);
        cumul2 += clock()-start2;

        double diff_re = fabs(GSL_REAL(el2_bf)-GSL_REAL(el2));
        double diff_im = fabs(GSL_IMAG(el2_bf)-GSL_IMAG(el2));
        if(diff_re > 1e-9 || diff_im > 1e-9)
        {
            printf("iteration: %d\n", i);
            printf("element brutf: (%.15lf, %.15lf)\n", GSL_REAL(el2_bf), GSL_IMAG(el2_bf));
            printf("element smart: (%.15lf, %.15lf)\n", GSL_REAL(el2), GSL_IMAG(el2));
            printf("difference: (%.15lf, %.15lf)\n", diff_re, diff_im);
        }

        gsl_matrix_complex_free(m1);
        gsl_matrix_complex_free(m2);
        gsl_matrix_complex_free(mm);
    }

    printf("time brutf %lf\n", (double)cumul1/(double)CLOCKS_PER_SEC);
    printf("time smart %lf\n", (double)cumul2/(double)CLOCKS_PER_SEC);
}





    
    

