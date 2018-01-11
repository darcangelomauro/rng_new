#define ACTION_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_double.h>
#include "update.h"
#include "matop.h"
#include "fileop.h"
#include "math.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include "global.h"
#include "macros.h"
#include "data.h"
#include "statistics.h"


double delta2(gsl_matrix_complex* dM, int uM)
{

    gsl_complex res = GSL_COMPLEX_ZERO;


    // CASE: s=1

    for(int i=0; i<nHL; i++)
    {
        // trace of the clifford part
        gsl_matrix_complex* m = gsl_matrix_complex_alloc(dimG, dimG);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, gamma[uM], gamma[i], GSL_COMPLEX_ZERO, m);
        gsl_complex A = trace(m);
        A = gsl_complex_mul_real(A, 4.);
        
        gsl_matrix_complex_free(m);

        if(GSL_REAL(A) != 0 || GSL_IMAG(A) != 0)
        {
            gsl_complex B;
            gsl_complex C;
            B = gsl_complex_mul(trace(dM), trace(MAT[i]));

            gsl_matrix_complex* m1 = gsl_matrix_complex_alloc(dim, dim);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, dM, MAT[i], GSL_COMPLEX_ZERO, m1);
            C = trace(m1);
            gsl_matrix_complex_free(m1);

            if(uM < nH && i < nH)
                C = gsl_complex_mul_real(C, (double)dim);
            
            else if(uM >= nH && i >= nH)
                C = gsl_complex_mul_real(C, (double)(-dim));

            B = gsl_complex_add(B, C);
            A = gsl_complex_mul(A, B);

            res = gsl_complex_add(res, A);
        }
    }


    // CASE: s=2

    gsl_matrix_complex* m2 = gsl_matrix_complex_alloc(dim, dim);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, dM, dM, GSL_COMPLEX_ZERO, m2);
    gsl_complex A = trace(m2);
    A = gsl_complex_mul_real(A, (double)dim);
    gsl_matrix_complex_free(m2);
    gsl_complex B = trace(dM);
    B = gsl_complex_mul(B, B);
    if(uM < nH)
        A = gsl_complex_add(A, B);
    else
        A = gsl_complex_sub(A, B);
    A = gsl_complex_mul_real(A, 2.*dimG);

    res = gsl_complex_add(res, A);

    //printf("%.15lf %.15lf\n", GSL_REAL(res), GSL_IMAG(res));

    return GSL_REAL(res);
}



