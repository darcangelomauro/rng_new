#define DIRAC_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include "matop.h"
#include "global.h"


// builds dirac operator using bruteforce (it's geometry independent)
void build_dirac()
{
    gsl_matrix_complex* a = gsl_matrix_complex_alloc(dim*dim, dim*dim);
    gsl_matrix_complex* b = gsl_matrix_complex_alloc(dimD, dimD);
    gsl_matrix_complex_set_zero(DIRAC);
    
    for(int i=0; i<nHL; i++)
    {
        if(i < nH)
        {
            anticommutator(MAT[i], a);
            tensor(gamma[i], a, b);
            gsl_matrix_complex_add(DIRAC, b);
        }
        else
        {
            commutator(MAT[i], a);
            tensor(gamma[i], a, b);
            gsl_matrix_complex_add(DIRAC, b);
        }
    }

    gsl_matrix_complex_free(a);
    gsl_matrix_complex_free(b);
}

double actionD2_bruteforce()
{
    build_dirac();
    
    gsl_matrix_complex* D2 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, DIRAC, DIRAC, GSL_COMPLEX_ZERO, D2);

    double trD2 = trace_herm(D2);

    gsl_matrix_complex_free(D2);

    return trD2;
}

double actionD4_bruteforce()
{
    build_dirac();
    
    gsl_matrix_complex* D2 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_matrix_complex* D4 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, DIRAC, DIRAC, GSL_COMPLEX_ZERO, D2);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, D2, D2, GSL_COMPLEX_ZERO, D4);

    double trD4 = trace_herm(D4);

    gsl_matrix_complex_free(D2);
    gsl_matrix_complex_free(D4);

    return trD4;
}

double actionD4D2_bruteforce()
{
    build_dirac();
    
    gsl_matrix_complex* D2 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_matrix_complex* D4 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, DIRAC, DIRAC, GSL_COMPLEX_ZERO, D2);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, D2, D2, GSL_COMPLEX_ZERO, D4);

    double trD2 = trace_herm(D2);
    double trD4 = trace_herm(D4);

    gsl_matrix_complex_free(D2);
    gsl_matrix_complex_free(D4);

    return trD2*G + trD4;
}


