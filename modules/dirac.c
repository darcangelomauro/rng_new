#define DIRAC_C

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include "update.h"
#include "matop.h"
#include "math.h"
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

gsl_complex actionD2_bruteforce()
{
    build_dirac();
    
    gsl_matrix_complex* D2 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, DIRAC, DIRAC, GSL_COMPLEX_ZERO, D2);

    gsl_complex trD2 = trace(D2);

    gsl_matrix_complex_free(D2);

    return trD2;
}

gsl_complex actionD4D2_bruteforce()
{
    build_dirac();
    
    gsl_matrix_complex* D2 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_matrix_complex* D4 = gsl_matrix_complex_calloc(dimD, dimD);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, DIRAC, DIRAC, GSL_COMPLEX_ZERO, D2);
    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, D2, D2, GSL_COMPLEX_ZERO, D4);

    gsl_complex trD2 = trace(D2);
    gsl_complex trD4 = trace(D4);

    gsl_matrix_complex_free(D2);
    gsl_matrix_complex_free(D4);

    return gsl_complex_add(gsl_complex_mul_real(trD2, G), trD4);
}


