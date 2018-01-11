#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include "matop.h"

#define NUM 4

int main()
{
    gsl_matrix_complex* A = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* B = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* C = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* D = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* res = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* res1 = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* res2 = gsl_matrix_complex_alloc(3, 3);
    

    gsl_matrix_complex_set(A, 0, 0, gsl_complex_rect(5., 2.));
    gsl_matrix_complex_set(A, 0, 1, gsl_complex_rect(6., 2.));
    gsl_matrix_complex_set(A, 0, 2, gsl_complex_rect(6., 2.));
    gsl_matrix_complex_set(A, 1, 0, gsl_complex_rect(8., 4.));
    gsl_matrix_complex_set(A, 1, 1, gsl_complex_rect(2., 6.));
    gsl_matrix_complex_set(A, 1, 2, gsl_complex_rect(4., 2.));
    gsl_matrix_complex_set(A, 2, 0, gsl_complex_rect(2, 2.));
    gsl_matrix_complex_set(A, 2, 1, gsl_complex_rect(7., 8.));
    gsl_matrix_complex_set(A, 2, 2, gsl_complex_rect(9., 5.));

    gsl_matrix_complex_set(B, 0, 0, gsl_complex_rect(3., 2.));
    gsl_matrix_complex_set(B, 0, 2, gsl_complex_rect(3., 2.));
    gsl_matrix_complex_set(B, 1, 0, gsl_complex_rect(4., 4.));
    gsl_matrix_complex_set(B, 1, 1, gsl_complex_rect(24., 5.));
    gsl_matrix_complex_set(B, 1, 2, gsl_complex_rect(4., 6.));
    gsl_matrix_complex_set(B, 2, 0, gsl_complex_rect(5, 6.));
    gsl_matrix_complex_set(B, 2, 1, gsl_complex_rect(4., 8.));
    gsl_matrix_complex_set(B, 2, 2, gsl_complex_rect(3., 9.));

    gsl_matrix_complex_set(C, 0, 0, gsl_complex_rect(1., 8.));
    gsl_matrix_complex_set(C, 0, 2, gsl_complex_rect(6., 7.));
    gsl_matrix_complex_set(C, 1, 0, gsl_complex_rect(2., 5.));
    gsl_matrix_complex_set(C, 1, 1, gsl_complex_rect(0., 4.));
    gsl_matrix_complex_set(C, 1, 2, gsl_complex_rect(5., 6.));
    gsl_matrix_complex_set(C, 2, 0, gsl_complex_rect(6, 6.));
    gsl_matrix_complex_set(C, 2, 1, gsl_complex_rect(6., 4.));
    gsl_matrix_complex_set(C, 2, 2, gsl_complex_rect(8., 3.));

    gsl_matrix_complex_set(D, 0, 0, gsl_complex_rect(1., 2.));
    gsl_matrix_complex_set(D, 0, 2, gsl_complex_rect(3., 2.));
    gsl_matrix_complex_set(D, 1, 0, gsl_complex_rect(3., 6.));
    gsl_matrix_complex_set(D, 1, 1, gsl_complex_rect(5., 6.));
    gsl_matrix_complex_set(D, 1, 2, gsl_complex_rect(5., 2.));
    gsl_matrix_complex_set(D, 2, 0, gsl_complex_rect(9., 5.));
    gsl_matrix_complex_set(D, 2, 1, gsl_complex_rect(6., 8.));
    gsl_matrix_complex_set(D, 2, 2, gsl_complex_rect(6., 9.));

    gsl_matrix_complex** group = malloc(NUM*sizeof(gsl_matrix_complex*));
    group[0] = A;
    group[1] = B;
    group[2] = C;
    group[3] = D;

    matrix_multiprod(group, NUM, res);
    printmat_complex(res, 3, 3);
    
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, A, B, GSL_COMPLEX_ZERO, res);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, C, D, GSL_COMPLEX_ZERO, res1);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, res, res1, GSL_COMPLEX_ZERO, res2);
    printmat_complex(res2, 3, 3);
    
}

