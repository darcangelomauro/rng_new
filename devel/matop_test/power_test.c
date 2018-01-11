#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include "matop.h"

#define POW 4
#define REP 10000000

int main()
{
    gsl_matrix_complex* A = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* B = gsl_matrix_complex_alloc(3, 3);

    gsl_matrix_complex_set(A, 0, 0, gsl_complex_rect(5., 2.));
    gsl_matrix_complex_set(A, 0, 1, gsl_complex_rect(6., 2.));
    gsl_matrix_complex_set(A, 0, 2, gsl_complex_rect(6., 2.));
    gsl_matrix_complex_set(A, 1, 0, gsl_complex_rect(8., 4.));
    gsl_matrix_complex_set(A, 1, 1, gsl_complex_rect(2., 6.));
    gsl_matrix_complex_set(A, 1, 2, gsl_complex_rect(4., 2.));
    gsl_matrix_complex_set(A, 2, 0, gsl_complex_rect(2, 2.));
    gsl_matrix_complex_set(A, 2, 1, gsl_complex_rect(7., 8.));
    gsl_matrix_complex_set(A, 2, 2, gsl_complex_rect(9., 5.));


    for(int i=0; i<REP; i++)
        matrix_power(A, POW, B);

    printmat_complex(B, 3, 3);
}

