#define GEOM21_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include "geom.h"
#include "global.h"

// TYPE (2,1) GEOMETRY

void geom_check21()
{
    if(nH != 3 || nL != 1)
    {
        printf("Error: geometry is not (2,1)\n");
        exit(EXIT_FAILURE);
    }
}

// build gamma matrices
void init_gamma21()
{
    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(dimG, dimG);

    // gamma0
    gsl_matrix_complex_set(gamma[0], 0, 0, gsl_complex_rect(1.0, 0.0));
    gsl_matrix_complex_set(gamma[0], 1, 1, gsl_complex_rect(-1.0, 0.0));
     
    // gamma1
    gsl_matrix_complex_set(gamma[1], 0, 1, gsl_complex_rect(1.0, 0.0));
    gsl_matrix_complex_set(gamma[1], 1, 0, gsl_complex_rect(1.0, 0.0));
     
    // gamma2
    gsl_matrix_complex_set(gamma[2], 0, 0, gsl_complex_rect(-1.0, 0.0));
    gsl_matrix_complex_set(gamma[2], 1, 1, gsl_complex_rect(-1.0, 0.0));
     
    // gamma3
    gsl_matrix_complex_set(gamma[3], 0, 1, gsl_complex_rect(0.0, 1.0));
    gsl_matrix_complex_set(gamma[3], 1, 0, gsl_complex_rect(0.0, -1.0));
}

