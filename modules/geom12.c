#define GEOM12_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include "geom.h"
#include "global.h"

// TYPE (1,2) GEOMETRY

void geom_check12()
{
    if(nH != 1 || nL != 3)
    {
        printf("Error: geometry is not (1,2)\n");
        exit(EXIT_FAILURE);
    }
}

// build gamma matrices
void init_gamma12()
{
    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(dimG, dimG);

    // gamma0
    gsl_matrix_complex_set(gamma[0], 0, 1, gsl_complex_rect(1.0, 0.0));
    gsl_matrix_complex_set(gamma[0], 1, 0, gsl_complex_rect(1.0, 0.0));
     
    // gamma1
    gsl_matrix_complex_set(gamma[1], 0, 0, gsl_complex_rect(-1.0, 0.0));
    gsl_matrix_complex_set(gamma[1], 1, 1, gsl_complex_rect(1.0, 0.0));
     
    // gamma2
    gsl_matrix_complex_set(gamma[2], 0, 1, gsl_complex_rect(0.0, 1.0));
    gsl_matrix_complex_set(gamma[2], 1, 0, gsl_complex_rect(0.0, -1.0));
     
    // gamma3
    gsl_matrix_complex_set(gamma[3], 0, 0, gsl_complex_rect(-1.0, 0.0));
    gsl_matrix_complex_set(gamma[3], 1, 1, gsl_complex_rect(-1.0, 0.0));
}
