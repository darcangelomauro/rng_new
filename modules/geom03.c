#define GEOM03_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_blas.h>
#include "geom.h"
#include "matop.h"
#include "global.h"

// TYPE (0,3) GEOMETRY

void geom_check03()
{
    if(nH != 1 || nL != 3)
    {
        printf("Error: geometry is not (0,3)\n");
        exit(EXIT_FAILURE);
    }
}

// build gamma matrices
void init_gamma03()
{

    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(2, 2);

    // gamma 0
    gsl_matrix_complex_set(gamma[0], 0, 0, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[0], 1, 1, gsl_complex_rect(1., 0.));

    // gamma 1
    gsl_matrix_complex_set(gamma[1], 0, 0, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gamma[1], 1, 1, gsl_complex_rect(0., -1.));

    // gamma 2
    gsl_matrix_complex_set(gamma[2], 0, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[2], 1, 0, gsl_complex_rect(-1., 0.));
    
    // gamma 3
    gsl_matrix_complex_set(gamma[3], 0, 1, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gamma[3], 1, 0, gsl_complex_rect(0., 1.));
}

