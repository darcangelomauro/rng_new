#define GEOM02_C

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



// TYPE (0,2) GEOMETRY

void geom_check02()
{
    if(nH != 0 || nL != 2)
    {
        printf("Error: geometry is not (0,2)\n");
        exit(EXIT_FAILURE);
    }
}

// build gamma matrices
void init_gamma02()
{

    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(2, 2);


    // gamma 1
    gsl_matrix_complex_set(gamma[0], 0, 0, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gamma[0], 1, 1, gsl_complex_rect(0., -1.));

    // gamma 2
    gsl_matrix_complex_set(gamma[1], 0, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[1], 1, 0, gsl_complex_rect(-1., 0.));
}

