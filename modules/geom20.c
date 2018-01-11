#define GEOM20_C

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



// TYPE (2,0) GEOMETRY

void geom_check20()
{
    if(nH != 2 || nL != 0)
    {
        printf("Error: geometry is not (2,0)\n");
        exit(EXIT_FAILURE);
    }
}


// build gamma matrices
void init_gamma20()
{

    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(2, 2);


    // gamma 1
    gsl_matrix_complex_set(gamma[0], 0, 0, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[0], 1, 1, gsl_complex_rect(-1., 0.));

    // gamma 2
    gsl_matrix_complex_set(gamma[1], 0, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[1], 1, 0, gsl_complex_rect(1., 0.));
}
    
