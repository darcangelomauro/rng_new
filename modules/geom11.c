#define GEOM11_C

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



// TYPE (1,1) GEOMETRY

void geom_check11()
{
    if(nH != 1 || nL != 1)
    {
        printf("Error: geometry is not (1,1)\n");
        exit(EXIT_FAILURE);
    }
}


// build gamma matrices
void init_gamma11()
{

    gamma[0] = gsl_matrix_complex_calloc(2, 2);
    gamma[0] = gsl_matrix_complex_calloc(2, 2);


    // gamma 1
    gsl_matrix_complex_set(gamma[0], 0, 0, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[0], 1, 1, gsl_complex_rect(-1., 0.));

    // gamma 2
    gsl_matrix_complex_set(gamma[1], 0, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[1], 1, 0, gsl_complex_rect(-1., 0.));
}

