#define GEOM10_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include "geom.h"
#include "global.h"

// TYPE (1,0) GEOMETRY

void geom_check10()
{
    if(nH != 1 || nL != 0)
    {
        printf("Error: geometry is not (1,0)\n");
        exit(EXIT_FAILURE);
    }
}

void init_gamma10()
{
    gamma[0] = gsl_matrix_complex_calloc(1, 1);
    gsl_matrix_complex_set(gamma[0], 0, 0, gsl_complex_rect(1., 0.));
}

