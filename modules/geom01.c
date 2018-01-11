#define GEOM01_C

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

// TYPE (0,1) GEOMETRY

void geom_check01()
{
    if(nH != 0 || nL != 1)
    {
        printf("Error: geometry is not (0,1)\n");
        exit(EXIT_FAILURE);
    }
}

void init_gamma01()
{
    gamma[0] = gsl_matrix_complex_calloc(1, 1);
    gsl_matrix_complex_set(gamma[0], 0, 0, gsl_complex_rect(0., -1.));
}

