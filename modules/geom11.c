#define GEOM11_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include "geom.h"
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
    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(dimG, dimG);

    // gamma0
    gsl_matrix_complex_set(gamma[0], 0, 0, gsl_complex_rect(1, 0));
    gsl_matrix_complex_set(gamma[0], 1, 1, gsl_complex_rect(-1, 0));
     
    // gamma1
    gsl_matrix_complex_set(gamma[1], 0, 1, gsl_complex_rect(0.0, 1.0));
    gsl_matrix_complex_set(gamma[1], 1, 0, gsl_complex_rect(0.0, -1.0));
     
}
