#define GEOM13_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_blas.h>
#include "geom.h"
#include "dirac.h"
#include "matop.h"
#include "global.h"

// TYPE (1,3) GEOMETRY

void geom_check13()
{
    if(nH != 2 || nL != 6)
    {
        printf("Error: geometry is not (1,3)\n");
        exit(EXIT_FAILURE);
    }
}

// build gamma matrices
void init_gamma13()
{
    
    for(int i=0; i<nHL; i++)
            gamma[i] = gsl_matrix_complex_calloc(dimG, dimG);


    // gamma H0
    gsl_matrix_complex_set(gamma[0], 0, 2, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[0], 1, 3, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[0], 2, 0, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[0], 3, 1, gsl_complex_rect(1., 0.));

    // gamma H1
    gsl_matrix_complex_set(gamma[1], 0, 2, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gamma[1], 1, 3, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gamma[1], 2, 0, gsl_complex_rect(0., -1.));
    gsl_matrix_complex_set(gamma[1], 3, 1, gsl_complex_rect(0., -1.));

    // gamma L0
    gsl_matrix_complex_set(gamma[2], 0, 3, gsl_complex_rect(-1., 0.));
    gsl_matrix_complex_set(gamma[2], 1, 2, gsl_complex_rect(-1., 0.));
    gsl_matrix_complex_set(gamma[2], 2, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[2], 3, 0, gsl_complex_rect(1., 0.));
    
    // gamma L1
    gsl_matrix_complex_set(gamma[3], 0, 3, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gamma[3], 1, 2, gsl_complex_rect(0., -1.));
    gsl_matrix_complex_set(gamma[3], 2, 1, gsl_complex_rect(0., -1.));
    gsl_matrix_complex_set(gamma[3], 3, 0, gsl_complex_rect(0., 1.));
    
    // gamma L2
    gsl_matrix_complex_set(gamma[4], 0, 2, gsl_complex_rect(-1., 0.));
    gsl_matrix_complex_set(gamma[4], 1, 3, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[4], 2, 0, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[4], 3, 1, gsl_complex_rect(-1., 0.));
    
    // gamma L3
    gsl_matrix_complex_set(gamma[5], 0, 2, gsl_complex_rect(0., -1.));
    gsl_matrix_complex_set(gamma[5], 1, 3, gsl_complex_rect(0., 1.));
    gsl_matrix_complex_set(gamma[5], 2, 0, gsl_complex_rect(0., -1.));
    gsl_matrix_complex_set(gamma[5], 3, 1, gsl_complex_rect(0., 1.));

    // gamma L4
    gsl_matrix_complex_set(gamma[6], 0, 3, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[6], 1, 2, gsl_complex_rect(-1., 0.));
    gsl_matrix_complex_set(gamma[6], 2, 1, gsl_complex_rect(1., 0.));
    gsl_matrix_complex_set(gamma[6], 3, 0, gsl_complex_rect(-1., 0.));

    // gamma L5
    gsl_matrix_complex_set(gamma[7], 0, 3, gsl_complex_rect(0., -1.));
    gsl_matrix_complex_set(gamma[7], 1, 2, gsl_complex_rect(0., -1.));
    gsl_matrix_complex_set(gamma[7], 2, 1, gsl_complex_rect(0., -1.));
    gsl_matrix_complex_set(gamma[7], 3, 0, gsl_complex_rect(0., -1.));
}


