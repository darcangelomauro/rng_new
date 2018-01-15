#ifndef ACTION_H
#define ACTION_H

#include "global.h"
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_rng.h>


extern double delta2(gsl_matrix_complex* dM, int uM);
extern double delta2_auto(gsl_matrix_complex* dM, int uM);
extern double delta4(gsl_matrix_complex* dM, int uM);
#endif
