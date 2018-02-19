#ifndef ACTION_H
#define ACTION_H

#include "global.h"
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_rng.h>


extern double delta2(int uM, int I, int J, gsl_complex z);
extern double delta4(int uM, int I, int J, gsl_complex z);
extern double delta4_BETA(int uM, int I, int J, gsl_complex z);
extern double dirac2();
extern double dirac4();
#endif
