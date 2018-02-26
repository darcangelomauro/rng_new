#ifndef ACTION_H
#define ACTION_H

#include <gsl/gsl_complex_math.h>

extern double delta2(int uM, int I, int J, gsl_complex z);
extern double delta4(int uM, int I, int J, gsl_complex z);
extern double delta42(int uM, int I, int J, gsl_complex z);
extern double dirac2();
extern double dirac4();
extern double dirac42();

#endif
