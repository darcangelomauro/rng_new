#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdlib.h>
#include <stdio.h>

extern double* binned_vector(double* v, int n, int dimbin);
extern int binomial_coeff(int n, int k);
extern int binned_size(int n, int dimbin);
#endif
