#define STATISTICS_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "statistics.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_blas.h>
#include "fileop.h"
#include "matop.h"

int binomial_coeff(n, k)
{
  if(k > n)
  {
      printf("Error: in binomial coefficient n must be > k\n");
      exit(EXIT_FAILURE);
  }

  if(k == 0)
    return 1;

  if(k > n/2)
    return binomial_coeff(n,n-k);
  
  return (n * binomial_coeff(n-1,k-1)) / k;
}

// returns size of binned vector
int binned_size(int n, int dimbin)
{
    int size = n/dimbin;
    if((double)n/dimbin != size)
        size++;

    return size;
}


// bins a n components vector with bins of size dimbin
// and returns the binned vector
double* binned_vector(double* v, int n, int dimbin)
{
    if(dimbin > n)
    {
        printf("Error, bin size cannot be greater than vector length\n");
        exit(EXIT_FAILURE);
    }

    int size = binned_size(n, dimbin);
    double* u = calloc(size, sizeof(double));
    int* count = calloc(size, sizeof(double));

    for(int i=0; i<n; i++)
    {
        int idx = i/dimbin;
        if(idx >= size)
        {
            printf("Error index out of range\n");
            exit(EXIT_FAILURE);
        }

        u[idx] += v[i];
        count[idx]++;
    }

    for(int i=0; i<size; i++)
        u[i] /= (double)count[i];

    return u;
}






