#ifndef GLOBAL_H
#define GLOBAL_H

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex_math.h>

#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

#define ABS2(a) gsl_complex_abs2(a)
#define CA(a,b) gsl_complex_add(a, b)
#define CAR(a,b) gsl_complex_add_real(a, b)
#define CM(a,b) gsl_complex_mul(a, b)
#define CMR(a,b) gsl_complex_mul_real(a, b)
#define CC(a) gsl_complex_conjugate(a)
#define MG(m,i,j) gsl_matrix_complex_get(m, i, j)


// CLIFFORD RELATED FUNCTIONS **********************************************
#define CLIFF_P 1
#define CLIFF_Q 3

// The following macro allows to simplify clifford-specific functions.
// For example, if CLIFF_P and CLIFF_Q are set to 2 and 0 respectively,
// then FUNCTION(init_gamma, CLIFF_P, CLIFF_Q) expands to init_gamma20.
#define PASTER(name, x, y) name ## x ## y
#define FUNCTION(name, x, y) PASTER(name, x, y)
// Note: PASTER is seemingly redundant, but it is necessary in order for
// the ## operator to expand correctly.


// GLOBAL VARIABLES ********************************************************
EXTERN int dim;
EXTERN int nH;
EXTERN int nL;
EXTERN double SCALE;        // metropolis scale factor
EXTERN double G;            // coupling constants
EXTERN double S;            // action
EXTERN int Nsw;             // number of simulation sweeps
#define GAP (100)
#define ADJ (1000)
#define nHL (nH+nL) 

// H and L matrices
EXTERN gsl_matrix_complex** MAT;
EXTERN int* e;

// dirac matrix
EXTERN int dimG;
EXTERN int dimD;
EXTERN gsl_matrix_complex* DIRAC;
EXTERN gsl_matrix_complex** gamma;
EXTERN gsl_complex** gamma_table;

#endif

