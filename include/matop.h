#ifndef MATOP_H
#define MATOP_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_rng.h>

extern void printmat_complex(gsl_matrix_complex *m, int n, int k);
extern void fprintmat_complex(gsl_matrix_complex *m, int n, int k, FILE* f);
extern void generate_HL(gsl_matrix_complex* m, int mode, int n, gsl_rng* r);
extern void matrix_power(gsl_matrix_complex* m, int n, gsl_matrix_complex* res);
extern void matrix_multiprod(gsl_matrix_complex** m, int n, gsl_matrix_complex* res);
extern void dispherm(gsl_matrix_complex* m, int n);
extern void diag(gsl_matrix_complex* m, gsl_vector* eval, int n);
extern gsl_complex trace(gsl_matrix_complex* m);
extern double trace_herm(gsl_matrix_complex* m);
extern double trace2(gsl_matrix_complex* m1, gsl_matrix_complex* m2);
extern double trace2_self(gsl_matrix_complex* m);
extern gsl_complex trace4(gsl_matrix_complex* m1, gsl_matrix_complex* m2, gsl_matrix_complex* m3, gsl_matrix_complex* m4);
extern gsl_complex trace3(gsl_matrix_complex* m1, gsl_matrix_complex* m2, gsl_matrix_complex* m3);
extern gsl_complex MG3(gsl_matrix_complex* m1, gsl_matrix_complex* m2, gsl_matrix_complex* m3, int i, int j);
extern gsl_complex MG2(gsl_matrix_complex* m1, gsl_matrix_complex* m2, int i, int j);
extern void traceless(gsl_matrix_complex* m, int n);
extern void tensor(gsl_matrix_complex* a, gsl_matrix_complex* b, gsl_matrix_complex* c);
extern void commutator(gsl_matrix_complex* m, gsl_matrix_complex* c);
extern void anticommutator(gsl_matrix_complex* m, gsl_matrix_complex* c);
extern void make_hermitian(gsl_matrix_complex* a);
#endif
