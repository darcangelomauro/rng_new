#ifndef UPDATE_H
#define UPDATE_H

#include "global.h"
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_rng.h>


extern void overwrite_init_file();
extern void init_data();
extern void init_data_analysis();
extern void init_minimal(void init_gamma());
extern void init_cold(gsl_complex Sfunc(), void init_gamma());
extern void init_custom(gsl_complex Sfunc(), void init_gamma(), char* filename);
extern void init_hot(gsl_complex Sfunc(), void init_gamma(), gsl_rng* r);
extern void simulation_free();
extern int measurement(int i, int step);
/*
extern int move(void Sfunc(double*, int*), int mode, gsl_rng* r);
extern double sweep(void Sfunc(double*, int*), int mode, gsl_rng* r);
extern void SCALE_autotune(double minTarget, double maxTarget, void Sfunc(double*, int*), gsl_rng* r);
extern char* simulation(void Sfunc(double*, int*), int mode, int renorm, void init_gamma(), gsl_rng* r);
extern void multicode_wrapper(void Sfunc(double*, int*), void init_gamma(), int renorm, double INCR_G, int REP_G, int INCR_dim, int REP_dim, gsl_rng* r);
*/
extern void hermitization();
extern void apply_renormalization(int b, FILE* fobsS, FILE* fobsHL);
#endif
