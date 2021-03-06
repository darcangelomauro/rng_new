#ifndef UPDATE_H
#define UPDATE_H

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_rng.h>


extern void overwrite_init_file(char* name_init);
extern void init_data(char* name_init);
extern void init_cold(double Sfunc());
extern void init_hot(double Sfunc(), gsl_rng* r);
extern void simulation_free();
extern int doit(int i, int step);
extern int move(double deltaS(int, int, int, gsl_complex), gsl_rng* r);
extern double sweep(double deltaS(int, int, int, gsl_complex), gsl_rng* r);
extern void SCALE_autotune(double minTarget, double maxTarget, double deltaS(int, int, int, gsl_complex), gsl_rng* r);
extern char* simulation(double Sfunc(), double deltaS(int, int, int, gsl_complex), int rank, char* varG_code, gsl_rng* r);
extern void multicode_wrapper(double Sfunc(), double deltaS(int, int, int, gsl_complex), double INCR_G, int REP_G, int rank, char* varG_code, gsl_rng* r);
extern void adjust();

#endif
