#define PRINTOP_C

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>
#include "global.h"


// PRINT SIMULATION DATA

void print_data(FILE* fdata)
{
    // print data for program
    fprintf(fdata, "%d %d %d %.15lf %.15lf %d %d\n", dim, nH, nL, SCALE, G, Nsw, GAP);

    // print data for human
    fprintf(fdata, "dim: %d\n", dim);
    fprintf(fdata, "nH: %d\n", nH);
    fprintf(fdata, "nL: %d\n", nL);
    fprintf(fdata, "SCALE: %lf\n", SCALE);
    fprintf(fdata, "G: %lf\n", G);
    fprintf(fdata, "Nsw: %d\n", Nsw);
    fprintf(fdata, "GAP: %d\n", GAP);
    fprintf(fdata, "p: %d\n", CLIFF_P);
    fprintf(fdata, "q: %d\n", CLIFF_Q);
}

void print_time(FILE* fdata, char* s)
{
    // format time
    time_t timer;
    char buffer[26];
    struct tm* tm_info;
    time(&timer);
    tm_info = localtime(&timer);
    strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);

    fprintf(fdata, "%s %s\n", s, buffer);
}

void print_simulation(FILE* fobsS, FILE* fobsHL)
{
    // print action
    fprintf(fobsS, "%.15lf\n", S);

    // print H matrices
    for(int l=0; l<nHL; l++)
    {
        for(int i=0; i<dim; i++)
        {
            for(int j=0; j<dim; j++)
            {
                gsl_complex z = gsl_matrix_complex_get(MAT[l], i, j);
                fprintf(fobsHL, "%.15lf ", GSL_REAL(z)); 
                fprintf(fobsHL, "%.15lf ", GSL_IMAG(z));
            }
        }
        fprintf(fobsHL, "\n");
    }
}
