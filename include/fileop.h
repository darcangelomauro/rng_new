#ifndef FILEOP_H
#define FILEOP_H

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

extern char* generate_code(int n, gsl_rng* r);
extern char* alloc_coded_filename(char* suffix, char* code);
extern char* alloc_folder_filename(char* suffix, char* folder);
extern char* alloc_rank_filename(int rank, char* prefix);
extern void avgfile(int row, int col, char *input, char *output);
extern void DfromH_t10(int row, int col, char *input, char *output);
extern void histo(char* inname, char* outname, int Nbin, int min, int max, int row, int col);
#endif
