#ifndef CODEOP_H
#define CODEOP_H

#include <gsl/gsl_rng.h>

extern char* generate_code(int n, gsl_rng* r);
extern char* alloc_coded_filename(char* suffix, char* code);
extern char* alloc_folder_filename(char* suffix, char* folder);
extern char* alloc_rank_filename(int rank, char* prefix);

#endif
