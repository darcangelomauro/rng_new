#ifndef GEOM_H
#define GEOM_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix_complex_double.h>

#ifndef GEOM10_C
extern void geom_check10();
extern void init_gamma10();
#endif

#ifndef GEOM01_C
extern void geom_check01();
extern void init_gamma01();
#endif

#ifndef GEOM20_C
extern void geom_check20();
extern void init_gamma20();
#endif

#ifndef GEOM02_C
extern void geom_check02();
extern void init_gamma02();
#endif

#ifndef GEOM11_C
extern void geom_check11();
extern void init_gamma11();
#endif

#ifndef GEOM03_C
extern void geom_check03();
extern void init_gamma03();
#endif

#ifndef GEOM13_C
extern void geom_check13();
extern void init_gamma13();
#endif

#ifndef GEOM04_C
extern void geom_check04();
extern void init_gamma04();
#endif

#ifndef GEOM40_C
extern void geom_check40();
extern void init_gamma40();
#endif

#ifndef GEOM22_C
extern void geom_check22();
extern void init_gamma22();
#endif

#ifndef GEOM31_C
extern void geom_check31();
extern void init_gamma31();
#endif



#endif
