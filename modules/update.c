#define UPDATE_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_double.h>
#include "update.h"
#include "matop.h"
#include "fileop.h"
#include "math.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector.h>
#include "global.h"
#include "macros.h"
#include "data.h"
#include "statistics.h"

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

#define NAME(vec, idx, alt, check)  (idx == check ? (alt) : (vec[idx]))


void overwrite_init_file()
{
    FILE* finit = fopen("init.txt", "w");
    fprintf(finit, "%d\n", dim);
    fprintf(finit, "%d\n", nH);
    fprintf(finit, "%d\n", nL);
    fprintf(finit, "%lf\n", SCALE);
    fprintf(finit, "%lf\n", G);
    fprintf(finit, "%d\n", Ntherm);
    fprintf(finit, "%d\n", Nsw);
    fclose(finit);
}


void init_data()
{
    // read input data
    FILE* finit = fopen("init.txt", "r");
    if(finit == NULL)
    {
        printf("Error: input.txt not found\n");
        exit(EXIT_FAILURE);
    }

    int narg;
    // initialize matrix dimension
    narg = fscanf(finit, "%d", &dim);
    // initialize #H matrices
    narg += fscanf(finit, "%d", &nH);
    // initialize #L matrices
    narg += fscanf(finit, "%d", &nL);
    // initialize SCALE
    narg += fscanf(finit, "%lf", &SCALE);
    // initialize coupling constant
    narg += fscanf(finit, "%lf", &G);
    // initialize number of thermalization sweeps
    narg += fscanf(finit, "%d", &Ntherm);
    // initialize number of simulation sweeps
    narg += fscanf(finit, "%d", &Nsw);

    if(narg < 7)
    {
        printf("Error: not enough data in init.txt\n");
        exit(EXIT_FAILURE);
    }

    // initialize dimD
    dimG = (int)pow(2., (cliff_p+cliff_q)/2);
    dimD = dimG*dim*dim;
}

void init_data_analysis()
{
    // read input data
    FILE* finit = fopen("init_analysis.txt", "r");
    if(finit == NULL)
    {
        printf("Error: input_analysis.txt not found\n");
        exit(EXIT_FAILURE);
    }

    int narg;
    // initialize matrix dimension
    narg = fscanf(finit, "%d", &dim);
    // initialize #H matrices
    narg += fscanf(finit, "%d", &nH);
    // initialize #L matrices
    narg += fscanf(finit, "%d", &nL);
    // initialize SCALE
    narg += fscanf(finit, "%lf", &SCALE);
    // initialize coupling constant
    narg += fscanf(finit, "%lf", &G);
    // initialize number of thermalization sweeps
    narg += fscanf(finit, "%d", &Ntherm);
    // initialize number of simulation sweeps
    narg += fscanf(finit, "%d", &Nsw);

    if(narg < 7)
    {
        printf("Error: not enough data in init_analysis.txt\n");
        exit(EXIT_FAILURE);
    }

    // initialize dimD
    if(nHL == 1)
        dimG = 1;
    else if(nHL<5)
        dimG = 2;
    else
        dimG = 4;
    dimD = dimG*dim*dim;
}

//initializes H matrices to unity and L to zero
void init_minimal(void init_gamma())
{
    /*
    //initialize H matrices to unity
    H = malloc(nH*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nH; i++)
        H[i] = gsl_matrix_complex_alloc(dim, dim);

    //initialize L matrices to zero
    L = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nL; i++)
        L[i] = gsl_matrix_complex_alloc(dim, dim);
    
    // allocate displacement matrix
    M = gsl_matrix_complex_alloc(dim, dim);
    
    // initialize gamma matrices and dirac
    gammaH = malloc(nH*sizeof(gsl_matrix_complex*));
    gammaL = malloc(nL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nH; i++)
        gammaH[i] = gsl_matrix_complex_calloc(dimG, dimG);
    for(int i=0; i<nL; i++)
        gammaL[i] = gsl_matrix_complex_calloc(dimG, dimG);
    init_gamma();
    
    DIRAC = gsl_matrix_complex_calloc(dimD, dimD);
    */
}
//initializes H matrices to unity and L to zero
void init_cold(gsl_complex Sfunc(), void init_gamma())
{
    //initialize H matrices to unity
    MAT = malloc(nHL*sizeof(gsl_matrix_complex*));
    tr = malloc(nH*sizeof(double));
    tr2 = malloc(nH*sizeof(double));
    for(int i=0; i<nHL; i++)
    {
        MAT[i] = gsl_matrix_complex_calloc(dim, dim);
        if(i < nH)
        {
            gsl_matrix_complex_set_identity(MAT[i]);
            tr[i] = dim;
            tr2[i] = dim;
        }
    }

    // allocate displacement matrix
    M = gsl_matrix_complex_calloc(dim, dim);
    

    // initialize action
    /*
    double var[14];
    int control[14] = {1,0,1,1,1,1,1,1,1,1,1,1,1,1};
    Sfunc(var, control);
    */


    // initialize gamma matrices and dirac
    gamma = malloc(nHL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(dimG, dimG);
    init_gamma();
    
    DIRAC = gsl_matrix_complex_calloc(dimD, dimD);
    
    S = GSL_REAL(Sfunc());
}

//initializes H and L with random matrices
void init_hot(gsl_complex Sfunc(), void init_gamma(), gsl_rng* r)
{
    //initialize H matrices randomly
    MAT = malloc(nHL*sizeof(gsl_matrix_complex*));
    tr = malloc(nH*sizeof(double));
    tr2 = malloc(nH*sizeof(double));
    for(int i=0; i<nHL; i++)
    {
        MAT[i] = gsl_matrix_complex_calloc(dim, dim);
        int mode = (i >= nH);
        generate_HL(MAT[i], mode, dim, r);
    }
    
    // allocate displacement matrix
    M = gsl_matrix_complex_calloc(dim, dim);
    

    // initialize action and traces
    /*
    double var[14];
    int control[14] = {1,0,1,1,1,1,1,1,1,1,1,1,1,1};
    Sfunc(var, control);

    S = var[0];
    for(int i=0; i<nH; i++)
    {
        tr[i] = var[nH+2];
        tr2[i] = var[nH+3];
    }
    */

    
    // initialize gamma matrices and dirac
    gamma = malloc(nHL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(dimG, dimG);
    init_gamma();
    
    DIRAC = gsl_matrix_complex_calloc(dimD, dimD);
    S = GSL_REAL(Sfunc());
}


void simulation_free()
{
    // free H matrices
    for(int i=0; i<nHL; i++)
        gsl_matrix_complex_free(MAT[i]);
    free(MAT);
    free(tr);
    free(tr2);

    // free displacement matrix
    gsl_matrix_complex_free(M);
    
    // free gamma matrices and dirac
    for(int i=0; i<nHL; i++)
        gsl_matrix_complex_free(gamma[i]);
    free(gamma);

    gsl_matrix_complex_free(DIRAC);
}



// returns 1 if ith sweep corresponds to
// a measurement (step is the distance between
// measurements)
int measurement(int i, int step)
{
    return !(i%step);
}

/*
int move(void Sfunc(double*, int*), int mode, gsl_rng* r)
{
    // check is what will be returned
    int accepted = 0;

    // MONTECARLO MOVE PROPOSAL
    double S1[14];
    int control[14] = {1,0,1,1,1,1,1,1,1,1,1,1,1,1};
    if(mode) control[1] = 1;

    
    // decide what matrix gets updated
    int uM = (int)(nHL*gsl_rng_uniform(r));
    while(uM == nHL)
        uM = (int)(nHL*gsl_rng_uniform(r));

    
    // buffer is needed to temporarily update the H or L matrices
    // but it never allocates new memory, so no need to free at the end
    gsl_matrix_complex* buffer;
    

    // generate displacement dM
    int mode = (uM >= nH);
    generate_HL(M, mode, dim, r);    //M is now the proposed displacement from MAT[uM]

    // compute displaced matrix
    gsl_matrix_complex_add(M, MAT[uM]);             //now M = MAT[uM] + displacement
    buffer = MAT[uM];
    MAT[uM] = M;
    Sfunc(S1, control);
    MAT[uM] = buffer;

    // now we evaluate the move
    if(S1[0] < S)
    {
        gsl_matrix_complex_memcpy(MAT[uM], M);
        if(!mode)
        {
            tr[uM] = S1[2*uM+2];
            tr2[uM] = S1[2*uM+3];
        }

        // update action
        S = S1[0];

        // move accepted
        accepted = 1;
    }
    else
    {
        double e = exp(S-S1[0]);
        double p = gsl_rng_uniform(r);

        if(e>p)
        {
            gsl_matrix_complex_memcpy(MAT[uM], M);
            if(!mode)
            {
                tr[uM] = S1[2*uM+2];
                tr2[uM] = S1[2*uM+3];
            }

            // update action
            S = S1[0];

            // move accepted
            accepted = 1;
        }
    }

    return accepted;

}

// a sweep is a collection of nHL proposed moves
// it returns the acceptance rate of that sweep
double sweep(void Sfunc(double*, int*), int mode, gsl_rng* r)
{
    double sum = 0.;
    for(int i=0; i<nHL; i++)
        sum += move(Sfunc, mode, r);


    return sum/(double)nHL;
}


// a routine to automatically set the SCALE factor to give acceptance rate
// between minTarget and maxTarget (very rudimental and ugly, seems to work reasonably well)
void SCALE_autotune(double minTarget, double maxTarget, void Sfunc(double*, int*), gsl_rng* r)
{
    int n1 = 50;
    int n2 = 50;
    int m = 50;
    double limit = 1e-5;
    double ar = 0.;

    // initialize ar
    for(int i=0; i<n1; i++)
        ar += sweep(Sfunc, 0, r);
    ar /= (double)n1;

    int count = 0;

    while((ar < minTarget || ar > maxTarget) && count < 2 )
    {
        count++;

        // tune delta
        double temp = SCALE;
        int order;
        for(order=0; order<1000; order++)
        {
            if((int)temp != 0)
                break;
            else
                temp *= 10.;
        }
        double delta = pow(10., -order);

        // tune SCALE
        while(delta > limit)
        {
            for(int j=0; j<m; j++)
            {
                if(ar > maxTarget)
                    SCALE += delta;
                else if(ar < minTarget)
                {
                    if((SCALE-delta) > 0)
                        SCALE -= delta;
                    else
                        break;
                }
                ar = 0.;
                for(int k=0; k<n2; k++)
                    ar += sweep(Sfunc, 0, r);
                ar /= (double)n2;
            }
            delta /= 10.;
        }
    }
}


// complete simulation routine
// first, it prints on file the simulation data
// the thermalization part outputs two files "thermalizationX.txt" (where X is 1 or 2) with the action value
// the simulation part outputs a file "simulation.txt" with the action and the H and L matrices
// returns acceptance rate
char* simulation(void Sfunc(double*, int*), int mode, int renorm, void init_gamma(), gsl_rng* r)
{
    init_data();
    
    GEOM_CHECK();
    
    init_cold(Sfunc, init_gamma);
    SCALE_autotune(0.21, 0.3, Sfunc, r);
    printf("Auto tuned SCALE: %lf\n", SCALE);

    // generate unique filename
    char* code = generate_code(5, r);
    printf("%s\n", code);
    char* name_data = alloc_coded_filename("data", code);
    char* name_therm1 = alloc_coded_filename("therm1", code);
    char* name_therm2 = alloc_coded_filename("therm2", code);
    char* name_simS = alloc_coded_filename("simS", code);
    char* name_simHL = alloc_coded_filename("simHL", code);
    char* name_simS_r = alloc_coded_filename("simS_r", code);
    char* name_simHL_r = alloc_coded_filename("simHL_r", code);

    FILE* fdata = fopen(name_data, "w");
    FILE* ftherm1 = fopen(name_therm1, "w");
    FILE* ftherm2 = fopen(name_therm2, "w");
    FILE* fsimS = fopen(name_simS, "w");
    FILE* fsimHL = fopen(name_simHL, "w");
    FILE* fsimS_r;
    FILE* fsimHL_r;
    if(renorm)
    {
        fsimS_r = fopen(name_simS_r, "w");
        fsimHL_r = fopen(name_simHL_r, "w");
    }

    free(name_data);
    free(name_therm1);
    free(name_therm2);
    free(name_simS);
    free(name_simHL);
    free(name_simS_r);
    free(name_simHL_r);

    // print simulation data
    print_data(fdata);

    // thermalization1
    init_cold(Sfunc, init_gamma);
    print_time(fdata, "start therm1:");
    for(int i=0; i<Ntherm; i++)
    {
        sweep(Sfunc, mode, r);
        print_thermalization_plus(ftherm1);
    }
    print_time(fdata, "end therm1:");
    fclose(ftherm1);
    simulation_free();
    
    // thermalization2
    init_hot(Sfunc, init_gamma, r);
    print_time(fdata, "start therm2:");
    for(int i=0; i<Ntherm; i++)
    {
        sweep(Sfunc, mode, r);
        print_thermalization_plus(ftherm2);
    }
    print_time(fdata, "end therm2:");
    fclose(ftherm2);

    
    // simulation (starts from second thermalization)
    double ar = 0.;
    print_time(fdata, "start simulation:");
    for(int i=0; i<Nsw; i++)
    {
        ar += sweep(Sfunc, mode, r);
        if(measurement(i, GAP))
        {
            print_simulation(fsimS, fsimHL);
            if(renorm)
                apply_renormalization(2, fsimS_r, fsimHL_r);
        }
    }
    print_time(fdata, "end simulation:");
    fprintf(fdata, "acceptance rate: %lf", ar/(double)Nsw);

    fclose(fsimS);
    fclose(fsimHL);
    if(renorm)
    {
        fclose(fsimS_r);
        fclose(fsimHL_r);
    }
    fclose(fdata);
    simulation_free();

    printf("acceptance rate: %lf\n", ar/(double)Nsw);

    return code;
}


// simulation routine for varying over matrix dimension and coupling constant G in a convenient way.
// for each matrix dimension generates a unique code XXXXX and three files:
// file "XXXXX_varG_data.txt": collects data on the range in which G varies, and total duration of the simulation
// file "XXXXX_varG_args.txt": collects list of codes needed to run multicode_analysis_main
// file "XXXXX_varG_G_args.txt": shows correspondence between G and code
// 
void multicode_wrapper(void Sfunc(double*, int*), void init_gamma(), int renorm, double INCR_G, int REP_G, int INCR_dim, int REP_dim, gsl_rng* r)
{
    // cycle over matrix dimension
    for(int j=0; j<REP_dim; j++)
    {
        // copy init.txt
        FILE* finit = fopen("init.txt", "r");
        if(finit == NULL)
        {
            printf("Error: input.txt not found\n");
            exit(EXIT_FAILURE);
        }
        int narg;
        int dim_, nH_, nL_, Ntherm_, Nsw_;
        double SCALE_, G_;
        // initialize matrix dimension
        narg = fscanf(finit, "%d", &dim_);
        // initialize #H matrices
        narg += fscanf(finit, "%d", &nH_);
        // initialize #L matrices
        narg += fscanf(finit, "%d", &nL_);
        // initialize SCALE
        narg += fscanf(finit, "%lf", &SCALE_);
        // initialize coupling constant
        narg += fscanf(finit, "%lf", &G_);
        // initialize number of thermalization sweeps
        narg += fscanf(finit, "%d", &Ntherm_);
        // initialize number of simulation sweeps
        narg += fscanf(finit, "%d", &Nsw_);
        if(narg < 7)
        {
            printf("Error: not enough data in init.txt\n");
            exit(EXIT_FAILURE);
        }
        fclose(finit);


        // generate code for varG simulation
        char* varG_code = generate_code(5, r);
        // generate output files
        char* name_varG_data = alloc_coded_filename("varG_data", varG_code);
        char* name_varG_args = alloc_coded_filename("varG_args", varG_code);
        char* name_varG_G_args = alloc_coded_filename("varG_G_args", varG_code);
        FILE* fvarG_data = fopen(name_varG_data, "w");
        FILE* fvarG_args = fopen(name_varG_args, "w");
        FILE* fvarG_G_args = fopen(name_varG_G_args, "w");
        
        // prompt some shenenigans
        printf("*********\n");
        printf("Starting variable G simulation with data:\n");
        printf("%d values of G uniformly distributed in range [%lf, %lf]\n", REP_G, G_, (G_ + (REP_G-1)*INCR_G));
        printf("dim: %d\n", dim_);
        printf("geometry: (%d, %d)\n", cliff_p, cliff_q);
        printf("Codename: %s\n", varG_code);
        printf("*********\n");

        // write same shenenigans
        fprintf(fvarG_data, "REP_G: %d\n", REP_G);
        fprintf(fvarG_data, "INCR_G: %lf\n", INCR_G);
        fprintf(fvarG_data, "G initial: %lf\n", G_);
        fprintf(fvarG_data, "G final: %lf\n", (G_ + (REP_G-1)*INCR_G));
        fprintf(fvarG_data, "dim: %d\n", dim_);
        fprintf(fvarG_data, "p: %d\n", cliff_p);
        fprintf(fvarG_data, "q: %d\n", cliff_q);

        // simulate with variable G
        print_time(fvarG_data, "start simulation:");
        for(int i=0; i<REP_G; i++)
        {
            char* code = simulation(Sfunc, 0, renorm, init_gamma, r);
            fprintf(fvarG_args, "%s ", code);
            fprintf(fvarG_G_args, "%lf %s\n", G, code);
            printf("dim: %d,    G: %lf\n", dim, G);
            G += INCR_G;
            overwrite_init_file();
            free(code);
        }
        print_time(fvarG_data, "end simulation:");

        // free memory
        free(varG_code);
        free(name_varG_data);
        free(name_varG_args);
        free(name_varG_G_args);
        fclose(fvarG_data);
        fclose(fvarG_args);
        fclose(fvarG_G_args);


        // restore init.txt, but update matrix dimension
        FILE* finit2 = fopen("init.txt", "w");
        if(finit2 == NULL)
        {
            printf("Error: unable to write on init.txt\n");
            exit(EXIT_FAILURE);
        }
        fprintf(finit2, "%d\n", dim_+INCR_dim);
        fprintf(finit2, "%d\n", nH_);
        fprintf(finit2, "%d\n", nL_);
        fprintf(finit2, "%lf\n", SCALE_);
        fprintf(finit2, "%lf\n", G_);
        fprintf(finit2, "%d\n", Ntherm_);
        fprintf(finit2, "%d\n", Nsw_);
        fclose(finit2);
    }
}
*/


void hermitization()
{
    for(int i=0; i<nHL; i++)
        make_hermitian(MAT[i]);
}


// the matrices get scaled by dim/b preserving hermitian behaviour
void apply_renormalization(int b, FILE* fobsS, FILE* fobsHL)
{
    /*
    if((double)dim/(double)b != dim/b)
    {
        printf("Error: dim is not a multiple of b\n");
        exit(EXIT_FAILURE);
    }

    int dim_r = dim/b;

    gsl_matrix_complex** MAT_r = malloc(nHL*sizeof(gsl_matrix_complex*));

    // build blocked matrices
    for(int i=0; i<nHL; i++)
    {
        MAT_r[i] = gsl_matrix_complex_alloc(dim_r, dim_r);
        renormalize(MAT[i], MAT_r[i]);
    }

    // print on file
    print_renormalized(fobsS, fobsHL, MAT_r, dim_r);

    // free memory
    for(int i=0; i<nHL; i++)
        gsl_matrix_complex_free(MAT_r[i]);
    free(MAT_r);
    */
}



