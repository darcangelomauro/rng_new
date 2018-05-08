#define UPDATE_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include "global.h"
#include "update.h"
#include "matop.h"
#include "codeop.h"
#include "printop.h"
#include "geom.h"


void overwrite_init_file(char* name_init)
{
    FILE* finit = fopen(name_init, "w");
    fprintf(finit, "%d\n", dim);
    fprintf(finit, "%d\n", nH);
    fprintf(finit, "%d\n", nL);
    fprintf(finit, "%lf\n", SCALE);
    fprintf(finit, "%lf\n", G);
    fprintf(finit, "%d\n", Nsw);
    fclose(finit);
}


void init_data(char* name_init)
{
    // read input data
    FILE* finit = fopen(name_init, "r");
    if(finit == NULL)
    {
        printf("Error: %s not found\n", name_init);
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
    // initialize number of simulation sweeps
    narg += fscanf(finit, "%d", &Nsw);

    if(narg < 6)
    {
        printf("Error: not enough data in %s\n", name_init);
        exit(EXIT_FAILURE);
    }

    // initialize dimD
    dimG = (int)pow(2., (int)(CLIFF_P + CLIFF_Q)/2);
    dimD = dimG*dim*dim;
}


void init_cold(double Sfunc())
{
    //initialize H matrices to unity
    MAT = malloc(nHL*sizeof(gsl_matrix_complex*));
    e = malloc(nHL*sizeof(int));
    for(int i=0; i<nHL; i++)
    {
        MAT[i] = gsl_matrix_complex_calloc(dim, dim);
        gsl_matrix_complex_set_identity(MAT[i]);
        if(i < nH)
            e[i] = 1;
        else
            e[i] = -1;
    }

    // initialize gamma matrices
    gamma = malloc(nHL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(dimG, dimG);
    FUNCTION(init_gamma, CLIFF_P, CLIFF_Q)();

    // initialize gamma_table
    gamma_table = malloc(5*sizeof(gsl_complex*));
    gamma_table[0] = malloc(1*sizeof(gsl_complex));
    gamma_table[1] = malloc(1*sizeof(gsl_complex));
    gamma_table[2] = malloc(nHL*nHL*sizeof(gsl_complex));
    gamma_table[3] = malloc(nHL*nHL*nHL*sizeof(gsl_complex));
    gamma_table[4] = malloc(nHL*nHL*nHL*nHL*sizeof(gsl_complex));

    for(int i=0; i<nHL; i++)
    {
        for(int j=0; j<nHL; j++)
        {
            gsl_matrix_complex* temp2 = gsl_matrix_complex_alloc(dimG, dimG);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, gamma[i], gamma[j], GSL_COMPLEX_ZERO, temp2);
            for(int k=0; k<nHL; k++)
            {
                gsl_matrix_complex* temp3 = gsl_matrix_complex_alloc(dimG, dimG);
                gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, temp2, gamma[k], GSL_COMPLEX_ZERO, temp3);
                for(int l=0; l<nHL; l++)
                {
                    gsl_matrix_complex* temp4 = gsl_matrix_complex_alloc(dimG, dimG);
                    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, temp3, gamma[l], GSL_COMPLEX_ZERO, temp4);

                    gamma_table[4][l + nHL*(k + nHL*(j + nHL*i))] = trace(temp4);
                    gsl_matrix_complex_free(temp4);
                }
                gamma_table[3][k + nHL*(j + nHL*i)] = trace(temp3);
                gsl_matrix_complex_free(temp3);
            }
            gamma_table[2][j + nHL*i] = trace(temp2);
            gsl_matrix_complex_free(temp2);
        }
    }
    
    // alloc dirac
    //DIRAC = gsl_matrix_complex_calloc(dimD, dimD);
    
    // initialize action
    S = Sfunc();
}


//initializes H and L with random matrices
void init_hot(double Sfunc(), gsl_rng* r)
{
    //initialize H matrices randomly
    MAT = malloc(nHL*sizeof(gsl_matrix_complex*));
    e = malloc(nHL*sizeof(int));
    for(int i=0; i<nHL; i++)
    {
        MAT[i] = gsl_matrix_complex_calloc(dim, dim);
        generate_HL(MAT[i], 0, dim, r);
        if(i < nH)
            e[i] = 1;
        else
            e[i] = -1;
    }
    
    // initialize gamma matrices
    gamma = malloc(nHL*sizeof(gsl_matrix_complex*));
    for(int i=0; i<nHL; i++)
        gamma[i] = gsl_matrix_complex_calloc(dimG, dimG);
    FUNCTION(init_gamma, CLIFF_P, CLIFF_Q)();

    // initialize gamma_table
    gamma_table = malloc(5*sizeof(gsl_complex*));
    gamma_table[0] = malloc(1*sizeof(gsl_complex));
    gamma_table[1] = malloc(1*sizeof(gsl_complex));
    gamma_table[2] = malloc(nHL*nHL*sizeof(gsl_complex));
    gamma_table[3] = malloc(nHL*nHL*nHL*sizeof(gsl_complex));
    gamma_table[4] = malloc(nHL*nHL*nHL*nHL*sizeof(gsl_complex));

    for(int i=0; i<nHL; i++)
    {
        for(int j=0; j<nHL; j++)
        {
            gsl_matrix_complex* temp2 = gsl_matrix_complex_alloc(dimG, dimG);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, gamma[i], gamma[j], GSL_COMPLEX_ZERO, temp2);
            for(int k=0; k<nHL; k++)
            {
                gsl_matrix_complex* temp3 = gsl_matrix_complex_alloc(dimG, dimG);
                gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, temp2, gamma[k], GSL_COMPLEX_ZERO, temp3);
                for(int l=0; l<nHL; l++)
                {
                    gsl_matrix_complex* temp4 = gsl_matrix_complex_alloc(dimG, dimG);
                    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, temp3, gamma[l], GSL_COMPLEX_ZERO, temp4);

                    gamma_table[4][l + nHL*(k + nHL*(j + nHL*i))] = trace(temp4);
                    gsl_matrix_complex_free(temp4);
                }
                gamma_table[3][k + nHL*(j + nHL*i)] = trace(temp3);
                gsl_matrix_complex_free(temp3);
            }
            gamma_table[2][j + nHL*i] = trace(temp2);
            gsl_matrix_complex_free(temp2);
        }
    }
    
    // alloc dirac
    //DIRAC = gsl_matrix_complex_calloc(dimD, dimD);

    // initialize action
    S = Sfunc();
}


void simulation_free()
{
    // free H matrices
    for(int i=0; i<nHL; i++)
        gsl_matrix_complex_free(MAT[i]);
    free(MAT);
    free(e);
    
    // free gamma matrices and dirac
    for(int i=0; i<nHL; i++)
        gsl_matrix_complex_free(gamma[i]);
    free(gamma);

    // free gamma table
    free(gamma_table[0]);
    free(gamma_table[1]);
    free(gamma_table[2]);
    free(gamma_table[3]);
    free(gamma_table[4]);
    free(gamma_table);

    //gsl_matrix_complex_free(DIRAC);
}



// returns 1 if i % step == 0
int doit(int i, int step)
{
    return !(i%step);
}


int move(double deltaS(int, int, int, gsl_complex), gsl_rng* r)
{
    // accepted is what will be returned
    int accepted = 0;

    // MONTECARLO MOVE PROPOSAL
    
    // decide what matrix gets updated
    int uM = (int)(nHL*gsl_rng_uniform(r));
    while(uM == nHL)
        uM = (int)(nHL*gsl_rng_uniform(r));
    
    // decide matrix entry that gets updated
    int i = (int)(dim*gsl_rng_uniform(r));
    while(i == dim)
        i = (int)(dim*gsl_rng_uniform(r));

    int j = (int)(dim*gsl_rng_uniform(r));
    while(j == dim)
        j = (int)(dim*gsl_rng_uniform(r));

    
    // generate random entry
    double x=0;
    double y=0;
    gsl_complex z;
    if(i==j)
    {
        // generate x uniformly between -1 and 1
        x = -1 + 2*gsl_rng_uniform(r);

        // metropolis scale
        x *= SCALE;

        z = gsl_complex_rect(x, 0.);
    }
    else
    {
        // generate x and y uniformly between -1 and 1
        x = -1 + 2*gsl_rng_uniform(r);
        y = -1 + 2*gsl_rng_uniform(r);

        // metropolis scale
        x *= SCALE;
        y *= SCALE;

        z = gsl_complex_rect(x, y);
    }
        
    
    // compute action difference
    double dS = deltaS(uM, i, j, z);


    // now evaluate the move
    if(dS < 0)
    {
        if(i != j)
        {
            gsl_complex w = MG(MAT[uM], i, j);
            gsl_matrix_complex_set(MAT[uM], i, j, CA(w, z));
            gsl_matrix_complex_set(MAT[uM], j, i, CA(CC(w), CC(z)));
        }
        else
        {
            double w = GSL_REAL(MG(MAT[uM], i, i));
            gsl_matrix_complex_set(MAT[uM], i, i, gsl_complex_rect(2.*x+w, 0.));
        }

        // update action
        S += dS;

        // move accepted
        accepted = 1;
    }
    else
    {
        double e = exp(-dS);
        double p = gsl_rng_uniform(r);

        if(e>p)
        {
            if(i != j)
            {
                gsl_complex w = MG(MAT[uM], i, j);
                gsl_matrix_complex_set(MAT[uM], i, j, CA(w, z));
                gsl_matrix_complex_set(MAT[uM], j, i, CA(CC(w), CC(z)));
            }
            else
            {
                double w = GSL_REAL(MG(MAT[uM], i, i));
                gsl_matrix_complex_set(MAT[uM], i, i, gsl_complex_rect(2.*x+w, 0.));
            }

            // update action
            S += dS;

            // move accepted
            accepted = 1;
        }
    }

    return accepted;
}

// a sweep is a collection of dim*dim*nHL proposed moves
// it returns the acceptance rate of that sweep
double sweep(double deltaS(int, int, int, gsl_complex), gsl_rng* r)
{
    double sum = 0.;
    double nmoves = ((double)dim*(dim-1.)/2. + dim)*nHL;
    for(int i=0; i<nmoves; i++)
        sum += move(deltaS, r);


    return sum/(double)(nmoves);
}


// a routine to automatically set the SCALE factor to give acceptance rate
// between minTarget and maxTarget (very rudimental and ugly, seems to work reasonably well)
void SCALE_autotune(double minTarget, double maxTarget, double deltaS(int, int, int, gsl_complex), gsl_rng* r)
{
    int n1 = 5;
    int n2 = 5;
    int m = 10;
    double limit = 1e-5;
    double ar = 0.;

    // initialize ar
    for(int i=0; i<n1; i++)
        ar += sweep(deltaS, r);
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
                    ar += sweep(deltaS, r);
                ar /= (double)n2;
            }
            delta /= 10.;
        }
    }
}


// complete simulation routine
// generates a 5 digit code that distinguishes simulations
// first, it prints on file the simulation data
// the thermalization part outputs two files "thermX.txt" (where X is 1 or 2) with the action value
// the simulation part outputs two files "simS.txt" "simHL.txt" with the action and the H and L matrices
// returns simulation code
char* simulation(double Sfunc(), double deltaS(int, int, int, gsl_complex), int rank, gsl_rng* r)
{
    // Initialize
    char* name_init = alloc_rank_filename(rank, "init");  
    init_data(name_init);
    free(name_init);
    FUNCTION(geom_check, CLIFF_P, CLIFF_Q)();
   
    // Generate unique filename
    char* code = generate_code(5, r);
    char* filename = alloc_rank_filename(rank, code);
    free(code);
    printf("%s\n", filename);
    char* name_data = alloc_coded_filename("data", filename);
    char* name_simS = alloc_coded_filename("simS", filename);
    char* name_simHL = alloc_coded_filename("simHL", filename);

    FILE* fdata = fopen(name_data, "w");
    FILE* fsimS = fopen(name_simS, "w");
    FILE* fsimHL = fopen(name_simHL, "w");

    free(name_data);
    free(name_simS);
    free(name_simHL);

    // Print simulation data
    print_data(fdata);


    // Simulation
    init_cold(Sfunc);
    SCALE_autotune(0.2, 0.3, deltaS, r);
    printf("Auto tuned SCALE: %lf\n", SCALE);
    print_time(fdata, "start simulation:");
    double ar = 0.;
    for(int i=0; i<Nsw; i++)
    {
        ar += sweep(deltaS, r);
        
        if(doit(i, ADJ))
            adjust();

        if(doit(i, GAP))
            print_simulation(fsimS, fsimHL);

    }
    print_time(fdata, "end simulation:");
    fprintf(fdata, "acceptance rate: %lf", ar/(double)Nsw);

    fclose(fsimS);
    fclose(fsimHL);
    fclose(fdata);
    simulation_free();

    printf("acceptance rate: %lf\n", ar/(double)Nsw);

    return filename;
}


// simulation routine for varying over matrix dimension and coupling constant G in a convenient way.
// for each matrix dimension generates a unique code XXXXX and three files:
// file "XXXXX_varG_data.txt": collects data on the range in which G varies, and total duration of the simulation
// file "XXXXX_varG_args.txt": collects list of codes needed to run multicode_analysis_main
// file "XXXXX_varG_G_args.txt": shows correspondence between G and code
// 
void multicode_wrapper(double Sfunc(), double deltaS(int, int, int, gsl_complex), double INCR_G, int REP_G, int rank, char* varG_code, gsl_rng* r)
{
    // copy init[rank].txt
    char* name_init = alloc_rank_filename(rank, "init");  
    FILE* finit = fopen(name_init, "r");
    if(finit == NULL)
    {
        printf("Error: %s not found\n", name_init);
        exit(EXIT_FAILURE);
    }
    
    int narg;
    int dim_, nH_, nL_, Nsw_;
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
    // initialize number of simulation sweeps
    narg += fscanf(finit, "%d", &Nsw_);
    
    if(narg < 6)
    {
        printf("Error: not enough data in %s\n", name_init);
        exit(EXIT_FAILURE);
    }
    fclose(finit);


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
    printf("geometry: (%d, %d)\n", CLIFF_P, CLIFF_Q);
    printf("Codename: %s\n", varG_code);
    printf("*********\n");

    // write same shenenigans
    fprintf(fvarG_data, "REP_G: %d\n", REP_G);
    fprintf(fvarG_data, "INCR_G: %lf\n", INCR_G);
    fprintf(fvarG_data, "G initial: %lf\n", G_);
    fprintf(fvarG_data, "G final: %lf\n", (G_ + (REP_G-1)*INCR_G));
    fprintf(fvarG_data, "dim: %d\n", dim_);
    fprintf(fvarG_data, "p: %d\n", CLIFF_P);
    fprintf(fvarG_data, "q: %d\n", CLIFF_Q);
    fprintf(fvarG_data, "Nsw: %d\n", Nsw_);
    fprintf(fvarG_data, "GAP: %d\n", GAP);
    fprintf(fvarG_data, "ADJ: %d\n", ADJ);

    // simulate with variable G
    print_time(fvarG_data, "start simulation:");
    for(int i=0; i<REP_G; i++)
    {
        char* code = simulation(Sfunc, deltaS, rank, r);
        fprintf(fvarG_args, "%s ", code);
        fprintf(fvarG_G_args, "%lf %s\n", G, code);
        printf("dim: %d,    G: %lf\n", dim, G);
        G += INCR_G;
        overwrite_init_file(name_init);
        free(code);
    }
    print_time(fvarG_data, "end simulation:");

    // free memory
    free(name_varG_data);
    free(name_varG_args);
    free(name_varG_G_args);
    fclose(fvarG_data);
    fclose(fvarG_args);
    fclose(fvarG_G_args);


    // restore init[rank].txt
    FILE* finit2 = fopen(name_init, "w");
    if(finit2 == NULL)
    {
        printf("Error: unable to write on %s\n", name_init);
        exit(EXIT_FAILURE);
    }
    fprintf(finit2, "%d\n", dim_);
    fprintf(finit2, "%d\n", nH_);
    fprintf(finit2, "%d\n", nL_);
    fprintf(finit2, "%lf\n", SCALE_);
    fprintf(finit2, "%lf\n", G_);
    fprintf(finit2, "%d\n", Nsw_);
    fclose(finit2);
    free(name_init);
}


void adjust()
{
    for(int i=0; i<nHL; i++)
    {
        make_hermitian(MAT[i]);

        if(i >= nH)
            traceless(MAT[i], dim);
    }
}
