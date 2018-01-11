#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include "fileop.h"
#include "statistics.h"
#include "global.h"

#define NUM 4
#define REP 26

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        printf("gimme dimension \n");
        exit(EXIT_FAILURE);
    }
    
    // tells me if a folder has been specified
    int folder = -2 + argc;

    // allocate code filename and output filename
    char* args;
    char* output;
    if(folder)
    {
        char* code = alloc_folder_filename(argv[1], argv[2]);
        args = alloc_coded_filename("args", code);
        output = alloc_coded_filename("error", code);
        free(code);
    }
    else
    {
        args = alloc_coded_filename("args", argv[1]);
        output = alloc_coded_filename("erorr", argv[1]);
    }

    // open code file and output file
    FILE* fargs = fopen(args, "r");
    if(fargs == NULL)
    {
        printf("Error: unable to open codes file: %s\n", args);
        exit(EXIT_FAILURE);
    }

    FILE* foutput = fopen(output, "w");
    if(foutput == NULL)
    {
        printf("Error: unable to open output file: %s\n", output);
        exit(EXIT_FAILURE);
    }


    char* read = malloc(200*sizeof(char));

    double G[NUM][REP];
    double var[NUM][REP];

    for(int i=0; i<NUM; i++) 
    {
        int t = fscanf(fargs, "%s", read);
        if(!t)
        {
            printf("Error: failed to read %s\n", read);
            exit(EXIT_FAILURE);
        }
        // open data file and simulation file
        char* avg; 
        if(folder)
        {
            char* code;
            code = alloc_folder_filename(read, argv[2]);
            avg = alloc_coded_filename("varG_average", code);
            free(code);
        }
        else
            avg = alloc_coded_filename("varG_average", read);
        
        FILE* favg = fopen(avg, "r");

        if(favg == NULL)
        {
            printf("Error: unable to read input file %s\n", read);
            exit(EXIT_FAILURE);
        }

        printf("Processing: %s\n", read); 

        for(int j=0; j<REP; j++)
        {
            int r;
            r = fscanf(favg, "%lf", &G[i][j]);
            if(r==0)
            {
                printf("Error: not able to read coupling const\n");
                exit(EXIT_FAILURE);
            }

            int r1 = 0;
            double temp;
            for(int k=0; k<5; k++)
            {
                r1 += fscanf(favg, "%lf", &temp);
            }
            if(r1 != 5)
            {
                printf("Error: not enough data in %s\n", avg);
                exit(EXIT_FAILURE);
            }

            r = fscanf(favg, "%lf", &var[i][j]);
            if(r==0)
            {
                printf("Error: not able to read error\n");
                exit(EXIT_FAILURE);
            }
            
            int r2 = 0;
            for(int k=0; k<12; k++)
            {
                r2 += fscanf(favg, "%lf", &temp);
            }
            if(r2 != 12)
            {
                printf("Error: not enough data in %s\n", avg);
                exit(EXIT_FAILURE);
            }
        }
        free(avg);
        fclose(favg);
    }

    double G_t[REP][NUM];
    double var_t[REP][NUM];

    for(int i=0; i<REP; i++)
    {
        for(int j=0; j<NUM; j++)
        {
            G_t[i][j] = G[j][i];
            var_t[i][j] = var[j][i];
        }
    }

    for(int i=0; i<REP; i++)
    {
        double meanG = gsl_stats_mean(G_t[i], 1, NUM);
        double meanvar = gsl_stats_mean(var_t[i], 1, NUM);
        double varvar = gsl_stats_variance(var_t[i], 1, NUM);
        fprintf(foutput, "%.15lf %.15lf %.15lf %.15lf\n", meanG, meanvar, sqrt(varvar/(double)NUM), varvar);
    }

    free(args);
    free(output);
    fclose(fargs);
    fclose(foutput);
    free(read);
}




