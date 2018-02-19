#define ACTION_C

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
#include <gsl/gsl_blas.h>
#include "global.h"
#include "macros.h"
#include "data.h"
#include "statistics.h"

#define ABS2(a) gsl_complex_abs2(a)
#define CA(a,b) gsl_complex_add(a, b)
#define CAR(a,b) gsl_complex_add_real(a, b)
#define CM(a,b) gsl_complex_mul(a, b)
#define CMR(a,b) gsl_complex_mul_real(a, b)
#define CC(a) gsl_complex_conjugate(a)
#define MG(m,i,j) gsl_matrix_complex_get(m, i, j)

double delta2(int uM, int I, int J, gsl_complex z)
{
    if(I != J)
        return 4.*dimG*dim*( 2.*GSL_REAL(CM(MG(MAT[uM], J, I), z) ) + ABS2(z) );
    else
    {
        double trM = trace_herm(MAT[uM]);
        return 8.*dimG*GSL_REAL(z)*( dim*(GSL_REAL(MG(MAT[uM], I, I)) + GSL_REAL(z)) + e[uM]*(trM + GSL_REAL(z)) );
    }
    
}

double delta4(int uM, int I, int J, gsl_complex z)
{
    double res = 0.;


    // D^3 dD part
    for(int i3=0; i3<nHL; i3++)
    {
        for(int i2=0; i2<nHL; i2++)
        {
            for(int i1=0; i1<=i3; i1++)
            {
                gsl_complex cliff = gamma_table[4][uM + nHL*(i3 + nHL*(i2 + nHL*i1))];

                if(GSL_REAL(cliff) != 0. || GSL_IMAG(cliff) != 0.)
                {
                    // alloc necessary matrices
                    gsl_matrix_complex* M1M2 = gsl_matrix_complex_alloc(dim, dim);
                    gsl_matrix_complex* M2M3 = gsl_matrix_complex_alloc(dim, dim);
                    gsl_matrix_complex* M1M3 = gsl_matrix_complex_alloc(dim, dim);
                    gsl_matrix_complex* M1M2M3 = gsl_matrix_complex_alloc(dim, dim);

                    // compute necessary matrix products
                    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i1], MAT[i2], GSL_COMPLEX_ZERO, M1M2);
                    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i2], MAT[i3], GSL_COMPLEX_ZERO, M2M3);
                    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i1], MAT[i3], GSL_COMPLEX_ZERO, M1M3);
                    gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i1], M2M3, GSL_COMPLEX_ZERO, M1M2M3);

                    // compute necessary traces
                    double trM1 = trace_herm(MAT[i1]);
                    double trM2 = trace_herm(MAT[i2]);
                    double trM3 = trace_herm(MAT[i3]);
                    double trM1M2 = GSL_REAL(trace(M1M2));
                    double trM2M3 = GSL_REAL(trace(M2M3));
                    double trM1M3 = GSL_REAL(trace(M1M3));
                    gsl_complex trM1M2M3 = trace(M1M2M3);
                    
                    // off-diagonal update
                    if(I != J)
                    {

                        // compute terms
                        // _______________________________________________________________________________________
                        gsl_complex T1 = CA(    CM(  MG(M1M2M3, J, I), z  ), CM(  MG(M1M2M3, I, J), CC(z)  )    );
                        T1 = CA(T1,  CMR(CC(T1), e[i1]*e[i2]*e[i3]*e[uM]));
                        T1 = CMR(T1, dim);

                        gsl_complex T2 = CA(    CM(  MG(M1M2, J, I), z  ), CM(  MG(M1M2, I, J), CC(z)  )    );
                        T2 = CA(CMR(T2, e[i3]),  CMR(CC(T2), e[i1]*e[i2]*e[uM]));
                        T2 = CMR(T2, trM3);
                        T1 = CA(T1, T2);

                        gsl_complex T3 = CA(    CM(  MG(M1M3, J, I), z  ), CM(  MG(M1M3, I, J), CC(z)  )    );
                        T3 = CA(CMR(T3, e[i2]),  CMR(CC(T3), e[i1]*e[i3]*e[uM]));
                        T3 = CMR(T3, trM2);
                        T1 = CA(T1, T3);

                        gsl_complex T4 = CA(    CM(  MG(M2M3, J, I), z  ), CM(  MG(M2M3, I, J), CC(z)  )    );
                        T4 = CA(CMR(T4, e[i1]),  CMR(CC(T4), e[i2]*e[i3]*e[uM]));
                        T4 = CMR(T4, trM1);
                        T1 = CA(T1, T4);

                        double T5 = trM1M2*(e[i1]*e[i2] + e[i3]*e[uM]);
                        T5 *= 2.*GSL_REAL(   CM(  MG(MAT[i3], J, I), z)  );
                        T1 = CAR(T1, T5);
                        
                        double T6 = trM2M3*(e[i2]*e[i3] + e[i1]*e[uM]);
                        T6 *= 2.*GSL_REAL(   CM(  MG(MAT[i1], J, I), z)  );
                        T1 = CAR(T1, T6);
                        
                        double T7 = trM1M3*(e[i1]*e[i3] + e[i2]*e[uM]);
                        T7 *= 2.*GSL_REAL(   CM(  MG(MAT[i2], J, I), z)  );
                        T1 = CAR(T1, T7);
                        //________________________________________________________________________________________
                        
                        

                        // add to total
                        if(i1 != i3)
                            res += 2.*GSL_REAL(CM(cliff, T1));
                        else
                            res += GSL_REAL(CM(cliff, T1));
                    }

                    
                    // diagonal update
                    else
                    {
                        
                        // compute terms
                        // _______________________________________________________________________________________
                        gsl_complex T1 = MG(M1M2M3, I, I);
                        T1 = CA(T1,  CMR(CC(T1), e[i1]*e[i2]*e[i3]*e[uM]));
                        T1 = CMR(T1, dim);

                        gsl_complex T2 = MG(M1M2, I, I);
                        T2 = CA(CMR(T2, e[i3]),  CMR(CC(T2), e[i1]*e[i2]*e[uM]));
                        T2 = CMR(T2, trM3);
                        T1 = CA(T1, T2);

                        gsl_complex T3 = MG(M1M3, I, I);
                        T3 = CA(CMR(T3, e[i2]),  CMR(CC(T3), e[i1]*e[i3]*e[uM]));
                        T3 = CMR(T3, trM2);
                        T1 = CA(T1, T3);

                        gsl_complex T4 = MG(M2M3, I, I);
                        T4 = CA(CMR(T4, e[i1]),  CMR(CC(T4), e[i2]*e[i3]*e[uM]));
                        T4 = CMR(T4, trM1);
                        T1 = CA(T1, T4);

                        double T5 = trM1M2*(e[i1]*e[i2] + e[i3]*e[uM]);
                        T5 *= GSL_REAL(MG(MAT[i3], I, I));
                        T1 = CAR(T1, T5);

                        double T6 = trM2M3*(e[i2]*e[i3] + e[i1]*e[uM]);
                        T6 *= GSL_REAL(MG(MAT[i1], I, I));
                        T1 = CAR(T1, T6);

                        double T7 = trM1M3*(e[i1]*e[i3] + e[i2]*e[uM]);
                        T7 *= GSL_REAL(MG(MAT[i2], I, I));
                        T1 = CAR(T1, T7);
                        
                        gsl_complex T8 = CA(CMR(CC(trM1M2M3), e[i1]*e[i2]*e[i3]),  CMR(trM1M2M3, e[uM]));
                        T1 = CA(T1, T8);
                        //________________________________________________________________________________________
                        
                        

                        // add to total
                        if(i1 != i3)
                            res += GSL_REAL(CM(cliff, T1))*4.*GSL_REAL(z);
                        else
                            res += GSL_REAL(CM(cliff, T1))*2.*GSL_REAL(z);
                    }

                    
                    // free memory
                    gsl_matrix_complex_free(M1M2);
                    gsl_matrix_complex_free(M2M3);
                    gsl_matrix_complex_free(M1M3);
                    gsl_matrix_complex_free(M1M2M3);
                }
            }
        }
    }

    res *= 4.;

   

    // D^2 dD^2 and D dD D dD term
    double temp = 0;
    for(int i=0; i<nHL; i++)
    {
        double cliff = GSL_REAL(gamma_table[4][uM + nHL*(i + nHL*(uM + nHL*i))]);
        
        // alloc necessary matrices
        gsl_matrix_complex* M1M1 = gsl_matrix_complex_alloc(dim, dim);

        // compute necessary matrix products
        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i], MAT[i], GSL_COMPLEX_ZERO, M1M1);

        // compute necessary traces
        double trM1 = trace_herm(MAT[i]);
        double trM1M1 = trace_herm(M1M1);

        // off-diagonal update
        if(I != J)
        {
            // compute terms D^2 dD^2
            // _______________________________________________________________________________________
            double T11 = 2.*dim*( GSL_REAL(MG(M1M1, I, I)) + GSL_REAL(MG(M1M1, J, J)) );

            double T21 = 4.*e[i]*trM1*( GSL_REAL(MG(MAT[i], I, I)) + GSL_REAL(MG(MAT[i], J, J)) );

            double T31 = GSL_REAL(  CM( MG(MAT[i], J, I), z)  );
            T31 *= T31*16.*e[i]*e[uM];
            //________________________________________________________________________________________
            
            // compute terms D dD D dD
            // _______________________________________________________________________________________
            
            double T12 = GSL_REAL(  CM( CM(MG(MAT[i], J, I), MG(MAT[i], J, I)), CM(z,z) )  );
            T12 += GSL_REAL( MG(MAT[i], I, I) ) * GSL_REAL( MG(MAT[i], J, J) ) * ABS2(z);
            T12 *= 4.*dim;
            
            double T22 = 4.*e[i]*trM1*( GSL_REAL(MG(MAT[i], I, I)) + GSL_REAL(MG(MAT[i], J, J)) );

            double T32 = 2.*GSL_REAL(  CM( MG(MAT[i], J, I), z)  );
            T32 *= T32*4.*e[i]*e[uM];
            //________________________________________________________________________________________
                    
                    
            
            // add to total
            temp += 2.*dimG*(ABS2(z)*(T11+T21+4.*trM1M1) + T31);
            temp += cliff*(T12 + ABS2(z)*(T22+4.*trM1M1) + T32);

        } 

        // diagonal update
        else
        {
            // compute terms D^2 dD^2
            // _______________________________________________________________________________________
            double T11 = 2.*dim*GSL_REAL(  MG(M1M1, I, I)  );

            double T21 = 4.*e[uM]*GSL_REAL(  MG(M1M1, I, I)  );

            double T31 = 4.*e[i]*trM1*GSL_REAL( MG(MAT[i], I, I) );

            double T41 = GSL_REAL(  MG(MAT[i], I, I)  );
            T41 *= T41*4.*e[i]*e[uM];
            //________________________________________________________________________________________
                    
            // compute terms D dD D dD
            // _______________________________________________________________________________________
            double T12 = GSL_REAL(  MG(MAT[i], I, I)  );
            T12 *= T12*2.*dim;

            double T22 = 4.*e[uM]*GSL_REAL(  MG(M1M1, I, I)  );

            double T32 = 4.*e[i]*trM1*GSL_REAL( MG(MAT[i], I, I) );

            double T42 = GSL_REAL(  MG(MAT[i], I, I)  );
            T42 *= T42*4.*e[i]*e[uM];
            //________________________________________________________________________________________
            
            // add to total
            temp += 8.*GSL_REAL(z)*GSL_REAL(z)*dimG*(T11+T21+T31+T41+2.*trM1M1);
            temp += 4.*GSL_REAL(z)*GSL_REAL(z)*cliff*(T12+T22+T32+T42+2.*trM1M1);
        }

        // free memory
        gsl_matrix_complex_free(M1M1);
    }

    res += 2.*temp;
    
                    

    // D dD^3 term

    // off-diangonal update
    if(I != J)
    {
        temp = 4.*dimG*(dim+6)*ABS2(z)*GSL_REAL(  CM( MG(MAT[uM], J, I), z )   );
        res += 4.*temp;
    }

    // diagonal update
    else
    {
        double trMum = trace_herm(MAT[uM]);
        double rez = 2.*GSL_REAL(z);
        temp = 2.*rez*rez*rez*dimG*(GSL_REAL(MG(MAT[uM], I, I))*(dim+3.*e[uM]+3.) + e[uM]*trMum);
        res += 4.*temp;
    }


    // dD^4 term

    // off-diangonal update
    if(I != J)
    {
        temp = dimG*4.*ABS2(z)*ABS2(z)*(dim+6.);
        res += temp;
    }
    
    // diagonal update
    else
    {
        double rez = GSL_REAL(z);
        temp = dimG*32.*(dim+3.+4*e[uM])*rez*rez*rez*rez;
        res += temp;
    }

    return res;
}


double delta4_BETA(int uM, int I, int J, gsl_complex z)
{
    double res = 0.;


    // D^3 dD part
    for(int i3=0; i3<nHL; i3++)
    {
        for(int i2=0; i2<nHL; i2++)
        {
            for(int i1=0; i1<=i3; i1++)
            {
                gsl_complex cliff = gamma_table[4][uM + nHL*(i3 + nHL*(i2 + nHL*i1))];

                if(GSL_REAL(cliff) != 0. || GSL_IMAG(cliff) != 0.)
                {
                    // compute necessary traces
                    double trM1 = trace_herm(MAT[i1]);
                    double trM2 = trace_herm(MAT[i2]);
                    double trM3 = trace_herm(MAT[i3]);
                    double trM1M2 = trace2(MAT[i1], MAT[i2]);
                    double trM2M3 = trace2(MAT[i2], MAT[i3]);
                    double trM1M3 = trace2(MAT[i1], MAT[i3]);
                    gsl_complex trM1M2M3 = trace3(MAT[i1], MAT[i2], MAT[i3]);
                    
                    // off-diagonal update
                    if(I != J)
                    {

                        // compute terms
                        // _______________________________________________________________________________________
                        gsl_complex T1 = CA(    CM(  MG3(MAT[i1], MAT[i2], MAT[i3], J, I), z  ), CM(  MG3(MAT[i1], MAT[i2], MAT[i3], I, J), CC(z)  )    );
                        T1 = CA(T1,  CMR(CC(T1), e[i1]*e[i2]*e[i3]*e[uM]));
                        T1 = CMR(T1, dim);

                        gsl_complex T2 = CA(    CM(  MG2(MAT[i1], MAT[i2], J, I), z  ), CM(  MG2(MAT[i1], MAT[i2], I, J), CC(z)  )    );
                        T2 = CA(CMR(T2, e[i3]),  CMR(CC(T2), e[i1]*e[i2]*e[uM]));
                        T2 = CMR(T2, trM3);
                        T1 = CA(T1, T2);

                        gsl_complex T3 = CA(    CM(  MG2(MAT[i1], MAT[i3], J, I), z  ), CM(  MG2(MAT[i1], MAT[i3], I, J), CC(z)  )    );
                        T3 = CA(CMR(T3, e[i2]),  CMR(CC(T3), e[i1]*e[i3]*e[uM]));
                        T3 = CMR(T3, trM2);
                        T1 = CA(T1, T3);

                        gsl_complex T4 = CA(    CM(  MG2(MAT[i2], MAT[i3], J, I), z  ), CM(  MG2(MAT[i2], MAT[i3], I, J), CC(z)  )    );
                        T4 = CA(CMR(T4, e[i1]),  CMR(CC(T4), e[i2]*e[i3]*e[uM]));
                        T4 = CMR(T4, trM1);
                        T1 = CA(T1, T4);

                        double T5 = trM1M2*(e[i1]*e[i2] + e[i3]*e[uM]);
                        T5 *= 2.*GSL_REAL(   CM(  MG(MAT[i3], J, I), z)  );
                        T1 = CAR(T1, T5);
                        
                        double T6 = trM2M3*(e[i2]*e[i3] + e[i1]*e[uM]);
                        T6 *= 2.*GSL_REAL(   CM(  MG(MAT[i1], J, I), z)  );
                        T1 = CAR(T1, T6);
                        
                        double T7 = trM1M3*(e[i1]*e[i3] + e[i2]*e[uM]);
                        T7 *= 2.*GSL_REAL(   CM(  MG(MAT[i2], J, I), z)  );
                        T1 = CAR(T1, T7);
                        //________________________________________________________________________________________
                        
                        

                        // add to total
                        if(i1 != i3)
                            res += 2.*GSL_REAL(CM(cliff, T1));
                        else
                            res += GSL_REAL(CM(cliff, T1));
                    }

                    
                    // diagonal update
                    else
                    {
                        
                        // compute terms
                        // _______________________________________________________________________________________
                        gsl_complex T1 = MG3(MAT[i1], MAT[i2], MAT[i3], I, I);
                        T1 = CA(T1,  CMR(CC(T1), e[i1]*e[i2]*e[i3]*e[uM]));
                        T1 = CMR(T1, dim);

                        gsl_complex T2 = MG2(MAT[i1], MAT[i2], I, I);
                        T2 = CA(CMR(T2, e[i3]),  CMR(CC(T2), e[i1]*e[i2]*e[uM]));
                        T2 = CMR(T2, trM3);
                        T1 = CA(T1, T2);

                        gsl_complex T3 = MG2(MAT[i1], MAT[i3], I, I);
                        T3 = CA(CMR(T3, e[i2]),  CMR(CC(T3), e[i1]*e[i3]*e[uM]));
                        T3 = CMR(T3, trM2);
                        T1 = CA(T1, T3);

                        gsl_complex T4 = MG2(MAT[i2], MAT[i3], I, I);
                        T4 = CA(CMR(T4, e[i1]),  CMR(CC(T4), e[i2]*e[i3]*e[uM]));
                        T4 = CMR(T4, trM1);
                        T1 = CA(T1, T4);

                        double T5 = trM1M2*(e[i1]*e[i2] + e[i3]*e[uM]);
                        T5 *= GSL_REAL(MG(MAT[i3], I, I));
                        T1 = CAR(T1, T5);

                        double T6 = trM2M3*(e[i2]*e[i3] + e[i1]*e[uM]);
                        T6 *= GSL_REAL(MG(MAT[i1], I, I));
                        T1 = CAR(T1, T6);

                        double T7 = trM1M3*(e[i1]*e[i3] + e[i2]*e[uM]);
                        T7 *= GSL_REAL(MG(MAT[i2], I, I));
                        T1 = CAR(T1, T7);
                        
                        gsl_complex T8 = CA(CMR(CC(trM1M2M3), e[i1]*e[i2]*e[i3]),  CMR(trM1M2M3, e[uM]));
                        T1 = CA(T1, T8);
                        //________________________________________________________________________________________
                        
                        

                        // add to total
                        if(i1 != i3)
                            res += GSL_REAL(CM(cliff, T1))*4.*GSL_REAL(z);
                        else
                            res += GSL_REAL(CM(cliff, T1))*2.*GSL_REAL(z);
                    }
                }
            }
        }
    }

    res *= 4.;

   

    // D^2 dD^2 and D dD D dD term
    double temp = 0;
    for(int i=0; i<nHL; i++)
    {
        double cliff = GSL_REAL(gamma_table[4][uM + nHL*(i + nHL*(uM + nHL*i))]);
        
        // compute necessary traces
        double trM1 = trace_herm(MAT[i]);
        double trM1M1 = trace2_self(MAT[i]);

        // off-diagonal update
        if(I != J)
        {
            // compute terms D^2 dD^2
            // _______________________________________________________________________________________
            double T11 = 2.*dim*( GSL_REAL(MG2(MAT[i], MAT[i], I, I)) + GSL_REAL(MG2(MAT[i], MAT[i], J, J)) );

            double T21 = 4.*e[i]*trM1*( GSL_REAL(MG(MAT[i], I, I)) + GSL_REAL(MG(MAT[i], J, J)) );

            double T31 = GSL_REAL(  CM( MG(MAT[i], J, I), z)  );
            T31 *= T31*16.*e[i]*e[uM];
            //________________________________________________________________________________________
            
            // compute terms D dD D dD
            // _______________________________________________________________________________________
            
            double T12 = GSL_REAL(  CM( CM(MG(MAT[i], J, I), MG(MAT[i], J, I)), CM(z,z) )  );
            T12 += GSL_REAL( MG(MAT[i], I, I) ) * GSL_REAL( MG(MAT[i], J, J) ) * ABS2(z);
            T12 *= 4.*dim;
            
            double T22 = 4.*e[i]*trM1*( GSL_REAL(MG(MAT[i], I, I)) + GSL_REAL(MG(MAT[i], J, J)) );

            double T32 = 2.*GSL_REAL(  CM( MG(MAT[i], J, I), z)  );
            T32 *= T32*4.*e[i]*e[uM];
            //________________________________________________________________________________________
                    
                    
            
            // add to total
            temp += 2.*dimG*(ABS2(z)*(T11+T21+4.*trM1M1) + T31);
            temp += cliff*(T12 + ABS2(z)*(T22+4.*trM1M1) + T32);

        } 

        // diagonal update
        else
        {
            // compute terms D^2 dD^2
            // _______________________________________________________________________________________
            double T11 = 2.*dim*GSL_REAL(  MG2(MAT[i], MAT[i], I, I)  );

            double T21 = 4.*e[uM]*GSL_REAL(  MG2(MAT[i], MAT[i], I, I)  );

            double T31 = 4.*e[i]*trM1*GSL_REAL( MG(MAT[i], I, I) );

            double T41 = GSL_REAL(  MG(MAT[i], I, I)  );
            T41 *= T41*4.*e[i]*e[uM];
            //________________________________________________________________________________________
                    
            // compute terms D dD D dD
            // _______________________________________________________________________________________
            double T12 = GSL_REAL(  MG(MAT[i], I, I)  );
            T12 *= T12*2.*dim;

            double T22 = 4.*e[uM]*GSL_REAL(  MG2(MAT[i], MAT[i], I, I)  );

            double T32 = 4.*e[i]*trM1*GSL_REAL( MG(MAT[i], I, I) );

            double T42 = GSL_REAL(  MG(MAT[i], I, I)  );
            T42 *= T42*4.*e[i]*e[uM];
            //________________________________________________________________________________________
            
            // add to total
            temp += 8.*GSL_REAL(z)*GSL_REAL(z)*dimG*(T11+T21+T31+T41+2.*trM1M1);
            temp += 4.*GSL_REAL(z)*GSL_REAL(z)*cliff*(T12+T22+T32+T42+2.*trM1M1);
        }

    }

    res += 2.*temp;
    
                    

    // D dD^3 term

    // off-diangonal update
    if(I != J)
    {
        temp = 4.*dimG*(dim+6)*ABS2(z)*GSL_REAL(  CM( MG(MAT[uM], J, I), z )   );
        res += 4.*temp;
    }

    // diagonal update
    else
    {
        double trMum = trace_herm(MAT[uM]);
        double rez = 2.*GSL_REAL(z);
        temp = 2.*rez*rez*rez*dimG*(GSL_REAL(MG(MAT[uM], I, I))*(dim+3.*e[uM]+3.) + e[uM]*trMum);
        res += 4.*temp;
    }


    // dD^4 term

    // off-diangonal update
    if(I != J)
    {
        temp = dimG*4.*ABS2(z)*ABS2(z)*(dim+6.);
        res += temp;
    }
    
    // diagonal update
    else
    {
        double rez = GSL_REAL(z);
        temp = dimG*32.*(dim+3.+4*e[uM])*rez*rez*rez*rez;
        res += temp;
    }

    return res;
}




double dirac2()
{
    double res = 0.;
    for(int i=0; i<nHL; i++)
    {
        gsl_matrix_complex* MM = gsl_matrix_complex_alloc(dim, dim);
        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i], MAT[i], GSL_COMPLEX_ZERO, MM);
        double trMM = trace_herm(MM);
        double trM = trace_herm(MAT[i]);
        gsl_matrix_complex_free(MM);

        res += (dim*trMM + e[i]*trM*trM);
    }

    return 2.*dimG*res;
}



double dirac4()
{
    double res = 0.;
    int* i = malloc(4*sizeof(int));


    for(i[3]=0; i[3]<nHL; i[3]++)
    {
        for(i[2]=0; i[2]<nHL; i[2]++)
        {
            for(i[1]=0; i[1]<=i[2]; i[1]++)
            {
                for(i[0]=0; i[0]<=i[3]; i[0]++)
                {
                    if(i[0]==i[1] && i[1]==i[2] && i[2]==i[3])
                    {
                        // alloc matrix products (trace needed)
                        gsl_matrix_complex* M0M0M0M0 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M0M0 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M0 = gsl_matrix_complex_alloc(dim, dim);
                            
                        // compute matrix products
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], MAT[i[0]], GSL_COMPLEX_ZERO, M0M0);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], M0M0, GSL_COMPLEX_ZERO, M0M0M0);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], M0M0M0, GSL_COMPLEX_ZERO, M0M0M0M0);

                        // compute traces
                        double trM0 = trace_herm(MAT[i[0]]);
                        double trM0M0 = trace_herm(M0M0);
                        double trM0M0M0 = trace_herm(M0M0M0);
                        double trM0M0M0M0 = trace_herm(M0M0M0M0);
                        
                        // free memory 
                        gsl_matrix_complex_free(M0M0M0M0);
                        gsl_matrix_complex_free(M0M0M0);
                        gsl_matrix_complex_free(M0M0);
                        
                        // add to total
                        double temp = 0;

                        // tr4 term
                        temp += dim*2.*trM0M0M0M0;

                        // tr3tr1 term
                        temp += 8.*e[i[0]]*trM0M0M0*trM0;

                        // tr2tr2 term
                        temp += 6.*trM0M0*trM0M0;
                        
                        res += dimG*temp;
                    }


                    
                    else if(i[0] == i[3] && i[1] == i[2])
                    {
                        // alloc matrix products (trace needed)
                        gsl_matrix_complex* M1M0M0M1 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M1M0 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M1M0M1 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M0 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M1M1 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M1 = gsl_matrix_complex_alloc(dim, dim);
                            
                        // compute matrix products
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], MAT[i[0]], GSL_COMPLEX_ZERO, M0M0);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[1]], MAT[i[1]], GSL_COMPLEX_ZERO, M1M1);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], MAT[i[1]], GSL_COMPLEX_ZERO, M0M1);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[1]], M0M1, GSL_COMPLEX_ZERO, M1M0M1);
                        gsl_blas_zhemm(CblasRight, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], M0M1, GSL_COMPLEX_ZERO, M0M1M0);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, M0M0, M1M1, GSL_COMPLEX_ZERO, M1M0M0M1);

                        // compute traces
                        double trM0 = trace_herm(MAT[i[0]]);
                        double trM1 = trace_herm(MAT[i[1]]);
                        double trM0M0 = trace_herm(M0M0);
                        double trM1M1 = trace_herm(M1M1);
                        double trM0M1 = GSL_REAL(trace(M0M1));
                        double trM0M1M0 = trace_herm(M0M1M0);
                        double trM1M0M1 = trace_herm(M1M0M1);
                        double trM1M0M0M1 = trace_herm(M1M0M0M1);
                        
                        // free memory 
                        gsl_matrix_complex_free(M1M0M0M1);
                        gsl_matrix_complex_free(M0M1M0);
                        gsl_matrix_complex_free(M1M0M1);
                        gsl_matrix_complex_free(M0M0);
                        gsl_matrix_complex_free(M1M1);
                        gsl_matrix_complex_free(M0M1);

                        // add to total
                        double temp = 0;

                        // tr4 term
                        temp += dim*2.*trM1M0M0M1;

                        // tr3tr1 term
                        temp += 4.*e[i[0]]*trM1M0M1*trM0;
                        temp += 4.*e[i[1]]*trM0M1M0*trM1;

                        // tr2tr2 term
                        temp += 4.*e[i[0]]*e[i[1]]*trM0M1*trM0M1;
                        temp += 2.*trM0M0*trM1M1;
                        
                        res += dimG*temp;
                    }
                    

                    else
                    {
                        gsl_complex cliff = gamma_table[4][i[3] + nHL*(i[2] + nHL*(i[1] + nHL*i[0]))];
                        
                        if(GSL_REAL(cliff) != 0. || GSL_IMAG(cliff) != 0.)
                        {
                            // alloc matrix products (trace needed)
                            gsl_matrix_complex* M0M1M2M3 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M0M1M2 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M0M1M3 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M0M2M3 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M1M2M3 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M0M1 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M0M2 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M0M3 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M1M2 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M1M3 = gsl_matrix_complex_alloc(dim, dim);
                            gsl_matrix_complex* M2M3 = gsl_matrix_complex_alloc(dim, dim);

                            // compute matrix products
                            gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], MAT[i[1]], GSL_COMPLEX_ZERO, M0M1);
                            gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], MAT[i[2]], GSL_COMPLEX_ZERO, M0M2);
                            gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], MAT[i[3]], GSL_COMPLEX_ZERO, M0M3);
                            gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[1]], MAT[i[2]], GSL_COMPLEX_ZERO, M1M2);
                            gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[1]], MAT[i[3]], GSL_COMPLEX_ZERO, M1M3);
                            gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[2]], MAT[i[3]], GSL_COMPLEX_ZERO, M2M3);
                            gsl_blas_zhemm(CblasRight, CblasUpper, GSL_COMPLEX_ONE, MAT[i[2]], M0M1, GSL_COMPLEX_ZERO, M0M1M2);
                            gsl_blas_zhemm(CblasRight, CblasUpper, GSL_COMPLEX_ONE, MAT[i[3]], M0M1, GSL_COMPLEX_ZERO, M0M1M3);
                            gsl_blas_zhemm(CblasRight, CblasUpper, GSL_COMPLEX_ONE, MAT[i[3]], M0M2, GSL_COMPLEX_ZERO, M0M2M3);
                            gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[1]], M2M3, GSL_COMPLEX_ZERO, M1M2M3);
                            gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], M1M2M3, GSL_COMPLEX_ZERO, M0M1M2M3);

                            // compute traces
                            gsl_complex trM0M1M2M3 = trace(M0M1M2M3);
                            gsl_complex trM0M1M2 = trace(M0M1M2);
                            gsl_complex trM0M1M3 = trace(M0M1M3);
                            gsl_complex trM0M2M3 = trace(M0M2M3);
                            gsl_complex trM1M2M3 = trace(M1M2M3);
                            double trM0M1 = GSL_REAL(trace(M0M1));
                            double trM0M2 = GSL_REAL(trace(M0M2));
                            double trM0M3 = GSL_REAL(trace(M0M3));
                            double trM1M2 = GSL_REAL(trace(M1M2));
                            double trM1M3 = GSL_REAL(trace(M1M3));
                            double trM2M3 = GSL_REAL(trace(M2M3));
                            double trM0 = trace_herm(MAT[i[0]]);
                            double trM1 = trace_herm(MAT[i[1]]);
                            double trM2 = trace_herm(MAT[i[2]]);
                            double trM3 = trace_herm(MAT[i[3]]);
                            
                            // free memory
                            gsl_matrix_complex_free(M0M1M2M3);
                            gsl_matrix_complex_free(M0M1M2);
                            gsl_matrix_complex_free(M0M1M3);
                            gsl_matrix_complex_free(M0M2M3);
                            gsl_matrix_complex_free(M1M2M3);
                            gsl_matrix_complex_free(M0M1);
                            gsl_matrix_complex_free(M0M2);
                            gsl_matrix_complex_free(M0M3);
                            gsl_matrix_complex_free(M1M2);
                            gsl_matrix_complex_free(M1M3);
                            gsl_matrix_complex_free(M2M3);


                            // add to total

                            // tr4 terms
                            gsl_complex T1 = CA(  trM0M1M2M3, CMR(CC(trM0M1M2M3), e[i[0]]*e[i[1]]*e[i[2]]*e[i[3]])  );
                            res += 2.*dim*GSL_REAL( CM(cliff, T1) );

                            // tr3tr1 terms
                            gsl_complex T3 = CA(  CMR(trM0M1M2, e[i[3]]), CMR(CC(trM0M1M2), e[i[0]]*e[i[1]]*e[i[2]])  );
                            T3 = CMR(T3, trM3);

                            gsl_complex T4 = CA(  CMR(trM0M1M3, e[i[2]]), CMR(CC(trM0M1M3), e[i[0]]*e[i[1]]*e[i[3]])  );
                            T4 = CMR(T4, trM2);
                            T3 = CA(T3, T4);

                            gsl_complex T5 = CA(  CMR(trM0M2M3, e[i[1]]), CMR(CC(trM0M2M3), e[i[0]]*e[i[2]]*e[i[3]])  );
                            T5 = CMR(T5, trM1);
                            T3 = CA(T3, T5);

                            gsl_complex T6 = CA(  CMR(trM1M2M3, e[i[0]]), CMR(CC(trM1M2M3), e[i[1]]*e[i[2]]*e[i[3]])  );
                            T6 = CMR(T6, trM0);
                            T3 = CA(T3, T6);

                            res += 2.*GSL_REAL(  CM(cliff, T3)  );

                            // tr2tr2 terms
                            double T7 = trM0M1*trM2M3*(e[i[0]]*e[i[1]] + e[i[2]]*e[i[3]]);
                            double T8 = trM0M2*trM1M3*(e[i[0]]*e[i[2]] + e[i[1]]*e[i[3]]);
                            double T9 = trM0M3*trM1M2*(e[i[0]]*e[i[3]] + e[i[1]]*e[i[2]]);

                            res += 2.*GSL_REAL(cliff)*(T7+T8+T9);
                        }
                    }
                }
            }
        }
    }
    
    for(i[3]=0; i[3]<nHL; i[3]++)
    {
        for(i[1]=0; i[1]<nHL; i[1]++)
        {
            for(i[2]=0; i[2]<i[1]; i[2]++)
            {
                for(i[0]=0; i[0]<i[3]; i[0]++)
                {
                    gsl_complex cliff = gamma_table[4][i[3] + nHL*(i[2] + nHL*(i[1] + nHL*i[0]))];
                    
                    if(GSL_REAL(cliff) != 0. || GSL_IMAG(cliff) != 0.)
                    {
                        // alloc matrix products (trace needed)
                        gsl_matrix_complex* M0M1M2M3 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M1M2 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M1M3 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M2M3 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M1M2M3 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M1 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M2 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M0M3 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M1M2 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M1M3 = gsl_matrix_complex_alloc(dim, dim);
                        gsl_matrix_complex* M2M3 = gsl_matrix_complex_alloc(dim, dim);

                        // compute matrix products
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], MAT[i[1]], GSL_COMPLEX_ZERO, M0M1);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], MAT[i[2]], GSL_COMPLEX_ZERO, M0M2);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], MAT[i[3]], GSL_COMPLEX_ZERO, M0M3);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[1]], MAT[i[2]], GSL_COMPLEX_ZERO, M1M2);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[1]], MAT[i[3]], GSL_COMPLEX_ZERO, M1M3);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[2]], MAT[i[3]], GSL_COMPLEX_ZERO, M2M3);
                        gsl_blas_zhemm(CblasRight, CblasUpper, GSL_COMPLEX_ONE, MAT[i[2]], M0M1, GSL_COMPLEX_ZERO, M0M1M2);
                        gsl_blas_zhemm(CblasRight, CblasUpper, GSL_COMPLEX_ONE, MAT[i[3]], M0M1, GSL_COMPLEX_ZERO, M0M1M3);
                        gsl_blas_zhemm(CblasRight, CblasUpper, GSL_COMPLEX_ONE, MAT[i[3]], M0M2, GSL_COMPLEX_ZERO, M0M2M3);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[1]], M2M3, GSL_COMPLEX_ZERO, M1M2M3);
                        gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, MAT[i[0]], M1M2M3, GSL_COMPLEX_ZERO, M0M1M2M3);

                        // compute traces
                        gsl_complex trM0M1M2M3 = trace(M0M1M2M3);
                        gsl_complex trM0M1M2 = trace(M0M1M2);
                        gsl_complex trM0M1M3 = trace(M0M1M3);
                        gsl_complex trM0M2M3 = trace(M0M2M3);
                        gsl_complex trM1M2M3 = trace(M1M2M3);
                        double trM0M1 = GSL_REAL(trace(M0M1));
                        double trM0M2 = GSL_REAL(trace(M0M2));
                        double trM0M3 = GSL_REAL(trace(M0M3));
                        double trM1M2 = GSL_REAL(trace(M1M2));
                        double trM1M3 = GSL_REAL(trace(M1M3));
                        double trM2M3 = GSL_REAL(trace(M2M3));
                        double trM0 = trace_herm(MAT[i[0]]);
                        double trM1 = trace_herm(MAT[i[1]]);
                        double trM2 = trace_herm(MAT[i[2]]);
                        double trM3 = trace_herm(MAT[i[3]]);
                        
                        // free memory
                        gsl_matrix_complex_free(M0M1M2M3);
                        gsl_matrix_complex_free(M0M1M2);
                        gsl_matrix_complex_free(M0M1M3);
                        gsl_matrix_complex_free(M0M2M3);
                        gsl_matrix_complex_free(M1M2M3);
                        gsl_matrix_complex_free(M0M1);
                        gsl_matrix_complex_free(M0M2);
                        gsl_matrix_complex_free(M0M3);
                        gsl_matrix_complex_free(M1M2);
                        gsl_matrix_complex_free(M1M3);
                        gsl_matrix_complex_free(M2M3);


                        // add to total

                        // tr4 terms
                        gsl_complex T1 = CA(  trM0M1M2M3, CMR(CC(trM0M1M2M3), e[i[0]]*e[i[1]]*e[i[2]]*e[i[3]])  );
                        res += dim*2.*GSL_REAL( CM(cliff, T1) );

                        // tr3tr1 terms
                        gsl_complex T3 = CA(  CMR(trM0M1M2, e[i[3]]), CMR(CC(trM0M1M2), e[i[0]]*e[i[1]]*e[i[2]])  );
                        T3 = CMR(T3, trM3);

                        gsl_complex T4 = CA(  CMR(trM0M1M3, e[i[2]]), CMR(CC(trM0M1M3), e[i[0]]*e[i[1]]*e[i[3]])  );
                        T4 = CMR(T4, trM2);
                        T3 = CA(T3, T4);

                        gsl_complex T5 = CA(  CMR(trM0M2M3, e[i[1]]), CMR(CC(trM0M2M3), e[i[0]]*e[i[2]]*e[i[3]])  );
                        T5 = CMR(T5, trM1);
                        T3 = CA(T3, T5);

                        gsl_complex T6 = CA(  CMR(trM1M2M3, e[i[0]]), CMR(CC(trM1M2M3), e[i[1]]*e[i[2]]*e[i[3]])  );
                        T6 = CMR(T6, trM0);
                        T3 = CA(T3, T6);

                        res += 2.*GSL_REAL(  CM(cliff, T3)  );

                        // tr2tr2 terms
                        double T7 = trM0M1*trM2M3*(e[i[0]]*e[i[1]] + e[i[2]]*e[i[3]]);
                        double T8 = trM0M2*trM1M3*(e[i[0]]*e[i[2]] + e[i[1]]*e[i[3]]);
                        double T9 = trM0M3*trM1M2*(e[i[0]]*e[i[3]] + e[i[1]]*e[i[2]]);

                        res += 2.*GSL_REAL(cliff)*(T7+T8+T9);
                    }
                }
            }
        }
    }

    free(i);
    return res;
}



















