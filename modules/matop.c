#define MATOP_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>
#include "matop.h"

#define REAL(a,i) (((double *) a)[2*(i)])
#define IMAG(a,i) (((double *) a)[2*(i)+1])
#define CONST_REAL(a,i) (((const double *) a)[2*(i)])
#define CONST_IMAG(a,i) (((const double *) a)[2*(i)+1])
#define mat_r(m, idx) (m->data[2*(idx)])
#define mat_i(m, idx) (m->data[2*(idx)+1])

void printmat_complex(gsl_matrix_complex *m, int n, int k)
{
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<k; j++)
        {
            gsl_complex z = gsl_matrix_complex_get(m, i, j);
            double x = GSL_REAL(z);
            double y = GSL_IMAG(z);
            int sgnx = x>=0;
            int sgny = y>=0;
            x = sqrt(x*x);
            y = sqrt(y*y);
            if(sgnx && sgny)
                printf("+%lf+i%lf  ", x, y);
            if(sgnx && !sgny)
                printf("+%lf-i%lf  ", x, y);
            if(!sgnx && sgny)
                printf("-%lf+i%lf  ", x, y);
            if(!sgnx && !sgny)
                printf("-%lf-i%lf  ", x, y);
        }
        
        printf("\n");
    }
}

void fprintmat_complex(gsl_matrix_complex *m, int n, int k, FILE* f)
{
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<k; j++)
        {
            gsl_complex z = gsl_matrix_complex_get(m, i, j);
            double x = GSL_REAL(z);
            double y = GSL_IMAG(z);
            fprintf(f, "%lf+i%lf  ", x, y);
        }
        
        fprintf(f, "\n");
    }
}

// generates nxn hermitian matrix H (mode 0) or traceless hermitian matrix L (mode 1)
// with entries between 0 and 1
void generate_HL(gsl_matrix_complex* m, int mode, int n, gsl_rng* r)
{
    // generate random hermitian matrix
    for(int i=0; i<n; i++)
    {

        // off diagonal part
        for(int j=0; j<i; j++)
        {
            double x, y;

            // generate x and y uniformly between -1 and 1
            x = -1 + 2*gsl_rng_uniform(r);
            y = -1 + 2*gsl_rng_uniform(r);

            gsl_matrix_complex_set(m, i, j, gsl_complex_rect(x,y));
            gsl_matrix_complex_set(m, j, i, gsl_complex_rect(x,-y));
        }


        // diagonal part
        double z;
        z = -1 + 2*gsl_rng_uniform(r);
        gsl_matrix_complex_set(m, i, i, gsl_complex_rect(z,0.));

    }



    // turn traceless if mode 1
    if(mode)
    {
        //antihermitian we don't need, the i can be included at any time
        //gsl_matrix_complex_scale(m, gsl_complex_rect(0.,1.));

        //render traceless
        traceless(m, n);
    }

    return;
}


void matrix_power(gsl_matrix_complex* m, int n, gsl_matrix_complex* res)
{
    int size = m->size1;
    gsl_matrix_complex* m_ = gsl_matrix_complex_alloc(size, size);
    gsl_matrix_complex_memcpy(m_, m);
    
    // res equals unity
    for(int i=0; i<size; i++)
    {
        for(int j=0; j<size; j++)
        {
            if(i==j)
                gsl_matrix_complex_set(res, i, j, GSL_COMPLEX_ONE);
            else
                gsl_matrix_complex_set(res, i, j, GSL_COMPLEX_ZERO);
        }
    }
    
    gsl_matrix_complex* r1 = gsl_matrix_complex_alloc(size, size);
    gsl_matrix_complex* mm = gsl_matrix_complex_alloc(size, size);
    while(n)
    {
        if(n % 2 == 1)
        {
            gsl_matrix_complex_memcpy(r1, res);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, r1, m_, GSL_COMPLEX_ZERO, res);
            n -= 1;
        }
        if(n)
        {
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, m_, m_, GSL_COMPLEX_ZERO, mm);
            gsl_matrix_complex_memcpy(m_, mm);

            n /= 2;
        }
    }
    gsl_matrix_complex_free(m_);
    gsl_matrix_complex_free(r1);
    gsl_matrix_complex_free(mm);
}

// calculates m[0]*m[1]*...*m[n-1] and outputs in res
void matrix_multiprod(gsl_matrix_complex** m, int n, gsl_matrix_complex* res)
{
    if(n<1)
    {
        printf("Error, give at least one matrix\n");
        exit(EXIT_FAILURE);
    }

    if(n==1)
        gsl_matrix_complex_memcpy(res, m[0]);

    else
    {
        int size = res->size1;
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, m[0], m[1], GSL_COMPLEX_ZERO, res);

        for(int i=2; i<n; i++)
        {
            gsl_matrix_complex* r1 = gsl_matrix_complex_alloc(size, size);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, res, m[i], GSL_COMPLEX_ZERO, r1);
            gsl_matrix_complex_memcpy(res, r1);
            gsl_matrix_complex_free(r1);
        }
    }
}








// dispherm gives a measure of how much
// the matrix m is distant from being hermitian
void dispherm(gsl_matrix_complex* m, int n)
{
    // diag stores the difference between 0 and
    // the imaginary part of the diagonal elements
    double *diag;
    diag = malloc(n*sizeof(double));
    for(int i=0; i<n; i++)
    {
        diag[i] = 0.0-GSL_IMAG(gsl_matrix_complex_get(m, i, i));
        printf("diag %d: %0.10lf\n", i, diag[i]);
    }

    // (re/im)offdiag stores the distance between
    // the real and imaginary part of opposite
    // off diagonal elements
    for(int i=0; i<n; i++)
    {
        for(int j=i+1; j<n; j++)
        {
            double temp1 = GSL_REAL(gsl_matrix_complex_get(m, i, j)); 
            double temp2 = GSL_REAL(gsl_matrix_complex_get(m, j, i)); 
            double temp3 = GSL_IMAG(gsl_matrix_complex_get(m, i, j)); 
            double temp4 = GSL_IMAG(gsl_matrix_complex_get(m, j, i)); 
            printf("re / im offdiag: %0.10lf / %0.10lf \n", temp1-temp2, temp3+temp4);
        }
    }


}



// outputs ordered eigenvalues of hermitiam matrix m into vector eval
// (eval must be already allocated)
void diag(gsl_matrix_complex* m, gsl_vector* eval, int n)
{

    // auxiliary matrix for diagonalization
    // (because gsl_eigen_herm, the diag algorithm, is destructive)
    gsl_matrix_complex* aux = gsl_matrix_complex_alloc(n, n); 
    gsl_matrix_complex_memcpy(aux, m);

    // allocation of workspace
    gsl_eigen_herm_workspace* w = gsl_eigen_herm_alloc(n);
    
    // diagonalization
    gsl_eigen_herm(aux, eval, w);
    gsl_matrix_complex_free(aux);
    gsl_eigen_herm_free(w);

    gsl_sort_vector(eval);

    return;
}


// outputs trace of matrix
gsl_complex trace(gsl_matrix_complex* m)
{
    int n1 = m->size1;
    int n2 = m->size2;
    int mtda = m->tda;
    if(n1 != n2)
    {
        printf("Error: cannot take trace of a non-square matrix\n");
        exit(EXIT_FAILURE);
    }

    double re=0;
    double im=0;
    for(int i=0; i<n1; i++)
    {
        re += m->data[2*(i*mtda + i)];
        im += m->data[2*(i*mtda + i)+1];
    }

    return gsl_complex_rect(re, im);
}

// outputs trace of hermitian matrix
double trace_herm(gsl_matrix_complex* m)
{
    int n1 = m->size1;
    int n2 = m->size2;
    int mtda = m->tda;
    if(n1 != n2)
    {
        printf("Error: cannot take trace of a non-square matrix\n");
        exit(EXIT_FAILURE);
    }

    double res=0;
    for(int i=0; i<n1; i++)
        res += m->data[2*(i*mtda + i)];

    return res;
}

// turn matrix traceless
void traceless(gsl_matrix_complex* m, int n)
{
    gsl_complex t = trace(m);
    t = gsl_complex_div_real(t, (double)n);
    for(int i=0; i<n; i++)
    {
        gsl_complex z = gsl_matrix_complex_get(m, i, i);
        gsl_matrix_complex_set(m, i, i, gsl_complex_sub(z, t));
    }
}


// tensor product, new and much faster version
// matrix a is (m x n), matrix b is (x x y)
// c must be (mx x ny)
// NOTE: ok so I don't understand when can tda be different from size2 (the number of columns),
// but apparently it can. But everytime I print it, tda=size2, so just think of tda
// as the number of columns.
// Anyway, to avoid potential errors, just rememebr this index identity:
// matrix element (i,j) corresponds to index (i*tda + j).
// If tda = size2, that formula is the usual row-major index identity; but since
// tda can be different from size2, the correct formula requires tda.
void tensor(gsl_matrix_complex* a, gsl_matrix_complex* b, gsl_matrix_complex* c) 
{
    int m = a->size1;
    int n = a->size2;
    int x = b->size1;
    int y = b->size2;
    int atda = a->tda;
    int btda = b->tda;
    int ctda = c->tda;

    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)
        {
            int idx1 = 2*(i*atda + j);
            double z_r = a->data[idx1];
            double z_i = a->data[idx1+1];
            for(int k=0; k<x; k++)
            {
                for(int l=0; l<y; l++)
                {
                    int idx2 = 2*(k*btda + l);
                    double w_r = b->data[idx2];
                    double w_i = b->data[idx2+1];

                    int idx3 = 2*((x*i + k)*ctda + y*j + l);
                    c->data[idx3] = z_r*w_r - z_i*w_i;
                    c->data[idx3+1] = z_i*w_r + z_r*w_i;
                }
            }
        }
    }
}

// faster implementation of commutator
void commutator(gsl_matrix_complex* m, gsl_matrix_complex* c)
{
    int n = m->size1;
    gsl_matrix_complex_set_zero(c);
    int mtda = m->tda;
    int ctda = c->tda;

    // build first term
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            int idx = 2*(i*mtda + j);
            double z_r = m->data[idx];
            double z_i = m->data[idx+1];
            for(int k=0; k<n; k++)
            {
                int idx1 = 2*((n*i+k)*ctda + n*j + k); 
                c->data[idx1] = z_r;
                c->data[idx1+1] = z_i;
            }
        }
    }

    // build second term
    for(int i=0; i<n; i++)
    {
        for(int k=0; k<n; k++)
        {
            for(int l=0; l<n; l++)
            {
                int idx = 2*(l*mtda + k);
                double z_r = m->data[idx];
                double z_i = m->data[idx+1];

                int idx1 = 2*((n*i+k)*ctda + n*i + l);
                c->data[idx1] -= z_r;
                c->data[idx1+1] -= z_i;
            }
        }
    }
}

// faster implementation of anticommutator
void anticommutator(gsl_matrix_complex* m, gsl_matrix_complex* c)
{
    int n = m->size1;
    gsl_matrix_complex_set_zero(c);
    int mtda = m->tda;
    int ctda = c->tda;

    // build first term
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            int idx = 2*(i*mtda + j);
            double z_r = m->data[idx];
            double z_i = m->data[idx+1];
            for(int k=0; k<n; k++)
            {
                int idx1 = 2*((n*i+k)*ctda + n*j + k); 
                c->data[idx1] = z_r;
                c->data[idx1+1] = z_i;
            }
        }
    }

    // build second term
    for(int i=0; i<n; i++)
    {
        for(int k=0; k<n; k++)
        {
            for(int l=0; l<n; l++)
            {
                int idx = 2*(l*mtda + k);
                double z_r = m->data[idx];
                double z_i = m->data[idx+1];

                int idx1 = 2*((n*i+k)*ctda + n*i + l);
                c->data[idx1] += z_r;
                c->data[idx1+1] += z_i;
            }
        }
    }
}


void make_hermitian(gsl_matrix_complex* a) 
{
    int m = a->size1;
    int n = a->size2;
    int atda = a->tda;
    if(m != n)
    {
        printf("Error: cannot hermitize a non-square matrix\n");
        return;
    }

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<=i; j++)
        {
            if(i != j)
            {
                int idx1 = 2*(i*atda + j);
                int idx2 = 2*(j*atda + i);
                double z1_r = a->data[idx1];
                double z1_i = a->data[idx1+1];
                double z2_r = a->data[idx2];
                double z2_i = -1.*(a->data[idx2+1]);

                double re = (z1_r + z2_r) / 2.;
                double im = (z1_i + z2_i) / 2.;
                a->data[idx1] = re;
                a->data[idx1+1] = im;
                a->data[idx2] = re;
                a->data[idx2+1] = -im;
            }
            else
            {
                int idx = 2*(i*atda + i);
                a->data[idx+1] = 0.;
            }
        }
    }
}


// ****** UNDERSTAND AND REWRITE THESE FUNCTIONS WHEN YOU HAVE TIME *******
// (they shoud make the computation faster, but they don't because the
// implementation is naive and BLAS is so fast)

double trace2(gsl_matrix_complex* m1, gsl_matrix_complex* m2)
{
    int n = m1->size1;
    int mtda = m1->tda;

    double res = 0.;
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<i; j++)
        {
            int idx = 2*(i*mtda + j);
            double re1 = m1->data[idx];
            double im1 = m1->data[idx+1];
            double re2 = m2->data[idx];
            double im2 = m2->data[idx+1];
            res += re1*re2 + im1*im2;
        }
    }

    res *= 2.;

    for(int i=0; i<n; i++)
    {
        int idx = 2*(i*mtda + i);
        double re1 = m1->data[idx];
        double re2 = m2->data[idx];
        res += re1*re2;
    }

    return res;
}

double trace2_self(gsl_matrix_complex* m)
{
    int n = m->size1;
    int mtda = m->tda;

    double res = 0.;
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<i; j++)
        {
            int idx = 2*(i*mtda + j);
            double re = m->data[idx];
            double im = m->data[idx+1];
            res += re*re + im*im;
        }
    }

    res *= 2.;

    for(int i=0; i<n; i++)
    {
        int idx = 2*(i*mtda + i);
        double re = m->data[idx];
        res += re*re;
    }

    return res;
}

gsl_complex trace3(gsl_matrix_complex* m1, gsl_matrix_complex* m2, gsl_matrix_complex* m3)
{
    int n = m1->size1;
    int mtda = m1->tda;

    double res_re = 0.;
    double res_im = 0.;

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            int idx1 = 2*(i*mtda + j);
            double a1 = m1->data[idx1];
            double b1 = m1->data[idx1+1];
            for(int k=0; k<n; k++)
            {
                //if(k != i)
                {
                    int idx2 = 2*(j*mtda + k);
                    int idx3 = 2*(k*mtda + i);
                    double a2 = m2->data[idx2];
                    double b2 = m2->data[idx2+1];
                    double a3 = m3->data[idx3];
                    double b3 = m3->data[idx3+1];
                    double A = (a1*a2-b1*b2);
                    double B = (a1*b2+a2*b1);

                    res_re += (A*a3 - B*b3);
                    res_im += (A*b3 + B*a3);
                }
            }
        }
    }

    return gsl_complex_rect(res_re, res_im);
}

gsl_complex trace4(gsl_matrix_complex* m1, gsl_matrix_complex* m2, gsl_matrix_complex* m3, gsl_matrix_complex* m4)
{
    int n = m1->size1;
    int mtda = m1->tda;

    double res_re = 0.;
    double res_im = 0.;

    for(int i=0; i<n; i++)
    {
        for(int k1=0; k1<n; k1++)
        {
            int idx1 = 2*(i*mtda + k1);
            double a1 = m1->data[idx1];
            double b1 = m1->data[idx1+1];
            for(int k2=0; k2<n; k2++)
            {
                int idx2 = 2*(k1*mtda + k2);
                double a2 = m2->data[idx2];
                double b2 = m2->data[idx2+1];
                double A = (a1*a2-b1*b2);
                double C = (a1*b2+a2*b1);
                for(int k3=0; k3<n; k3++)
                {
                    int idx3 = 2*(k2*mtda + k3);
                    int idx4 = 2*(k3*mtda + i);

                    double a3 = m3->data[idx3];
                    double b3 = m3->data[idx3+1];
                    double a4 = m4->data[idx4];
                    double b4 = m4->data[idx4+1];
                    
                    double B = (a3*a4-b3*b4);
                    double D = (a3*b4+a4*b3);
                    res_re += A*B - C*D;
                    res_im += B*C + A*D;
                }
            }
        }
    }

    return gsl_complex_rect(res_re, res_im);
}

gsl_complex MG3(gsl_matrix_complex* m1, gsl_matrix_complex* m2, gsl_matrix_complex* m3, int i, int j)
{
    double res_re = 0.;
    double res_im = 0.;
    
    int n = m1->size1;
    int mtda = m1->tda;

    for(int k1=0; k1<n; k1++)
    {
        int idx1 = 2*(i*mtda + k1);
        double a1 = m1->data[idx1];
        double b1 = m1->data[idx1+1];
        for(int k2=0; k2<n; k2++)
        {
            int idx2 = 2*(k1*mtda + k2);
            int idx3 = 2*(k2*mtda + j);
            double a2 = m2->data[idx2];
            double b2 = m2->data[idx2+1];
            double a3 = m3->data[idx3];
            double b3 = m3->data[idx3+1];
            double A = (a1*a2-b1*b2);
            double B = (a1*b2+a2*b1);

            res_re += (A*a3 - B*b3);
            res_im += (A*b3 + B*a3);
        }
    }

    return gsl_complex_rect(res_re, res_im);
}

gsl_complex MG2(gsl_matrix_complex* m1, gsl_matrix_complex* m2, int i, int j)
{
    double res_re = 0.;
    double res_im = 0.;
    
    int n = m1->size1;
    int mtda = m1->tda;

    for(int k=0; k<n; k++)
    {
        int idx1 = 2*(i*mtda + k);
        int idx2 = 2*(j*mtda + k);
        double a1 = m1->data[idx1];
        double b1 = m1->data[idx1+1];
        double a2 = m2->data[idx2];
        double b2 = m2->data[idx2+1];

        res_re += a1*a2 + b1*b2;
        res_im += a2*b1 - a1*b2;
    }

    return gsl_complex_rect(res_re, res_im);
}



