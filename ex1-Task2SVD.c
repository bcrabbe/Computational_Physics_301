#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
/*** SVD ***/


int main (void)
{
    FILE *FP;
    FP =fopen("/Users/ben/Documents/Phys-yr-3/computational_physics_2/ex1/SVD.txt","w");

    // clock_t t_before = time(NULL);
    int n=3;//size of matrix to be inverted

    // double k, step=0.1e-14;
    //for(n=2;n<=300;n++)//for looking at time scaling
    //for(k=0.00000011;k>=0;k-=step)//for looking at singular behavior
    //{

    /** fill this with the values on the matrix you want to invert **/
  /*  double a_values[] = { 1, 1, 1,
                          1, 2, -1,
                          2, 3, k, };

*/
    /** or fill it with random values **/
    int i;
    double a_values[(n*n)];
    srand ( time(NULL) );// this will change the random numbers each time
    for(i=1;i<=(n*n);i++)
    {
        a_values[i] = rand()%500;//will generate a random integer between 0 and 500
    }

    gsl_matrix_view A = gsl_matrix_view_array (a_values, n, n);

    /*** SVDecomp ***/
    gsl_vector * S = gsl_vector_alloc(n);
    gsl_vector * work = gsl_vector_alloc(n);
    gsl_matrix * V = gsl_matrix_alloc(n,n);

    gsl_linalg_SV_decomp(&A.matrix, V, S, work); // A becomes U,

    /*** makes a diagonal matrix out of 1/S, and checks its condition number ***/
    gsl_matrix *D = gsl_matrix_calloc(n,n);
    int j;
    double wmax=0, wmin=DBL_MAX;
    for(i=0; i<n; i++)
    {
        if(gsl_vector_get(S,i)>=wmax){wmax=gsl_vector_get(S,i);}
        if(gsl_vector_get(S,i)<=wmin){wmin=gsl_vector_get(S,i);}
        gsl_matrix_set(D,i,i, (1/gsl_vector_get(S, i)) );//produces diagonal matrix of 1/(values of S)
    }
    double condition_number = wmax/wmin;
    if(condition_number != condition_number){printf("\n!!!WARNING: Matrix is singular\n");}
    if((1/condition_number)<1e-12){printf("\n!!!WARNING: Matrix is close to singular ---> rounding errors occuring\n");}

    /*** computes A^-1 = V*D*U^T  ***/
    gsl_matrix * inverse = gsl_matrix_alloc(n,n);
    gsl_matrix *step = gsl_matrix_alloc(n,n);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, D, &A.matrix, 0.0, step);//Matrix multiplication...computes  step = D * U(transposed)
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0,V, step, 0.0, inverse);//computes A^-1 = V * step

    /*** prints inverse matrix to screen ***/
    printf("\n A^-1:\n");
    for( i=0;i<n;i++)
    {
        for( j=0; j<n;j++)
        {
            printf("%g\t",gsl_matrix_get(inverse,i,j));
        }
        printf("\n");
    }

    /*** calculates the ave. magnitude of elements of inv A ***/
 /*   double tot=0;
    for( i=0;i<n;i++)
    {
        for( j=0; j<n;j++)
        {
            tot+=fabs(gsl_matrix_get(inverse,i,j));
        }
    }
    double ave = tot/(n*n);

    fprintf(FP,"%0.14lf\t%lf\n",k,ave);
*/
    gsl_matrix_free (step);
    gsl_matrix_free (inverse);
    gsl_matrix_free (V);
    gsl_matrix_free (D);
    gsl_vector_free (S);
    gsl_vector_free (work);

    //clock_t t_after = time(NULL);
  //double dt = t_after - t_before;
 // fprintf(FP,"%d\t%lf\n",n,dt);
 // printf("%d\t%lf\n",n,dt);

    fclose(FP);
    return 0;
}
