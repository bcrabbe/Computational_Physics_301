#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
/*** LUD ***/
int
main (void)
{
    FILE *FP;
    FP =fopen("/Users/ben/Documents/Phys-yr-3/computational_physics_2/ex1/LUD.txt","w");
    //clock_t t_before = time(NULL);
    int n=3;//size of matrix to be inverted
  //  double k, step=1e-14;

   // for(n=2;n<=500;n++)//for looking at time scaling
  //  for(k=0.000000011;k>=0;k-=step)//for looking at singular behaviour
    {


    /** fill this with the values on the matrix you want to invert **/
  /*  double a_values[] = { 1, 1, 1,
                          1, 2, -1,
                          2, 3, k, };
*/
    /** or fill it with random values **/
    int i;
    double a_values[(n*n)];
    srand ( time(NULL) );// this will change the random numbers each time
    for(i=0;i<(n*n)-1;i++)
    {
        a_values[i] = rand()%500;//will generate a random integer between 0 and 500
    }

    /** LU Decomp **/
    int s;
    gsl_matrix_view A = gsl_matrix_view_array (a_values, n, n);
    gsl_permutation * p = gsl_permutation_alloc (n);
    gsl_linalg_LU_decomp (&A.matrix, p, &s);
    /** inverse calc **/
    gsl_matrix * inverse = gsl_matrix_alloc(n,n);
    gsl_linalg_LU_invert(&A.matrix, p,inverse);
    /** prints to screen**/
    printf ("A^-1 = \n");
    int j;
    for( i=0;i<n;i++)
    {
        for( j=0; j<n;j++)
        {
           printf("%f\t",gsl_matrix_get(inverse,i,j));
        }
        printf("\n");
    }

    /**  for singular test- calculates the ave. magnitude of elements of inv A **/
  /*  double tot=0;
    for( i=0;i<n;i++)
    {
        for( j=0; j<n;j++)
        {
            tot+=fabs(gsl_matrix_get(inverse,i,j));
        }
    }
    double ave = tot/(n*n);
    fprintf(FP,"%lf\t%lf\n",k,ave);
*/
    gsl_matrix_free (inverse);
    gsl_permutation_free (p);
    /*** for time test ***/
  //clock_t t_after = time(NULL);
  //double dt = t_after - t_before;
  //fprintf(FP,"%d\t%lf\n",n,dt);
   // printf("%d\t%lf\n",n,dt);
    }
    fclose(FP);
  return 0;
}
















