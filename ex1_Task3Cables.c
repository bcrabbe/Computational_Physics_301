#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>



int main()
{
    const int n=3;//number of cables
    const double Z=-10.0, mass=50.0; // the camera is 10m below the hanging points throughout...and has a mass of 50kg
    const double B=121.2435565, H=105.0; // the dimensions (base, height) of the triangle formed by the 3 hanging points
    const double g=9.80665;//m/s-2

    double x, y, step=0.5;//the camera's x and y position and the spacing of the points we calculate the tension at
    double theta1,theta2,theta3,phi1,phi2,phi3;
    double maxtension=0, mtx,mty;

    FILE * FP =fopen("/Users/ben/Documents/Phys-yr-3/computational_physics_2/ex1/camera-cables/output/tension.txt","w");
    /*** A*T=c ***/
    gsl_vector * T = gsl_vector_alloc(n);//vector for the tensions
    gsl_vector * c = gsl_vector_calloc(n);//set vector for the RHS to ...
    gsl_vector_set(c, 2, (mass*g)); // .. 0, 0, mg

    int s;//for LUD
    gsl_permutation * p = gsl_permutation_alloc (3);//for LUD

    for(y=0;y<=H;y+=step)
    {
        double xmin=(y-H)/1.732050808; //half the length of the base at of a triangle height H-y (the denom is sqrt(3) )
        double xmax=(H-y)/1.732050808;
        for(x=xmin;x<=xmax;x+=step)
        {
            /*** Computes the elements and then forms the coefficient matrix A ***/
            theta1=atan(-x/(H-y));
            theta2=atan(y/((B/2)+x));
            theta3=atan(y/((B/2)-x));
            phi1=atan(Z/sqrt(pow(x,2)+pow((H-y),2)));//denom is the x-y disp form point 1
            phi2=atan(Z/sqrt(pow((x+(B/2)),2)+pow(y,2)));//denom is the x-y disp form point 2
            phi3=atan(Z/sqrt(pow(((B/2)-x),2)+pow(y,2)));//denom is the x-y disp form point 3
            double a_values[]={sin(theta1), -cos(theta2), cos(theta3),
                        cos(theta1), -sin(theta2), -sin(theta3),
                        sin(phi1),    sin(phi2),    sin(phi3) };
            gsl_matrix_view A = gsl_matrix_view_array (a_values, 3, 3);

            /*** does LU decomp ***/

            gsl_linalg_LU_decomp (&A.matrix, p, &s);

            /*** Computes tensions ***/
            gsl_linalg_LU_solve(&A.matrix, p, c, T);

            /*** prints tensions at x,y to file ***/
            fprintf(FP,"%lf %lf %lf\n",x,y,fabs(gsl_vector_get(T,0)));

            /*** finds the max tension ***/
            if(fabs(gsl_vector_get(T,0))>maxtension)
            {
                maxtension=fabs(gsl_vector_get(T,0));
                mtx=x;
                mty=y;
            }
        }
    }
    printf("max tension = %lf at %lf, %lf", maxtension, mtx,mty);
    gsl_permutation_free(p);
    gsl_vector_free (c);
    gsl_vector_free (T);

    return 0;
}
