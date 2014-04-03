#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <math.h>
# define VN 1
# define DIRICHLET 2
//how to i pass the gsl_vector struct ??
//fprintVector(gsl_vector * phiT);

static int L = 50;//length of rod is 50 cm
static int N = 26;//number of nodes
static int Case = VN;//if 1 there is no heat loss from the far end of the poker i.e. Von Neumann BCs
//and
//if 2 the far end of the poker is immersed in a block of ice at 0◦C. i.e. Dirichlet boundary conds
int main ()
{
    double h = L/(N-1);//∆x (m)
    double dt = 0.01;
    double t=0;

    /*** Set material constants ******************/
    double thermalConductivity = 0.059;// W/cm/K
    double specificHeatCapacity = 450;// J/kg/K
    double density = 7900e-6;// kg/cm^3
    double alpha = thermalConductivity/(specificHeatCapacity*density);//the thermal diffusivity

    /*** form the equations ******************/
    gsl_vector * phiT = gsl_vector_alloc(N+2);//contains an element for every node of the rod and 2 to comunicate bcs
    gsl_vector_set(phiT,0,1000);//sets 1000 boundary condition
    int row;
    for(row=1; row<=N; row++)
    {
        gsl_vector_set(phiT,row,20);//rest is set to 20
    }
    if(Case==DIRICHLET){gsl_vector_set(phiT,(N+1),0);}//sets boundary condition of other end to 0

    gsl_vector * phiTdt = gsl_vector_calloc(N+2);//vector of forward time phi values
    gsl_matrix * A = gsl_matrix_calloc(N+2,N+2);//matrix initialised to zero.

    /*** Here we fill coefficient matrix *************/
    double a = -(dt*alpha)/pow(h,2);
    double b = (1+(2*dt*alpha)/pow(h,2));


    gsl_matrix_set(A,0,0,1);//first boundary node
    for(row=1; row<(N+1);row++)
    {
        gsl_matrix_set(A,row,(row-1),a);
        gsl_matrix_set(A,row,row,b);
        gsl_matrix_set(A,row,(row+1),a);
    }
    gsl_matrix_set(A,N+1,N+1,1);//2nd boundary node

    /*** Here we do LUD on the coefficient matrix ***************/
    gsl_permutation * p = gsl_permutation_alloc (N+2);//for LUD
    int s;//for LUD
    gsl_linalg_LU_decomp (A, p, &s);

    /*** iterating forwards in time *********/
    while( (t/dt)<=1e9)
    {
        gsl_linalg_LU_solve (A, p, phiT, phiTdt);

        if(Case == VN)
        { //sets phiTdt_N-1 = phiTdt_N-2
            //Von Neumann boundary conditions so we want dphi/dx = 0 at righthand edge
            //using backwards difference approximation since we are on the end of the rod
            //gives that last element should equal the one to its left with an O(h) error
            gsl_vector_set(phiTdt,(N) ,gsl_vector_get(phiTdt,(N-1)));
        }

        for(row=1; row<=N; row++)
        {//sets phiT equal to phiTdt
            gsl_vector_set(phiT,row,gsl_vector_get(phiTdt,row));
        }
        t+=dt;
        //fprints out every power of 10*dt
        if( fabs(t-1*dt)<dt || fabs(t-10*dt)<(dt/2) || fabs(t-1e2*dt)<(dt/2) || fabs(t-1e3*dt)<(dt/2)
        || fabs(t-1e4*dt)<(dt/2) || fabs(t-1e5*dt)<(dt/2) || fabs(t-1e6*dt)<(dt/2) || fabs(t-1e7*dt)<(dt/2)
           || fabs(t-1e8*dt)<(dt/2) || fabs(t-1e9*dt)<(dt/2) )
        {
            FILE * FP;
            char fileName[256];
            sprintf(fileName,
            "t=%.0lf.txt",t);
                FP = fopen(fileName,"w");
            if ((FP == NULL))
            {
                perror ("Error opening out put file:");//error message if file doesnt open
            }
            else
            {
                printf("t=%lf\n",t);
                for(row = 1; row <= N; row++)
                {
                    float x = (row-1)*h;
                    printf("%lf\n",gsl_vector_get(phiT,row));
                    fprintf(FP,"%lf\t%lf\n",x,gsl_vector_get(phiT,row));
                }
            }
            printf("\n");
            fclose(FP);
        }
    }

    gsl_permutation_free (p);
    gsl_vector_free (phiT);
    gsl_vector_free (phiTdt);
    gsl_matrix_free (A);

    return 0;
}

//fprintVector
