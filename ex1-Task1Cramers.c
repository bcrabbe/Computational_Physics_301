#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>



float **matrix(long nrl, long nrh, long ncl, long nch); //allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
float **minormatrix(int i, int j, int n, float **m);//will return the minor matrix corresponding to element ij of the matrix m, also need to supply the size (n) of the matrix
float det(float **a,int n); //calculates a determinant of the matrix pointed to by a, using lapace's method.
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);

int main()
{
    FILE *FP;
    FP =fopen("cramers.txt","w");
    int n=3;//size of matrix to be inverted

    //clock_t t_before = time(NULL);//time before calculation
    //float dt = 0;
    //double k, step=1e-14;//for singularity test
    //for(k=0.000000011;k>=0;k-=step)
   // {


/*** Creates matrix ***/

    float **a;
    a = matrix(1,n,1,n); //creates a n x n matrix referenced by a[1..n][1...n]


    int i,j;
/*** Initialises randomly ***/

    srand ( time(NULL) );// this will change the random numbers each time
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=n;j++)
        {
            float aij = rand()%500;//will generate a random integer between 0 and 500
            a[i][j] = aij;
        }
    }


/*** Prints the Matrix to screen  ***/
    printf("\nA=\n");
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
        {
            printf("%f\t",a[i][j]);
        }
    printf("\n");
    }

/*** Calculates the determinant ***/
     float determinant;
    determinant = det(a,n);


/*** Makes the inverse matrix ***/
    float **inva, **Mij, Cij, detMij;
    int sign ,n1=(n-1);
    inva = matrix(1,n,1,n);
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
        {
            Mij = minormatrix(i,j,n,a);
            detMij = det(Mij,n1);
            sign = pow(-1,(i+j));// + - + - \n - + - +
            Cij = sign*detMij;
            inva[j][i] = (Cij/determinant); //note idices swapped since its transpose
            free_matrix(Mij,1,n1,1,n1);

        }
    }

/*** Prints inverse matrix ***/
    printf("\nThe inverse:\n");
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
        {
            printf("%f\t",inva[i][j]);
        }
    printf("\n");
    }

    /*** calculates the ave. magnitude of elements of inv A for singularity test ***/
/*    double tot=0;
    for( i=1;i<=n;i++)
    {
        for( j=1; j<=n;j++)
        {
            tot+=fabs(inva[i][j]);
        }
    }
    double ave = tot/(n*n);
    fprintf(FP,"%lf\t%lf\n",k,ave);
*/
    free_matrix(a,1,n,1,n);
    free_matrix(inva,1,n,1,n);

    /*** time test ***/
    //clock_t t_after = time(NULL);
    //double dt = t_after - t_before;
    //fprintf(FP,"%d\t%lf\n",n,dt);
    //printf("%d\t%lf\n",n,dt);



    fclose(FP);
    return 0;
}



/*pass it the row column numbers i j and the size of the matrix and the matrix itself and it will return
the minor matrix corresponding to this element */
float **minormatrix(int i, int j, int n, float **m)
{
    float **minor;
    int u=1,r=1;
    int u1, r1, n1;
    n1 = n -1;
    minor = matrix(1, n1, 1, n1);
    for(u=1; u<=n; u++)
    {
        if(u == i)
        {
            //Do nothing
        }
        else
        {
            for(r=1; r<=n; r++)
            {
                u1 = u-1;
                r1= r-1;
                if(r==j)
                {
                   //Do nothing
                }
                if( (u>i) && (r>j) )
                {
                    minor[u1][r1] = m[u][r];
                }
                if( (u>i) && (r<j) )
                {
                    minor[u1][r] = m[u][r];
                }
                if( (u<i) && (r>j) )
                {
                    minor[u][r1] = m[u][r];
                }
                if( (u<i) && (r<j) )
                {
                    minor[u][r] = m[u][r];
                }
            }
        }
    }
    return minor;
}



float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+1)*sizeof(float*)));
	m += 1;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+1)*sizeof(float)));
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++){
	m[i]=m[i-1]+ncol;
	}

	/* return pointer to array of pointers to rows */
	return m;
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((char*) (m[nrl]+ncl-1));
	free((char*) (m+nrl-1));
}


float det(float **m, int n)
{
    float determinant=0;
    if(n==2)
    {
        determinant = (m[1][1]*m[2][2]) - (m[2][1]*m[1][2]);
        if(fabs(determinant)<(1e-6))
        {
            printf("\n!!!!!!!!!!!!!!!\n det = %f -->Rounding errors \n!!!!!!!!!!!!!!!\n",determinant);
        }
        return determinant;
    }
    else
    {
        float Cij,detMij;//cofactor terms
        float **Mij;//minormatrix ij
        int i=n,j;//will always expand by row nth so i=n -- j(column no.) will move across the row
        int n1=n-1;
        int sign;
        for(j=1; j<=n; j++)
        {
            Mij = minormatrix(i,j,n,m);//returns the minor matrix for element ij of the matrix m (size n)- the produced matrix will be size n-1
            detMij = det(Mij,n1);
            sign = pow(-1,(i+j));// + - + - + - +
            Cij = sign*detMij;
            determinant += m[i][j]*Cij;//det m = m11*C11 + m12*C12 + ... + m1n*C1n
            free_matrix(Mij,1,n1,1,n1);
        }
        return determinant;
    }
}


