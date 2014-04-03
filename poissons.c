#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <errno.h>

# define DISPLAY 1
# define PLOT 2
# define VECTOR 3

# define POINT 1
# define SQUARE 2
# define CAPACITOR 3

double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void boundaryInit(double **grid, long numberOfNodesX1,long numberOfNodesY1);
void gridInit(double **grid, long numberOfNodesX,long numberOfNodesY);
void iterate(double **grid, long numberOfNodesX,long numberOfNodesY);
void printGrid(double **grid, long numberOfNodesX1,long numberOfNodesY1);
void fprintGrid(double **grid, long numberOfNodesX1,long numberOfNodesY1, int type);
void createCapacitor(double **grid, long numberOfNodesX, long numberOfNodesY,int a, int d);

double maxChange=10;
double hX,hY;//node spacing in X and Y directions (m)
double X=1, Y=1; //the dimensions of the area we wish to model (m)
double Q=10;//potential of the configurations
int source = CAPACITOR; // choose POINT, CAPACITOR or SQUARE.
double convergence=1e-6;
int main()
{
    /*** Creates Grid **************************************/
    long numberOfNodesX=101, numberOfNodesX1=numberOfNodesX+1;//adjust to change spacing
    long numberOfNodesY=101, numberOfNodesY1=numberOfNodesY+1;
    hX = X/(numberOfNodesX-1);//node spacing X diiection (m)
    hY = Y/(numberOfNodesY-1);//nodes spacing Y (m)
    printf("grid spacings:\n x:%lf\ny:%lf\n",hX,hY);
    double **grid;
    grid = dmatrix(0,numberOfNodesX1,0,numberOfNodesY1);
    //forms a grid with the required number of nodes: [1->numberOfNodesY][1->numberOfNodesX]
    //aswell as an additional layer on each side to store boundary conditions.

    /***Initializes Grid ************************************/
    boundaryInit(grid,numberOfNodesX1,numberOfNodesY1);
    gridInit(grid,numberOfNodesX,numberOfNodesY);
    a=10, d=10; //length and separation of capacitor plates in number of nodes
    if( (source==CAPACITOR)  && (d>numberOfNodesY || a>numberOfNodesX) )
    {
        printf("ERROR: Capacitor too large.\n");
        return -1;
    }

    createCapacitor(grid,numberOfNodesX,numberOfNodesY,a,d);
    //fprintGrid(grid,numberOfNodesX1,numberOfNodesY1, DISPLAY);//prints initial grid to display file

    /*** Iterates to solution*****************************/
    int numOfIterations=0;
    while(maxChange>convergence)
    {
        if(numOfIterations>0){createCapacitor(grid,numberOfNodesX,numberOfNodesY,a,d);}
        iterate(grid,numberOfNodesX,numberOfNodesY);

        printf("maxChange=%lf\n",maxChange);
        ++numOfIterations;

    }
    createCapacitor(grid,numberOfNodesX,numberOfNodesY,a,d);

    /*** Prints final grid *******************************/
    //fprintGrid(grid,numberOfNodesX1,numberOfNodesY1,DISPLAY);
    fprintGrid(grid,numberOfNodesX1,numberOfNodesY1,PLOT);
    printf("\nnumber of iterations required: %d\n",numOfIterations);

    /***  E field *****************************/
    if(source == CAPACITOR){fprintGrid(grid,numberOfNodesX1, numberOfNodesY1, VECTOR);}
    free_matrix(grid,0,numberOfNodesX1, 0, numberOfNodesY1);
    return 0;
}






void iterate(double **grid, long numberOfNodesX,long numberOfNodesY)
{//iterates using Gauss-Siedel method.
    int i,j,maxChangeNodeI, maxChangeNodeJ;
    double change,newValue;
    maxChange=0;
    for(i = 1; i <= numberOfNodesY; i++)
    {
        for(j = 1; j <= numberOfNodesX; j++)
        {
            newValue = (grid[i-1][j] + grid[i+1][j] + grid[i][j-1] + grid[i][j+1])/4;
            change = fabs(newValue - grid[i][j] );
            if( change>maxChange && fabs(fabs(grid[i][j])-Q)>1e-6 )//checks that i j is not part of the conductor
            {
                maxChange=change;
                maxChangeNodeI=i;
                maxChangeNodeJ=j;
            }
            grid[i][j] = newValue;
        }
    }
    printf("maxChange in grid[%d][%d]",maxChangeNodeI,maxChangeNodeJ);
}

void printGrid(double **grid, long numberOfNodesX1,long numberOfNodesY1)
{
    int i,j;
    for(i = 0; i <= numberOfNodesY1; i++)
    {
        for(j = 0; j <= numberOfNodesX1; j++)
        {
            printf("%lf\t",grid[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void fprintGrid(double **grid, long numberOfNodesX1,long numberOfNodesY1, int type)
{
    static int numberOfCalls=0;
    ++numberOfCalls;
    char fileName[256];
    if(type == PLOT)//printing for gnuplot:
    {
        FILE *FP;
        sprintf(fileName,
            "Potential_d=%d.txt",d);
        FP = fopen(fileName,"w");
        if ((FP == NULL))
        {
            perror ("Error opening out put file:");//error message if file doesnt open
        }

        else
        {
            int i,j;
            for(i = 0; i <= numberOfNodesY1; i++)
            {
                for(j = 0; j <= numberOfNodesX1; j++)
                {
                    double y = (i-1)*hY;
                    double x = (j-1)*hX;//the x,y coordinates (m) of node i,j
                    fprintf(FP,"%lf\t%lf\t%lf\n",x,y,grid[i][j]);// prints x,y,z
                }
                fprintf(FP,"\n");
            }
            fprintf(FP,"\n");
        }
        fclose(FP);
    }
    else if(type == DISPLAY)
    {//and prints in matrix form
        FILE *FP2;
        sprintf(fileName,
            "point_matrix_h-%lf.txt",hX);
        FP2 = fopen(fileName,"w");
        if ((FP2 == NULL))
        {
            perror ("Error opening out put file:");
        }
        else
        {
            int i,j;
            for(i = 0; i <= numberOfNodesY1; i++)
            {
                for(j = 0; j <= numberOfNodesX1; j++)
                {



                    fprintf(FP2,"%lf\t",grid[i][j]);
                }
                fprintf(FP2,"\n");
            }
            fprintf(FP2,"\n");
        }
    fclose(FP2);
    }
    else if(type == VECTOR)
    {
        FILE *FP3;
        FILE *FP4;
        sprintf(fileName,
            "Evector_d=%d.txt",d);
        FP3 = fopen(fileName,"w");
        if ((FP3 == NULL))
        {
            perror ("Error opening out put file:");
        }
        sprintf(fileName,
            "Emag_ad=%f.txt",(float)a/(float)d);
        FP4 = fopen(fileName,"w");
        if ((FP4 == NULL))
        {
            perror ("Error opening out put file:");
        }
        else
        {   int i,j,numberOfNodesX = numberOfNodesX1-1,numberOfNodesY=numberOfNodesY1-1;
            for(i = 1; i <= numberOfNodesY ; ++i/*+=1/hY*/)
            {
                for(j = 1; j <= numberOfNodesX; ++j/*+=1/hX*/)
                {//here we use the symetric 3 point formula.
                    //this minus| is specific to EM since: E = -grad(V)
                    double Ey = (grid[i-1][j]-grid[i+1][j])/(2*hY);//vertical E component
                    double Ex = -(grid[i][j+1]-grid[i][j-1])/(2*hX);//horizontal E component
                    double y = (i-1)*hY;
                    double x = (j-1)*hX;//the x,y coordinates (m) of node i,j
                    double Emag = sqrt(pow(Ex,2) +pow(Ey,2));
                    if( fabs((fabs(grid[i][j]) -Q)) < 0.1 )
                    {
                        //this stops it prinint a huge arrow at the centre.
                        Ex = 0;
                        Ey = 0;
                    }
                    if( (i%(numberOfNodesY/10)==0) && (j%(numberOfNodesX/10)==0) )
                    {
                        fprintf(FP3,"%lf\t%lf\t0\t%lf\t%lf\t0\n",x,y,(Ex/(12*70)),(Ey/(12*70)));
                    }
                    fprintf(FP4,"%lf\t%lf\t%lf\n",x,y,Emag);
                }
                 fprintf(FP4,"\n");
            }
        }
    }
    else{printf("ERROR: fprintf type undefined\n");}
}

void createCapacitor(double **grid, long numberOfNodesX, long numberOfNodesY, int a, int d)
{

    /**** forms a Capacitor *********************************************/
    if(source == CAPACITOR)
    {//attempts to make it central
        float topPlateRow = ((float)numberOfNodesY/2) - ((float)d/2) ;
        float bottomPlateRow = 1 + ((float)numberOfNodesY/2) + ((float)d/2) ;
        float plateColMin = 1 + ((float)numberOfNodesX/2) - ((float)a/2) ;
        float plateColMax = 1 + ((float)numberOfNodesX/2) + ((float)a/2) ;
        int i;
        for(i=(int)plateColMin; i<(int)plateColMax; i++ )
        {
            grid[(int)topPlateRow][i]= Q;
            grid[(int)bottomPlateRow][i]= -Q;
        }

    }


    /*** forms square conductors of various sizes ****************/
    // Grid spacing must be selected to one of these values
    if(source == SQUARE)
    {

        if(numberOfNodesX==501)
        {
            int i, j;

            for(i=200; i<=300; i++)
            {
                for(j=200; j<=300; j++)
                {
                    grid[i][j] = Q;
                }
            }
        }
        else{printf("ERROR: make n comply with square source\n");}
    }
    else{printf("ERROR: make number of nodes a value specified in function createCapacitor \n");}
    }

    /**** Forms a Point source***************************/
    //note: number of nodes must be odd for there be a central point.
    if(source == POINT)
    {
        int centreX = (numberOfNodesX+1)/2;
        int centreY = (numberOfNodesY+1)/2;
        grid[centreY][centreX] = Q;
    }


}

void gridInit(double **grid, long numberOfNodesX,long numberOfNodesY)
{
    int i,j;
    for(i = 1; i <= numberOfNodesY; i++)
    {
        for(j = 1; j <= numberOfNodesX; j++)
        {
           // srand ( time(NULL) );
            double gridInitValue =0;// rand()%10;
            grid[i][j] = gridInitValue;
        }
    }
}

void boundaryInit(double **grid, long numberOfNodesX1,long numberOfNodesY1)
{
    int i;
    //top and bottom edges
    for(i = 0; i <= numberOfNodesX1; i++)
    {
        grid[0][i] = 0;
        grid[numberOfNodesY1][i] = 0;
    }
    //left and right edges
    for(i = 0 ;i <= numberOfNodesY1; i++)
    {
        grid[i][0] = 0;
        grid[i][numberOfNodesX1] = 0;
    }

}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
	if (!m) printf("allocation failure 1 in matrix()");
	m += 1;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+1)*sizeof(double)));
	if (!m[nrl]) printf("allocation failure 2 in matrix()");
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((char*) (m[nrl]+ncl-1));
	free((char*) (m+nrl-1));
}
