#include <stdio.h>
#include <stdlib.h>
#include <gsl_rng.h>
#include <math.h>
#include <gsl_histogram.h>
#include <time.h>

int main()
{
    ///initialise random number generator:
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r, 5489);//sets seed

    size_t numberOfDataPoints = 1e5;

    ///gsl implementation for histogramming data:
    size_t numberOfBins = cbrt((double)numberOfDataPoints);
    gsl_histogram * h = gsl_histogram_alloc(numberOfBins);
    gsl_histogram_set_ranges_uniform(h, 0.00, M_PI);

    printf("generating data....\n");
    size_t i;
    for(i=0; i<=numberOfDataPoints; ++i)
    {
        double rand01 = gsl_rng_uniform(r);             //generated random number between 0 and 1
        double randSin0Pi = acos(1-2*rand01);           //transformed by quantile function to sin(x) distribution
        gsl_histogram_increment(h, randSin0Pi);         //bins data
    }
    printf("done.\n");

    FILE * FP2;
    FP2 = fopen("histogram.txt","w");

    double mean = gsl_histogram_mean(h);
    //1st moment is mean:
    double moments[6] = {mean,0,0,0,0,0};
    double sinMoments[6] = {M_PI_2,0,0,0,0,0};
    double Y=0, sinY=0;
    for(i=0; i<numberOfBins; ++i)
    {
        ///gets x and y coord and fprints data
        double y = 2*gsl_histogram_get (h, i)/gsl_histogram_sum(h);
        double binLowerLimit,  binUpperLimit;
        gsl_histogram_get_range(h, i, &binLowerLimit, &binUpperLimit);
        double x = (binLowerLimit + binUpperLimit)/2;
        fprintf(FP2,"%lf %lf %lf\n",x,y,sin(x));

        ///moments calculations
        Y += y;
        sinY+=sin(x);
        int s;
        double delta =(x-mean),sinDelta=(x-M_PI_2);
        for(s=1; s<6; ++s)
        {
            moments[s] += (pow(delta,(s+1)) - moments[s])*(y/Y);
            sinMoments[s] += (pow(sinDelta,(s+1)) - sinMoments[s])*(sin(x)/sinY);
        }
    }
    ///prints moments
    printf("moments:");
    printf("P(x')\t sin(x)\n");
    int s;
    for(s=0; s<6; ++s)
    {
       printf("%lf\t%lf\n",moments[s],sinMoments[s]);
    }


    gsl_histogram_free (h);
    gsl_rng_free (r);

    return 0;
}
