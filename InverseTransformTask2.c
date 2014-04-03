#include <stdio.h>
#include <stdlib.h>
#include <gsl_rng.h>
#include <time.h>
#include <math.h>
#include <gsl_histogram2d.h>
#include <gsl_randist.h>

double VectorMag(double *v);

int main()
{

    double v=2e-3;            //meters per microsecond
    double meanLifetime=520;//microseconds
    double zDetector=2;     //meters
    double xDetectorSmear=0.1, yDetectorSmear=0.3;//meters
    double xDetectorSize=10,yDetectorSize=10;//meters
    size_t numberOfEvents=1e6;
    ///initialise random number generator:
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r, 5489);       //sets seed to recomended value

    ///gsl implementation for histogramming data:
    size_t xNumberOfBins=200,yNumberOfBins=200;
    gsl_histogram2d * h = gsl_histogram2d_alloc(xNumberOfBins,yNumberOfBins);
    gsl_histogram2d_set_ranges_uniform (h, -xDetectorSize/2, xDetectorSize/2, -yDetectorSize/2, yDetectorSize/2);

    double gammaVelocityDirection[3], decayVertex[3], detectorPositionTruth[2],detectorPositionMeasured[2];
    size_t i;
    printf("generating data....\n");
    for(i=0; i<numberOfEvents; ++i)
    {   //first find the z coord of decay:
        decayVertex[2] = (v*gsl_ran_poisson(r, meanLifetime));

         //we require that it decayed in the decay volume, otherwise it misses the detector
        //the chance of this happening is vanishingly small
        if(decayVertex[2]<0 || decayVertex[2]>zDetector )
        {
            ++i;
        }
        else
        {
            //then the z direction velocity vector of the produced gamma:
            double z = gsl_rng_uniform(r)*2 - 1;   //produces a uniform deviate between -1 and 1

            if(z>0)       //if the z direction of gamma is +ve we can extrapolate to detector
            {                                     //if not then the event misses detector so we dont continue
                //now we generate random x and y velocity directions aswell
                double phi = gsl_rng_uniform(r)*M_PI*2;
                double theta = asin(z);

                gammaVelocityDirection[0]=cos(theta)*cos(phi);
                gammaVelocityDirection[1]=cos(theta)*sin(phi);
                gammaVelocityDirection[2]=z;
                //works out point gamma hits detector:
                double timeTravelledToDetector = (zDetector-decayVertex[2])/gammaVelocityDirection[2];
                detectorPositionTruth[0] = timeTravelledToDetector*gammaVelocityDirection[0];
                detectorPositionTruth[1] = timeTravelledToDetector*gammaVelocityDirection[1];
                //applys detector resolution smearing:
                detectorPositionMeasured[0] = detectorPositionTruth[0] + gsl_ran_gaussian_ziggurat( r, xDetectorSmear);
                detectorPositionMeasured[1] = detectorPositionTruth[1] + gsl_ran_gaussian_ziggurat( r, yDetectorSmear);
                //bin data:
                gsl_histogram2d_increment(h, detectorPositionMeasured[0], detectorPositionMeasured[1]);
            }
        }
    }
    printf("number of events: %zu\nnumber of data: %lf\n",i,gsl_histogram2d_sum (h));
    printf("std deviation: x = %lf, y = %lf, average = %lf\n",gsl_histogram2d_xsigma(h),gsl_histogram2d_ysigma(h),
                (gsl_histogram2d_xsigma(h)+gsl_histogram2d_ysigma(h))/2 );
    printf("FWHM: x = %lf, y = %lf, average = %lf\n",2.355*gsl_histogram2d_xsigma(h), 2.355*gsl_histogram2d_ysigma(h),
                (2.355*gsl_histogram2d_xsigma(h)+2.355*gsl_histogram2d_ysigma(h))/2 );
    ///outputs histogrammed data:
    FILE * FP1 = fopen("detectorData.txt","w");
    FILE * FP2 = fopen("2dDetectorData.txt","w");
    size_t j;
    for(i=0; i<xNumberOfBins; ++i)
    {
        double binLowerLimit, binUpperLimit;
        gsl_histogram2d_get_xrange(h,i,&binLowerLimit,&binUpperLimit);
        double x = (binLowerLimit+binUpperLimit)/2;
        fprintf(FP2,"%lf %lf\n",x,gsl_histogram2d_get(h, i, yNumberOfBins/2));
        for(j=0; j<yNumberOfBins; j++)
        {
            gsl_histogram2d_get_yrange(h,j,&binLowerLimit,&binUpperLimit);
            double y = (binLowerLimit+binUpperLimit)/2;
            fprintf(FP1,"%lf %lf %f\n",x,y,gsl_histogram2d_get(h, i, j));
        }
        fprintf(FP1,"\n");
    }

    gsl_histogram2d_free (h);
    return 0;
}



double VectorMag(double *v)
{
    return sqrt( pow(v[0],2) + pow(v[1],2) + pow(v[2],2) );
}
