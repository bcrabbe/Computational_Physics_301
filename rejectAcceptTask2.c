#include <stdio.h>
#include <stdlib.h>
#include <gsl_rng.h>
#include <time.h>
#include <math.h>
#include <gsl_histogram2d.h>
#include <gsl_randist.h>

#define FALSE 0
#define TRUE 1

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
    size_t i=0;
    printf("generating data....\n");
    while(i<numberOfEvents)
    {
         //first find the z coord of decay:
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
            gammaVelocityDirection[2] = gsl_rng_uniform(r)*2 - 1;   //produces a uniform deviate between -1 and 1

            if(gammaVelocityDirection[2]<0) //emitted gamma moves in opposite direction so wont hit the detector
            {                               //so dont need to continue
                ++i;
            }
            else                            //if the z direction of gamma is +ve we can extrapolate to detector
            {
                int rejected=FALSE;
rejectionSample://program returns to this point if velocity is rejected, if that happens we must resample z veolcity
                //however we required it to be +ve to get to here so now can only be 0>zVel>=1
                if(rejected==TRUE){gammaVelocityDirection[2] = gsl_rng_uniform(r);}

                //now we generate random x and y velocity directions aswell
                //this gives us points randomly generated in a unit cube
                gammaVelocityDirection[0] = gsl_rng_uniform(r)*2 - 1;
                gammaVelocityDirection[1] = gsl_rng_uniform(r)*2 - 1;

                //find its radial distance:
                double gammaVelocityDirectionMag = VectorMag(gammaVelocityDirection);

                if(gammaVelocityDirectionMag>1)
                {                                  //if the point lies outside the unit sphere then it is rejected.
                    rejected=TRUE;
                    goto rejectionSample;           //and we try again
                }
                else
                {//otherwise it is accepted and then normalised to become a point on the surface of the unit sphere

                    gammaVelocityDirection[0] = gammaVelocityDirection[0]/gammaVelocityDirectionMag;
                    gammaVelocityDirection[1] = gammaVelocityDirection[1]/gammaVelocityDirectionMag;
                    gammaVelocityDirection[2] = gammaVelocityDirection[2]/gammaVelocityDirectionMag;
                    //works out point gamma hits detector:
                    double timeTravelledToDetector = (zDetector-decayVertex[2])/gammaVelocityDirection[2];
                    detectorPositionTruth[0] = timeTravelledToDetector*gammaVelocityDirection[0];
                    detectorPositionTruth[1] = timeTravelledToDetector*gammaVelocityDirection[1];
                    //applys detector resolution smearing:
                    detectorPositionMeasured[0] = detectorPositionTruth[0] + gsl_ran_gaussian_ziggurat( r, xDetectorSmear);
                    detectorPositionMeasured[1] = detectorPositionTruth[1] + gsl_ran_gaussian_ziggurat( r, yDetectorSmear);
                    //bin data:
                    gsl_histogram2d_increment(h, detectorPositionMeasured[0], detectorPositionMeasured[1]);
                    ++i;
                }

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
    size_t j;
    for(i=0; i<xNumberOfBins; ++i)
    {
        double binLowerLimit, binUpperLimit;
        gsl_histogram2d_get_xrange(h,i,&binLowerLimit,&binUpperLimit);
        double x = (binLowerLimit+binUpperLimit)/2;
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
