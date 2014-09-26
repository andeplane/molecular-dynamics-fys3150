#include <math/random.h>
#include <cmath>

Random::Random(long seed_) :
    iy(0),
    seed(seed_)
{

}

double Random::nextGauss(double mean, double standardDeviation) {
    double randomNumberFromStandardNormalDistribution = sqrt( -2.0*log(1.0 - nextDouble()) )* cos( 6.283185307 * nextDouble() );

    return standardDeviation*randomNumberFromStandardNormalDistribution + mean;
}

double Random::nextDouble()
{
   int             j;
   long            k;
   double          temp;

   if (seed <= 0 || !iy) {
      if (-(seed) < 1) seed=1;
      else seed = -(seed);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (seed)/IQ;
         seed = IA*(seed - k*IQ) - IR*k;
         if(seed < 0) seed += IM;
         if(j < NTAB) iv[j] = seed;
      }
      iy = iv[0];
   }
   k     = (seed)/IQ;
   seed = IA*(seed - k*IQ) - IR*k;
   if(seed < 0) seed += IM;
   j     = iy/NDIV;
   iy    = iv[j];
   iv[j] = seed;
   if((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}
