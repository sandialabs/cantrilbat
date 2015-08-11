

#include <cstdio>
#include <cmath>

 double calcV( double relExtent)
{
   double volts;
   double xLi = 1.0 - relExtent;
/*
 *   Line 4174 of dualfoil 5.1

   MCMB 2510 carbon (Bellcore)
 c      c1=-0.160d0
 c      c2=1.32d0
 c      c3=-3.0d0
 c      g0=c1+c2*expg(c3*xx(2+mpa,j)/ct1)
 c      g0=g0+10.d0*expg(-2000.d0*xx(2+mpa,j)/ct1)
 c      g1=c2*c3*expg(c3*xx(2+mpa,j)/ct1)/ct1
 c      g1=g1-10.d0*2000.d0/ct1*expg(-2000.d0*xx(2+mpa,j)/ct1)
 c     MCMB 2528 graphite measured by Chris Bogatu 2000, Telcordia and PolyStor materials
 c     for 0.01 < x < 0.9
 *
 *   (note this agrees exactly with dualfoil)
 */
  
   volts = ( 0.194 + 1.5 * exp(-120.0 * xLi)
	     + 0.0351 * tanh( (xLi - 0.286)   / 0.083)
	     - 0.0045 * tanh( (xLi - 0.849)   /0.119)
	     - 0.035 *  tanh( (xLi - 0.9233)  /0.05)
	     - 0.0147 * tanh( (xLi - 0.5)     /0.034)
	     - 0.102 *  tanh( (xLi - 0.194)   /0.142)
	     - 0.022 *  tanh( (xLi - 0.9)     /0.0164)
	     - 0.011 *  tanh( (xLi - 0.124)   /0.0226)
	     + 0.0155 * tanh( (xLi - 0.105)   /0.029)) ;



return volts;
} 

int main()
{
   int numP = 51;
 
   printf("DUALFOIL Negative electrode OCV data \n\n");
   
   printf("   relExtent           xLi          OCV \n");
   for (int i = 0; i < numP; i++) {
      double relE = 0.0 + (double (i) / (numP - 1.0));
      double volts = calcV(relE);
      printf("%15.5E %15.5E %20.13f \n", relE, 1.0-relE, volts);
    }
    return 0;
}
