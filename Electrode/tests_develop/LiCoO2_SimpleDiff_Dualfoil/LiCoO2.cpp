

#include <cstdio>
#include <cmath>

 double calcV( double relExtent)
{
   double volts;
   double xLi = relExtent;


        /*
         *   Line 4580 of dualfoil 5.1
	 * 
	 *	  Measured by Oscar Garcia 2001 using Quallion electrodes for
	 *  c     0.5 < y < 0.99.  Fit revised by Karen Thomas in May 2003 to
         *  c     match Doyle's fit for y < 0.4 and Garcia's data at larger y.
	 *  c     Valid for 0 < y < 0.99. Note that capacity fade is found to
	 *  c     occur experimentally if y goes below 0.5; this is not included
	 *  c     in the model
	 */

   volts = ( 2.16216 + 0.07645 * tanh(30.834 - 54.4806 * xLi)
		  + 2.1581  * tanh( 52.294  - 50.294  * xLi)
		  - 0.14169 * tanh( 11.0923 - 19.8543 * xLi)
		  + 0.2051 *  tanh( 1.4684  - 5.4888  * xLi)
		  + 0.2531 *  tanh( (-xLi + 0.56478) / 0.1316)
		  - 0.02167 * tanh( ( xLi - 0.525)   / 0.006)      );


return volts;
} 

int main()
{
   int numP = 51;
 
   printf("DUALFOIL Positive electrode OCV data \n\n");
   
   printf("   relExtent           xLi    xV       OCV \n");
   for (int i = 0; i < numP; i++) {
      double relE = 0.0 + (double (i) / (numP - 1.0));
      double volts = calcV(relE);
      printf("%15.5E %15.5E %15.5E %20.13f \n", relE, relE, 1.0 - relE, volts);
    }
    return 0;
}
