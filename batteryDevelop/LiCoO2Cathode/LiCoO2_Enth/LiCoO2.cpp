

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

double calcV_mao( double relExtent)
{
   double volts;
   double xLi = relExtent;


        /*
         *   Line 4580 of dualfoil 5.1
         * 
         *        Measured by Oscar Garcia 2001 using Quallion electrodes for
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

    return   volts;
}

double calcV_dUdT_mao( double relExtent)
{
   double xLi = relExtent;
   double xLi2 = xLi * xLi;
   double xLi3 = xLi2 * xLi;
   double xLi4 = xLi3 * xLi;


    double A = -0.19952 + 0.92837 * xLi - 1.36455 * xLi2 + 0.61154 * xLi3;
    double B = 1.0 - 5.66148 * xLi + 11.47636 * xLi2 - 9.82431 * xLi3 + 3.04876 * xLi4;

    double dUdT  = A / B;
    return 1.0E-3 * dUdT;
}

int main()
{
   int numP = 51;
 
   printf("DUALFOIL Mao et al (2014) (they are the same) Positive electrode OCV data \n\n");
   
   double dd = calcV_dUdT_mao(1.0);

   printf("  relExtent(xLi)       xV           OCV           dVdT         DeltaS         DeltaH       DeltaG\n");
   for (int i = 0; i < numP; i++) {
      double relE = 0.0 + (double (i) / (numP - 1.0));
      double voltsM = calcV_mao(relE);

      double dvoltsdT = calcV_dUdT_mao(relE);

      double F = 96485.3E3;   // Joules/kmol
      double DeltaS = - F * dvoltsdT;
      double DeltaG = voltsM * F;
      double T = 300.;
      double DeltaH = DeltaG + T * DeltaS;

    
      printf("%15.5E %15.5E %10.5f      % 10.5E  % 10.5E   % 10.5E  % 10.5E\n",
             relE, 1.0-relE, voltsM,  dvoltsdT, DeltaS , DeltaH, DeltaG);


    }
    return 0;
}
