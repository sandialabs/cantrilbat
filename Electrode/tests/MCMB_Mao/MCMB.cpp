/*
 *  Parameterizations of the OCV for MCMB contained within the paper.
 */

#include <cstdio>
#include <cmath>

double calcV_dold( double relExtent)
{
   double volts;
   double xLi = 1.0 - relExtent;
   volts = (  0.124   + 1.5 * exp(-150.0 * xLi) 
                    + 0.0351 * tanh( (xLi - 0.286) / 0.083)
                    - 0.0045 * tanh( (xLi - 0.90)  / 0.119)
                    - 0.035  * tanh( (xLi - 0.99)  / 0.05 )
                    - 0.0147 * tanh( (xLi - 0.50)  / 0.034)
                    - 0.102  * tanh( (xLi - 0.194) / 0.142)
                    - 0.022  * tanh( (xLi - 0.98 ) / 0.0164)
                    - 0.011  * tanh( (xLi - 0.124) / 0.0226)
                    + 0.0155 * tanh( (xLi - 0.105) / 0.029));

    return volts;
} 
//
//  Mao et al. expression for the open circuit voltage for the MCMB electrode
//
double calcV_mao( double relExtent)
{
   double volts;
   double xLi = 1.0 - relExtent;


   volts = (  0.194   + 1.5 * exp(-120.0 * xLi) 
                    + 0.0351 * tanh( (xLi - 0.286) / 0.083)
                    - 0.0045 * tanh( (xLi - 0.849)  / 0.119)
                    - 0.035  * tanh( (xLi - 0.9233)  / 0.05 )
                    - 0.0147 * tanh( (xLi - 0.50)  / 0.034)
                    - 0.102  * tanh( (xLi - 0.194) / 0.142)
                    - 0.022  * tanh( (xLi - 0.90 ) / 0.0164)
                    - 0.011  * tanh( (xLi - 0.124) / 0.0226)
                    + 0.0155 * tanh( (xLi - 0.105) / 0.029));

    return volts;
} 

double calcdVdT_mao( double relExtent)
{
   double DvoltsDT;
   double xLi = 1.0 - relExtent;
   double xLi2 = xLi * xLi;
   double xLi3 = xLi * xLi * xLi;
   double xLi4 = xLi3 * xLi;
   double xLi5 = xLi4 * xLi;
   double xLi6 = xLi5 * xLi;
   double xLi7 = xLi6 * xLi;
   double xLi8 = xLi7 * xLi;

   double A = (0.00527 + 3.29927 * xLi - 91.79326 * xLi2 + 1004.91101 * xLi3 - 5812.27813 * xLi4 + 19329.75490 * xLi5
               - 37147.89470 * xLi6 + 38379.18127 * xLi7 - 16515.05308 * xLi8);

   double B = (1.0 - 48.09287 * xLi + 1017.23480 * xLi2 - 10481.80419 * xLi3 + 59431.30001 * xLi4
              - 195881.64880 * xLi5 + 374577.31520 * xLi6 - 385821.16070 * xLi7 + 165705.85970 * xLi8);

   DvoltsDT = A / B * 1.0E-3;

   return DvoltsDT;
}


int main()
{
   int numP = 51;
    
   printf("   relExtent           xLi          OCV(Dold)   OCV(Mao et al.)    dVdT       DeltaS       DeltaH       DeltaG\n");
   for (int i = 0; i < numP; i++) {
      double relE = 0.0 + (double (i) / (numP - 1.0));
      double volts_A = calcV_dold(relE);
      double volts_B = calcV_mao(relE);
      double dvoltsdt = calcdVdT_mao(relE);
      double F = 96485.3E3;   // Joules/kmol
      double DeltaS = - F * dvoltsdt;
      double DeltaG = volts_B * F;
      double T = 298.15;
      double DeltaH = DeltaG + T * DeltaS;
        
      printf("%15.5E %15.5E %10.5f  %10.5f       % 10.5E  % 10.5E   % 10.5E  % 10.5E\n", 
             relE, 1.0-relE, volts_A, volts_B, dvoltsdt, DeltaS , DeltaH, DeltaG);
    }
    return 0;
}
