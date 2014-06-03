

#include <cstdio>
#include <cmath>

 double calcV( double relExtent)
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

int main()
{
   int numP = 51;
    
   printf("   relExtent           xLi          OCV \n");
   for (int i = 0; i < numP; i++) {
      double relE = 0.0 + (double (i) / (numP - 1.0));
      double volts = calcV(relE);
      printf("%15.5E %15.5E %10.5f \n", relE, 1.0-relE, volts);
    }
    return 0;
}
