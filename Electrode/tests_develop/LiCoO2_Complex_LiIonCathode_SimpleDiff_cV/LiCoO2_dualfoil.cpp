

#include <cstdio>
#include <cmath>

/*
c     Marc Doyle's fit
c     r1=4.825510d0
c     r2=0.950237d0
c     r3=0.913511d0
c     r4=0.600492d0
c     g0=r1-r2*expg(-((xx(2+mpa,j)/ct3-r3)/r4)**2)
c     g1=2.0d0*r2*(xx(2+mpa,j)/ct3-r3)*expg(-((xx(2+mpa,j)/ct3
c    1-r3)/r4)**2)/r4/r4/ct3
c
c     Measured by Oscar Garcia 2001 using Quallion electrodes for
c     0.5 < y < 0.99.  Fit revised by Karen Thomas in May 2003 to
c     match Doyle's fit for y < 0.4 and Garcia's data at larger y.
c     Valid for 0 < y < 0.99. Note that capacity fade is found to
c     occur experimentally if y goes below 0.5; this is not included
c     in the model.

  g0 = 2.16216d0+0.07645d0*dtanh(30.834d0-54.4806d0*sto)
    & + 2.1581d0*dtanh(52.294d0-50.294d0*sto)
     & - 0.14169d0*dtanh(11.0923d0-19.8543d0*sto)
     & + 0.2051d0*dtanh(1.4684d0-5.4888d0*sto)
     & + 0.2531d0*dtanh((-sto+0.56478d0)/0.1316d0)
     & - 0.02167d0*dtanh((sto-0.525d0)/0.006d0)

*/

 double calcV( double relExtent)
{
   double volts;

   relExtent += 0.04;
   if (relExtent > 1.0 ) {
      // relExtent = 1.0; 
   }
   double sto = relExtent;
   volts = (    2.16216 
              + 0.07645 * tanh( 30.834  - 54.4806 * sto)
              + 2.1581  * tanh( 52.294  - 50.294  * sto)
              - 0.14169 * tanh( 11.0923 - 19.8543 * sto)
              + 0.2051  * tanh( 1.4684  - 5.4888  * sto)
              + 0.2531  * tanh((-sto + 0.56478) / 0.1316)
              - 0.02167 * tanh((sto  -0.525)    / 0.006)
           );

return volts;
} 

int main()
{
   int numP = 101;
    
   printf("   relExtent           xLi          OCV \n");
   for (int i = 0; i < numP; i++) {
      double relE = 0.0 + (double (i) / (numP - 1.0));
      double volts = calcV(relE);
      printf("%15.5E %15.5E %10.5f \n", relE, relE, volts);
    }
    return 0;
}
