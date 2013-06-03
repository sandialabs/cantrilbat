/**
 *
 *  @file 
*
*   Small program to reproduce Eqn. 36 Silvester and Pitzer
 *  J. Phys. Chem. , 81, 1822 (1977).
 *
 *  Delta G0 for NaCl(solid) -> Na+ + Cl-
 *  in water
 */
                                   
/*
 *  $Author: hkmoffa $
 *  $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
 *  $Revision: 5 $
 */
#include <stdio.h>
#include <float.h>
#include <math.h>

using namespace std;
double q16 = 41587.11;
double q17 = -315.90;
double q18 = 0.8514;
double q19 = -8.3637E-4;
double q20 = -1729.93;
double q20delta = -0.001669;

double delG(double T) {
    double tt = T * T;
    double tmp = q16 - q17 * T * log(T) - q18 * tt - 0.5 * q19 * tt * T
	+ q20 * T + q20delta * T;
    return tmp;
}

double delH(double T) {
    double tt = T * T;
    double tmp = q16 + q17 * T + q18 * tt + q19 * tt * T;
    return tmp;
}

double delCp(double T) {
    double tt = T * T;
    double tmp =  q17  + 2.0 * q18 * T + 3.0 * q19 * tt;
    return tmp;
}

double delS(double T) {
    double tt = T * T;
    double tmp =  q17 * (log(T) + 1.0) + 2.0 * q18 * T + 3. / 2. * q19 * tt -(q20+q20delta);
    return tmp;
}


int main(int argc, char **argv)
{

    double T[100];

    int itop = 30;
    T[0] = 293.15;
    T[1] = 298.15;
    double tt = 303.15;
    double tmp;
    int i;
    for (i = 2; i < itop; i++) {
	T[i] = tt;
	tt += 10.;
    }

    printf(" Temperature  DeltaG0       DeltaG0        DeltaH0         DeltaCp0          DeltaS0\n");
    printf("              (cal/gmol)  (kJ/gmol)        (kJ/gmol)       (J/gmolK)         (J/gmolK) \n");    
    for (i = 0; i < itop; i++) {
      tmp = delG(T[i]);
      double tmp2 = tmp * 4.184 / 1.0E3;
      double dh   = delH(T[i]);
      double dcp  = delCp(T[i]);
      double ds   = delS(T[i]);
      dh *= 4.184 / 1.0E3;
      dcp *= 4.184;
      printf("%g        %g          %g        %g      %g         %g\n", T[i], tmp, tmp2, dh, dcp, ds); 
    }

}
