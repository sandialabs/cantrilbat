/*
 * $Id: electrodeCell_findI.cpp 544 2013-03-01 16:12:29Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#include "zuzax/equil/vcs_internal.h"

#include "electrodeCell_prep.h"
#include "electrodeCell_kin.h"

#include <stdio.h>

using namespace std;
using namespace Zuzax;


double findV(Electrode *electrode, double Itarget,
			      double Elow, double Ehigh,
			      int printLvl, double err,
			      int maxsteps) {
 
  double deltaVoltage;
  double dT, dTa, dTmax, Enew;
  double Inow;
  double Ilow = Undef;
  double Ihigh = Undef;
  double Ierr, IConvErr;
  double Enow = 0.0;
  double dIdV = 0.0;

  InterfaceKinetics *iK = electrode->reactingSurface(0);

  deltaVoltage = electrode->voltage();
  double Epast = deltaVoltage - 0.01;
  //electrode->setDeltaVoltage(Epast);
  electrode->setVoltages(Epast, 0.0);
  double Ipast = processGERCurrent(electrode->m_rmcEGR[0], electrode, 0, *iK,  *(electrode->m_egr[0]));
  //electrode->setDeltaVoltage(Epast + 0.01);
  electrode->setVoltages(Epast + 0.01, 0.0);

  for (int n = 0; n < maxsteps; n++) {
    // start with a loose error tolerance, but tighten it as we get
    // close to the final temperature
 
    try {
      deltaVoltage = electrode->voltage();
      Enow = deltaVoltage;
   
      // May have to determine the sign of the slope here. We will see.
      Inow = processGERCurrent(electrode->m_rmcEGR[0], electrode, 0, *iK,  *(electrode->m_egr[0]));
      
      if (n == 0) {
	dIdV = (Inow - Ipast)/(Enow - Epast);
      }
   
      if (printLvl > 0) {
	plogf("Inow = %g Enow = %g\n", Inow, Enow);
	plogendl();
      }
      if (Inow < Itarget) {
	if (dIdV > 1.0E-7) {
	  if (Enow > Elow) {
	    Elow = Enow;
	    Ilow = Inow;
	  }
	} else if (dIdV < -1.0E-7) {
	  if (Enow < Ehigh) {
	    Ehigh = Enow;
	    Ihigh = Inow;
	  }
	}
      } else {
	if (dIdV > 1.0E-7) {
	  if (Enow < Ehigh) {
	    Ehigh = Enow;
	    Ihigh = Inow;
	  }
	} else if (dIdV < -1.0E-7) {
	  if (Enow > Elow) {
	    Elow = Enow;
	    Ilow = Inow;
	  }
	}
      }
      if (Ilow != Undef && Ihigh != Undef) {
	dIdV = (Ihigh - Ilow)/(Ehigh - Elow);
	dT = (Itarget - Inow)/dIdV;
	dTa = fabs(dT);
	dTmax = 0.5*fabs(Ehigh - Elow);
	if (dTa > dTmax) dT *= dTmax/dTa;
      }
      else {
	double Enew = 0.5*(Elow + Ehigh);
	dT = Enew - Enow;
	if (dT < -0.1) dT = -0.1;
	if (dT >  0.1) dT =  0.1;
      }
      double acpb = std::max(fabs(dIdV), 1.0E-6);
      double denom = std::max(fabs(Itarget), acpb);
      Ierr = Itarget - Inow;
      IConvErr = fabs((Ierr)/denom);
  
      if (printLvl > 0) {
	plogf("   findV: It = %d, Ecurr  = %g Icurr = %g, Itarget = %g\n",
	      n, Enow, Inow, Itarget);
	plogf("                   I rel error = %g, dI/dE = %g, IConvErr = %g\n",
	      Ierr, dIdV, IConvErr);
      }

      if (IConvErr < err) { // || dTa < 1.0e-4) {

	if (printLvl > 0) {
	  plogf("   findV: CONVERGENCE: Ifinal  = %g Efinal = %g, Its = %d \n",
		Inow, Enow, n);
	  plogf("                   I rel error = %g, dIdE = %g, IConvErr = %g\n",
		Ierr, dIdV, IConvErr);
	}
	goto done;
      }
      Enew = Enow + dT;
      // electrode->setDeltaVoltage(Enew);
      electrode->setVoltages(Enew, 0.0);

    }
    catch (ZuzaxError err) {
   
      Enew = 0.5*(Enow + Ehigh);
      if (fabs(Enew - Enow) < 1.0) Enew = Enow + 1.0;
      //electrode->setDeltaVoltage(Enew);   
      electrode->setVoltages(Enew, 0.0);
    }

  }

  throw ZuzaxError("ELECTRODE_MODEL:findV","No convergence for I");
 done:;
  deltaVoltage = electrode->voltage();
  return deltaVoltage;
}
