/* ====================================================================== */
/* $RCSfile: water_SSThermo.cpp,v $ */
/* $Author: hkmoffa $ */
/* $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $ */
/* $Revision: 5 $ */
/* ====================================================================== */

#include <iostream>

#include "ct_defs.h"
#include "TTPhase.h"
#include "PhaseList.h"

#include "IdealGasMix.h"

#include "TPX_SSThermo.h"



using namespace Cantera;
using namespace std;

#ifndef TRUE
# define TRUE 1
#endif
#ifndef FALSE
# define FALSE 0
#endif

#ifndef BOOLEAN
# define BOOLEAN int
#endif

#ifndef MAX
# define MAX(x,y) (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
# define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

#ifndef SWAP
# define SWAP(x1, x2, temp) ((temp) = (x1), (x1) = (x2), (x2) = (temp))
#endif

#ifndef SQUARE
# define SQUARE(x) ((x) * (x))
#endif

#ifndef DSIGN
# define DSIGN(x) (( (x) == (0.0) ) ? (0.0) : ( ((x) > 0.0) ? 1.0 : -1.0 ))
#endif

/**************************************************************************/
/*
*    If we have a second function from another program that is compiled
* in, then we won't use vcs_second function.
*
*/
#ifdef HAVE_MD_TIMER
#  define vcs_second second
extern double second(void);
#endif


int main()
{
  try {
    int subID;
    phase_t *phase = new phase_t();
    bool stoichp = true;
    int phaseLimitation = LIQUID;
    int spID = TPX_SSThermo::preparePhase(Pure_Water, phase,
					  phaseLimitation, stoichp);
    TPX_SSThermo w(Pure_Water, phase, spID, 
		   stoichp, phaseLimitation);

    w.setState_TP(500., 1.013E7);
    double result = w.satPressure(500.);
    printf("result = %g\n", result);

    printf("Critical Temperature = %g Kelvin -> CRC value = %g \n",
	   w.critTemperature(), 373.99 + 273.15);
    printf("Critical Pressure    = %g Pa -> CRC value = %g\n",  
	   w.critPressure(), 22.064 * 1.0E6);
    printf("Critical Volume = %g m^3/kg -> CRC value = %g\n",
	   1.0 /w.critDensity(), 1.0/ (0.322) / 1.0E6 * 1.0E3);
    printf("Molecular Weight = %g gm/mol > CRC value = %g\n",
	 phase->molecularWeight(0), 18.01528);

  w.setState_TP(273.161, 1.013E5);
  printf("PSat(0 C) = %g Pa -> CRC value = %g\n",
	 w.satPressure(273.161), 0.6113 * 1.0E3);

  w.setState_TP(273.15 + 20, 1.013E5);
  printf("PSat(20C) = %g Pa -> CRC value = %g\n",
	 w.satPressure(273.15 + 20), 2.3388 * 1.0E3);

  w.setState_TP(273.15 + 25, 1.013E5);
  printf("PSat(25C) = %g Pa\n", w.satPressure(273.15 + 25.0));

  w.setState_TP(273.15 + 40, 1.013E5);
  printf("PSat(40C) = %g Pa -> CRC value = %g\n",
	 w.satPressure(273.15 + 40.0), 7.3814 * 1.0E3);

  w.setState_TP(273.15 + 60, 1.013E5);
  printf("PSat(60C) = %g Pa -> CRC value = %g\n",
	 w.satPressure(273.15 + 60.), 19.932 * 1.0E3);

  w.setState_TP(273.15 + 100, 1.013E5);
  printf("PSat(100C) = %g Pa -> CRC value = %g\n",
	 w.satPressure(373.15), 101.325 * 1.0E3);

  w.setState_TP(298.15, 1.013E5);
  printf("s = %g J/(kg K)\n", w.entropy_mass());
  printf("h = %g J/kg\n", w.enthalpy_mass());
  printf("s = %g J/(K gmol) (CRC = 70.0) \n", w.entropy_mole() * 1.0E-3);
  printf("h = %g kJ/gmol (CRC = -285.8)\n",
	 w.enthalpy_mole() * 1.0E-6);
  printf("Cp = %g J/(K gmol) (CRC = 75.3)\n", 
	 w.cp_mole() * 1.0E-3);
  printf("G = %g kJ/gmol\n",
	 w.gibbs_mole() * 1.0E-6);
    
    delete phase;
  } catch(CanteraError &CE) {
    cout << "ERROR: Cantera threw an error message:" << endl;
    cout << CE.errorMessage() << endl;
  } 
}
 
