/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 508 $
 *
 * 
 */

#ifndef DEBUG_HKM
#define DEBUG_HKM
#endif

#include "cantera/IdealGasMix.h"
#include <cantera/thermo/IonsFromNeutralVPSSTP.h>
#include <cantera/thermo/MargulesVPSSTP.h>
#include <cantera/thermo/MolalityVPSSTP.h>
#include <cantera/thermo/PureFluidPhase.h>
#include "cantera/equilibrium.h"
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <fstream>

using namespace std;
using namespace Cantera;

int main(int argc, char **argv) {
#ifdef DEBUG_CHEMEQUIL
  Cantera::ChemEquil_print_lvl = 0;
#endif
  try {

    //ofstream textFile("gasData.txt");
    //ofstream csvFile("gasData.csv");

    IdealGasPhase* gas = new IdealGasMix("air_below6000K.cti","air_below6000K");
    /*
    //GibbsExcessVPSSTP* test1;
    IonsFromNeutralVPSSTP* test2 = new IonsFromNeutralVPSSTP( "LiKCl_recipmoltenSalt_trans.xml" );
    MargulesVPSSTP* test3 = new MargulesVPSSTP( "LiKCl_Margules_trans.xml" );
    //MolalityVPSSTP* test4;
    //PureFluidPhase* test5;
    //ThermoPhase* test6;

    gas->reportCSV(textFile,csvFile);
    //test1->reportCSV(textFile,csvFile);
    test2->reportCSV(textFile,csvFile);
    test3->reportCSV(textFile,csvFile);
    //test4->reportCSV(textFile,csvFile);
    //test5->reportCSV(textFile,csvFile);
    //test6->reportCSV(textFile,csvFile);

    */

    vector_fp IndVar2(6, 0.0);
    IndVar2[0] = 1.5E5;
    IndVar2[1] = 3.0E5;
    IndVar2[2] = 9.0E5;
    IndVar2[3] = 2.7E6;
    IndVar2[4] = 6.7E6;
    IndVar2[5] = 1.0E7;

    vector_fp IndVar1(7, 0.0);
    IndVar1[0] = 1.0E-8;
    IndVar1[1] = 1.0E-7;
    IndVar1[2] = 1.0E-6;
    IndVar1[3] = 1.0E-5;
    IndVar1[4] = 1.0E-4;
    IndVar1[5] = 1.0E-3;
    IndVar1[6] = 1.0E-2;
    int nj = 6;
    int ni = 7;

    for (int j=0; j<nj; j++) {
      for (int i=0; i<ni; i++) {
	stringstream fileName;
	fileName << "gas_" << i << "_" << j <<"_Data.csv";
	ofstream csvFile(fileName.str().c_str());
	double offset = -301471.39;
	gas->setState_UV(IndVar2[j]+offset,1.0/IndVar1[i]);
	double tkelvin = gas->temperature();
        if (tkelvin > 6000.) {
          tkelvin = 5900.;
        }
	double pres = gas->pressure();
	printf("Initial T = %g, pres = %g atm \n", tkelvin, pres/OneAtm);
	//textFile << setw(20) << "Initial T (K)," << setw(20) << "Pressure (atm)\n"; 
	//beginLogGroup("topEquil", -1);
	equilibrate(*gas,"UV", -1);
	//endLogGroup("topEquil");
	cout << gas->report() << endl;
	//gas->reportCSV(textFile,csvFile);
#if defined (CANTERA_VERSION_18_LTD)
       	gas->reportCSV(csvFile);
#endif
	tkelvin = gas->temperature();
	pres = gas->pressure();
	printf("Final T = %g, pres = %g atm\n", tkelvin, pres/OneAtm);
	//textFile << setw(20) << "Initial T (K)," << setw(20) << "Pressure (atm)\n"; 
   
      }
    }
    delete gas;
  }

  catch (CanteraError) {
    showErrors();
  }
  return 0;
}
