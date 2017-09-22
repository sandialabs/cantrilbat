/*
 * $Id: LiCoO2_Cathode_3_cc.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "cantera/numerics/NonlinearSolver.h"
#include "cantera/equilibrium.h"
#include "Electrode_Factory.h"
#include "importPL.h"
#include "BE_BlockEntry.h"
#include <iomanip>
#include "EState_XML.h"

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
using namespace esmodel;

// a lvl of one prints out the .csv file
int mpequil_debug_print_lvl = 1;
int VCS_Debug_Print_Lvl = 3;

void printUsage() {
    cout << "usage: electrodeCell [-h] [-help_cmdfile] [-d #] [mpequil.inp]"
         <<  endl;
    cout << "    -h           help" << endl;
    cout << "    -d           #   : level of debug printing" << endl;
    cout << "  electrodeCell.inp    : command file" << endl;
    cout << "                     : (if missing, assume mpequil.inp)" 
	 << endl;
    cout << endl;
}
//=====================================================================================================
int main(int argc, char **argv)
{
  int retn = 0;
  //bool doCathode = false;
  string commandFileNet = "cell.inp";

  string commandFileA = "anode.inp";
  // printed usage

  //VCSnonideal::vcs_timing_print_lvl = 0;
  NonlinearSolver::s_TurnOffTiming = true;
  NonlinearSolver::s_print_NumJac = true;


  ZZCantera::XML_Node* xEout =  getElectrodeOutputFile("solnSave_0_0.xml", 1);
  
  EState_ID_struct e_id;
  e_id.readIdentificationFromXML(*xEout);

  EState* es = newEStateObject(e_id.EState_Type_String);
  es->readIdentificationFromXML(*xEout); 

  int globalTimeStepNum = 0;
  ZZCantera::XML_Node* x = selectLastGlobalTimeStepInterval(xEout, globalTimeStepNum);

  double timeVal;
  ZZCantera::XML_Node* xSt =  locateTimeLast_GlobalTimeStepIntervalFromXML(*x, timeVal, 1);
  es->readStateFromXML(*xSt);
  //es->readStateFromXML(*x);

  printf("read global time step num %d\n", globalTimeStepNum);

  return retn;

} 
