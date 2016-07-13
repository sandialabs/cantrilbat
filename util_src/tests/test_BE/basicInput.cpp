/**
 *  @file basicInput.cpp
 *     This is a simple test program that reads 5 basic types of block
 *     entries, 3 different ways, and makes sure they all work.
 */
/*
 * $Author: hkmoffa $
 * $Revision: 17 $
 * $Date: 2012-03-23 11:42:20 -0600 (Fri, 23 Mar 2012) $
 */
/*
 * Copywrite 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BE_BlockEntry.h"
#include "LE_PickList.h"
#include "BE_MoleComp.h"
#include "BE_UnitConversionPressure.h"
#include "BE_UnitConversionLength.h"
#include "BE_MultiBlock.h"
#include "LE_OneDblUnits.h"
#include "LE_OneStr.h"
#include "LE_OneBoolInt.h"
#include "LE_OneBool.h"
#include "LE_OneDbl.h"
#include "LE_OneInt.h"
#include "BI_DepIntMaxMin.h"
#include "md_timer.h"
#include "mdp_allo.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <cstdio>
#include <climits>
#include <cfloat>
#include <cstdlib>

using namespace std;
using namespace BEInput;
using namespace mdpUtil;

static void printCW() {
    printf("\n");
    printf("Copywrite 2005 Sandia Corporation. Under the terms of Contract\n");
    printf("DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government\n");
    printf("retains certain rights in this software.\n");
    printf("\n");
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/



/*
 * Setup options for the program. -> These get initialized
 * via interaction with the user through the input file.
 */
class ProgramOptions {
public:
  /*
   * TDcads Model Parameters Block
   */
  bool IncludeGasSpecies;
  int IncludeTemperature;
  int BoundaryConditions;
  int GasID;

  /*
   * TDcads Time Step Parameters Block
   */
  int    Debugging;
  double TFinal;
  int PrintFlag;
  int printSolnSteps;
  int printSolnInterval;

  /*
   * TDcads Initial Conditions Block
   */  
  double Temperature;
  double Pressure;
  double VolInit;
  double *mfInit;  
  string GasFile;
  bool   InitPartSmallPositive;
  string ParticlePhaseFile;
  string SolnFile;

  /*
   * Specification of the inflow conditions
   * when psr option is used
   */
  double InflowTemperature;
  double InflowMolarFlowRate;
  double *InflowGasMoleFractions;
   
   
  /*
   * Default values for the options
   */
  ProgramOptions () :
    IncludeGasSpecies (false),
    IncludeTemperature(0),
    BoundaryConditions(0),
    GasID(0),
    Debugging(0),
    TFinal(5000.0),  // Time in seconds
    PrintFlag(1),
    printSolnSteps(1),
    printSolnInterval(0),
    Temperature(298.15), // Temperature in Kelvin
    Pressure(101325.),  // Pressure in pascals
    VolInit(1.0E-6),    //1 cm**3 as an initial volume
    mfInit(0),
    InitPartSmallPositive(true),
    InflowTemperature(300.),
    InflowMolarFlowRate(0.),
    InflowGasMoleFractions(0)  {
  }

  ~ProgramOptions() {
    mdp_safe_free((void **) &mfInit);
    mdp_safe_free((void **) &InflowGasMoleFractions);

  }
} PO;


/******************************************************************************
 *
 *
 */
void printUsage() {
  cout << "usage: TDcads [-h] [-help_cmdfile] cmdfile.txt"
       <<  endl;
  cout << "    -h           help" << endl;
  cout << "    -help_cmdfile Provides a listing of "
    "cmdfile.txt keylines" << endl;
  cout << "                  "
    "(with the cmdfile, it will print out species names)" 
       << endl;
}


/******************************************************************************
 *
 *
 */
void setup_input_pass(BlockEntry *topbe, bool printInputFormat) {
 
  /*
   * Create a list of species to be used in the input file.
   */
  int nSpecies = 4;
  char **GSList = mdp_alloc_VecFixedStrings(nSpecies,
					    MAX_INPUT_STR_LN+1);
  strcpy( GSList[0] , "H2" );
  strcpy( GSList[1] , "N2" );
  strcpy( GSList[2] , "O2" );
  strcpy( GSList[3] , "CH4" );


  /*
   *  create a section method description block and start writing 
   *  line elements in it. 
   *       BLOCK TDCADS MODEL DEFINITION
   */
  BlockEntry *bemd = new BlockEntry("TDcads Model Definition");
  topbe->addSubBlock(bemd);

  /* --------------------------------------------------------------
   * Include Gas Species in Solution Vector = [boolean]
   * - Defaults to no.
   */
  LE_OneBool *bincg =
    new LE_OneBool("Include Gas Species in Solution Vector",
		   &PO.IncludeGasSpecies, 1);
  bincg->set_default(false);
  bemd->addLineEntry(bincg);

  /* --------------------------------------------------------------
   * Include Temperature in Solution Vector = [boolean]
   * - Defaults to no.
   */
  LE_OneBoolInt *binct =
    new LE_OneBoolInt("Include Temperature in Solution Vector",
		      &PO.IncludeTemperature, 1);
  binct->set_default(0);
  bemd->addLineEntry(binct);

  /* --------------------------------------------------------------
   * Boundary Conditions = ["constant pressure"]
   *                       ["constant volume"]
   *                       ["psr"]
   * This specifies how the gas conditions change as a function
   * of the reactions that are occurring.
   * psr means that the time dependent psr equations are solved.
   */
  const char *cbc[3] = {"constant pressure", "constant volume", "psr" };
  LE_PickList *lebc =
    new LE_PickList("Boundary Conditions",
		    &PO.BoundaryConditions, cbc, 3, 1,
		    "BoundaryConditions");
  lebc->set_default("constant pressure");
  bemd->addLineEntry(lebc);

  /* --------------------------------------------------------------
   * Volume of Reactor = [double] [vol units]
   *    (optional) (default = 1 cm**3)
   *    This is used as the initial or constant volume of the
   *    reactor. Lengths are converted to the mks units.
   */
  BE_UnitConversion *ucL = new BE_UnitConversionLength();
  LE_OneDblUnits *levol = 
    new LE_OneDblUnits("Volume of Reactor", &(PO.VolInit), 0,
		       "PO.VolInit", ucL);
  levol->set_default(1.0E-6);
  levol->set_limits(1.0E100, 1.0E-100);
  bemd->addLineEntry(levol);



  /*****************************************************************
   *  create a subblock
   *      
   */
  BlockEntry *betic = new BlockEntry("TDcads Inflow Conditions");
  bemd->addSubBlock(betic);
  /*
   * If the boundary condition is "psr", i.e., 2, we require the
   * TDcads inflow Conditions block.
   * We require that the "Boundary Condition entry be previously
   * processed and have a value of 2.
   */
  BI_DepIntMaxMin* deptic =
    new BI_DepIntMaxMin(lebc, BIDT_INTMAXMIN, 2, 2,
			BIDRT_ONENUMTR);
  betic->declareDependency(deptic);

  /* --------------------------------------------------------------
   *  Inflow Molar Flow Rate = [ double ]
   *     (required)
   *  Temperature of the inflow stream.
   *  This is used to establish the enthalpy of the input stream.
   */
  LE_OneDbl *ledimfr = new LE_OneDbl("Inflow Molar Flow Rate", 
				     &(PO.InflowMolarFlowRate), 1,
				     "PO.InflowMolarFlowRate");
  betic->addLineEntry(ledimfr);



  /* --------------------------------------------------------------
   * Inflow Gas Mole Fractions
   *
   *   START BLOCK Inflow Gas Mole Fraction
   *      H2 = 0.1
   *      O2 = 0.54
   *      ...
   *   END BLOCK Inflow Gas Mole Fraction
   */
  PO.InflowGasMoleFractions = mdp_alloc_dbl_1(nSpecies, 0.0);
  PO.InflowGasMoleFractions[0] = 1.0;
  BE_MoleComp *bigmf = 
    new BE_MoleComp("Inflow Gas Mole Fractions",
		    &(PO.InflowGasMoleFractions), 1,
		    GSList, nSpecies, 1, 
		    "InflowGasMoleFractions");
  bigmf->generateDefLE();
  betic->addSubBlock(bigmf);
 

  mdp_safe_free((void **) &GSList);
  BaseEntry::set_SkipUnknownEntries(false);
}



/******************************************************************************
 *  Read the input file
 *
 *  printFlag 0 = no output
 *            1 = error message output
 *            2 = output
 */
bool process_input(BlockEntry *cf, string fileName, int printFlag) {
  static int pass = 0;
  pass++;
  cf->ZeroLineCount();
  const TKInput::TOKEN tok_in;
  TKInput::TOKEN tok_out;
  FILE *ifp = fopen(fileName.c_str(), "r");
  if (!ifp) {
    if (printFlag) {
      cout << "ERROR can't open file " << fileName<< endl;
    }
    return false;
  }
  if (printFlag > 1) {
    printf("==========================================================\n");
    printf(" STARTING PROCESSING COMMAND FILE %s, PASS # %d\n",
	   fileName.c_str(), pass);
    printf("==========================================================\n");
  }
  try {
    /*
     * Call the block read function at the main level
     */
    cf->ZeroLineCount();
    cf->read_block(ifp, &tok_out, &tok_in, 0);
  } catch (BI_InputError &bi) {
    /*
     * This catches error messages
     */
    cout << bi.errorMessage() << endl;
    return false;
  }
  fclose(ifp);
  if (printFlag > 1) {
    printf("=========================================================\n");
    printf(" FINISHED PROCESSING COMMAND FILE %s, PASS # %d\n",
	   fileName.c_str(), pass);
    printf("=========================================================\n");
  }
  return true;
}

void read_info(ProgramOptions& P1, BlockEntry *topbe) {

  /***************************************/
  // One Bool check
  LineEntry *le = topbe->searchLineEntry("Include Gas Species in Solution Vector");
  const bool *bptr = static_cast<const bool *>(le->currentValueAsVoidP());
  P1.IncludeGasSpecies = *bptr;
 
  LE_OneBool *le_bool = dynamic_cast<LE_OneBool *>(le);
  bool bvalue = le_bool->currentTypedValue();
  if (bvalue != P1.IncludeGasSpecies) {
    printf ("we have an inconsistency: bvalue != P1.IncludeGasSpecies\n");
    exit(-1);
  }

  if (PO.IncludeGasSpecies !=  P1.IncludeGasSpecies) {
   printf ("we have an inconsistency: PO.IncludeGasSpecies !=  P1.IncludeGasSpecies\n");
   exit(-1);
  }

  /***************************************/
  // One BoolInt check
  LineEntry *le1 = topbe->searchLineEntry("Include Temperature in Solution Vector");
  const int *iptr = static_cast<const int *>(le1->currentValueAsVoidP());
  P1.IncludeTemperature = *iptr;
 
  LE_OneBoolInt *le_boolInt = dynamic_cast<LE_OneBoolInt *>(le1);
  int ivalue = le_boolInt->currentTypedValue();
  if (ivalue != P1.IncludeTemperature) {
    printf ("we have an inconsistency: ivalue != P1.IncludeTemperature\n");
    exit(-1);
  }

  if (PO.IncludeTemperature != P1.IncludeTemperature) {
   printf ("we have an inconsistency: PO.IncludeTemperature != P1.IncludeTemperature \n");
   exit(-1);
  }

  /***************************************/
  // PickList check
  LineEntry *le2 = topbe->searchLineEntry("Boundary Conditions");
  const int *i2ptr = static_cast<const int *>(le2->currentValueAsVoidP());
  P1.BoundaryConditions = *i2ptr;
 
  LE_PickList *le_pl = dynamic_cast<LE_PickList *>(le2);
  int i2value = le_pl->currentTypedValue();
  if (i2value != P1.BoundaryConditions) {
    printf ("we have an inconsistency: i2value != P1.BoundaryConditions\n");
    exit(-1);
  }

  if (PO.BoundaryConditions != P1.BoundaryConditions) {
   printf ("we have an inconsistency: PO.BoundaryConditions != P1.BoundaryConditions \n");
   exit(-1);
  }

  /***************************************/
  // OneDblUnits check
  le = topbe->searchLineEntry("Volume of Reactor");
  const double *dptr = static_cast<const double *>(le->currentValueAsVoidP());
  P1.VolInit = *dptr;
 
  LE_OneDblUnits *le_du = dynamic_cast<LE_OneDblUnits *>(le);
  double dvalue = le_du->currentTypedValue();
  if (dvalue != P1.VolInit) {
    printf ("we have an inconsistency: dvalue != P1.VolInit\n");
    exit(-1);
  }

  if (PO.VolInit != P1.VolInit) {
   printf ("we have an inconsistency: PO.VolInit != P1.VolInit \n");
   exit(-1);
  }


  /***************************************/
  // OneDbl check
  le = topbe->searchLineEntry("Inflow molar flow rate");
  const double *d2ptr = static_cast<const double *>(le->currentValueAsVoidP());
  P1.InflowMolarFlowRate = *d2ptr;
 
  LE_OneDbl *le_du2 = dynamic_cast<LE_OneDbl *>(le);
  double dvalue2 = le_du2->currentTypedValue();
  if (dvalue2 != P1.InflowMolarFlowRate) {
    printf ("we have an inconsistency: dvalue2 != P1.InflowMolarFlowRate\n");
    exit(-1);
  }

  if (PO.InflowMolarFlowRate != P1.InflowMolarFlowRate) {
   printf ("we have an inconsistency: PO.InflowMolarFlowRate != P1.InflowMolarFlowRate \n");
   exit(-1);
  }

  /***********************************************/
  // Mole Frac check
  int nSpecies = 4;
  mdp_safe_free((void **) &P1.InflowGasMoleFractions);
  P1.InflowGasMoleFractions = mdp_alloc_dbl_1(nSpecies, 0.0);

  
  BlockEntry *be = topbe->searchBlockEntry("Inflow Gas Mole Fractions");
  const double *moleFractions = 
    static_cast<const double *>(be->currentValueAsVoidP());
  for (int i = 0; i < nSpecies; i++) {
    (P1.InflowGasMoleFractions)[i] = moleFractions[i];
  }
 
  
  BE_MoleComp *be_dbl = dynamic_cast<BE_MoleComp *>(be);
  moleFractions = be_dbl->currentTypedValue();

  for (int i = 0; i < nSpecies; i++) {
    if ((P1.InflowGasMoleFractions)[i] != moleFractions[i]) {
      printf ("we have an inconsistency: P1.InflowGasMoleFractions != moleFractions \n");
      exit(-1);
    }
  }

  for (int i = 0; i < nSpecies; i++) {
    if ((P1.InflowGasMoleFractions)[i] != (PO.InflowGasMoleFractions)[i]) {
      printf ("we have an inconsistency: "
	      "P1.InflowGasMoleFractions != PO.InflowGasMoleFractions \n");
      exit(-1);
    }
  }

}

/******************************************************************************
 *
 *
 */
int main(int argc, char** argv) { 

  bool ok = true;
  bool printInputFormat = false; // print cmdfile.txt format
  bool printedUsage = false; // bool indicated that we have already
  // printed usage
  bool noProperInit = false; // Sometimes whn you want to print
  // the cmdfilfe format, you don't have
  // a gas file or anything.
  string cmdfile;
  int printBIProclevel = 2;
  /**
   *
   *      Look for command-line options
   */
  if (argc > 1) {
    string tok;
    for (int j = 1; j < argc; j++) {
      tok = string(argv[j]);
      if (tok[0] == '-') {
	int nopt = static_cast<int>(tok.size());
	for (int n = 1; n < nopt; n++) {
	  if (!strcmp(tok.c_str() + 1, "help_cmdfile")) {
	    printInputFormat = true;
	  } else if (tok[n] == 'h') {
	    printUsage();
	    printedUsage = true;
	  } else {
	    printUsage();
	    printedUsage = true;
	    exit(1);
	  }
	}
      } else if (cmdfile == "") {
	cmdfile = tok;

      } else {
	printUsage();
	printedUsage = true;
	exit(1);
      }
    }
  }
  if (!printInputFormat) {
    if (cmdfile == "") {
      if (!printedUsage) {
	cout << "ERROR: command file not found" << endl;
	printUsage();
      }
      exit(1);
    }
  } else {
    if (cmdfile == "") {
      printBIProclevel = 0;
      BEInput::BI_SetPrintLevel(0);
      noProperInit = true;
    }
  }

  printCW();

  /*
   * General Catch block to trap Errors and print them
   */
  try {

    /**
     * Initialize a block input structure for the command filfe
     */
    BlockEntry *cf = new BlockEntry("command_file");
    /**
     * Setup and process the input deck from standard input
     * -> Might have to print out the input and quit as well.
     */
    setup_input_pass(cf, printBIProclevel);
   
  

    ok = process_input(cf, cmdfile, printBIProclevel);
    if (!ok) {
      if (!printInputFormat) {
	return -1;
      }
    }

    /*
     * OK, now go get all of the information by direct query
     */

    ProgramOptions * P1_ptr = new ProgramOptions();

    read_info(*P1_ptr, cf);

    /*
     * we are suppose to do multiple passes on the input file.
     * Let's do it, and see if we get the correct result again,
     * after zeroing the counters.
     */
    cf->ZeroLineCount();
    ok = process_input(cf, cmdfile, printBIProclevel);
    if (!ok) {
      if (!printInputFormat) {
	return -1;
      }
    }
    ProgramOptions * P2_ptr = new ProgramOptions();

    read_info(*P2_ptr, cf);
    read_info(*P1_ptr, cf);

    delete cf;

    delete P1_ptr;
    delete P2_ptr;

  }
  /*
   * This block print out all unidentified errors and exits
   * with an error flag.
   */
  catch (BI_InputError &bie) {
    cout << bie.what() << endl;
    return 1;
  }

  return 0;
}

/*****************************************************************************/
