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


#include "BlockEntryGlobal.h"
#include "mdp_allo.h"
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <fstream>
#include <cstdlib>

#include <climits>
#include <cfloat>


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

class InputPipe {
public:
  double InflowMolarFlowRate;
  double *InflowGasMoleFractions;

  InputPipe() :
    InflowMolarFlowRate(0.0),
    InflowGasMoleFractions(0)
  {
    InflowGasMoleFractions = mdp_alloc_dbl_1(4, 0.0);
  }

  ~InputPipe() {
    mdp_safe_free((void **) &InflowGasMoleFractions);
  }
};


/*
 * Setup options for the program. -> These get initialized
 * via interaction with the user through the input file.
 */
class ProgramOptions {
public:
  /*
   * TDcads Model Parameters Block
   */
  int BoundaryConditions;
  int GasID;

  /*
   * TDcads Time Step Parameters Block
   */
  int    Debugging;
  double TFinal;
  int PrintFlag;
  int printSolnSteps;


  /*
   * TDcads Initial Conditions Block
   */  
  double Temperature;
  double Pressure;
  double VolInit;
  InputPipe **ipl;
  int numIPL;

  /*
   * Specification of the inflow conditions
   * when psr option is used
   */
  double InflowTemperature;
  double InflowMolarFlowRate;
   
   
  /*
   * Default values for the options
   */
  ProgramOptions () :
    BoundaryConditions(0),
    GasID(0),
    Debugging(0),
    TFinal(5000.0),  // Time in seconds
    PrintFlag(1),
    printSolnSteps(1),
    Temperature(298.15), // Temperature in Kelvin
    Pressure(101325.),  // Pressure in pascals
    VolInit(1.0E-6),    //1 cm**3 as an initial volume
    ipl(0),
    numIPL(0),
    InflowTemperature(300.),
    InflowMolarFlowRate(0.)    {
  }

  ~ProgramOptions() {
    if (ipl) {
      // Note the extra empty structure
      for (int i = 0; i < numIPL+1; i++) {
	delete ipl[i];
      }
      mdp_safe_free((void **) &ipl);
    }
  }

} PO;



/******************************************************************************
 *
 *
 */
void printUsage() {
  cout << "usage: basicMultiBlock [-h] [-help_cmdfile] cmdfile.txt"
       <<  endl;
  cout << "    -h           help" << endl;
  cout << "    -help_cmdfile Provides a listing of "
    "cmdfile.txt keylines" << endl;
  cout << "                  "
    "(with the cmdfile, it will print out species names)" 
       << endl;
}

/***********************************************************************
 *
 * generic function Wrapper around new, in order to create a function
 * pointer for BE_MultiBlock
 */
void *getNewIPInput(void *data_loc) {
  InputPipe *ptr = new InputPipe();
  return (void *)ptr;
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
  
  strcpy( GSList[0] , "H2");
  strcpy( GSList[1] , "N2");
  strcpy( GSList[2] , "O2");
  strcpy( GSList[3] , "CH4");

  PO.ipl = (InputPipe **) mdp_alloc_ptr_1(2);
  PO.ipl[0] = new InputPipe();
  InputPipe *ifl_ptr =  PO.ipl[0];
  /*
   *
   */
  BE_MultiBlock *bemd = new BE_MultiBlock("Inflow Pipe Definition",
					  &(PO.numIPL),
					  (void ***) &(PO.ipl),
					  getNewIPInput, 0, 0);

  topbe->addSubBlock(bemd);

  /* --------------------------------------------------------------
   *  Inflow Molar Flow Rate = [ double ]
   *     (required)
   *  Temperature of the inflow stream.
   *  This is used to establish the enthalpy of the input stream.
   */
  LE_OneDbl *ledimfr = new LE_OneDbl("Inflow Molar Flow Rate", 
				     &(ifl_ptr->InflowMolarFlowRate), 1,
				     "ifl_ptr->InflowMolarFlowRate");
  bemd->addLineEntry(ledimfr);


  /* --------------------------------------------------------------
   * Inflow Gas Mole Fractions
   *
   *   START BLOCK Inflow Gas Mole Fraction
   *      H2 = 0.1
   *      O2 = 0.54
   *      ...
   *   END BLOCK Inflow Gas Mole Fraction
   */

  BE_MoleComp *bigmf = 
    new BE_MoleComp("Inflow Gas Mole Fractions",
		    &(ifl_ptr->InflowGasMoleFractions), 1,
		    GSList, nSpecies, 1, 
		    "ifl_ptr->InflowGasMoleFractions");
  bigmf->generateDefLE();
  bemd->addSubBlock(bigmf);
 

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


/********************************************************************/
void read_info(ProgramOptions& P1, BlockEntry *topbe) {

  InputPipe *ipl_ptr = 0;
  int mc = 0;
  if (!P1.ipl) {
    P1.ipl = (InputPipe **) mdp_alloc_ptr_1(2);
    P1.numIPL = 0;
  }
  for ( ; ; ) {
    BlockEntry *mbbe = topbe->match_block("Inflow pipe Definition", mc);
    if (!mbbe) {
      break;
    }
    int numTimesProc = mbbe->get_NumTimesProcessed();
    if (numTimesProc <= 0) {
      break;
    }
    ipl_ptr = P1.ipl[mc];
    if (P1.numIPL < mc +1 && !ipl_ptr) {
      mdp_realloc_ptr_1((void ***) &(P1.ipl), P1.numIPL+2, P1.numIPL);
      P1.ipl[mc] = new InputPipe();
    }
    P1.numIPL = mc + 1;
    InputPipe *ipl_ptr = P1.ipl[mc];
 
   /***************************************/
    // OneDbl check
    LineEntry *le = mbbe->searchLineEntry("Inflow molar flow rate");
    const double *d2ptr = static_cast<const double *>(le->currentValueAsVoidP());
    ipl_ptr->InflowMolarFlowRate = *d2ptr;

    LE_OneDbl *le_du2 = dynamic_cast<LE_OneDbl *>(le);
    double dvalue2 = le_du2->currentTypedValue();
    if (dvalue2 != ipl_ptr->InflowMolarFlowRate) {
      printf ("we have an inconsistency: dvalue2 != ipl_ptr->InflowMolarFlowRate\n");
      exit(-1);
    }
    InputPipe *ipl0_ptr =  PO.ipl[mc];
    if (ipl0_ptr->InflowMolarFlowRate != ipl_ptr->InflowMolarFlowRate) {
      printf ("we have an inconsistency: ipl0_ptr->InflowMolarFlowRate "
	      "!= ipl_ptr->InflowMolarFlowRate \n");
      exit(-1);
    }


    /***********************************************/
    // Mole Frac check
    int nSpecies = 4;
    double *inflowGasMoleFractions = ipl_ptr->InflowGasMoleFractions;
 

    BlockEntry *be = mbbe->searchBlockEntry("Inflow Gas Mole Fractions");
    const double *moleFractions =
      static_cast<const double *>(be->currentValueAsVoidP());
    for (int i = 0; i < nSpecies; i++) {
      (inflowGasMoleFractions)[i] = moleFractions[i];
    }


    BE_MoleComp *be_dbl = dynamic_cast<BE_MoleComp *>(be);
    moleFractions = be_dbl->currentTypedValue();

    for (int i = 0; i < nSpecies; i++) {
      if ((inflowGasMoleFractions)[i] != moleFractions[i]) {
	printf ("we have an inconsistency: inflowGasMoleFractions != moleFractions \n");
	exit(-1);
      }
    }

    double *POinflowGasMoleFractions = ipl0_ptr->InflowGasMoleFractions;
    for (int i = 0; i < nSpecies; i++) {
      if ((inflowGasMoleFractions)[i] != (POinflowGasMoleFractions)[i]) {
	printf ("we have an inconsistency: "
		"inflowGasMoleFractions != POinflowGasMoleFractions \n");
	exit(-1);
      }
    }
    mc++;
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

    for (int num = 0; num < PO.numIPL; num++) {
      InputPipe *ip = PO.ipl[num];
      
      printf("Pipe Number %d\n", num);
      printf("\t Flow rate = %g\n", ip->InflowMolarFlowRate);
      for (int k = 0; k < 4; k++) {
	printf("\t\tMF[%d] = %g\n", k, ip->InflowGasMoleFractions[k]);
      }

    }


    for (int num = 0; num < P1_ptr->numIPL; num++) {
      InputPipe *ip = P1_ptr->ipl[num];
      
      printf("Pipe Number %d\n", num);
      printf("\t Flow rate = %g\n", ip->InflowMolarFlowRate);
      for (int k = 0; k < 4; k++) {
	printf("\t\tMF[%d] = %g\n", k, ip->InflowGasMoleFractions[k]);
      }

    }

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

    read_info(*P1_ptr, cf);

    delete cf;

    for (int num = 0; num < PO.numIPL; num++) {
      InputPipe *ip = PO.ipl[num];
      
      printf("Pipe Number %d\n", num);
      printf("\t Flow rate = %g\n", ip->InflowMolarFlowRate);
      for (int k = 0; k < 4; k++) {
	printf("\t\tMF[%d] = %g\n", k, ip->InflowGasMoleFractions[k]);
      }

    }


    for (int num = 0; num < P1_ptr->numIPL; num++) {
      InputPipe *ip = P1_ptr->ipl[num];
      
      printf("Pipe Number %d\n", num);
      printf("\t Flow rate = %g\n", ip->InflowMolarFlowRate);
      for (int k = 0; k < 4; k++) {
	printf("\t\tMF[%d] = %g\n", k, ip->InflowGasMoleFractions[k]);
      }

    }


    delete P1_ptr;

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
