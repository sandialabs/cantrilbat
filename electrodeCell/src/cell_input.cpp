/*
 * $Id: cell_input.cpp 504 2013-01-07 22:32:48Z hkmoffa $
 */

#include "tok_input_util.h"



#include "cantera/equilibrium.h"
#include "cell_input.h"


#include "BlockEntryGlobal.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace Cantera;
using namespace std;
using namespace BEInput;
using namespace TKInput;

/*************************************************************************
 *
 */

CELL_KEY_INPUT CellO;

CELL_KEY_INPUT::CELL_KEY_INPUT () :
  anodeElectrodeFileName_("anode.inp"),
  cathodeElectrodeFileName_("cathode.inp"),
  doCathode_(false),
  netCurrent_(0.0),
  AnodeVoltage_(0.0),
  CathodeVoltage_(0.0),
  NetVoltage_(0.0),
  Title_("")
{
    Title_ = "cell Problem";
}
/****************************************************************************
 *
 */
CELL_KEY_INPUT::~CELL_KEY_INPUT () {
 
}
/****************************************************************************
 *
 */
void CELL_KEY_INPUT::InitForInput() {
 
}
   
/****************************************************************************
 *
 */
static void setup_input_pass1(BlockEntry *cf)
{

  /* ---------------------------------------------------------------
   *
   */
  LE_OneStr *sA = new LE_OneStr("Anode File Name",
				&CellO.anodeElectrodeFileName_,
				1, 1,  1, "anodeFileName");
  cf->addLineEntry(sA);
  
  /* ---------------------------------------------------------------
   *
   */
  LE_OneStr *sC = new LE_OneStr("Cathode File Name",
				&CellO.cathodeElectrodeFileName_,
				1, 1,  0, "cathodeFileName");
  cf->addLineEntry(sC);
  
  /* --------------------------------------------------------------
   * Current
   */
  LE_OneDbl *d1 = new LE_OneDbl("Net Current",
				&CellO.netCurrent_, 1, "netCurrent");
  d1->set_default(0.0);
  d1->set_limits(3000., -3000.);
  cf->addLineEntry(d1);

  /* --------------------------------------------------------------
   * Anode Voltage
   */
  LE_OneDbl *dAV = new LE_OneDbl("Anode Voltage",
				&CellO.AnodeVoltage_, 1, "anodeVoltage");
  dAV->set_default(0.0);
  dAV->set_limits(3000., -3000.);
  cf->addLineEntry(dAV);
  
  /* ---------------------------------------------------------------
   *
   */
  LE_OneStr *s2 = new LE_OneStr("Title",
				&CellO.Title_, 100000, 1,  0, "ElectrodeCellTitle");
  cf->addLineEntry(s2);


  BaseEntry::set_SkipUnknownEntries(false);
}


   

/******************************************************************************
 *  Read the input file
 *
 *  printFlag 0 = no output
 *            1 = error message output
 *            2 = output
 */
bool process_cell_input(BlockEntry *cf, string fileName, int printFlag) {
  static int pass = 0;
  pass++;
  cf->ZeroLineCount();
  const TOKEN tok_in;
  TOKEN tok_out;
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


/*****************************************************/
static void print_char(const char letter, const int num)
{
   int i;
   for (i = 0; i < num; i++) printf("%c", letter);
}

/*****************************************************************************
 *
 * cell_input():
 *
 *    INput for vcs_Cantera. This routine will combine a text input file
 *    with a cantera xml or cti file to create an equilibrium problem
 *    to solve.
 */
int
cell_input(std::string commandFile)
{
   int  retn = MPEQUIL_SUCCESS;

   
   printf("\n"); print_char('=', 80); printf("\n");
   print_char('=', 20); printf(" cell_input: START OF PROBLEM STATEMENT ");
   print_char('=', 21); printf("\n");
   print_char('=', 80); printf("\n\n");  

   int printBIProclevel = 2;

   /**
     * Initialize a block input structure for the command file
     */
    BlockEntry *cf = new BlockEntry("command_file");
    /**
     * Setup and process the input deck for first time.
     * -> Might have to print out the input and quit as well.
     */
    setup_input_pass1(cf);
    bool ok = process_cell_input(cf, commandFile, printBIProclevel);
    if (!ok) {
      return -1;
    }
   


    if (CellO.anodeElectrodeFileName_.size() == 0) {
      printf("Cantera anode file name must be specified\n");
      exit (-1);
    }


    LineEntry *le = cf->searchLineEntry("Cathode File Name");
    int cathodeSpecified = le->get_NumTimesProcessed();
    if (cathodeSpecified) {
      CellO.doCathode_ = true;
    }
     

   /*
    *   Printout of the input information
    */
   printf("\n"); print_char('-', 80); printf("\n");
   printf("             Information about problem\n");
   printf("  Anode   file name = %s\n", CellO.anodeElectrodeFileName_.c_str());
   printf("  Cathode file name = %s\n", CellO.cathodeElectrodeFileName_.c_str());
   printf("   docathode = %d\n", CellO.doCathode_);
   printf("  Net Current = %g\n", CellO.netCurrent_);
   printf("  Anode Voltage  = %g\n", CellO.AnodeVoltage_);
   printf("  Cathode Voltage = %g\n", CellO.CathodeVoltage_);
   printf("  Title = %s\n", CellO.Title_.c_str());

   print_char('=', 80); printf("\n\n");

   delete cf;
   return retn;
}
/**************************************************************************/
