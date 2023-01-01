/**
 *  @file example2.cpp
 *
 */


#include "IdealReactingGas.h"
                                  
#include "BlockInput.h"
#include "LE_PickList.h"
#include "BE_MoleComp.h"
#include "BE_UnitConversionPressure.h"
#include "LE_OneDblUnits.h"
#include "LE_OneStr.h"
#include "LE_OneBoolInt.h"
#include "LE_OneDbl.h"
#include "LE_OneInt.h"
#include "LE_StrListDbl.h"
#include "md_timer.h"
#include "mdp_allo.h"


#include "stdio.h"
#include <string>

//  Example 2
//
//  Read a mechanism, and print to the standard output stream a
//  well-formatted Chemkin ELEMENT section.
//

#include "mdp_allo.h"
using namespace Cantera;
using namespace std;
#include <iostream>

int Sint = 34;
int Debug = 0;
int tracking = 0;
double xmole[30];
double CBathMol[30];
double XBathMol[30];
char **CList = 0;
int BG_ID = 0;
string XMLInputFile;

/*****************************************************************/
void setup_lists()
{
  CList =  mdp_alloc_VecFixedStrings(30, 40);
  for (int i = 0; i < 30; i++) {
    sprintf(CList[i], "species_%d", i);
  }
}

/*****************************************************************/
void setup_input(BlockEntry *cf)
{


    LE_OneBoolInt *b1 = new LE_OneBoolInt("Debugging", &Debug);
    b1->set_default(0);
    cf->addLineEntry(b1);

    b1 = new LE_OneBoolInt("set tracking", &Debug, 1);
    b1->set_default(true);
    cf->addLineEntry(b1);


    LE_OneInt *a1 = new LE_OneInt("This is a one int entry", &Sint);
    a1->set_default(45);
    a1->set_limits(400, 0);
    a1->set_PrintString("GlobalIntVal");
    cf->addLineEntry(a1);

    LE_PickList *pl1 = new LE_PickList("Bath Gas ID", &BG_ID,
				       CList, 30, 1, "BG_ID");
    pl1->set_default(1);
    cf->addLineEntry(pl1); 

    LE_OneStr *s1 = new LE_OneStr("XML Input File", &XMLInputFile,
				  1, 1, 1);
    s1->set_default("chem.xml");
    cf->addLineEntry(s1);

    LE_StrListDbl *sl1 = new LE_StrListDbl("Xmole",
					   (double **) &(xmole) ,
					   CList, 30, 1, "Xmol");
    sl1->set_limits(1.0, 0.0);
    sl1->set_default(0.0);
    cf->addLineEntry(sl1);

    BE_StrDbl *besb = new BE_StrDbl("Bath Gas Composition",
				    (double **) &(CBathMol),
				    0, CList, 30, 0, "CBathmol");
    besb->set_default(0.0);
    besb->set_limits(1.0, 0.0);
    besb->generateDefLE();
    cf->addSubBlock(besb);



    BE_MoleComp * bmc = new BE_MoleComp("Bath Gas Mole Fraction",
					(double **) &(XBathMol), 1,
					CList, 30, 0, "XBathmol");
  bmc->generateDefLE();
  cf->addSubBlock(bmc);
}

/*****************************************************************
 *
 * process_input():
 *
 *  This routine catches error messages
 */
void process_input(BlockEntry *cf)
{
  const TOKEN tok_in;
  TOKEN tok_out;
  try {
    /*
     * Call the block read function at the main level
     */
    cf->read_block(stdin, &tok_out, &tok_in, 0);
  } catch (BI_InputError &bi) {
    /*
     * This catches error messages
     */
    cout << bi.errorMessage() << endl;
  }
}

/*****************************************************************/


int main(int argc, char** argv) { 

  BlockEntry *cf = new BlockEntry("command_file");
  
  setup_lists();
  setup_input(cf);

  cf->print_usage();

  process_input(cf);
  delete cf;

  cout << "input file = " << XMLInputFile << endl;
}
/***********************************************************/
