/**
 *  @file example2.cpp
 *
 *  $Id: chInput.cpp 508 2013-01-07 22:54:04Z hkmoffa $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "chInput.h"
#include "changeHf.h"



#include "LE_PickList.h"
#include "BE_MoleComp.h"
#include "BE_UnitConversionEnergy.h"
#include "LE_OneDblUnits.h"
#include "LE_OneStr.h"
#include "LE_OneBool.h"
#include "LE_OneBoolInt.h"
#include "LE_OneDbl.h"
#include "LE_OneInt.h"
#include "LE_VecDbl.h"
#include "LE_VecDblVarLength.h"
#include "md_timer.h"
#include "mdp_allo.h"
#include "BE_MolalityComp.h"
#include "zuzax/thermo/ThermoPhase.h"

//#include "zuzax/kernel/thermo/MolalityVPSSTP.h"

using namespace BEInput;
using namespace TKInput;

IOoptions IOO;
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 *  This routine should work with g_ptr = 0 for 
 *  the sake of documenting the input file.
 */
void setup_input(BlockEntry *cf)
{

    /* ---------------------------------------------------------------
     * IOO.FileName - 
     *
     *  Sets the name of the file to be modified
     */
    LE_OneStr *fN = new LE_OneStr("File Name", &(IOO.FileName),
				  1, 1, 1, "FileName");
    cf->addLineEntry(fN);

    /* ---------------------------------------------------------------
     * IOO.FileName - 
     *
     *  Sets the name of the file to be modified
     */
    LE_OneStr *dfN = new LE_OneStr("Destination File Name",
				   &(IOO.DestFileName),
				   1, 1, 0, "DestFileName");
    cf->addLineEntry(dfN);

    /* --------------------------------------------------------------
     * Energy Units = ["Kcal_cgs", "KJoule"]
     *    (optional)
     *    default = Kcal_cgs
     *
     *    Select the units for specification of the enthalpy
     */
    const char *cunits[2] = {"Kcal_cgs", "KJoule"};
    LE_PickList *lepunits =
	new LE_PickList("Energy Units", &(IOO.OutputUnits),
			cunits, 2, 0, "EnergyUnits");
    lepunits->set_default(0);
    cf->addLineEntry(lepunits);

    /* ---------------------------------------------------------------
     * IOO.SpeciesName - 
     *
     *  Sets the name of the species to be modified
     */
    LE_OneStr *sN = new LE_OneStr("Species Name", &(IOO.SpeciesName),
				  1, 1, 1, "SpeciesName");
    cf->addLineEntry(sN);

    /* ---------------------------------------------------------------
     * IOO.DeltaVal - 
     *
     *  Sets the delta value of the Hf
     */
    BE_UnitConversionEnergy *ucEn = new BE_UnitConversionEnergy();
    LE_OneDblUnits *dDV = 
	new LE_OneDblUnits("Delta Hf", &(IOO.DeltaValue), 0, 
			   "IOO.DeltaVal", ucEn);
    dDV->set_default(0.);
    dDV->set_limits(1.0E10, -1.0E10);
    cf->addLineEntry(dDV);

 
    /* ---------------------------------------------------------------
     * IOO.HFValue - 
     *
     *  Sets the absolute value of the Hf
     */
   BE_UnitConversionEnergy *ucEn2 = new BE_UnitConversionEnergy();
    LE_OneDblUnits *dHV = 
	new LE_OneDblUnits("Absolute Hf", &(IOO.HfValue), 0, 
		      "IOO.HfValue", ucEn2);
    dHV->set_default(HFVALUE_UNSET);
    dHV->set_limits(1.0E10, -1.0E10);
    cf->addLineEntry(dHV);
 
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 *  Read the input file
 */
void process_input(BlockEntry *cf, FILE *cmdptr)
{
    const TOKEN tok_in;
    TOKEN tok_out;
    try {
      /*
       * Set the print level by calling a static function
       * of the BaseEntry baseclass.
       */
      BaseEntry::set_printProcessedLine(DebugPrinting);

      /*
       * Set whether the input file is mimed back to stdout
       */
      set_tok_input_print_flag(DebugPrinting);
      /*
       * Call the block read function at the main level
       */
      cf->read_block(cmdptr, &tok_out, &tok_in, 0);
    } catch (BI_InputError &bi) {
      /*
       * This catches error messages
       */
      cout << bi.errorMessage() << endl;
      exit(-1);
    }

 
    /*
     * Setup the strings and the modifiers for the output
     */
    UIO.setup(IOO.OutputUnits);

    if (IOO.DestFileName == "") {
      IOO.DestFileName = IOO.FileName;
    }

    /*
     * Check requirements on the specification of the change
     * in the Heat of formation
     */
    if (IOO.DeltaValue != 0.0 && 
	(IOO.HfValue != HFVALUE_UNSET)) {
      printf("changeHf ERROR: Both Delta Hf and Absolute Hf were set\n");
      exit(-1);
    }
    
}
/*****************************************************************/


bool doubleEqual(double d1, double d2, double atol) {
    double sum = fabs(d1) + fabs(d2); 
    if (sum > 1.0E-275) {
      double denom = sum + fabs(atol);
      if (fabs(d1-d2) < denom * 1.0E-10) return true;
      return false;
    } 
    if (fabs(d1-d2) < atol) return true;
    return false;
}
/*****************************************************************/
/*****************************************************************/
