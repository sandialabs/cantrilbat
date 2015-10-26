/**
 *  @file example2.cpp
 *
 *  $Id: cttInput.cpp 497 2013-01-07 21:17:04Z hkmoffa $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cttables.h"
#include "cttInput.h"


#include "LE_PickList.h"
#include "BE_MoleComp.h"
#include "BE_UnitConversionPressure.h"
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
#include "cantera/thermo/ThermoPhase.h"
#include "LE_MultiCStr.h"
#include "BE_MultiBlock.h"

#include "cantera/thermo/MolalityVPSSTP.h"

using namespace std;
using namespace Cantera;
using namespace BEInput;
using namespace mdpUtil;

IOoptions IOO;

/************************************************************************
 *
 *
 * generic function Wrapper around new, in order to create a function
 * pointer for BE_MultiBlock
 */
void *getNewEGRInput(void *data_loc) {
  void *ptr = new EGRInput();
  return ptr;
}

void *getNewERSSpec(void *data_loc) {
  void *ptr = new ERSSpec();
  return ptr;
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 *  This routine should work with g_kin_ptr = 0 and g_ptr = 0 for 
 *  the sake of documenting the input file.
 */
void setup_input_pass1(BlockEntry *cf, Kinetics *g_kin_ptr, ThermoPhase *g_ptr)
{
  int numSpecies;
  char **CList = 0;
  BaseEntry::set_SkipUnknownEntries(3);
  /*
   * Create a list of species to be used in the input file.
   */
  if (!g_ptr) {
    numSpecies = 3;
    CList = mdp_alloc_VecFixedStrings(numSpecies,
				      MAX_INPUT_STR_LN+1);
    strncpy(CList[0], "Species_0", MAX_INPUT_STR_LN);
    strncpy(CList[1], "Species_1", MAX_INPUT_STR_LN);
    strncpy(CList[2], ".........", MAX_INPUT_STR_LN);
  } else {
    numSpecies = g_ptr->nSpecies();
    CList = mdp_alloc_VecFixedStrings(numSpecies,
				      MAX_INPUT_STR_LN+1);
    const vector<string>& sn = g_ptr->speciesNames();

    for (int i = 0; i < numSpecies; i++) {
      strncpy(CList[i], sn[i].c_str(), MAX_INPUT_STR_LN);
    }
  }

  /* --------------------------------------------------------------
   * NumberCanteraFiles - 
   *
   * Configure the number of cantera files
   */
  LE_OneInt *in1 = new LE_OneInt("Number of Cantera Files",
				 &IOO.NumberCanteraFiles);
  in1->set_default(1);
  cf->addLineEntry(in1);

  /* --------------------------------------------------------------
   * DebugPrinting - 
   *
   * Configure Input of debugging variable. This can also
   * be done on the command line.
   */
  LE_OneBoolInt *b1 = new LE_OneBoolInt("DebugPrinting",
					&DebugPrinting);
  b1->set_default(DebugPrinting);
  cf->addLineEntry(b1);

  mdp_safe_free((void **) &CList);
}

void setup_input_pass2(BlockEntry *cf, Kinetics *g_kin_ptr, 
		       ThermoPhase *g_ptr) {

  BaseEntry::set_SkipUnknownEntries(3);
  /*
   *  Get the input deck for
   *  Cantera description of the model.
   */
  LE_MultiCStr *s1 =
    new LE_MultiCStr("Cantera File Name", &IOO.CanteraFileNames,
		     1, 1,  0, "CanteraFileNames");
  s1->set_default("gas.cti");

  /*
   * Set up a dependency on the input from the Number of Cantera
   * Files card
   */
  LineEntry *sle1 = cf->searchLineEntry("Number of Cantera Files");
  int numF = 1;
  (void) sle1->ansDepCheckOneInt(numF);
  s1->set_NumTimesRequired(numF);

  cf->addLineEntry(s1);

}

// Pass 3
/*
 * In this pass, we can assume that the PhaseList is fully formed.
 */
void setup_input_pass3(BlockEntry *cf,
		       Kinetics *g_kin_ptr, ThermoPhase *g_ptr,
		       PhaseList *pl)
{

  int iph;
  BaseEntry::set_SkipUnknownEntries(0);
    
  /* --------------------------------------------------------------
   * Output Units = ["Kcal_cgs", "KJoule"]
   *    (optional)
   *    default = Kcal_cgs
   *
   *    Select the units for output.
   */
  const char * cunits[2] = {"Kcal_cgs", "KJoule"};
  LE_PickList *lepunits =
    new LE_PickList("Output Units", &(IOO.OutputUnits),
		    cunits, 2, 0, "OutputUnits");
  lepunits->set_default(0);
  cf->addLineEntry(lepunits);

  /* --------------------------------------------------------------
   * Add Chemical Potential Column = [ true, false ]
   *
   * Specify whether you want an extra column containing
   * the raw value of the chemical potential.
   */
  LE_OneBool *bCPC =
    new LE_OneBool("Add Chemical Potential Column", 
		   &IOO.ChemPotColumn, 0, "ChemPotColumn");
  bCPC->set_default(false);
  cf->addLineEntry(bCPC);

  /* --------------------------------------------------------------
   * Add Internal Energy Column = [ true, false ]
   *
   * Specify whether you want an extra column containing
   * the relative value of the internal energy.
   */
  LE_OneBool *bIEC =
    new LE_OneBool("Add Internal Energy Column", 
		   &IOO.IntEnergyColumn, 0, "IntEnergyColumn");
  bIEC->set_default(false);
  cf->addLineEntry(bIEC);
  bIEC = 0;


  /* --------------------------------------------------------------
   * Skip Transport = [ true, false ]
   *
   * Specify whether you want to skip the inclusion of transport numbers into the thermodynamics
   * tables for each species. Default = false , i.e., the transport is included
   */
  LE_OneBool *bITR = new LE_OneBool("Skip Transport", &IOO.SkipTransport, 0, "SkipTransport");
  bITR->set_default(false);
  cf->addLineEntry(bITR);
  bITR = 0;
    
  /* --------------------------------------------------------------
   * BG.Temperature - 
   *
   * Configure Input of Bath Temperature
   */
  LE_OneDbl *b2 = new LE_OneDbl("Bath Temperature", 
				&(BG.Temperature), 0, 
				"BG.Temperature");
  b2->set_default(298.15);
  b2->set_limits(20000., 0.0);
  cf->addLineEntry(b2);
  b2 = 0;

  /* --------------------------------------------------------------
   * BG.Pressure - 
   * 
   * Configure Input of Bath Pressure
   */
  BE_UnitConversion *ucPres = new BE_UnitConversionPressure();
  LE_OneDblUnits *b3 = new LE_OneDblUnits("Bath Pressure",
					  &(BG.Pressure), 0, 
					  "BG.Pressure", ucPres);
  b3->set_default(101325.);
  b3->set_limits(1.E20, 0.0);
  cf->addLineEntry(b3);

  int nVolPhases = pl->nVolPhases();
  int nSurPhases = pl->nSurPhases();

  BG.Xmol = mdp_alloc_dbl_1(IOO.nTotSpecies + 2, 0.0);

  mdp_safe_free((void **) &BG.Molalities);
  BG.Molalities = mdp_alloc_dbl_1(IOO.nTotSpecies + 2, 0.0);

  BG.XmolPLSpecVec = mdp_alloc_dbl_1(IOO.nTotSpecies + 2, 0.0);
  BG.MolalitiesPLSpecVec = mdp_alloc_dbl_1(IOO.nTotSpecies + 2, 0.0);
  BG.XmolPLPhases = (double **) mdp_alloc_ptr_1(nVolPhases + nSurPhases);
  BG.MolalitiesPLPhases = (double **)mdp_alloc_ptr_1(nVolPhases + nSurPhases);
  BG.BathSpeciesIDVec = mdp_alloc_int_1(nVolPhases + nSurPhases, 0);
  BG.PotentialPLPhases =  mdp_alloc_dbl_1(nVolPhases + nSurPhases, 0.0);
  BG.PhaseMoles.resize(nVolPhases + nSurPhases, 0.0);

  for (iph = 0; iph < pl->nVolPhases(); iph++) {
    int kstart =  pl->getGlobalSpeciesIndexVolPhaseIndex(iph);
    BG.XmolPLPhases[iph] =   BG.XmolPLSpecVec + kstart;
    BG.MolalitiesPLPhases[iph] =  BG.MolalitiesPLSpecVec + kstart;
  }
  for (iph = 0; iph < pl->nSurPhases(); iph++) {
   int tph = iph + pl->nVolPhases();
   int kstart =  pl->getGlobalSpeciesIndexSurPhaseIndex(iph);
   BG.XmolPLPhases[tph] =   BG.XmolPLSpecVec + kstart;
   BG.MolalitiesPLPhases[tph] =  BG.MolalitiesPLSpecVec + kstart;
  }

  for (iph = 0; iph < pl->nVolPhases(); iph++) {
    string phaseBath = "Bath Specification for Phase ";
    ThermoPhase *tp = &(pl->volPhase(iph));
    string phaseNm = tp->name();
    int nSpecies = tp->nSpecies();
    phaseBath += phaseNm; 
    /*
     *  create a section method description block and start writing
     *  line elements in it.
     */
    BlockEntry *bbathphase = new BlockEntry(phaseBath.c_str());
    cf->addSubBlock(bbathphase);
    int kstart =  pl->getGlobalSpeciesIndexVolPhaseIndex(iph);

    /* --------------------------------------------------------------
     * BG.BathSpeciesID - 
     *
     * Create a PickList Line Element made out of the list of 
     * species
     */
    LE_PickList *plp1 =
      new LE_PickList("Bath Species ID", &(BG.BathSpeciesIDVec[iph]),
		      (const char **)IOO.SpeciesNames+kstart, nSpecies, 1, "BG_ID");
    plp1->set_default(0);
    bbathphase->addLineEntry(plp1);

    /* --------------------------------------------------------------
     * BG.PotentialPLPhases
     *  Input the bath voltage for the phase
     */
    LE_OneDbl *iVolt = 
      new LE_OneDbl("Voltage", &(BG.PotentialPLPhases[iph]), 0, 
		    "Voltage");
    iVolt->set_default(0.0);
    iVolt->set_limits(10., -10.);
    bbathphase->addLineEntry(iVolt);
    iVolt = 0;

    /* --------------------------------------------------------------
     * BG.PhaseMoles
     *  Input the number of moles for the phase
     */
    LE_OneDbl *iTMoles = 
      new LE_OneDbl("Phase Moles", &(BG.PhaseMoles[iph]), 0, 
		    "PhaseMoles");
    iTMoles->set_default(0.0);
    iTMoles->set_limits(1.0E9, 0.0);
    bbathphase->addLineEntry(iTMoles);

    /* --------------------------------------------------------------
     * BG.BathSpeciesMoleFractions - 
     *
     * Create a PickList Line Element made out of the list of 
     * species
     */
    // double *xstart = BG.XmolPLSpecVec + kstart;
    BE_MoleComp *bpmc = new BE_MoleComp("Bath Species Mole Fraction",
					&(BG.XmolPLPhases[iph]), 0,
					IOO.SpeciesNames+kstart, nSpecies, 0, "XBathmol");
    bpmc->generateDefLE();
    bbathphase->addSubBlock(bpmc);

    if (tp->activityConvention() == cAC_CONVENTION_MOLALITY) {
      MolalityVPSSTP * m_ptr = dynamic_cast<MolalityVPSSTP *>(tp);
      if (m_ptr == 0) {
	printf("Dynamic cast failed for some reason\n");
	exit(-1);
      }
      int indS = m_ptr->solventIndex();
      double mwS = m_ptr->molecularWeight(indS);
      BE_MolalityComp *bmolal = 
	new BE_MolalityComp("Bath Species Molalities",
			    &(BG.MolalitiesPLPhases[iph]), 0,
			    IOO.SpeciesNames+kstart, nSpecies, 
			    indS, mwS,
			    "MolalitiesBath");
      //bmolal->generateDefLE();
      bbathphase->addSubBlock(bmolal);

    }

  }

  /*
   * Surface Phases
   */

  for (int iphS = 0; iphS < pl->nSurPhases(); iphS++) {
    string phaseBath = "Bath Specification for Phase ";
    iph = pl->nVolPhases() + iphS;
    ThermoPhase *tp = &(pl->surPhase(iphS));
    string phaseNm = tp->name();
    int nSpecies = tp->nSpecies();
    phaseBath += phaseNm; 
    /*
     *  create a section method description block and start writing
     *  line elements in it.
     */
    BlockEntry *bbathphase = new BlockEntry(phaseBath.c_str());
    cf->addSubBlock(bbathphase);
    int kstart =  pl->getGlobalSpeciesIndexSurPhaseIndex(iphS);

    /* --------------------------------------------------------------
     * BG.BathSpeciesID - 
     *
     * Create a PickList Line Element made out of the list of 
     * species
     */
    LE_PickList *plp1 =
      new LE_PickList("Bath Species ID", &(BG.BathSpeciesIDVec[iph]),
		      (const char **)IOO.SpeciesNames+kstart, nSpecies, 1, "BG_ID");
    plp1->set_default(0);
    bbathphase->addLineEntry(plp1);


    /* --------------------------------------------------------------
     * BG.BathSpeciesMoleFractions - 
     *
     * Create a PickList Line Element made out of the list of 
     * species
     */
    BE_MoleComp *bpmc = new BE_MoleComp("Bath Species Mole Fraction",
					&(BG.XmolPLPhases[iph]), 0,
					IOO.SpeciesNames+kstart, nSpecies, 0, "XBathmol");
    bpmc->generateDefLE();
    bbathphase->addSubBlock(bpmc);

  }

 

  /* --------------------------------------------------------------
   * Use Reference Pressure in Thermo Tables
   *
   * Specify whether you want to use the reference pressure or 
   * the bath gas pressure in the Thermo Tables for each
   * species
   *       default = true
   *       optional
   */
  LE_OneBool *bRPTT =
    new LE_OneBool("Use Reference Pressure in Thermo Tables", 
		   &IOO.UseRefPressureInThermoTables,
		   0, "UseRefPressureInThermoTables");
  bRPTT->set_default(true);
  cf->addLineEntry(bRPTT);
    
  /** --------------------------------------------------------------
   *
   * Subblock for specification of the temperature table format
   */
  BlockEntry* beTT = new BlockEntry("Temperature Table Format", 0);
  cf->addSubBlock(beTT);

  /* ---------------------------------------------------------------
   * IOO.m_TTnpts - 
   *
   *  Sets the number of points in the temperature tables
   */
  LE_OneInt *iTTnpts = 
    new LE_OneInt("Number of Points", &(IOO.m_TTnpts), 0, "TTpts");
  iTTnpts->set_default(14);
  iTTnpts->set_limits(10000, 2);
  beTT->addLineEntry(iTTnpts);

 
  /* ---------------------------------------------------------------
   * IOO.m_TTDeltaT - 
   *
   *  Sets the deltaT between points in the temperature table.
   */
  LE_OneDbl *iTTDeltaT = 
    new LE_OneDbl("Delta Temperature", &(IOO.m_TTDeltaT), 0, 
		  "TTDeltaT");
  iTTDeltaT->set_default(100.);
  iTTDeltaT->set_limits(10000., 1.0E-10);
  beTT->addLineEntry(iTTDeltaT);


 
  /* ---------------------------------------------------------------
   * IOO.m_TTTlow - 
   *
   *  Sets the low temperature in the temperature table.
   */
  LE_OneDbl *iTTTlow = 
    new LE_OneDbl("Low Temperature", &(IOO.m_TTTlow), 0, 
		  "TTTlow");
  iTTTlow->set_default(300.);
  iTTTlow->set_limits(10000., 1.0E-10);
  beTT->addLineEntry(iTTTlow);
 
  /* ----------------------------------------------------------------
   *  Input a vector of doubles on a single line
   *  LE_VecDbl(const char *lineName, double **hndlAAA, int vecLength,
   *	      int numRL = 0, const char *varName = 0);
   *
   */
  LE_VecDblVarLength *v1 =
    new LE_VecDblVarLength("Added Temperatures",
			   &(IOO.AddedTemperatures),
			   1, 0,"AddedTemperatures");
  beTT->addLineEntry(v1);

  /* ------------------------------------------------------------------*/

   /** --------------------------------------------------------------
   *
   * Subblock for specification of the voltage table format
   */
  BlockEntry* beVV = new BlockEntry("Voltage Table Format", 0);
  cf->addSubBlock(beVV);

  /* ---------------------------------------------------------------
   * IOO.m_VVnpts - 
   *
   *  Sets the number of points in the voltage tables
   */
  LE_OneInt *iVVnpts = 
    new LE_OneInt("Number of Points", &(IOO.m_VVnpts), 0, "VVpts");
  iVVnpts->set_default(11);
  iVVnpts->set_limits(10000, 2);
  beVV->addLineEntry(iVVnpts);
  iVVnpts = 0;
 
  /* ---------------------------------------------------------------
   * IOO.m_VVDeltaV - 
   *
   *  Sets the deltaV between points in the voltage table.
   */
  LE_OneDbl *iVVDeltaV = 
    new LE_OneDbl("Delta Voltage", &(IOO.m_VVDeltaV), 0, 
		  "VVDeltaV");
  iVVDeltaV->set_default(0.2);
  iVVDeltaV->set_limits(100., 1.0E-10);
  beVV->addLineEntry(iVVDeltaV);
 
  /* ---------------------------------------------------------------
   * IOO.m_VVVlow - 
   *
   *  Sets the low voltage in the voltage table.
   */
  LE_OneDbl *iVVVlow = 
    new LE_OneDbl("Low Voltage", &(IOO.m_VVVlow), 0, 
		  "VVVlow");
  iVVVlow->set_default(-1.0);
  iVVVlow->set_limits(100., -100.);
  beVV->addLineEntry(iVVVlow);
  iVVVlow = 0;

  /* ---------------------------------------------------------------
   * IOO.VVincZero  
   *
   *  Sets whether the V = 0 point should be added to the table
   */
  LE_OneBool *iVVb1 = 
    new LE_OneBool("Include Zero Voltage", &(IOO.VVincZero), 0, "VVincZero");
  iVVb1->set_default(false);
  beVV->addLineEntry(iVVb1);
  iVVb1 = 0;

 /* ---------------------------------------------------------------
   * IOO.VVincEzero  
   *
   *  Sets whether the E_naught = 0 point should be added to the table
   */
  LE_OneBool *iVVb2 = 
    new LE_OneBool("Include Ezero Voltage", &(IOO.VVincEzero), 0, "VVincEzero");
  iVVb2->set_default(false);
  beVV->addLineEntry(iVVb2);
  iVVb2 = 0;

  /* ---------------------------------------------------------------
   * IOO.VVincEeq  
   *
   *  Sets whether the E_eq = 0 point should be added to the table
   */
  LE_OneBool *iVVb3 = 
    new LE_OneBool("Include Eeq Voltage", &(IOO.VVincEeq), 0, "VVincEeq");
  iVVb3->set_default(true);
  beVV->addLineEntry(iVVb3);
  iVVb3 = 0;
 
  /* ----------------------------------------------------------------
   *  Input a vector of doubles on a single line
   *  LE_VecDbl(const char *lineName, double **hndlAAA, int vecLength,
   *	      int numRL = 0, const char *varName = 0);
   *
   */
  LE_VecDblVarLength *vv =
    new LE_VecDblVarLength("Added Voltages",
			   &(IOO.AddedVoltages),
			   1, 0,"AddedVoltages");
  beVV->addLineEntry(vv);
  vv = 0;
  
  /* --------------------------------------------------------------- */

  /*
   * Delete the lists created in this section. Note the Line
   * Element Objects's ownership has been transfered to the
   * BlockInput object. They will be deleted when that object
   * is deleted.
   */
  IOO.m_EGRList = (struct EGRInput **) mdp_alloc_ptr_1(2);
  IOO.m_EGRList[0] = new EGRInput();

  BlockEntry *sbEGR = new BE_MultiBlock("Extra Global Reaction",
				      &(IOO.numExtraGlobalRxns),
				      (void ***) &(IOO.m_EGRList),
				      getNewEGRInput, (void *) 0,
				      0);
  cf->addSubBlock(sbEGR);

  /*
   * OK now define items in the block
   */
  struct EGRInput *egr_ptr = IOO.m_EGRList[0];
  /* --------------------------------------------------------------
   * EGR -> Special Species
   *
   * Create a PickList Line Element made out of the list of 
   * species
   */
  LE_PickList *egrSS =
    new LE_PickList("Special Species", &(egr_ptr->m_SS_KinSpeciesKindex),
		    (const char **)IOO.SpeciesNames, IOO.nTotSpecies, 1, "Special_Species");
  egrSS->set_default(0);
  sbEGR->addLineEntry(egrSS);
  

  egr_ptr->m_ERSList =(struct ERSSpec **) mdp_alloc_ptr_1(2);
  egr_ptr->m_ERSList[0] = new ERSSpec();

  BlockEntry *sbERS = new BE_MultiBlock("Elementary Reaction Specification",
					&(egr_ptr->m_numElemReactions),
					(void ***) &(egr_ptr->m_ERSList),
					getNewERSSpec, (void *) 0,
					0);
  sbEGR->addSubBlock(sbERS);

  /*
   * OK now define items in the block
   */
  struct ERSSpec *ers_ptr = egr_ptr->m_ERSList[0];
  

  /* ---------------------------------------------------------------
   * 
   *
   *  Sets the number of points in the temperature tables
   */
  LE_OneInt *iRxnIndex = 
    new LE_OneInt("Reaction Index", &(ers_ptr->m_reactionIndex), 1, "ReactionIndex");
  iRxnIndex->set_default(-1);
  sbERS->addLineEntry(iRxnIndex);

  LE_OneDbl *iRxnMult = 
    new LE_OneDbl("Reaction Multiplier", &(ers_ptr->m_reactionMultiplier), 1, 
		  "ReactionMultiplier");
  iRxnMult->set_default(0.0);
  sbERS->addLineEntry(iRxnMult);


}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 *  Read the input file
 */
int process_input(BlockEntry *cf,
		  string commandFile,
		  Kinetics *gKinetics,
		  ThermoPhase *gThermo,
		  PhaseList *pl)
{
  const TK_TOKEN tok_in;
  TK_TOKEN tok_out;
  static int pass = 0;
  pass++;
  IOO.reprep();
  cf->ZeroLineCount();
  FILE *ifp = fopen(commandFile.c_str(), "r");
  if (!ifp) {
    cout << "ERROR can't open file " << commandFile<< endl;
    return false;
  }
 
  try {
    /*
     * Set the print level by setting the static member
     * of the LineEntry baseclass.
     */
    BaseEntry::set_printProcessedLine(DebugPrinting);

    /*
     * Set whether the input file is mimed back to stdout
     */
    TKInput::set_tok_input_print_flag(DebugPrinting);
    /*
     * Call the block read function at the main level
     */
    cf->read_block(ifp, &tok_out, &tok_in, 0);
  } catch (BI_InputError &bi) {
    /*
     * This catches error messages
     */
    cout << bi.errorMessage() << endl;
    exit(-1);
  }
  fclose(ifp);

  
  if (pass >= 3) {

    LineEntry *v1 =
      cf->searchLineEntry("Added Temperatures");
    IOO.NumAddedTemps = v1->get_NumTimesProcessed();
  }
  /*
   * Setup the strings and the modifiers for the output
   */
  UIO.setup(IOO.OutputUnits);

  return 1;
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
/*
 *
 */
void IOoptions::InitForInput(PhaseList *pl) {
    nTotPhases = pl->nPhases();
    nTotSpecies = pl->nSpecies();
    nTotElements = pl->nElements();

    /*
     * Include all Phases by default
     */
    PhaseInclude = mdp_alloc_int_1(nTotPhases, 1);
    MoleNumber   = mdp_alloc_dbl_1(nTotSpecies, 0.0);
    ElementAbundances = mdp_alloc_dbl_1(nTotElements, 0.0);


    SpeciesNames = mdp_alloc_VecFixedStrings(nTotSpecies,
                                             132);
    PhaseNames = mdp_alloc_VecFixedStrings(nTotPhases,
                                           132);
    int kT = 0;
    for (int iphase = 0; iphase < nTotPhases; iphase++) {
      ThermoPhase *tPhase = &(pl->thermo(iphase));
      string id = tPhase->id();
      strncpy(PhaseNames[iphase], id.c_str(), 132);
      int nspecies = tPhase->nSpecies();
      for (int k = 0; k < nspecies; k++) {
        string sname = tPhase->speciesName(k);
        strncpy(SpeciesNames[kT], sname.c_str(), 132);
        kT++;
      }
    }

    /*
     * Create a list of elements for matching wrt the input.
     */
    ElementNames = mdp_alloc_VecFixedStrings(nTotElements,
                                             132);
    const Elements *eObj = pl->getGlobalElements();
    for (int e = 0; e < nTotElements; e++) {
      string eName = eObj->elementName(e);
      strncpy(ElementNames[e], eName.c_str(), 132);
    }
}

IOoptions::IOoptions() :
  NumberCanteraFiles(1),
  CanteraFileNames(0),
  ProcessAll(1),
  PrintThermoTable(0),
  OutputUnits(UNITS_KCAL_CGS),
  ChemPotColumn(false),
  IntEnergyColumn(false),
  SkipTransport(false),
  m_TTnpts(14),
  m_TTDeltaT(100),
  TTinc298(true),
  m_TTTlow(300.),
  NumAddedTemps(0),
  AddedTemperatures(0),
  m_VVnpts(11),
  m_VVDeltaV(0.2),
  VVincZero(false),
  VVincEzero(false),
  VVincEeq(true),
  m_VVVlow(-1.0),
  NumAddedVoltages(0),
  AddedVoltages(0),
  UseRefPressureInThermoTables(true),
  nTotPhases(0),
  nTotSpecies(0),
  nTotElements(0),
  numExtraGlobalRxns(0),
  m_EGRList(0)
{
}
/*
 * Destructor():
 */
IOoptions::~IOoptions() {
  if (CanteraFileNames) {
    for (int i = 0; i <  NumberCanteraFiles; i++) {
      mdp_safe_free((void **) &(CanteraFileNames[i]));
    }
    mdp_safe_free((void **) &CanteraFileNames);
  }
  mdp_safe_free((void **) &PrintThermoTable);
  mdp_safe_free((void **) &PhaseInclude);
  mdp_safe_free((void **) &MoleNumber);
  mdp_safe_free((void **) &PhaseNames);
  mdp_safe_free((void **) &ElementNames);
  mdp_safe_free((void **) &ElementAbundances);
  mdp_safe_free((void **) &SpeciesNames);
  mdp_safe_free((void **) &AddedTemperatures);

  struct EGRInput **ptr;
  if (m_EGRList != 0) {
    for (ptr = m_EGRList; *ptr != 0; ptr++) {
      delete *ptr;
    }
    mdp_safe_free((void **) &m_EGRList);
  }
}

void IOoptions::reprep() {
  if (CanteraFileNames) {
    for (int i = 0; i <  NumberCanteraFiles; i++) {
      if (CanteraFileNames[i]) {
	mdp_safe_free((void **) &(CanteraFileNames[i]));
      } else {
	break;
      }
    }
    mdp_safe_free((void **) &CanteraFileNames);
  }
}
