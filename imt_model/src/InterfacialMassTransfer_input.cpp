/*
 * $Id: InterfacialMassTransfer_input.cpp 507 2013-01-07 22:48:29Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"



#include "cantera/equilibrium.h"
#include "cantera/thermo/MolalityVPSSTP.h"

#include "PhaseList.h"

#include "BlockEntryGlobal.h"


#include "importAllCTML.h"
#include "ApplBase_print.h"
#include "InterfacialMassTransfer_input.h"
#include "InterfacialMassTransfer.h"
#include "InterfacialMassTransfer_Factory.h"

#include "importPL.h"

using namespace Cantera;
using namespace std;
using namespace BEInput;
using namespace TKInput;
using namespace mdpUtil;
using namespace ca_ab;

#include <string>

namespace Cantera {
/*************************************************************************
 *
 * MPEQUIL_KEY_INPUT(): constructor
 */
IMT_KEY_INPUT::IMT_KEY_INPUT () :
  CanteraFNSurface(""),
  NumberCanteraFiles(1),
  CanteraFileNames(0),
  Temperature(300.),
  PressureA(OneAtm),
  m_BG(0),
  MoleNumber(0),
  MoleFraction(0),
  PhaseInclude(0),
  ProblemType(TP),
  SpeciesNames(0),
  PhaseNames(0),
  ElementNames(0),
  ElementAbundances(0),
  specifiedElementAbundances(false),
  specifiedBlockKmolSpecies(0),
  PhaseAName(""),
  PhaseBName(""),
  solnAIndex_(-1),
  solnBIndex_(-1),
  Title(""),
  SurfaceArea(1.0),
  BLThickness_A(1.0),
  BLThickness_B(1.0),
  nTotPhases(0),
  nTotSpecies(0),
  nTotElements(0),
  m_pl(0)
{
  m_BG = new ElectrodeBath();
  m_pl = new PhaseList();
}
//======================================================================================================================
  IMT_KEY_INPUT::IMT_KEY_INPUT (const IMT_KEY_INPUT & right) :
  CanteraFNSurface(""),
  NumberCanteraFiles(1),
  CanteraFileNames(0),
  Temperature(300.),
  PressureA(OneAtm),
  m_BG(0),
  MoleNumber(0),
  MoleFraction(0),
  PhaseInclude(0),
  ProblemType(TP),
  SpeciesNames(0),
  PhaseNames(0),
  ElementNames(0),
  ElementAbundances(0),
  specifiedElementAbundances(false),
  specifiedBlockKmolSpecies(0),
  PhaseAName(""),
  PhaseBName(""),
  solnAIndex_(-1),
  solnBIndex_(-1),
  Title(""),
  SurfaceArea(1.0),
  BLThickness_A(1.0),
  BLThickness_B(1.0),
  nTotPhases(0),
  nTotSpecies(0),
  nTotElements(0),
  m_pl(0)
  {
    operator=(right);
  }
//======================================================================================================================

IMT_KEY_INPUT & IMT_KEY_INPUT::operator=(const IMT_KEY_INPUT & right)
{
  /*
   * Check for self assignment.
   */
  if (this == &right) return *this;
  
  NumberCanteraFiles = right.NumberCanteraFiles;

  if (CanteraFileNames) {
    for (int i = 0; CanteraFileNames[i] != 0; i++) {
      mdp_safe_free((void **) &(CanteraFileNames[i]));
    }
    mdp_safe_free((void **) &(CanteraFileNames));
  }
  if (right.CanteraFileNames) {
    CanteraFileNames = (char **) mdp_alloc_ptr_1(right.NumberCanteraFiles);
    
    for (int i = 0; i < right.NumberCanteraFiles; i++) {
      CanteraFileNames[i] = mdp_copy_string(CanteraFileNames[i]);
    }

  }

  Temperature = right.Temperature;
  PressureA = right.PressureA;

  if (m_BG) {
    delete m_BG;
  }
  m_BG = new ElectrodeBath(*(right.m_BG));

  nTotPhases = right.nTotPhases;
  nTotSpecies = right.nTotSpecies;
  nTotElements = right.nTotElements;

  delete MoleFraction;
  MoleFraction = mdp_alloc_dbl_1(nTotSpecies, 0.0);
  mdpUtil::mdp_copy_dbl_1(MoleFraction , right.MoleFraction, nTotSpecies);

  delete PhaseInclude;
  PhaseInclude = mdp_alloc_int_1(nTotPhases, 1);
  mdpUtil::mdp_copy_int_1(PhaseInclude, right.PhaseInclude, nTotPhases);

  delete MoleNumber;
  MoleNumber = mdp_alloc_dbl_1(nTotSpecies, 0.0);
  mdpUtil::mdp_copy_dbl_1(MoleNumber, right.MoleNumber, nTotSpecies);

  delete ElementAbundances;
  ElementAbundances = mdp_alloc_dbl_1(nTotElements, 0.0);
  mdpUtil::mdp_copy_dbl_1(ElementAbundances, right.ElementAbundances, nTotElements);

  delete SpeciesNames;
  SpeciesNames = mdp_alloc_VecFixedStrings(nTotSpecies,  MPEQUIL_MAX_NAME_LEN_P1);
  mdpUtil::mdp_copy_VecFixedStrings(SpeciesNames, (const char **) right.SpeciesNames, nTotSpecies, MPEQUIL_MAX_NAME_LEN_P1);

  delete PhaseNames;
  PhaseNames = mdp_alloc_VecFixedStrings(nTotPhases, MPEQUIL_MAX_NAME_LEN_P1);
  mdpUtil::mdp_copy_VecFixedStrings(PhaseNames, (const char **) right.PhaseNames, nTotPhases,MPEQUIL_MAX_NAME_LEN_P1);

  ProblemType = right.ProblemType;
  specifiedElementAbundances = right.specifiedElementAbundances;
  specifiedBlockKmolSpecies = right.specifiedBlockKmolSpecies;
  PhaseAName              = right.PhaseAName;
  PhaseBName              = right.PhaseBName;
  solnAIndex_             = right.solnAIndex_;
  solnBIndex_             = right.solnBIndex_;
  Title                   = right.Title;
  SurfaceArea             = right.SurfaceArea;
  BLThickness_A           = right.BLThickness_A;
  BLThickness_B           = right.BLThickness_B;
  nTotPhases              = right.nTotPhases;
  nTotSpecies             = right.nTotSpecies;
  nTotElements            = right.nTotElements;

  /*
   * Return the reference to the current object
   */
  return *this;
}
//======================================================================================================================
/****************************************************************************
 *
 */
IMT_KEY_INPUT::~IMT_KEY_INPUT () {
  if (CanteraFileNames) {
    for (int i = 0; CanteraFileNames[i] != 0; i++) {
      mdp_safe_free((void **) &(CanteraFileNames[i]));
    }
    mdp_safe_free((void **) &(CanteraFileNames));
  }
  delete(m_BG); m_BG=0;

  mdp_safe_free((void **) &PhaseInclude);
  mdp_safe_free((void **) &MoleNumber);
  mdp_safe_free((void **) &MoleFraction);
  mdp_safe_free((void **) &SpeciesNames);
  mdp_safe_free((void **) &PhaseNames);
  mdp_safe_free((void **) &ElementNames);
  mdp_safe_free((void **) &ElementAbundances);


  delete m_pl; m_pl = 0;
}

/****************************************************************************
 *
 *  In this routine, we fill in the fields in the IMT_KEY_INPUT 
 *  structure that will be needed to parse the keyword input.
 *  We also size the vectors in the structure appropriately, now
 *  we know the extent of the problem.
 *  The information to do this has already been gathered into the PhaseList 
 *  structure
 *
 *  These are:
 *         nTotPhases 
 *         nTotSpecies
 *         nTotElements
 *         SpeciesNames
 *         PhaseNames
 *         ElementNames 
 *         CanteraFNSurface
 *
 *  The sized fields include:
 *         PhaseInclude
 *         MoleNumber 
 *         MoleNumberIG 
 *         ElementAbundances 
 */
void IMT_KEY_INPUT::InitForInput(const Cantera::PhaseList * const pl) {
  nTotPhases  = pl->nPhases();
  nTotSpecies = pl->nSpecies();
  nTotElements = pl->nElements();
    
  /*
   * Include all Phases by default
   */
  PhaseInclude = mdp_alloc_int_1(nTotPhases, 1);
  MoleNumber   = mdp_alloc_dbl_1(nTotSpecies, 0.0);
  MoleFraction = mdp_alloc_dbl_1(nTotSpecies, 0.0);

  ElementAbundances = mdp_alloc_dbl_1(nTotElements, 0.0);
    

  SpeciesNames = mdp_alloc_VecFixedStrings(nTotSpecies,
					   MPEQUIL_MAX_NAME_LEN_P1);
  PhaseNames = mdp_alloc_VecFixedStrings(nTotPhases,
					 MPEQUIL_MAX_NAME_LEN_P1);
  int kT = 0;
  for (int iphase = 0; iphase < nTotPhases; iphase++) {
    ThermoPhase *tPhase = &(pl->thermo(iphase));

    string id = tPhase->id();
    strncpy(PhaseNames[iphase], id.c_str(), MPEQUIL_MAX_NAME_LEN);
    int nspecies = tPhase->nSpecies();
    tPhase->getMoleFractions(MoleFraction + kT);
    for (int k = 0; k < nspecies; k++) {
      string sname = tPhase->speciesName(k);
      strncpy(SpeciesNames[kT], sname.c_str(), MPEQUIL_MAX_NAME_LEN);
      kT++;
    }
 
  }

  if (pl->nSurPhases() > 0) {
    CanteraFNSurface = pl->firstSurfaceFile();
  }

  ElementNames = mdp_alloc_VecFixedStrings(nTotElements,
					   MPEQUIL_MAX_NAME_LEN_P1);
  const Elements *eObj = pl->getGlobalElements();

  for (int e = 0; e < nTotElements; e++) {
    std::string eName = eObj->elementName(e);
    strncpy(ElementNames[e], eName.c_str(), MPEQUIL_MAX_NAME_LEN);
  }
}


/************************************************************************
 *
 */
static void setup_input_pass1(BlockEntry *cf, IMT_KEY_INPUT *ei)
{
  /* ---------------------------------------------------------------
   *
   */
  LE_OneStr *mn = new LE_OneStr("Interfacial Mass Transport Model Name",
				&(ei->IMT_ModelName), 100000, 1,  1, "IMTModelName");
  cf->addLineEntry(mn);

  /* ---------------------------------------------------------------
   *
   * Obtain the number of cantera files to be read
   */
  LE_OneInt *s1 = new LE_OneInt("Number of Cantera Files",
				&(ei->NumberCanteraFiles), 0,
				"NumCanteraFiles");
  s1->set_default(1);
  cf->addLineEntry(s1);
  BaseEntry::set_SkipUnknownEntries(true);

  /* ---------------------------------------------------------------
   *
   */
  LE_OneStr *sA = new LE_OneStr("Phase A Name",
				&(ei->PhaseAName), 100000, 1,  1, "PhaseAName");
  cf->addLineEntry(sA);

  /* ---------------------------------------------------------------
   *
   */
  LE_OneStr *sB = new LE_OneStr("Phase B Name",
				&(ei->PhaseBName), 100000, 1,  1, "PhaseBName");
  cf->addLineEntry(sB);

  /* --------------------------------------------------------------
   * Gross Surface Area
   */
  LE_OneDbl *ds = new LE_OneDbl("Gross Surface Area of Interface",
				&(ei->SurfaceArea), 1, "surfaceArea");
  ds->set_default(0.01);
  ds->set_limits(3000., 1.0E-15);
  cf->addLineEntry(ds);


  /* --------------------------------------------------------------
   * Phase A Boundary Layer Thickness
   */
  LE_OneDbl *dbl_a = new LE_OneDbl("Phase A Boundary Layer Thickness",
				&(ei->BLThickness_A), 1, "blthick_solnA");
  dbl_a->set_default(0.01);
  dbl_a->set_limits(3000., 1.0E-15);
  cf->addLineEntry(dbl_a);


  /* --------------------------------------------------------------
   * Phase B Boundary Layer Thickness
   */
  LE_OneDbl *dbl_b = new LE_OneDbl("Phase B Boundary Layer Thickness",
				&(ei->BLThickness_B), 1, "blthick_solnB");
  dbl_b->set_default(0.01);
  dbl_b->set_limits(3000., 1.0E-15);
  cf->addLineEntry(dbl_b);

}
//======================================================================================================================
/****************************************************************************
 *
 */
static void setup_input_pass2(BlockEntry *cf, IMT_KEY_INPUT *ei)
{
  LineEntry *sle1 = 0;
  /*
   *  Get the input deck for
   *  Cantera description of the model.
   */
  LE_MultiCStr *s1 =
    new LE_MultiCStr("Cantera File Name", &(ei->CanteraFileNames),
		     1, 1,  0, "CanteraFileNames");
  s1->set_default("gas.cti");

  /*
   * Set up a dependency on the input from the Number of cantera
   * Files card
   */
  sle1 = cf->searchLineEntry("Number of Cantera Files");
  int numF = 1; 
  (void) sle1->ansDepCheckOneInt(numF);

  s1->set_NumTimesRequired(numF);
 
    
  cf->addLineEntry(s1);
  BaseEntry::set_SkipUnknownEntries(true);
}


//======================================================================================================================             
/****************************************************************************
 *
 */
static void setup_input_pass3(BlockEntry *cf, 
			      IMT_KEY_INPUT *ei)
{

  PhaseList *pl = ei->m_pl;
  int iph;
  /* ---------------------------------------------------------------
   *
   */
  LE_OneStr *s2 = new LE_OneStr("Interface Mass Transport Title",
				&(ei->Title), 100000, 1,  0, "ElectrodeName");
  cf->addLineEntry(s2);


  /* --------------------------------------------------------------
   * Temperature
   */
  LE_OneDbl *d1 = new LE_OneDbl("Temperature",
				&(ei->Temperature), 0, "Temperature");
  d1->set_default(300.);
  d1->set_limits(3000., 0.0);
  cf->addLineEntry(d1);


  /* --------------------------------------------------------------
   * Pressure -
   *
   * Configure the application Pressure
   */
  BE_UnitConversion *ucPres = new BE_UnitConversionPressure();
  LE_OneDblUnits *b5 = new LE_OneDblUnits("Phase A Pressure",
					  &(ei->PressureA), 0,
					  "PO.PressureA", ucPres);
  b5->set_default(OneAtm);
  b5->set_limits(1.E20, 0.0);
  cf->addLineEntry(b5);

  for (iph = 0; iph < pl->nVolPhases(); iph++) {
    ThermoPhase *tp = &(pl->volPhase(iph));
    string phaseNm = tp->name();
    if (ei->PhaseAName == phaseNm) {
      ei->solnAIndex_ = iph;
    }
    if (ei->PhaseBName == phaseNm) {
      ei->solnBIndex_ = iph;
    }
  }
  if (ei->solnAIndex_ < 0) {
    throw CanteraError("InterfacialMassTransfer_input::setup_input_pass3()",
		       "Phase A name, " + ei->PhaseAName + ", was not found in the list of Volume phases");
  }
  if (ei->solnBIndex_ < 0) {
    throw CanteraError("InterfacialMassTransfer_input::setup_input_pass3()",
		       "Phase B name, " + ei->PhaseAName + ", was not found in the list of Volume phases");
  }



  /*
   *  Set up the Bath Gas BG object to receive input
   */
  int nVolPhases = pl->nVolPhases();

  ElectrodeBath &BG = *(ei->m_BG);


  BG.XmolPLSpecVec = mdp_alloc_dbl_1(ei->nTotSpecies + 2, 0.0);
  BG.MolalitiesPLSpecVec = mdp_alloc_dbl_1(ei->nTotSpecies + 2, 0.0);
  BG.XmolPLPhases = (double **) mdp_alloc_ptr_1(nVolPhases + pl->nSurPhases());
  BG.MolalitiesPLPhases = (double **) mdp_alloc_ptr_1(nVolPhases + pl->nSurPhases());

  BG.PhaseMoles.resize(nVolPhases + pl->nSurPhases(), 0.0);
  BG.PhaseMass.resize(nVolPhases + pl->nSurPhases(), 0.0);

  for (iph = 0; iph < nVolPhases; iph++) {
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


  /* ---------------------------------------------------------------------------------------------------------
   *  Specify the phase A Boundary Conditions
   */
  string phaseBath = "Boundary Condition Specification for Phase ";
  ThermoPhase *tpA = &(pl->volPhase(ei->solnAIndex_));
  int nSpeciesA = tpA->nSpecies();
  phaseBath += ei->PhaseAName; 
  /*
   *  create a section method description block and start writing
   *  line elements in it.
   */
  BlockEntry *bphaseA = new BlockEntry(phaseBath.c_str());
  cf->addSubBlock(bphaseA);
  int kstart =  pl->getGlobalSpeciesIndexVolPhaseIndex(ei->solnAIndex_);

  /* --------------------------------------------------------------
   *     Specify the phase A mole fractions
   *
   * Create a PickList Line Element made out of the list of  species
   */
  string phaseMF = "Phase " + ei->PhaseAName + " Mole Fraction";
  BE_MoleComp_VecDbl *bmfA = new BE_MoleComp_VecDbl(phaseMF.c_str(),
						    &(ei->XmfPhaseA_), 1,
						    ei->SpeciesNames+kstart,
						    nSpeciesA, 0, "XMoleFractionA");
  bmfA->generateDefLE();
  bphaseA->addSubBlock(bmfA);
  /* 
   * --------------------------------------------------------------------------------------------------------- 
   */


  /* ---------------------------------------------------------------------------------------------------------
   *  Specify the phase B Boundary Conditions
   */
  phaseBath = "Boundary Condition Specification for Phase ";
  ThermoPhase *tpB = &(pl->volPhase(ei->solnBIndex_));
  int nSpeciesB = tpB->nSpecies();
  phaseBath += ei->PhaseBName; 
  /*
   *  create a section method description block and start writing
   *  line elements in it.
   */
  BlockEntry *bphaseB = new BlockEntry(phaseBath.c_str());
  cf->addSubBlock(bphaseB);
  kstart =  pl->getGlobalSpeciesIndexVolPhaseIndex(ei->solnBIndex_);

  /* --------------------------------------------------------------
   *     Specify the phase A mole fractions
   *
   * Create a PickList Line Element made out of the list of  species
   */ 
  phaseMF = "Phase " + ei->PhaseBName + " Mole Fraction";  
  BE_MoleComp_VecDbl *bmfB = new BE_MoleComp_VecDbl(phaseMF.c_str(),
						    &(ei->XmfPhaseB_), 1,
						    ei->SpeciesNames+kstart,
						    nSpeciesB, 0, "XMoleFractionB");
  bmfB->generateDefLE();
  bphaseB->addSubBlock(bmfB);
  /* 
   * --------------------------------------------------------------------------------------------------------- 
   */

  

  /*
   *  Specify a block for each Phase to receive inputs on composition
   *  and voltage
   */

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
     * BG.PhaseMoles
     *  Input the number of moles for the phase
     */
    LE_OneDbl *iTMoles = new LE_OneDbl("Phase Moles", &(BG.PhaseMoles[iph]), 0, "PhaseMoles");
    iTMoles->set_default(0.0);
 
    bbathphase->addLineEntry(iTMoles);
    // BI_Dependency * depimm_ig = new BI_Dependency(dpNum, BIDT_ENTRYPROCESSED, BIDRT_ZERONUMTIMESREQUIRED);
    // iTMoles->declareDependency(depimm_ig);
    //BI_Dependency *dep_phaseMoles_Porosity = new BI_Dependency(ePorosity, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    //iTMoles->declareDependency(dep_phaseMoles_Porosity);

    /* --------------------------------------------------------------
     * BG.PhaseMass
     *  Input the mass for the phase
     */
    LE_OneDbl *iTMass =   new LE_OneDbl("Phase Mass", &(BG.PhaseMass[iph]), 0, "PhaseMass");
    iTMass->set_default(0.0);

    bbathphase->addLineEntry(iTMass);
    //   If the 'Particle Number to Follow' card is in the input deck, then the 'Electrode Porosity'
    //   card is mandatory.
    //   QUESTION -- CAN THIS APPLY TO JUST ONE OF THE PHASES OR MUST IT APPLY TO ALL
    BI_Dependency *dep_phaseMass_phaseMoles = new BI_Dependency(iTMoles, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    iTMass->declareDependency(dep_phaseMass_phaseMoles);
    //BI_Dependency *dep_phaseMass_Porosity = new BI_Dependency(ePorosity, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    //iTMoles->declareDependency(dep_phaseMass_Porosity);

    /* --------------------------------------------------------------
     * BG.BathSpeciesMoleFractions - 
     *
     * Create a PickList Line Element made out of the list of 
     * species
     */
    // double *xstart = BG.XmolPLSpecVec + kstart;
    BE_MoleComp *bpmc = new BE_MoleComp("Bath Species Mole Fraction",
					&(BG.XmolPLPhases[iph]), 0,
					ei->SpeciesNames+kstart,
					nSpecies, 0, "XBathmol");
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
			    ei->SpeciesNames+kstart, nSpecies, 
			    indS, mwS,
			    "MolalitiesBath");
      //bmolal->generateDefLE();
      bbathphase->addSubBlock(bmolal);

    }

  }



  /*
   *  Specify a block for each surface Phase to receive inputs on composition
   *  and voltage
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
     * BG.BathSpeciesMoleFractions - 
     *
     * Create a PickList Line Element made out of the list of 
     * species
     */
    BE_MoleComp *bpmc = new BE_MoleComp("Bath Species Mole Fraction",
					&(BG.XmolPLPhases[iph]), 0,
					ei->SpeciesNames+kstart, nSpecies, 0, "XBathmol");
    bpmc->generateDefLE();
    bbathphase->addSubBlock(bpmc);

  }

  /*
   *   SPECIFY OTHER WAYS TO INPUT THE CONCENTRATION, without
   *   breaking it down via phases.
   */


  /* ------------------------------------------------------------------
   * Block Input For initial number of moles of species
   *
   */
  BlockEntry * BESIM = 
    new BE_StrDbl("Species Initial KMoles", &ei->MoleNumber,
		  0, 0, ei->SpeciesNames, ei->nTotSpecies, 1,
		  "MoleNumber");
  cf->addSubBlock(BESIM);

 
  /* ------------------------------------------------------------------
   * Define a block that may occur multiples times
   *  The number of occurances of the block is .... PO.numExtraGlobalRxns
   *  All of the data defined in the block is
   *  in the pointer to the vector of structures .. PO.m_EGRList
   *
   */ 

  
  BaseEntry::set_SkipUnknownEntries(false);
}
//======================================================================================================================
/******************************************************************************
 *  Read the input file
 *
 *  printFlag 0 = no output
 *            1 = error message output
 *            2 = output
 */
bool process_electrode_input(BlockEntry *cf, string fileName, int printFlag) {
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



/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/*
 *  Set the bath gas for the one thermophase. Note, we hardcode
 *  what we want here. An input deck might be the more appropriate
 *  place to set these properties.
 */
void
setElectrodeBathSpeciesConditions(ThermoPhase& g, IMT_KEY_INPUT &EI, ElectrodeBath &BG, int iph, int printLvl) {

  int nsp = g.nSpecies();


  /*
   * We need to set the state here. Note, mass fractions
   * don't matter, since we are only querying the 
   * standard state of the species.
   */
  g.setState_TPX(EI.Temperature, EI.PressureA, BG.XmolPLPhases[iph]);



  /*
   * Print out table summarizing bath gas conditions
   */ 
  if (printLvl) {
    double *C = new double [nsp];
    double *act = new double [nsp];
    g.getConcentrations(C);
    g.getActivities(act);
    print_char('=', 100);
    cout << endl;
    dnt(0);
    cout << " SUMMARY OF SPECIES IN THE MECHANISM WITH A DESCRIPTION "
	 << " OF BATH COMPOSITION: " << endl; 
    cout << endl;
    dnt(1);
    cout << "Total pressure = " << EI.PressureA * 760. / OneAtm;
    cout << " torr " << endl;
    dnt(1); cout << "Temperature (where needed) = " << EI.Temperature
		 << " Kelvin" << endl;



    dnt(4); 
    cout << "Number       Name      Mole_fraction  Concentration (gmol/cm**3)   Activities";
    cout << endl;
    dnt(4);
    cout << "-----------------------------------------------------------------------------";
    cout << endl;
    string spN;
    for (size_t k = 0; k < g.nSpecies(); k++) {
      dnt(4); 
      pr_if(k+1, 5);
      spN = g.speciesName(k);
      pr_sf(spN, 16);  
      pr_df(BG.XmolPLPhases[iph][k], 16, 4); 
      pr_de(C[k] * 1.0E-3, 16, 3);
      pr_de(act[k], 16, 3);
      cout << endl;
    }
    dnt(4);
    cout << "------------------------------------------------------------------------------";
    cout << endl;
    delete [] C;
    delete [] act;
 
  }
}
 



/*****************************************************************************
 *
 * imt_input():
 *
 *    INput for vcs_Cantera. This routine will combine a text input file
 *    with a cantera xml or cti file to create an equilibrium problem
 *    to solve.
 */
int
imt_input(IMT_KEY_INPUT *ei, string commandFile, BlockEntry *cf)
{
  int  retn = MPEQUIL_SUCCESS;
   
  printf("\n"); print_char('=', 80); printf("\n");
  print_char('=', 20); printf(" electrode_input: START OF PROBLEM STATEMENT ");
  print_char('=', 21); printf("\n");
  print_char('=', 80); printf("\n\n");  

  int printBIProclevel = 2;

  /*
   * Setup and process the input deck for first time.
   * -> Might have to print out the input and quit as well.
   */
  setup_input_pass1(cf, ei);
  bool ok = process_electrode_input(cf, commandFile, printBIProclevel);
  if (!ok) {
    return -1;
  }


  /*
   * Possibly change the print level for the last
   */
  IMT_Types_Enum ieos =  string_to_IMT_Types_Enum(ei->IMT_ModelName);
  if (ieos == UNKNOWN_IMT) {
    printf("unknown imt model\n");
  }

  /**
   * Setup and process the input deck for second time.
   * -> Might have to print out the input and quit as well.
   */
  setup_input_pass2(cf, ei);
  ok = process_electrode_input(cf, commandFile, printBIProclevel);
  if (!ok) {
    return -1;
  }

  int ifiles = 0;
  for (; ei->CanteraFileNames[ifiles] != 0; ifiles++) {
  }
  if (ifiles != ei->NumberCanteraFiles) {
    printf("Number of requested files differ\n");
    exit(-1);
  }

  /*
   * Read in all of the phase specifications from the cantera
   * input files into the PhaseList structure.
   */
  PhaseList *pl = ei->m_pl;
  std::string fn;
  bool surNotFound = true;
  for (int i = 0; i < ei->NumberCanteraFiles; i++) {
    fn = ei->CanteraFileNames[i];
    importAllCTMLIntoPhaseList(pl, fn);
    if (surNotFound && (pl->nSurPhases() > 0)) {
      surNotFound = false;
      //   pl->CanteraFNSurface = fn;
    }
  }
  /*
   * Setup internally for next pass through the input file.
   */
  ei->InitForInput(pl);

  /**
   * Setup and process the input deck for third time
   * -> Might have to print out the input and quit as well.
   */
  setup_input_pass3(cf, ei);
     
  /*
   * Process the third pass of the input file ->
   *   We read everything this time.
   */
  ok = process_electrode_input(cf, commandFile, printBIProclevel);
  if (!ok) {
    return -1;
  }


  /*
   *     OK, We have Everything Read in
   *
   *   Query the input file for some more information
   */

  bool molesSpecified = false;


  /*
   * Determine the total number of kmols for each species
   */
  // Note that as of 11/17/10 the variable molesSpecified 
  // and member data ei->specifiedBlockKmolSpecies are not 
  // used anywhere but these few lines.

  BlockEntry *be = cf->searchBlockEntry("Species Initial KMoles");
  ei->specifiedBlockKmolSpecies = be->get_NumTimesProcessed();
  if (ei->specifiedBlockKmolSpecies > 0) {
    molesSpecified = true;
  }

  ElectrodeBath *BG = ei->m_BG;


  for (int iph = 0; iph < pl->nVolPhases(); iph++) {
    string phaseBath = "Bath Specification for Phase ";
    ThermoPhase *tp = &(pl->volPhase(iph));
    int nSpecies = tp->nSpecies();
    string phaseNm = tp->name();
    //int nSpecies = tp->nSpecies();
    phaseBath += phaseNm; 
    BlockEntry *be = cf->searchBlockEntry(phaseBath.c_str());
    int numTimes = be->get_NumTimesProcessed();
    double *molF = BG->XmolPLPhases[iph];
    if (numTimes > 0) {
      double kmol = BG->PhaseMoles[iph]; //if not kmol given, compute from mass
      int kstart =  pl->getGlobalSpeciesIndexVolPhaseIndex(iph);
      for (int k = 0; k < nSpecies; k++) {
	ei->MoleFraction[kstart + k] = molF[k];
      }
      tp->setMoleFractions(molF);
  

      // if the mass and not the number of moles was set...
      if (!(kmol > 0.0)) {
	kmol = BG->PhaseMass[iph] / tp->meanMolecularWeight();
	BG->PhaseMoles[iph] = kmol;
      } else if ( BG->PhaseMass[iph] > 0.0 ) {
	throw CanteraError("electrode_input()",
			   "both number of moles and mass of phase specified");
      }
 
      //update number of moles (gets done again in electrode_model_init())
      for (int k = 0; k < nSpecies; k++) {
	ei->MoleNumber[kstart + k] += molF[k] * kmol;
      }
    }
  }

  imt_model_init(ei, cf);

  return retn;
}




/**************************************************************************/

//=========================================================================================
 int imt_model_init(IMT_KEY_INPUT *ei,  BEInput::BlockEntry *cf)
{

  int iph;

 
  ElectrodeBath *BG_ptr = ei->m_BG;

  PhaseList *pl = ei->m_pl;
  /*
   *  Loop Over all phases in the PhaseList, adding these
   *  formally to the electrodeCell object.
   */
  for (iph = 0; iph < pl->nPhases(); iph++) {
  
    ThermoPhase *tphase = &(pl->thermo(iph));
    int nSpecies = tphase->nSpecies();

    // Find the name of the input block
    string phaseBath = "Bath Specification for Phase "; 
    string phaseNm = tphase->name();
    phaseBath += phaseNm; 

    
    bool molVecSpecified = false;
    BEInput::BlockEntry *pblock = cf->searchBlockEntry(phaseBath.c_str());
    if (pblock) {
      BEInput::BlockEntry *pbsmf =
	pblock->searchBlockEntry("Bath Species Mole Fraction");
      if (pbsmf) {
	if (pbsmf->get_NumTimesProcessed() > 0) {
	  molVecSpecified = true;
	}
      }
    }

    double *molesSpecies = new double[nSpecies]; 

    bool molalVecSpecified = false;
    if (pblock) {
      BEInput::BlockEntry *pbsmm =
	pblock->searchBlockEntry("Bath Species Molalities");
      if (pbsmm) {
	if (pbsmm->get_NumTimesProcessed() > 0) {
	  molalVecSpecified = true;
	}
      }
    }
    if (molalVecSpecified) {
      MolalityVPSSTP * m_ptr = dynamic_cast<MolalityVPSSTP *>(tphase);
      if (m_ptr == 0) {
	printf("Dynamic cast failed for some reason\n");
	exit(-1);
      }
      m_ptr->setState_TPM(ei->Temperature, ei->PressureA,
			  BG_ptr->MolalitiesPLPhases[iph]);
      m_ptr->getMoleFractions(BG_ptr->XmolPLPhases[iph]);
    }
    //    electrodeA->phaseMoles_[iph] = BG_ptr->PhaseMoles[iph];
    double totalMoles = BG_ptr->PhaseMoles[iph];

    /*
     *  Setup the global arrays in the electrode object
     */
    //    int kstart = electrodeA->getGlobalSpeciesIndex(iph, 0);
    int  istart = pl->getGlobalSpeciesIndex(iph, 0);
    for (int k = 0; k < nSpecies; k++) {
      molesSpecies[k] =  totalMoles * BG_ptr->XmolPLPhases[iph][k];
      // electrodeA->spMf[kstart+k] = BG_ptr->XmolPLPhases[iph][k];
      //      electrodeA->spMoles[kstart+k] = totalMoles * BG_ptr->XmolPLPhases[iph][k];
      ei->MoleNumber[istart + k] = molesSpecies[k];

    }
    //electrodeA->setPhaseMoleNumbers(iph, molesSpecies);
    delete [] molesSpecies;

  }

  /*
   * Set up the MultiPhase object. Right now it will contain all of the
   * volume phases only.
   */
  //  electrodeA->downloadMP(); 



  return 0;
}
//=========================================================================================
int imt_model_print(Cantera::InterfacialMassTransfer *electrodeA,  IMT_KEY_INPUT *ei,
		    BEInput::BlockEntry *cf)
{

  int i, iph;

 
  ElectrodeBath *BG_ptr = ei->m_BG;


  /*
   * Load up the temperature and pressure
   */
  //electrodeA->T = ei->Temperature;
  //  electrodeA->Pres = ei->Pressure;

  electrodeA->setState_TP(ei->Temperature, ei->PressureA);
  /*
   *  Loop Over all phases in the PhaseList, adding these
   *  formally to the electrodeCell object.
   */
  for (iph = 0; iph < electrodeA->nPhases(); iph++) {
  
    ThermoPhase *tphase = &(electrodeA->thermo(iph));
    int nSpecies = tphase->nSpecies();

    // Find the name of the input block
    string phaseBath = "Bath Specification for Phase "; 
    string phaseNm = tphase->name();
    phaseBath += phaseNm; 

    /*
     * Set up the Bath Gas Conditions in each of
     * the ThermoPhase objects
     */

    
    bool molVecSpecified = false;
    BEInput::BlockEntry *pblock = cf->searchBlockEntry(phaseBath.c_str());
    if (pblock) {
      BEInput::BlockEntry *pbsmf =
	pblock->searchBlockEntry("Bath Species Mole Fraction");
      if (pbsmf) {
	if (pbsmf->get_NumTimesProcessed() > 0) {
	  molVecSpecified = true;
	}
      }
    }


    double *molesSpecies = new double[nSpecies]; 

    bool molalVecSpecified = false;
    if (pblock) {
      BEInput::BlockEntry *pbsmm =
	pblock->searchBlockEntry("Bath Species Molalities");
      if (pbsmm) {
	if (pbsmm->get_NumTimesProcessed() > 0) {
	  molalVecSpecified = true;
	}
      }
    }
    if (molalVecSpecified) {
      MolalityVPSSTP * m_ptr = dynamic_cast<MolalityVPSSTP *>(tphase);
      if (m_ptr == 0) {
	printf("Dynamic cast failed for some reason\n");
	exit(-1);
      }
      m_ptr->setState_TPM(ei->Temperature, ei->PressureA, BG_ptr->MolalitiesPLPhases[iph]);
      m_ptr->getMoleFractions(BG_ptr->XmolPLPhases[iph]);
    }
    //    electrodeA->phaseMoles_[iph] = BG_ptr->PhaseMoles[iph];
    double totalMoles = BG_ptr->PhaseMoles[iph];

    setElectrodeBathSpeciesConditions(*tphase, *ei, *BG_ptr, iph, 0); 

    /*
     *  Setup the global arrays in the electrode object
     */
    //    int kstart = electrodeA->getGlobalSpeciesIndex(iph, 0);
    for (int k = 0; k < nSpecies; k++) {
      molesSpecies[k] =  totalMoles * BG_ptr->XmolPLPhases[iph][k];
      // electrodeA->spMf[kstart+k] = BG_ptr->XmolPLPhases[iph][k];
      //      electrodeA->spMoles[kstart+k] = totalMoles * BG_ptr->XmolPLPhases[iph][k];
    }
    electrodeA->setPhaseMoleNumbers(iph, molesSpecies);
    delete [] molesSpecies;

    //    tphase->getPartialMolarVolumes(&(electrodeA->VolPM[kstart]));
    //   tphase->getElectrochemPotentials(&(electrodeA->spElectroChemPot[kstart]));


  }



  /*
   * Query whether an initial estimate has been made and then set iest.
   * Copy guess into vprobin
   */


  /*
   *          Printout the species information: PhaseID's and mole nums
   */
  printf("\n"); print_char('-', 80); printf("\n");
  printf("             Phase IDs of species\n");
  printf("            species     phaseID        phaseName   ");
  printf(" Initial_Estimated_KMols\n");

  for (iph = 0; iph < electrodeA->nPhases(); iph++) {
    ThermoPhase *tphase = &(electrodeA->thermo(iph));
    int nspeciesP = tphase->nSpecies();
    string pName = tphase->id();
    int kstart =  electrodeA->getGlobalSpeciesIndex(iph, 0);
    for (i = 0; i < nspeciesP; i++) {
      int kT = kstart + i;
      string spName = tphase->speciesName(i);
      printf("%16s      %5d   %16s",
	     spName.c_str(), iph, pName.c_str()); 
      printf("             %-10.5g", electrodeA->moleNumSpecies(kT));
      printf("\n");
   
    }
  }

  /*
   *   Printout of the Phase structure information
   */
  printf("\n"); print_char('-', 80); printf("\n");
  printf("             Information about phases\n");
  printf("  PhaseName    PhaseNum SingSpec GasPhase NumSpec");
  printf("  TMolesInert       TKmols\n");
   
  for (iph = 0; iph < electrodeA->nPhases(); iph++) {
    ThermoPhase *tphase = &(electrodeA->thermo(iph));
    int nspeciesP = tphase->nSpecies();
    string pName = tphase->id();
    printf("%16s %5d ", pName.c_str(), iph);
    if (nspeciesP > 1) {
      printf("  no  ");
    } else {
      printf("  yes ");
    }
    printf("%8d  ", nspeciesP);
    printf(" %11g  ", 0.0);
    printf("%16e\n", electrodeA->phaseMoles(iph));
  }

  /*
   *   Printout the Element information
   */
  printf("\n"); print_char('-', 80); printf("\n");
  printf("             Information about Elements\n");
  printf("     ElementName  Abundance_Kmols\n");
  for (i = 0; i < electrodeA->nElements(); ++i) {
    string eName = electrodeA->elementName(i);
    printf("%12s ", eName.c_str());
    printf("  %11g \n", electrodeA->elementMoles(i));
  }


  printf("\n"); print_char('=', 80); printf("\n");
  print_char('=', 20);
  printf(" mpequil: END OF PROBLEM STATEMENT ");
  print_char('=', 23); printf("\n");
  print_char('=', 80); printf("\n\n");

  
  return 0;
}


//=========================================================================================

}
//=========================================================================================
