/*
 * $Id: epequil_input.cpp 511 2013-01-07 23:32:30Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"



#include "cantera/equilibrium.h"
#include "epequil_input.h"


#include "BlockEntryGlobal.h"
#include "importAllPhases.h"
#include "mdp_allo.h"

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
using namespace std;
using namespace BEInput;
using namespace mdpUtil;

/*************************************************************************
 *
 */

EPEQUIL_KEY_INPUT PO;
/*************************************************************************
 * EPEQUIL_KEY_INPUT(): constructor
 */
EPEQUIL_KEY_INPUT::EPEQUIL_KEY_INPUT () :
    CanteraFN1("gas.cti"),
    NumberCanteraFiles(1),
    CanteraFileNames(0),
    Temperature(300.),
    Pressure(OneAtm),
    Vol(1.0),
    MoleNumber(0),
    MoleNumberIG(0),
    PhaseInclude(0),
    ProblemType(0),
    SpeciesNames(0),
    PhaseNames(0),
    ElementNames(0),
    ElementAbundances(0),
    specifiedElementAbundances(false)
{
    Title = "epequil Cantera Problem";
}
/****************************************************************************
 *
 */
EPEQUIL_KEY_INPUT::~EPEQUIL_KEY_INPUT () {
  if (CanteraFileNames) {
    for (int i = 0; CanteraFileNames[i] != 0; i++) {
      mdp_safe_free((void **) &(CanteraFileNames[i]));
    }
    mdp_safe_free((void **) &(CanteraFileNames));
  }
  mdp_safe_free((void **) &PhaseInclude);
  mdp_safe_free((void **) &MoleNumber);
  mdp_safe_free((void **) &MoleNumberIG);
  mdp_safe_free((void **) &SpeciesNames);
  mdp_safe_free((void **) &PhaseNames);
  mdp_safe_free((void **) &ElementNames);
  mdp_safe_free((void **) &ElementAbundances);
}
/****************************************************************************
 *
 */
void EPEQUIL_KEY_INPUT::InitForInput(MultiPhase *mp) {
    int nTotPhases = mp->nPhases();
    int nTotSpecies = mp->nSpecies();
    int nTotElements = mp->nElements();
    
    /*
     * Include all Phases by default
     */
    PhaseInclude = mdp_alloc_int_1(nTotPhases, 1);
    MoleNumber   = mdp_alloc_dbl_1(nTotSpecies, 0.0);
    MoleNumberIG = mdp_alloc_dbl_1(nTotSpecies, 0.0);
    ElementAbundances = mdp_alloc_dbl_1(nTotElements, 0.0);


    SpeciesNames = mdp_alloc_VecFixedStrings(nTotSpecies,
                                             EPEQUIL_MAX_NAME_LEN_P1);
    PhaseNames = mdp_alloc_VecFixedStrings(nTotPhases,
					   EPEQUIL_MAX_NAME_LEN_P1);
    int kT = 0;
    for (int iphase = 0; iphase < nTotPhases; iphase++) {
      ThermoPhase *tPhase = &(mp->phase(iphase));
      string id = tPhase->id();
      strncpy(PhaseNames[iphase], id.c_str(), EPEQUIL_MAX_NAME_LEN);
      int nspecies = tPhase->nSpecies();
      for (int k = 0; k < nspecies; k++) {
	string sname = tPhase->speciesName(k);
	strncpy(SpeciesNames[kT], sname.c_str(), EPEQUIL_MAX_NAME_LEN);
	kT++;
      }
    }

    ElementNames = mdp_alloc_VecFixedStrings(nTotElements,
					     EPEQUIL_MAX_NAME_LEN_P1);
    for (int e = 0; e < nTotElements; e++) {
      string eName = mp->elementName(e);
      strncpy(ElementNames[e], eName.c_str(), EPEQUIL_MAX_NAME_LEN);
    }
}
   
/****************************************************************************
 *
 */
static void setup_input_pass1(BlockEntry *cf)
{
    /*
     * Obtain the number of cantera files to be read
     */
    LE_OneInt *s1 = new LE_OneInt("Number of Cantera Files",
				  &PO.NumberCanteraFiles, 0,
				  "NumCanteraFiles");
    s1->set_default(1);
    cf->addLineEntry(s1);
    BaseEntry::set_SkipUnknownEntries(3);
}

/****************************************************************************
 *
 */
static void setup_input_pass2(BlockEntry *cf)
{
    LineEntry *sle1 = 0;
    /*
     *  Get the input deck for
     *  Cantera description of the model.
     */
    LE_MultiCStr *s1 =
	new LE_MultiCStr("Cantera File Name", &PO.CanteraFileNames,
			 1, 1,  0, "CanteraFileNames");
    s1->set_default("gas.cti");

    /*
     * Set up a dependency on the input from the Number of cantera
     * Files card
     */
    sle1 = cf->searchLineEntry("Number of Cantera Files");
    int numF = 1; 
    bool okbefore = sle1->ansDepCheckOneInt(numF);
    if (okbefore) {
      printf("Found it before\n");
      s1->set_NumTimesRequired(numF);
    } else {
      printf("Num Lines not in input deck -> no dependency\n");
      // Note -> this is not right -> should be one or more dependencies.
    }
    

    cf->addLineEntry(s1);
    BaseEntry::set_SkipUnknownEntries(3);
}
                   
/****************************************************************************
 *
 */
static void setup_input_pass3(BlockEntry *cf, EPEQUIL_INPUT *pi)
{
    MultiPhase *mp = pi->m_mp;
    /* ---------------------------------------------------------------
     *
     */
    LE_OneStr *s2 = new LE_OneStr("Title",
				  &PO.Title, 100000, 1,  0, "Title");
    cf->addLineEntry(s2);

    /* --------------------------------------------------------------
     * Temperature
     */
    LE_OneDbl *d1 = new LE_OneDbl("Temperature",
				  &PO.Temperature, 0, "Temperature");
    d1->set_default(300.);
    d1->set_limits(3000., 0.0);
    cf->addLineEntry(d1);


    /* --------------------------------------------------------------
     * Pressure -
     *
     * Configure the application Pressure
     */
    BE_UnitConversion *ucPres = new BE_UnitConversionPressure();
    LE_OneDblUnits *b5 = new LE_OneDblUnits("Pressure",
                                            &(PO.Pressure), 0,
                                            "PO.Pressure", ucPres);
    b5->set_default(OneAtm);
    b5->set_limits(1.E20, 0.0);
    cf->addLineEntry(b5);
    
    /* ------------------------------------------------------------------
     * Block Input For Element abundances
     *
     */
    BlockEntry * BEA = 
	new BE_StrDbl("KMolar Element Abundances", &PO.ElementAbundances,
		      0, 0, PO.ElementNames, mp->nElements(), 1,
		      "ElementAbundances");
    cf->addSubBlock(BEA);
    /* ------------------------------------------------------------------
     * Block Input For the initial guess for species
     *
     */
    BlockEntry * BESIMG = 
	new BE_StrDbl("Species Initial KMoles Guess", &PO.MoleNumberIG,
		      0, 0, PO.SpeciesNames, mp->nSpecies(), 1,
		      "MoleNumberIG");
    cf->addSubBlock(BESIMG);
    /* ------------------------------------------------------------------
     * Block Input For initial number of moles of species
     *
     */
    BlockEntry * BESIM = 
	new BE_StrDbl("Species Initial KMoles", &PO.MoleNumber,
		      0, 0, PO.SpeciesNames, mp->nSpecies(), 1,
		      "MoleNumber");
    cf->addSubBlock(BESIM);

    BaseEntry::set_SkipUnknownEntries(0);
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
    const TK_TOKEN tok_in;
    TK_TOKEN tok_out;
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


/**************************************************************************
 *
 *  EPEQUIL_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.  
 */

EPEQUIL_INPUT::EPEQUIL_INPUT() :
    tplist(0),
    m_mp(0),
    prob_type(0),
    nspecies(0),
    ne(0),
    nphase(0),
    T(298.15),
    Pres(1.0E5),
    Vol(0),
    VolPM(0),
    spMoles(0),
    spMf(0),
    spChemPot(0),
    elementMoles(0),
    phaseMoles(0),
    Title(0),
    iest(-1),
    specifiedElementAbundances(false)
{
}

/**************************************************************************
 *
 *  EPEQUIL_CANTERAPROB_INPUT:destructor
 *
 * We need to manually free all of the arrays.
 */
EPEQUIL_INPUT::~EPEQUIL_INPUT() 
{
    for (int i = 0; i < nphase; i++) {
      delete tplist[i];
    }
    mdp_safe_free((void **) &tplist);
    delete m_mp; m_mp = 0;
    mdp_safe_free((void **) &VolPM);
    mdp_safe_free((void **) &spMoles);
    mdp_safe_free((void **) &spMf);
    mdp_safe_free((void **) &spChemPot);
    mdp_safe_free((void **) &elementMoles);
    mdp_safe_free((void **) &phaseMoles);
    mdp_safe_free((void **) &Title);
}

/*****************************************************/
static void print_char(const char letter, const int num)
{
   int i;
   for (i = 0; i < num; i++) printf("%c", letter);
}

/*****************************************************************************
 *
 * epequil_input():
 *
 *    Input for vcs_Cantera. This routine will combine a text input file
 *    with a cantera xml or cti file to create an equilibrium problem
 *    to solve.
 */
int
epequil_input(EPEQUIL_INPUT *pi, string commandFile)
{
   int   i, k, e;
   int  retn = EPEQUIL_SUCCESS;
   int   iphase, kT;
   ThermoPhase *tphase = 0;
   
   
   printf("\n"); print_char('=', 80); printf("\n");
   print_char('=', 20); printf(" vcs_input: START OF PROBLEM STATEMENT ");
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
    bool ok = process_input(cf, commandFile, printBIProclevel);
    if (!ok) {
      return -1;
    }
    /**
     * Setup and process the input deck for second time.
     * -> Might have to print out the input and quit as well.
     */
    setup_input_pass2(cf);
    ok = process_input(cf, commandFile, printBIProclevel);
    if (!ok) {
      return -1;
    }

    int ifiles = 0;
    for (; PO.CanteraFileNames[ifiles] != 0; ifiles++) {
    }
    if (ifiles != PO.NumberCanteraFiles) {
      printf("Number of requested files differ\n");
      exit(-1);
    }

    if (PO.CanteraFN1.size() == 0) {
      printf("Cantera file name must be specified\n");
      exit (-1);
    }
    /**
     * Read in all of the phase specifications from the cantera
     * input files into Cantera's MultiPhase structure.
     */
    pi->m_mp = new MultiPhase();
    MultiPhase *mp = pi->m_mp;
    for (int i = 0; i < PO.NumberCanteraFiles; i++) {
      PO.CanteraFN1 = PO.CanteraFileNames[i];
      importAllCTML(pi, PO.CanteraFN1);
    }

    /*
     * Setup internally for next pass through the input file.
     */
    PO.InitForInput(mp);

    /**
     * Setup and process the input deck for second time
     * -> Might have to print out the input and quit as well.
     */
    setup_input_pass3(cf, pi);
     
    /*
     * Process the first pass of the input file ->
     *   We are just after the information needed to initialize
     *   the Cantera structures and size the problem
     */
    ok = process_input(cf, commandFile, printBIProclevel);
    if (!ok) {
      return -1;
    }
    mp->init();
   
    pi->nspecies= mp->nSpecies();
    pi->ne = mp->nElements();
    pi->nphase = mp->nPhases();
    pi->spMoles   = mdp_alloc_dbl_1(pi->nspecies, 0.0);
    pi->spMf      = mdp_alloc_dbl_1(pi->nspecies, 0.0);
    pi->VolPM     = mdp_alloc_dbl_1(pi->nspecies, 0.0);
    pi->spChemPot = mdp_alloc_dbl_1(pi->nspecies, 0.0);
    pi->elementMoles= mdp_alloc_dbl_1(pi->ne, 0.0);
    pi->phaseMoles  = mdp_alloc_dbl_1(pi->nphase, 0.0);
 
    /*
     * Query whether an initial estimate has been made and then set iest.
     * Copy guess into vprobin
     */
    BlockEntry *be = cf->searchBlockEntry("Species Initial KMoles Guess");
    int guessSpecified= be->get_NumTimesProcessed();
    if (guessSpecified > 0) {
      pi->iest = 0;
      mdp_copy_dbl_1(pi->spMoles, PO.MoleNumberIG, pi->nspecies);
    } else {
      pi->iest = -1;
    }

   /*
    * Store the temperature and pressure
    */
   pi->T = PO.Temperature;
   pi->Pres = PO.Pressure; /* note resulting pres will have mks units */
   pi->Vol = PO.Vol;
   pi->Title = mdp_copy_string(PO.Title.c_str());

   /*
    * Determine the total number of kmols for each species
    */
   be = cf->searchBlockEntry("KMolar Element Abundances");
   PO.specifiedElementAbundances = be->get_NumTimesProcessed();

   be = cf->searchBlockEntry("Species Initial KMoles");
   int specifiedSpeciesMoles = be->get_NumTimesProcessed();

   if (PO.specifiedElementAbundances && specifiedSpeciesMoles) {
     printf("Can't specify initial conditions 2 ways\n");
     exit(-1);
   }

   if (specifiedSpeciesMoles) {
     if (!guessSpecified) {
       mdp_copy_dbl_1(pi->spMoles, PO.MoleNumber, pi->nspecies);
     }
     kT = 0;
     for (iphase = 0; iphase < pi->nphase; iphase++) {
       ThermoPhase *tPhase = &(mp->phase(iphase));
       pi->phaseMoles[iphase] = 0.0;
       for (k = 0; k < (int) tPhase->nSpecies(); k++) {
	 pi->phaseMoles[iphase] += PO.MoleNumber[kT];
	 for (e = 0; e < pi->ne; e++) {
	   double natoms = mp->nAtoms(kT, e);
	   pi->elementMoles[e] += natoms * PO.MoleNumber[kT];
	 }
	 kT++;
       }
     }
   } else {
     if (PO.specifiedElementAbundances) {
       /*
	* If we have specified the element abundances, we just copy
	* keeping the units as kmols.
	*/
       pi->specifiedElementAbundances = true;
       for (i = 0; i < pi->ne; i++) {
	 pi->elementMoles[i] = PO.ElementAbundances[i];
       }
     } else {
       printf("Either specify species moles or elemental "
	      "abundances\n");
       exit(-1);
     }
   }


   /*
    *          Printout the species information: PhaseID's and mole nums
    */
   printf("\n"); print_char('-', 80); printf("\n");
   printf("             Phase IDs of species\n");
   printf("            species     phaseID        phaseName   ");
   printf(" Initial_Estimated_KMols\n");
   kT = 0;
   for (iphase = 0; iphase < pi->nphase; iphase++) {
     tphase = &(mp->phase(iphase));
     int nspeciesP = tphase->nSpecies();
     string pName = tphase->id();
     for (i = 0; i < nspeciesP; i++) {
       string spName = mp->speciesName(kT);
       printf("%16s      %5d   %16s",
	      spName.c_str(), iphase, pName.c_str()); 
       if (! pi->specifiedElementAbundances) {
	 printf("             %-10.5g\n", pi->spMoles[kT]);
       } else {     
	 printf("                N/A\n");
       }
       kT++;
     }
   }

   /*
    *   Printout of the Phase structure information
    */
   printf("\n"); print_char('-', 80); printf("\n");
   printf("             Information about phases\n");
   printf("  PhaseName    PhaseNum SingSpec GasPhase NumSpec");
   printf("  TMolesInert       TKmols\n");
   
   for (iphase = 0; iphase < pi->nphase; iphase++) {
     tphase = &(mp->phase(iphase));
     int nspeciesP = tphase->nSpecies();
     string pName = tphase->id();
     printf("%16s %5d ", pName.c_str(), iphase);
     if (nspeciesP > 1) {
       printf("  no  ");
     } else {
       printf("  yes ");
     }
     printf("%8d  ", nspeciesP);
     printf(" %11g  ", 0.0);
     if (! pi->specifiedElementAbundances) {
       printf("%16e\n", pi->phaseMoles[iphase]);
     } else {
       printf("        N/A\n");	
     }
   }

   /*
    *   Printout the Element information
    */
   printf("\n"); print_char('-', 80); printf("\n");
   printf("             Information about Elements\n");
   printf("     ElementName  Abundance_Kmols\n");
   for (i = 0; i < pi->ne; ++i) {
     string eName = mp->elementName(i);
     printf("%12s ", eName.c_str());
     printf("  %11g \n", pi->elementMoles[i]);
   }


   printf("\n"); print_char('=', 80); printf("\n");
   print_char('=', 20);
   printf(" epequil: END OF PROBLEM STATEMENT ");
   print_char('=', 23); printf("\n");
   print_char('=', 80); printf("\n\n");

   delete cf;
   return retn;
}
/**************************************************************************/
