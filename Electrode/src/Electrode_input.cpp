/*
 * $Id: Electrode_input.cpp 576 2013-03-27 23:13:53Z hkmoffa $
 */



#include "Electrode_input.h"
#include "Electrode_Factory.h"
#include "importPL.h"

#include "BlockEntryGlobal.h"

using namespace std;
using namespace BEInput;
using namespace TKInput;
using namespace ca_ab;
using namespace mdpUtil;

#include <string>
//================================================================================================
ElectrodeBath::ElectrodeBath() :
    XmolPLSpecVec(0),
    XmolPLPhases(0),
    MolalitiesPLSpecVec(0),
    MolalitiesPLPhases(0),
    CapLeftCoeffPhases(0),
    CapLeftCoeffSpecVec(0),
    CapZeroDoDCoeffPhases(0),
    CapZeroDoDCoeffSpecVec(0)
{
}
//================================================================================================
ElectrodeBath::ElectrodeBath(const ElectrodeBath &right) :
    XmolPLSpecVec(0),
    XmolPLPhases(0),
    MolalitiesPLSpecVec(0),
    MolalitiesPLPhases(0),
    CapLeftCoeffPhases(0),
    CapLeftCoeffSpecVec(0),
    CapZeroDoDCoeffPhases(0),
    CapZeroDoDCoeffSpecVec(0)
{
    operator=(right);
}
//================================================================================================
 ElectrodeBath& ElectrodeBath::operator=(const ElectrodeBath& right)
 {
    if (this == &right) {
        return *this;
     }

    // Shallow pointer representation -> this is wrong and must be fixed up in parent routine.
    XmolPLSpecVec                = right.XmolPLSpecVec;
    XmolPLPhases                 = right.XmolPLPhases;
    MolalitiesPLSpecVec          = right.MolalitiesPLSpecVec;
    MolalitiesPLPhases           = right.MolalitiesPLPhases;
    CapLeftCoeffPhases           = right.CapLeftCoeffPhases;
    CapLeftCoeffSpecVec          = right.CapLeftCoeffSpecVec;
    CapZeroDoDCoeffPhases        = right.CapZeroDoDCoeffPhases;
    CapZeroDoDCoeffSpecVec       = right.CapZeroDoDCoeffSpecVec;
 
    PhaseMoles                   = right.PhaseMoles;
    PhaseMass                    = right.PhaseMass;
 
    return *this;
 }
//================================================================================================
ElectrodeBath::~ElectrodeBath() 
{
    mdpUtil::mdp_safe_free((void**) &XmolPLSpecVec);
    mdpUtil::mdp_safe_free((void**) &XmolPLPhases);
    mdpUtil::mdp_safe_free((void**) &MolalitiesPLSpecVec);
    mdpUtil::mdp_safe_free((void**) &MolalitiesPLPhases);
    mdpUtil::mdp_safe_free((void**) &CapLeftCoeffPhases);
    mdpUtil::mdp_safe_free((void**) &CapLeftCoeffSpecVec);
    mdpUtil::mdp_safe_free((void**) &CapZeroDoDCoeffPhases);
    mdpUtil::mdp_safe_free((void**) &CapZeroDoDCoeffSpecVec);
}
/*************************************************************************
 *
 * MPEQUIL_KEY_INPUT(): constructor
 */
ELECTRODE_KEY_INPUT::ELECTRODE_KEY_INPUT(int printLvl) :
    printLvl_(printLvl),
    commandFile_(""),
    CanteraFNSurface(""),
    NumberCanteraFiles(1),
    CanteraFileNames(0),
    Temperature(300.),
    Pressure(OneAtm),
    Vol(1.0),
    m_BG(0),
    MoleNumber(0),
    MoleFraction(0),
    PotentialPLPhases(0),
    PhaseInclude(0),
    ProblemType(TP),
    SpeciesNames(0),
    PhaseNames(0),
    ElementNames(0),
    ElementAbundances(0),
    specifiedElementAbundances(false),
    specifiedBlockKmolSpecies(0),
    xmlStateInfoLevel(0),
    electrodeModelName(""),
    electrodeCapacityType(0),
    electrodeName(""),
    particleDiameter(1.0E-6),
    particleNumberToFollow(-1.0),
    electrodeGrossArea(-1.0),
    electrodeGrossDiameter(-1.0),
    electrodeGrossThickness(-1.0),
    porosity(-1.0),
    RxnExtTopLimit(-1.0),
    RxnExtBotLimit(-1.0),
    numExtraGlobalRxns(0),
    m_EGRList(0),
    nTotPhases(0),
    nTotSpecies(0),
    nTotElements(0),
    RelativeCapacityDischargedPerMole(-1.0),
    m_pl(0)
{
    m_BG = new ElectrodeBath();
    m_pl = new PhaseList();

    //m_EGRList = new EGRInput*[2];
    m_EGRList = (EGRInput**) mdpUtil::mdp_alloc_ptr_1(2);
    m_EGRList[0] = new EGRInput();
    m_EGRList[1] = 0;
}
//===============================================================================================
ELECTRODE_KEY_INPUT::ELECTRODE_KEY_INPUT(const ELECTRODE_KEY_INPUT &right) :
    printLvl_(right.printLvl_),
    commandFile_(""),
    CanteraFNSurface(""),
    NumberCanteraFiles(1),
    CanteraFileNames(0),
    Temperature(300.),
    Pressure(OneAtm),
    Vol(1.0),
    m_BG(0),
    MoleNumber(0),
    MoleFraction(0),
    PotentialPLPhases(0),
    PhaseInclude(0),
    ProblemType(TP),
    SpeciesNames(0),
    PhaseNames(0),
    ElementNames(0),
    ElementAbundances(0),
    specifiedElementAbundances(false),
    specifiedBlockKmolSpecies(0),
    xmlStateInfoLevel(0),
    electrodeModelName(""),
    electrodeCapacityType(0),
    electrodeName(""),
    particleDiameter(1.0E-6),
    particleNumberToFollow(-1.0),
    electrodeGrossArea(-1.0),
    electrodeGrossDiameter(-1.0),
    electrodeGrossThickness(-1.0),
    porosity(-1.0),
    RxnExtTopLimit(-1.0),
    RxnExtBotLimit(-1.0),
    numExtraGlobalRxns(0),
    m_EGRList(0),
    nTotPhases(0),
    nTotSpecies(0),
    nTotElements(0),
    RelativeCapacityDischargedPerMole(-1.0),
    m_pl(0)
{   
    ELECTRODE_KEY_INPUT::operator=(right);
}
//===============================================================================================
 ELECTRODE_KEY_INPUT&  ELECTRODE_KEY_INPUT::operator=(const ELECTRODE_KEY_INPUT& right)
 {
     if (this == &right) {
        return *this;
     }

     printLvl_                       = right.printLvl_;
     commandFile_                    = right.commandFile_;
     lastBlockEntryPtr_              = right.lastBlockEntryPtr_;
     CanteraFNSurface                = right.CanteraFNSurface;
     NumberCanteraFiles              = right.NumberCanteraFiles;

     if (CanteraFileNames) {
        for (int i = 0; CanteraFileNames[i] != 0; i++) {
            free(CanteraFileNames[i]);
        }
        free(CanteraFileNames);
     }
     if (!right.CanteraFileNames) {
       CanteraFileNames = 0;
     } else {
       CanteraFileNames = (char**) mdpUtil::mdp_alloc_ptr_1(NumberCanteraFiles +1);
       for (int i = 0; i <  NumberCanteraFiles     ; i++) {
          CanteraFileNames[i] = copy_string(right.CanteraFileNames[i]);
       }
     }

     Temperature                      = right.Temperature;
     Pressure                         = right.Pressure;
     Vol                              = right.Vol;
     nTotSpecies                      = right.nTotSpecies;
     nTotPhases                       = right.nTotPhases;
     nTotElements                     = right.nTotElements;
     if (m_BG) {
	 delete m_BG;
     } 
    
     m_BG = new ElectrodeBath(*right.m_BG);

     mdpUtil::mdp_realloc_dbl_1(&(m_BG->XmolPLSpecVec), nTotSpecies+2, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(m_BG->XmolPLSpecVec, (right.m_BG)->XmolPLSpecVec, nTotSpecies);

     mdpUtil::mdp_realloc_ptr_1((void ***) &(m_BG->XmolPLPhases), nTotPhases, 0);

     mdpUtil::mdp_realloc_dbl_1(&(m_BG->MolalitiesPLSpecVec), nTotSpecies+2, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(m_BG->MolalitiesPLSpecVec, (right.m_BG)->MolalitiesPLSpecVec, nTotSpecies);

     mdpUtil::mdp_realloc_ptr_1((void ***) &(m_BG->MolalitiesPLPhases), nTotPhases, 0);

     mdpUtil::mdp_realloc_dbl_1(&(m_BG->CapLeftCoeffSpecVec), nTotSpecies+2, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(m_BG->CapLeftCoeffSpecVec, (right.m_BG)->CapLeftCoeffSpecVec, nTotSpecies);

     mdpUtil::mdp_realloc_ptr_1((void ***) &(m_BG->CapLeftCoeffPhases), nTotPhases, 0);
     
     mdpUtil::mdp_realloc_dbl_1(&(m_BG->CapZeroDoDCoeffSpecVec), nTotSpecies+2, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(m_BG->CapZeroDoDCoeffSpecVec, (right.m_BG)->CapZeroDoDCoeffSpecVec, nTotSpecies);
   
     mdpUtil::mdp_realloc_ptr_1((void ***) &(m_BG->CapZeroDoDCoeffPhases), nTotPhases, 0);

     int nVolPhases = right.m_pl->nVolPhases();
     for (int iph = 0; iph < nVolPhases; iph++) {
	 int kstart =  right.m_pl->getGlobalSpeciesIndexVolPhaseIndex(iph);
	 m_BG->XmolPLPhases[iph]          = m_BG->XmolPLSpecVec + kstart;
	 m_BG->MolalitiesPLPhases[iph]    = m_BG->MolalitiesPLSpecVec + kstart;
	 m_BG->CapLeftCoeffPhases[iph]    = m_BG->CapLeftCoeffSpecVec + kstart;
	 m_BG->CapZeroDoDCoeffPhases[iph] = m_BG->CapZeroDoDCoeffSpecVec + kstart;
     }

     //---
   
     mdpUtil::mdp_realloc_dbl_1(&MoleNumber, nTotSpecies, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(MoleNumber, right.MoleNumber, nTotSpecies);

     mdpUtil::mdp_realloc_dbl_1(&MoleFraction, nTotSpecies, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(MoleFraction, right.MoleFraction, nTotSpecies);

     mdpUtil::mdp_realloc_dbl_1(&PotentialPLPhases, nTotPhases, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(PotentialPLPhases, right.PotentialPLPhases, nTotPhases);

     mdpUtil::mdp_realloc_int_1(&PhaseInclude, nTotPhases, 0);
     mdpUtil::mdp_copy_int_1(PhaseInclude, right.PhaseInclude, nTotPhases);
     
     ProblemType                       = right.ProblemType;

     if (SpeciesNames) {
        for (int i = 0; i < nTotSpecies; i++) {
            free(SpeciesNames[i]);
        }
        free(SpeciesNames);
     }
     SpeciesNames = (char**) mdp_alloc_ptr_1(nTotSpecies);
     for (int i = 0; i < nTotSpecies; i++) {
	SpeciesNames[i] = copy_string(right.SpeciesNames[i]);
     }

     if (PhaseNames) {
        for (int i = 0; i < nTotPhases; i++) {
            free(PhaseNames[i]);
        }
        free(PhaseNames);
     }
     PhaseNames = (char**) mdp_alloc_ptr_1(nTotPhases);
     for (int i = 0; i < nTotPhases; i++) {
	PhaseNames[i] = copy_string(right.PhaseNames[i]);
     }

     if (ElementNames) {
        for (int i = 0; i < nTotElements; i++) {
            free(ElementNames[i]);
        }
        free(ElementNames);
     }
     ElementNames = (char**) mdp_alloc_ptr_1(nTotElements);
     for (int i = 0; i < nTotElements; i++) {
	ElementNames[i] = copy_string(right.ElementNames[i]);
     }

     mdp_realloc_dbl_1(&PotentialPLPhases, nTotPhases, 0, 0.0);
     mdp_copy_dbl_1(PotentialPLPhases, right.PotentialPLPhases, nTotPhases);

     mdp_realloc_int_1(&PhaseInclude, nTotPhases, 0);
     mdp_copy_int_1(PhaseInclude, right.PhaseInclude, nTotPhases);
     
     ProblemType                       = right.ProblemType;

     if (SpeciesNames) {
        for (int i = 0; i < nTotSpecies; i++) {
            free(SpeciesNames[i]);
        }
        free(SpeciesNames);
     }
     SpeciesNames = (char**) mdp_alloc_ptr_1(nTotSpecies);
     for (int i = 0; i < nTotSpecies; i++) {
	SpeciesNames[i] = copy_string(right.SpeciesNames[i]);
     }

     if (PhaseNames) {
        for (int i = 0; i < nTotPhases; i++) {
            free(PhaseNames[i]);
        }
        free(PhaseNames);
     }
     PhaseNames = (char**) mdp_alloc_ptr_1(nTotPhases);
     for (int i = 0; i < nTotPhases; i++) {
	PhaseNames[i] = copy_string(right.PhaseNames[i]);
     }

     if (ElementNames) {
        for (int i = 0; i < nTotElements; i++) {
            free(ElementNames[i]);
        }
        free(ElementNames);
     }
     ElementNames = (char**) mdp_alloc_ptr_1(nTotElements);
     for (int i = 0; i < nTotElements; i++) {
	ElementNames[i] = copy_string(right.ElementNames[i]);
     }

     mdp_realloc_dbl_1(&ElementAbundances, nTotElements, 0, 0.0);
     mdp_copy_dbl_1(ElementAbundances, right.ElementAbundances, nTotElements);

     specifiedElementAbundances          = right.specifiedElementAbundances;
     specifiedBlockKmolSpecies           = right.specifiedBlockKmolSpecies;
     xmlStateInfoLevel                   = right.xmlStateInfoLevel;
     electrodeModelName                  = right.electrodeModelName;
     electrodeCapacityType               = right.electrodeCapacityType;
     electrodeName                       = right.electrodeName;
     particleDiameter                    = right.particleDiameter;
     particleNumberToFollow              = right.particleNumberToFollow;
     electrodeGrossArea                  = right.electrodeGrossArea;
     electrodeGrossDiameter              = right.electrodeGrossDiameter;
     electrodeGrossThickness             = right.electrodeGrossThickness;
     porosity                            = right.porosity;
     RxnExtTopLimit                      = right.RxnExtTopLimit;
     RxnExtBotLimit                      = right.RxnExtBotLimit;
     numExtraGlobalRxns                  = right.numExtraGlobalRxns;

     if (m_EGRList) {
        EGRInput** ptr;
        for (ptr = m_EGRList; *ptr != 0; ptr++) {
            delete *ptr;
        }
        free(m_EGRList);
    }
     m_EGRList = (EGRInput**) mdp_alloc_ptr_1(numExtraGlobalRxns+1);
     for (int i = 0; i < numExtraGlobalRxns  ; i++) {
	 m_EGRList[i] = new EGRInput(* right.m_EGRList[i]);
     }

     RelativeCapacityDischargedPerMole   = right.RelativeCapacityDischargedPerMole;

     if (m_pl) {
	 delete m_pl;
     }
     m_pl = new Cantera::PhaseList(*(right.m_pl));

   
     return *this;
 }
//===============================================================================================

ELECTRODE_KEY_INPUT::~ELECTRODE_KEY_INPUT()
{
    if (CanteraFileNames) {
        for (int i = 0; CanteraFileNames[i] != 0; i++) {
            free(CanteraFileNames[i]);
        }
        free(CanteraFileNames);
    }
    if (m_BG) {
      delete(m_BG);
    }
    m_BG=0;

    delete [] PhaseInclude;
    delete [] MoleNumber;
    delete [] MoleFraction;
    delete [] PotentialPLPhases;
    free(SpeciesNames);
    free(PhaseNames);
    free(ElementNames);
    delete [] ElementAbundances;

    if (m_EGRList) {

      
        EGRInput** ptr;
        for (ptr = m_EGRList; *ptr != 0; ptr++) {
            delete *ptr;
        }

        free(m_EGRList);
    }

    delete m_pl;
    m_pl = 0;
}

/****************************************************************************
 *
 *  In this routine, we fill in the fields in the ELECTRODE_KEY_INPUT
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
void ELECTRODE_KEY_INPUT::InitForInput(const Cantera::PhaseList* const pl)
{
    nTotPhases  = pl->nPhases();
    nTotSpecies = pl->nSpecies();
    nTotElements = pl->nElements();

    /*
     * Include all Phases by default
     */
    PhaseInclude = new int[nTotPhases];
    std::fill(PhaseInclude, PhaseInclude+nTotPhases, 1);
    MoleNumber = new double[nTotSpecies];
    std::fill(MoleNumber, MoleNumber+nTotSpecies, 0.0);
    MoleFraction = new double[nTotSpecies];
    std::fill(MoleFraction, MoleFraction+nTotSpecies, 0.0);
    PotentialPLPhases = new double[nTotPhases];
    std::fill(PotentialPLPhases, PotentialPLPhases+nTotPhases, 0.0);
    ElementAbundances = new double[nTotElements];
    std::fill(ElementAbundances, ElementAbundances+nTotElements, 0.0);


    SpeciesNames = mdp_alloc_VecFixedStrings(nTotSpecies,
                   MPEQUIL_MAX_NAME_LEN_P1);
    PhaseNames = mdp_alloc_VecFixedStrings(nTotPhases,
                                           MPEQUIL_MAX_NAME_LEN_P1);
    int kT = 0;
    for (int iphase = 0; iphase < nTotPhases; iphase++) {
        ThermoPhase* tPhase = &(pl->thermo(iphase));

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
    const Elements* eObj = pl->getGlobalElements();

    for (int e = 0; e < nTotElements; e++) {
        std::string eName = eObj->elementName(e);
        strncpy(ElementNames[e], eName.c_str(), MPEQUIL_MAX_NAME_LEN);
    }
}
//======================================================================================================================
EGRInput::EGRInput() :
    m_SS_KinSpeciesKindex(0),
    m_numElemReactions(0),
    m_ERSList(0)
{
    m_ERSList = (ERSSpec**) mdpUtil::mdp_alloc_ptr_1(2);
    m_ERSList[0] = new ERSSpec();
}
//======================================================================================================================
EGRInput::EGRInput(const EGRInput& right) :
    m_SS_KinSpeciesKindex(0),
    m_numElemReactions(0),
    m_ERSList(0)
{
    m_SS_KinSpeciesKindex = right.m_SS_KinSpeciesKindex;
    m_numElemReactions = right.m_numElemReactions;
    m_ERSList = (ERSSpec**) mdpUtil::mdp_alloc_ptr_1(m_numElemReactions + 2);
    for (int i = 0; i < m_numElemReactions; i++) {
        m_ERSList[i] = new ERSSpec(*(right.m_ERSList[i]));
    }
}
//======================================================================================================================
EGRInput& EGRInput::operator=(const EGRInput& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    if (m_ERSList) {
        ERSSpec** ptr;
        for (ptr = m_ERSList; *ptr != 0; ptr++) {
            delete *ptr;
        }
        mdpUtil::mdp_safe_free((void**) &m_ERSList);
    }

    m_SS_KinSpeciesKindex = right.m_SS_KinSpeciesKindex;
    m_numElemReactions = right.m_numElemReactions;
    m_ERSList = (ERSSpec**) mdpUtil::mdp_alloc_ptr_1(m_numElemReactions + 2);
    for (int i = 0; i < m_numElemReactions; i++) {
        m_ERSList[i] = new ERSSpec(*(right.m_ERSList[i]));
    }
    /*
     * Return the reference to the current object
     */
    return *this;

}
//======================================================================================================================
EGRInput::~EGRInput()
{
    if (m_ERSList) {
        ERSSpec** ptr;
        for (ptr = m_ERSList; *ptr != 0; ptr++) {
            delete *ptr;
        }
    }
    mdpUtil::mdp_safe_free((void**) &m_ERSList);
}
//======================================================================================================================
void ELECTRODE_KEY_INPUT::setup_input_pass1(BlockEntry* cf)
{

    /*
     *  Store the pointer to the BlockEntry structure that is used to parse the command file
     */
    lastBlockEntryPtr_ = cf;

    /* --------------------------------------------------------------
     * Electrode Name  = Identifying name to be used on printouts
     *    (required)
     *    default = ""
     */
    LE_OneStr* eName = new LE_OneStr("Electrode Name",  &(electrodeName),
                                     15, 1, 0, "ElectrodeName");
    cf->addLineEntry(eName);

    /* --------------------------------------------------------------
     * Electrode Type = ["anode", "cathode"]
     *    (required)
     *    default = anode
     */
    const char* etype[2] = {"anode", "cathode"};
    LE_PickList* leptype =
        new LE_PickList("Electrode Type", &(electrodeCapacityType), etype, 2, 1, "ElectrodeType");
    leptype->set_default(0);
    cf->addLineEntry(leptype);

    /*
     * Obtain the number of cantera files to be read
     */
    LE_OneInt* s1 = new LE_OneInt("Number of Cantera Files", &(NumberCanteraFiles), 0, "NumCanteraFiles");
    s1->set_default(1);
    cf->addLineEntry(s1);


    /* --------------------------------------------------------------
     * Electrode Model Name  = ["BaseType", "InfCapacity", "MP_RxnExtent",
     *			      "MultiPlateau_NoDiff", "SimpleDiff", "DiffTALE"
     *     		      "SimplePhaseChangeDiffusion"};
     *
     *  This must correspond to the name used by the factory method for instanteating
     *  Electrode objects.
     *
     *    (required)
     *    default = ""
     */
    LE_OneStr* eModelName = new LE_OneStr("Electrode Model Name",  &(electrodeModelName),
                                          1, 1, 1, "ElectrodeModelName");
    cf->addLineEntry(eModelName);

    BaseEntry::set_SkipUnknownEntries(true);
}
//========================================================================================================================
/****************************************************************************
 *
 */
void  ELECTRODE_KEY_INPUT::setup_input_pass2(BlockEntry* cf)
{
    /*
     *  Store the pointer to the BlockEntry structure that is used to parse the command file
     */
    lastBlockEntryPtr_ = cf;

    LineEntry* sle1 = 0;
    /*
     *  Get the input deck for
     *  Cantera description of the model.
     */
    LE_MultiCStr* s1 =
        new LE_MultiCStr("Cantera File Name", &(CanteraFileNames),
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

/************************************************************************
 *
 *
 * generic function Wrapper around new, in order to create a function
 * pointer for BE_MultiBlock
 */
void* getNewEGRInput(void* data_loc)
{
    void* ptr = new EGRInput();
    return ptr;
}

void* getNewERSSpec(void* data_loc)
{
    void* ptr = new ERSSpec();
    return ptr;
}


/****************************************************************************
 *
 */
void  ELECTRODE_KEY_INPUT::setup_input_pass3(BlockEntry* cf)
{
    /*
     *  Store the pointer to the BlockEntry structure that is used to parse the command file
     */
    lastBlockEntryPtr_ = cf;

    PhaseList* pl = m_pl;
    int iph;

    /* --------------------------------------------------------------
     * Temperature -  double  [optional]  [default = 300.0K]
     *
     *   The temperature is initialized to this value at the electrode level.
     *   This means that all ThermoPhases are initialized to this level as an initial condition.
     *   The temperature may be overrided by programs that use the Electrode object.
     */
    LE_OneDbl* d1 = new LE_OneDbl("Temperature", &(Temperature), 0, "Temperature");
    d1->set_default(300.);
    d1->set_limits(3000., 0.0);
    cf->addLineEntry(d1);

    /* --------------------------------------------------------------
     * Pressure -
     *
     * Configure the application Pressure
     */
    BE_UnitConversion* ucPres = new BE_UnitConversionPressure();
    LE_OneDblUnits* b5 = new LE_OneDblUnits("Pressure", &(Pressure), 0,
                                            "PO.Pressure", ucPres);
    b5->set_default(OneAtm);
    b5->set_limits(1.E20, 0.0);
    cf->addLineEntry(b5);

    /* ---------------------------------------------------------------------------
     * Particle diameter -
     *
     * Configure the characteristic particle size
     */
    BE_UnitConversion* ucLength = new BE_UnitConversionLength();
    LE_OneDblUnits* dpSize = new LE_OneDblUnits("Particle Diameter", &(particleDiameter),
            1, "particleDiameter", ucLength);
    dpSize->set_default(1.0E-6);
    dpSize->set_limits(1.0, 0.0);
    cf->addLineEntry(dpSize);

    /* ---------------------------------------------------------------------------
     * Particle Number to Follow -
     *
     * Configure the number of particles to follow -> we are setting up an
     * extrinsic measure for the size of the system
     */
    LE_OneDbl* dpNum = new LE_OneDbl("Particle Number to Follow", &(particleNumberToFollow),
                                     0, "particleNumbertoFollow");
    dpNum->set_default(-1.0);
    dpNum->set_limits(1.0E28, 1.0);
    cf->addLineEntry(dpNum);

    /* ---------------------------------------------------------------------------
     * Electrode Gross Area -
     *
     * Configure the area of the electrode -> we are setting up an
     * extrinsic measure for the size of the system
     */
    BE_UnitConversion* ucLength2 = new BE_UnitConversionLength();
    LE_OneDblUnits* eArea = new LE_OneDblUnits("Electrode Gross Area", &(electrodeGrossArea),
            0, "electrodeGrossArea", ucLength2);
    eArea->set_default(-1.0);
    eArea->set_limits(1.0E28, 1.0E-20);
    cf->addLineEntry(eArea);
    // If we specify particle number, we do not specify electrode area
    BI_Dependency* dep_eArea_dpNum = new BI_Dependency(dpNum, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    eArea->declareDependency(dep_eArea_dpNum);

    /* ---------------------------------------------------------------------------
     * Electrode Diameter -
     *
     * Configure the area of the electrode -> we are setting up an
     * extrinsic measure for the size of the system
     */
    BE_UnitConversion* ucLength2d = new BE_UnitConversionLength();
    LE_OneDblUnits* eDiameter = new LE_OneDblUnits("Electrode Gross Diameter", &(electrodeGrossDiameter),
            0, "electrodeGrossDiameter", ucLength2d);
    eDiameter->set_default(-1.0);
    eDiameter->set_limits(1.0E28, 1.0E-20);
    cf->addLineEntry(eDiameter);
    // If we specify particle number, we do not specify electrode diameter
    BI_Dependency* dep_eDia_dpNum = new BI_Dependency(dpNum, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    eDiameter->declareDependency(dep_eDia_dpNum);
    // If we specify electrode area, do not also specify electrode diameter
    BI_Dependency* dep_eDia_eArea = new BI_Dependency(eArea, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    eDiameter->declareDependency(dep_eDia_eArea);

    /* ---------------------------------------------------------------------------
     * Electrode Thickness -
     *
     * Configure the width of the electrode -> we are setting up an
     * extrinsic measure for the size of the system
     */
    BE_UnitConversion* ucLength3 = new BE_UnitConversionLength();
    LE_OneDblUnits* eThickness = new LE_OneDblUnits("Electrode Gross Thickness", &(electrodeGrossThickness),
            0, "electrodeGrossThickness", ucLength3);
    eThickness->set_default(-1.0);
    eThickness->set_limits(1.0E28, 1.0E-20);
    cf->addLineEntry(eThickness);

    //   If the 'Electrode Area' card is in the input deck, then the 'Electrode Thickness' card is mandatory.
    BI_Dependency* dep_eThickness_eArea = new BI_Dependency(eArea, BIDT_ENTRYPROCESSED, BIDRT_ONENUMTR);
    eThickness->declareDependency(dep_eThickness_eArea);

    //    If the "Particle Number to Follow" card is in the input deck, then the "Electrode Thickness" card
    //    appearance is an error
    BI_Dependency* depw = new BI_Dependency(dpNum, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    eThickness->declareDependency(depw);

    //    If the "Electrode Thickness" card is in the input deck, then the "Particle Number to Follow"
    //    card appearance is an error.
    BI_Dependency* depw_a = new BI_Dependency(eThickness, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    dpNum->declareDependency(depw_a);

    /* ---------------------------------------------------------------------------
     * Electrode porosity
     *
     * Configure the fraction of free area within the electrode.
     */
    LE_OneDbl* ePorosity = new LE_OneDbl("Electrode Porosity", &(porosity), 0, "electrodePorosity");
    ePorosity->set_default(-1.0);
    ePorosity->set_limits(1.0E28, 1.0E-20);
    cf->addLineEntry(ePorosity);

    //   If the 'Particle Number to Follow' card is in the input deck, then the 'Electrode Porosity'
    //   card is mandatory.
    //   Note that we don't really need porosity IF you bother to compute the volume of the electrode properly
    BI_Dependency* dep_ePorosity_dpNum = new BI_Dependency(dpNum, BIDT_ENTRYPROCESSED, BIDRT_ONENUMTR);
    ePorosity->declareDependency(dep_ePorosity_dpNum);



    //   If the 'Electrode Area' card is in the input deck, then the 'Electrode Porosity' card is prohibited.
    //BI_Dependency *dep_ePorosity_eArea = new BI_Dependency(eArea, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    //ePorosity->declareDependency(dep_ePorosity_eArea);

    //   If the 'Electrode Diameter' card is in the input deck, then the 'Electrode Porosity' card is prohibited.
    //BI_Dependency *dep_ePorosity_eDiameter = new BI_Dependency(eDiameter, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
    //ePorosity->declareDependency(dep_ePorosity_eDiameter);


    const char* cxml[4] = {"none", "minimal", "globalSteps", "intermediateSteps"};
    LE_PickList* lepXML = new LE_PickList("Level of XML State Information Created", &(xmlStateInfoLevel),
                                          cxml, 4, 0, "XMLLevel");
    lepXML->set_default(0);
    cf->addLineEntry(lepXML);



    /*
     *  Set up the Bath Gas BG object to receive input
     */
    int nVolPhases = pl->nVolPhases();

    ElectrodeBath& BG = *(m_BG);


    BG.XmolPLSpecVec = mdpUtil::mdp_alloc_dbl_1(nTotSpecies + 2, 0.0);
    BG.MolalitiesPLSpecVec = mdpUtil::mdp_alloc_dbl_1(nTotSpecies + 2, 0.0);
    BG.XmolPLPhases = (double**) mdpUtil::mdp_alloc_ptr_1(nVolPhases + pl->nSurPhases());
    BG.MolalitiesPLPhases = (double**) mdpUtil::mdp_alloc_ptr_1(nVolPhases + pl->nSurPhases());

    BG.CapLeftCoeffSpecVec = mdpUtil::mdp_alloc_dbl_1(nTotSpecies + 2, 0.0);
    BG.CapLeftCoeffPhases = (double**) mdpUtil::mdp_alloc_ptr_1(nVolPhases + pl->nSurPhases());
    BG.CapZeroDoDCoeffSpecVec = mdp_alloc_dbl_1(nTotSpecies + 2, 0.0);
    BG.CapZeroDoDCoeffPhases = (double**) mdp_alloc_ptr_1(nVolPhases + pl->nSurPhases());

    BG.PhaseMoles.resize(nVolPhases + pl->nSurPhases(), 0.0);
    BG.PhaseMass.resize(nVolPhases + pl->nSurPhases(), 0.0);

    for (iph = 0; iph < nVolPhases; iph++) {
        int kstart =  pl->getGlobalSpeciesIndexVolPhaseIndex(iph);
        BG.XmolPLPhases[iph] =   BG.XmolPLSpecVec + kstart;
        BG.MolalitiesPLPhases[iph] =  BG.MolalitiesPLSpecVec + kstart;
        BG.CapLeftCoeffPhases[iph] = BG.CapLeftCoeffSpecVec + kstart;
        BG.CapZeroDoDCoeffPhases[iph] = BG.CapZeroDoDCoeffSpecVec + kstart;
    }
    for (iph = 0; iph < pl->nSurPhases(); iph++) {
        int tph = iph + pl->nVolPhases();
        int kstart =  pl->getGlobalSpeciesIndexSurPhaseIndex(iph);
        BG.XmolPLPhases[tph] =   BG.XmolPLSpecVec + kstart;
        BG.MolalitiesPLPhases[tph] =  BG.MolalitiesPLSpecVec + kstart;
        BG.CapLeftCoeffPhases[tph] = BG.CapLeftCoeffSpecVec + kstart;
        BG.CapZeroDoDCoeffPhases[tph] = BG.CapZeroDoDCoeffSpecVec + kstart;
    }

    /*
     *  Specify a block for each Phase to receive inputs on composition
     *  and voltage
     */

    for (iph = 0; iph < pl->nVolPhases(); iph++) {
        string phaseBath = "Bath Specification for Phase ";
        ThermoPhase* tp = &(pl->volPhase(iph));
        string phaseNm = tp->name();
        int nSpecies = tp->nSpecies();
        phaseBath += phaseNm;
        /*
         *  create a section method description block and start writing
         *  line elements in it.
         */
        BlockEntry* bbathphase = new BlockEntry(phaseBath.c_str());
        cf->addSubBlock(bbathphase);
        int kstart = pl->getGlobalSpeciesIndexVolPhaseIndex(iph);

        /* --------------------------------------------------------------
         * BG.PotentialPLPhases
         *  Input the voltage for the phase
         */
        LE_OneDbl* iVolt = new LE_OneDbl("Voltage", &(PotentialPLPhases[iph]), 0,  "Voltage");
        iVolt->set_default(0.0);
        iVolt->set_limits(10., -10.);
        bbathphase->addLineEntry(iVolt);

        /* --------------------------------------------------------------
         * BG.PhaseMoles
         *  Input the number of moles for the phase
         */
        LE_OneDbl* iTMoles = new LE_OneDbl("Phase Moles", &(BG.PhaseMoles[iph]), 0, "PhaseMoles");
        iTMoles->set_default(0.0);
        iVolt->set_limits(0.0, 1.0E9);
        bbathphase->addLineEntry(iTMoles);
        // BI_Dependency * depimm_ig = new BI_Dependency(dpNum, BIDT_ENTRYPROCESSED, BIDRT_ZERONUMTIMESREQUIRED);
        // iTMoles->declareDependency(depimm_ig);
        //BI_Dependency *dep_phaseMoles_Porosity = new BI_Dependency(ePorosity, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
        //iTMoles->declareDependency(dep_phaseMoles_Porosity);

        /* --------------------------------------------------------------
         * BG.PhaseMass
         *  Input the mass for the phase
         */
        LE_OneDbl* iTMass =   new LE_OneDbl("Phase Mass", &(BG.PhaseMass[iph]), 0, "PhaseMass");
        iTMass->set_default(0.0);
        iVolt->set_limits(0.0, 1.0E9);
        bbathphase->addLineEntry(iTMass);
        //   If the 'Particle Number to Follow' card is in the input deck, then the 'Electrode Porosity'
        //   card is mandatory.
        //   QUESTION -- CAN THIS APPLY TO JUST ONE OF THE PHASES OR MUST IT APPLY TO ALL
        BI_Dependency* dep_phaseMass_phaseMoles = new BI_Dependency(iTMoles, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
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
        BE_MoleComp* bpmc = new BE_MoleComp("Bath Species Mole Fraction", &(BG.XmolPLPhases[iph]), 0,
                                            SpeciesNames+kstart, nSpecies, 0, "XBathmol");
        bpmc->generateDefLE();
        bbathphase->addSubBlock(bpmc);

        /* --------------------------------------------------------------
         * BG.CapLeftCoeffPhases -
         *
         *  Create a StrDbl block vector made out of the list of species.
	 *  We store the capacity coefficients for the species here
         */
        BE_StrDbl* pclc= new BE_StrDbl("Capacity Left Coefficients", &(BG.CapLeftCoeffPhases[iph]), 0, 0,
                                       SpeciesNames+kstart, nSpecies, 0, "CapLeftCoeff");
        pclc->generateDefLE();
        bbathphase->addSubBlock(pclc);

        /* --------------------------------------------------------------
         * BG.CapZeroDoDCoeffPhases -
         *
         * Create a PickList Line Element made out of the list of  species
         */
        BE_StrDbl* pczc= new BE_StrDbl("Capacity At Zero DoD Coefficients",
                                       &(BG.CapZeroDoDCoeffPhases[iph]), 0, 0, SpeciesNames+kstart,
                                       nSpecies, 0, "CapZeroDoDCoeff");
        pczc->generateDefLE();
        bbathphase->addSubBlock(pczc);



        /* --------------------------------------------------------------
         *
         */

        if (tp->activityConvention() == cAC_CONVENTION_MOLALITY) {
            MolalityVPSSTP* m_ptr = dynamic_cast<MolalityVPSSTP*>(tp);
            if (m_ptr == 0) {
                printf("Dynamic cast failed for some reason\n");
                exit(-1);
            }
            int indS = m_ptr->solventIndex();
            double mwS = m_ptr->molecularWeight(indS);
            BE_MolalityComp* bmolal =
                new BE_MolalityComp("Bath Species Molalities", &(BG.MolalitiesPLPhases[iph]), 0,
                                    SpeciesNames+kstart, nSpecies, indS, mwS,
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
        ThermoPhase* tp = &(pl->surPhase(iphS));
        string phaseNm = tp->name();
        int nSpecies = tp->nSpecies();
        phaseBath += phaseNm;
        /*
         *  create a section method description block and start writing
         *  line elements in it.
         */
        BlockEntry* bbathphase = new BlockEntry(phaseBath.c_str());
        cf->addSubBlock(bbathphase);
        int kstart =  pl->getGlobalSpeciesIndexSurPhaseIndex(iphS);



        /* --------------------------------------------------------------
         * BG.BathSpeciesMoleFractions -
         *
         * Create a PickList Line Element made out of the list of
         * species
         */
        BE_MoleComp* bpmc = new BE_MoleComp("Bath Species Mole Fraction",
                                            &(BG.XmolPLPhases[iph]), 0,
                                            SpeciesNames+kstart, nSpecies, 0, "XBathmol");
        bpmc->generateDefLE();
        bbathphase->addSubBlock(bpmc);

    }

    BlockEntry* rb = new BlockEntry("Reaction Extent Limits", 0);
    cf->addSubBlock(rb);
    LE_OneDbl* iRETop =  new LE_OneDbl("Reaction Extent Top Limit", &(RxnExtTopLimit), 1,
                                       "RxnExtTopLimit");
    iRETop->set_default(-1.0);
    rb->addLineEntry(iRETop);

    LE_OneDbl* iREBot =  new LE_OneDbl("Reaction Extent Bottom Limit", &(RxnExtBotLimit), 1,
                                       "RxnExtBotLimit");
    iREBot->set_default(-1.0);
    rb->addLineEntry(iREBot);



    /*
     *   SPECIFY OTHER WAYS TO INPUT THE CONCENTRATION, without
     *   breaking it down via phases.
     */


    /* ------------------------------------------------------------------
     * Block Input For initial number of moles of species
     *
     */
    BlockEntry* BESIM =
        new BE_StrDbl("Species Initial KMoles", &MoleNumber,
                      0, 0, SpeciesNames, nTotSpecies, 1,
                      "MoleNumber");
    cf->addSubBlock(BESIM);

    /* ------------------------------------------------------------------
     * Specify the initial relative amount of capacity discharged
     *
     */
    LE_OneDbl* iCapDis =  new LE_OneDbl("Capacity Discharged Per Initial Mole",
                                        &(RelativeCapacityDischargedPerMole), 0,
                                        "RelCapacityDischarged");
    iCapDis->set_default(-1.0);
    cf->addLineEntry(iCapDis);

    /* ------------------------------------------------------------------
     * Define a block that may occur multiples times
     *  The number of occurances of the block is .... PO.numExtraGlobalRxns
     *  All of the data defined in the block is
     *  in the pointer to the vector of structures .. PO.m_EGRList
     *
     */


    BlockEntry* sbEGR = new BE_MultiBlock("Extra Global Reaction",
                                          &(numExtraGlobalRxns),
                                          (void***) &(m_EGRList),
                                          getNewEGRInput, (void*) 0,
                                          0);
    cf->addSubBlock(sbEGR);

    /*
     * OK now define items in the block
     */
    EGRInput* egr_ptr = m_EGRList[0];

    /* --------------------------------------------------------------
     * EGR -> Special Species
     *
     * Create a PickList Line Element made out of the list of
     * species
     */
    LE_PickList* egrSS =
        new LE_PickList("Special Species", &(egr_ptr->m_SS_KinSpeciesKindex),
                        (const char**) SpeciesNames, nTotSpecies,
                        1, "Special_Species");
    egrSS->set_default(0);
    sbEGR->addLineEntry(egrSS);

    BlockEntry* sbERS = new BE_MultiBlock("Elementary Reaction Specification",
                                          &(egr_ptr->m_numElemReactions),
                                          (void***) &(egr_ptr->m_ERSList),
                                          getNewERSSpec, (void*) 0, 0);
    sbEGR->addSubBlock(sbERS);
    /*
     * OK now define items in the block
     */
    ERSSpec* ers_ptr = egr_ptr->m_ERSList[0];


    /* ------------------------------------------------------------------
     *  Define a reaction index
     */
    LE_OneInt* iRxnIndex =
        new LE_OneInt("Reaction Index", &(ers_ptr->m_reactionIndex), 1, "ReactionIndex");
    iRxnIndex->set_default(-1);
    sbERS->addLineEntry(iRxnIndex);

    LE_OneDbl* iRxnMult =
        new LE_OneDbl("Reaction Multiplier", &(ers_ptr->m_reactionMultiplier), 1,
                      "ReactionMultiplier");
    iRxnMult->set_default(0.0);
    sbERS->addLineEntry(iRxnMult);


    BaseEntry::set_SkipUnknownEntries(true);
}
//======================================================================================================================
/******************************************************************************
 *  Read the input file
 *
 *  printFlag 0 = no output
 *            1 = error message output
 *            2 = output -> processed lines
 *            3 = output -> original line and processed lines
 */
bool process_electrode_input(BlockEntry* cf, string fileName, int printFlag, int pass)
{

    if (printFlag > 2) {
        TKInput::set_tok_input_print_flag(1);
        BaseEntry::set_printProcessedLine(true);
    } else {
        TKInput::set_tok_input_print_flag(0);
        if (printFlag == 2) {
            BaseEntry::set_printProcessedLine(true);
        } else {
            BaseEntry::set_printProcessedLine(false);
        }
    }
    cf->ZeroLineCount();
    const TOKEN tok_in;
    TOKEN tok_out;
    FILE* ifp = fopen(fileName.c_str(), "r");
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
    } catch (BI_InputError& bi) {
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
setElectrodeBathSpeciesConditions(ThermoPhase& g, ELECTRODE_KEY_INPUT& EI, ElectrodeBath& BG, int iph, int printLvl)
{

    int nsp = g.nSpecies();


    /*
     * We need to set the state here. Note, mass fractions
     * don't matter, since we are only querying the
     * standard state of the species.
     */
    g.setState_TPX(EI.Temperature, EI.Pressure, BG.XmolPLPhases[iph]);

    /*
     * Set the electric potential
     */
    g.setElectricPotential(EI.PotentialPLPhases[iph]);

    /*
     * Print out table summarizing bath gas conditions
     */
    if (printLvl) {
        double* C = new double [nsp];
        double* act = new double [nsp];
        g.getConcentrations(C);
        g.getActivities(act);
        print_char('=', 100);
        cout << endl;
        dnt(0);
        cout << " SUMMARY OF SPECIES IN THE MECHANISM WITH A DESCRIPTION "
             << " OF BATH COMPOSITION: " << endl;
        cout << endl;
        dnt(1);
        cout << "Total pressure = " << EI.Pressure * 760. / OneAtm;
        cout << " torr " << endl;
        dnt(1);
        cout << "Temperature (where needed) = " << EI.Temperature
             << " Kelvin" << endl;
        dnt(1);
        cout << "Voltage (where needed) = " << EI.PotentialPLPhases[iph]
             << " Volts" << endl;


        dnt(4);
        cout << "Number       Name      Mole_fraction  Concentration (gmol/cm**3)   Activities";
        cout << endl;
        dnt(4);
        cout << "-----------------------------------------------------------------------------";
        cout << endl;
        string spN;
        for (int k = 0; k < (int) g.nSpecies(); k++) {
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



//======================================================================================================================
int
ELECTRODE_KEY_INPUT::electrode_input(std::string commandFile, BlockEntry* cf)
{
    int  retn = MPEQUIL_SUCCESS;

    if (printLvl_ > 1) {
        printf("\n");
        print_char('=', 80);
        printf("\n");
        print_char('=', 20);
        printf(" electrode_input: START OF PROBLEM STATEMENT ");
        print_char('=', 21);
        printf("\n");
        print_char('=', 80);
        printf("\n\n");
    }

    int printBIProclevel = 1;
    if (printLvl_ == 2 || printLvl_ == 4) {
        printBIProclevel = 2;
    }
    if (printLvl_ >= 5) {
        printBIProclevel = 3;
    }


    /*
     *  Save the name of the input file
     */
    commandFile_ = commandFile;

    /*
     *  Store the pointer to the BlockEntry structure that is used to parse the command file
     */
    lastBlockEntryPtr_ = cf;

    /*
     * Setup and process the input deck for first time.
     * -> Might have to print out the input and quit as well.
     */
    setup_input_pass1(cf);
    bool ok = process_electrode_input(cf, commandFile, printBIProclevel, 1);
    if (!ok) {
        return -1;
    }
    /**
     * Setup and process the input deck for second time.
     * -> Might have to print out the input and quit as well.
     */
    setup_input_pass2(cf);
    ok = process_electrode_input(cf, commandFile, printBIProclevel, 2);
    if (!ok) {
        return -1;
    }

    int ifiles = 0;
    for (; CanteraFileNames[ifiles] != 0; ifiles++) {
    }
    if (ifiles != NumberCanteraFiles) {
        printf("Number of requested files differ\n");
        exit(-1);
    }

    /*
     * Read in all of the phase specifications from the cantera
     * input files into the PhaseList structure.
     */
    PhaseList* pl = m_pl;
    std::string fn;
    bool surNotFound = true;
    for (int i = 0; i < NumberCanteraFiles; i++) {
        fn = CanteraFileNames[i];
        importAllCTMLIntoPhaseList(pl, fn);
        if (surNotFound && (pl->nSurPhases() > 0)) {
            surNotFound = false;
            //   pl->CanteraFNSurface = fn;
        }
    }
    /*
     * Setup internally for next pass through the input file.
     */
    InitForInput(pl);

    /*
     * Setup and process the input deck for third time
     * -> Might have to print out the input and quit as well.
     */
    setup_input_pass3(cf);
    //cf->print_usage();
    //exit(-1);

    /*
     * Possibly change the print level for the last
     */
    Electrode_Types_Enum ieos =  string_to_Electrode_Types_Enum(electrodeModelName);
    if (ieos != MP_RXNEXTENT_ET) {
        if (printLvl_ >= 3) {
            printBIProclevel = 3;
        }
    }
    /*
     * Process the third pass of the input file ->
     *   We read everything this time.
     */
    ok = process_electrode_input(cf, commandFile, printBIProclevel, 3);
    if (!ok) {
        return -1;
    }

    /*
     *  Post process the results
     *   -> here we figure out what the MoleNumber[] and MoleFraction[] vectors should be
     */
    post_input_pass3(cf);
    //cf->print_usage();
    //exit(-1);

    return retn;
}
//======================================================================================================================
int ELECTRODE_KEY_INPUT::electrode_input_child(std::string commandFile, BlockEntry* cf)
{
    cf->clear();
    /*
     *  Redo the base input
     */
    electrode_input(commandFile, cf);

    //cf->print_usage();
    /*
     *  Extend the input options for the child
     */
    setup_input_child1(cf);

    int printBIProclevel = 1;
    if (printLvl_ == 2 || printLvl_ == 4) {
        printBIProclevel = 2;
    } else if (printLvl_ >= 5) {
        printBIProclevel = 3;
    }
    bool ok = process_electrode_input(cf, commandFile, printBIProclevel, 4);
    if (!ok) {
        return -1;
    }
    post_input_child1(cf);

    //cf->print_usage();
    setup_input_child2(cf);
    //cf->print_usage();

    printBIProclevel = 1;
    if (printLvl_ == 2) {
        printBIProclevel = 2;
    } else if (printLvl_ >= 3) {
        printBIProclevel = 3;
    }
    ok = process_electrode_input(cf, commandFile, printBIProclevel, 5);
    if (!ok) {
        return -1;
    }
    post_input_child2(cf);

    setup_input_child3(cf);

    ok = process_electrode_input(cf, commandFile, printBIProclevel, 5);
    if (!ok) {
        return -1;
    }
    post_input_child3(cf);


    return 0;
}
//=========================================================================================
//   Initialize some of the fields in the ELECTRODE_KEY_INPUT structure by
//   post processing some of the data entries
/*
 *     Fields that are filled in and changed by this routine:
 *         MoleNumber[]
 *         MoleFraction[]
 *
 *  @param cf  Block entry used to parse the input file.
 *
 *   @return   Returns 0
 */
int ELECTRODE_KEY_INPUT::post_input_pass3(const BEInput::BlockEntry* cf)
{
    int nt, iph;
    /*
     *  Loop Over all phases in the PhaseList, locating any default specifications
     *  of the mole numbers and mole fractions that exist in the input file.
     */
    for (iph = 0; iph < m_pl->nPhases(); iph++) {
        int  kstart = m_pl->getGlobalSpeciesIndex(iph, 0);
        ThermoPhase* tphase = &(m_pl->thermo(iph));
        int nSpecies = tphase->nSpecies();
        // Find the name of the input block
        string phaseBath = "Bath Specification for Phase ";
        string phaseNm = tphase->name();
        phaseBath += phaseNm;
        double* molF = m_BG->XmolPLPhases[iph];
        bool molVecSpecified = false;
        BEInput::BlockEntry* pblock = cf->searchBlockEntry(phaseBath.c_str());
        int numTimes = 0;
        if (pblock) {
            numTimes = pblock->get_NumTimesProcessed();
        }
        /*
         *  If we have information about the specification of the phase, then process it
         */
        if (numTimes > 0) {
            BEInput::LineEntry*  pbpmoles = pblock->searchLineEntry("Phase Moles");
            int nt_pmoles = pbpmoles->get_NumTimesProcessed();
            double kmol = m_BG->PhaseMoles[iph];

            BEInput::BlockEntry* pbsmf = pblock->searchBlockEntry("Bath Species Mole Fraction");
            if (pbsmf) {
                if (pbsmf->get_NumTimesProcessed() > 0) {
                    molVecSpecified = true;
                }
            }
            if (molVecSpecified) {
                for (int k = 0; k < nSpecies; k++) {
                    MoleFraction[kstart + k] = molF[k];
                }
                tphase->setMoleFractions(molF);
            } else {
                tphase->getMoleFractions(MoleFraction + kstart);
            }
	    /*
	     *  Process the Phase Mass entry
	     */
            if (nt_pmoles == 0) {
                BEInput::LineEntry*  pbpm = pblock->searchLineEntry("Phase Mass");
                if (pbpm) {
                    nt = pbpm->get_NumTimesProcessed();
                    if (nt > 0) {
                        kmol = m_BG->PhaseMass[iph] / tphase->meanMolecularWeight();
                        m_BG->PhaseMoles[iph] = kmol;
                    }
                }
            }

            if (!molVecSpecified) {
                bool molalVecSpecified = false;
                if (pblock) {
                    BEInput::BlockEntry* pbsmm =
                        pblock->searchBlockEntry("Bath Species Molalities");
                    if (pbsmm) {
                        if (pbsmm->get_NumTimesProcessed() > 0) {
                            molalVecSpecified = true;
                        }
                    }
                }
                if (molalVecSpecified) {
                    MolalityVPSSTP* m_ptr = dynamic_cast<MolalityVPSSTP*>(tphase);
                    if (m_ptr == 0) {
                        printf("Dynamic cast failed for some reason\n");
                        exit(-1);
                    }
                    m_ptr->setState_TPM(Temperature, Pressure,
                                        m_BG->MolalitiesPLPhases[iph]);
                    m_ptr->getMoleFractions(m_BG->XmolPLPhases[iph]);
                    m_ptr->getMoleFractions(MoleFraction + kstart);
                }
            }
            /*
             * From the bath specification, get the total number of moles of the phase
             */
            double totalMoles = m_BG->PhaseMoles[iph];
	    /*
	     *  Set the electric potential of the phase 
	     */
            tphase->setElectricPotential(PotentialPLPhases[iph]);
            /*
             *  Setup the global MoleNumber array in the Electrode object
             */
            for (int k = 0; k < nSpecies; k++) {
                MoleNumber[kstart + k] = totalMoles * molF[k];
            }
        } else {
            /*
             *  We have no information about the specification of the phase.
             *  All we can do is set the mole numbers of all of the species to zero
             *  We also read in the mole fraction vector from the ThermoPhase object
             */
            for (int k = 0; k < nSpecies; k++) {
                MoleNumber[kstart + k] = 0.0;
            }
            tphase->getMoleFractions(MoleFraction + kstart);
        }

    }
    return 0;
}
//======================================================================================================================
// setup for child1
/*
 *  Virtual function that may get added onto
 *
 *  @param cf Pointer to the BlockEntry to be added onto by this routine
 */
void ELECTRODE_KEY_INPUT::setup_input_child1(BEInput::BlockEntry* cf)
{
}
//======================================================================================================================
// setup for child2
/*
 *  Virtual function that may get added onto
 *
 *  @param cf Pointer to the BlockEntry to be added onto by this routine
 */
void ELECTRODE_KEY_INPUT::setup_input_child2(BEInput::BlockEntry* cf)
{
}
//======================================================================================================================
// setup for child3
/*
 *  Virtual function that may get added onto
 *
 *  @param cf Pointer to the BlockEntry to be added onto by this routine
 */
void ELECTRODE_KEY_INPUT::setup_input_child3(BEInput::BlockEntry* cf)
{
}
//======================================================================================================================
// Post processing done after the input
/*
 *  Virtual function that may get added onto
 *
 *  @param cf Pointer to the BlockEntry to be added onto by this routine
 */
void ELECTRODE_KEY_INPUT::post_input_child1(BEInput::BlockEntry* cf)
{
}
//======================================================================================================================
// Post processing done after the input
/*
 *  Virtual function that may get added onto
 *
 *  @param cf Pointer to the BlockEntry to be added onto by this routine
 */
void ELECTRODE_KEY_INPUT::post_input_child2(BEInput::BlockEntry* cf)
{
}
//======================================================================================================================
// Post processing done after the input
/*
 *  Virtual function that may get added onto
 *
 *  @param cf Pointer to the BlockEntry to be added onto by this routine
 */
void ELECTRODE_KEY_INPUT::post_input_child3(BEInput::BlockEntry* cf)
{
}
//======================================================================================================================
int electrode_model_print(Cantera::Electrode* electrodeA,  ELECTRODE_KEY_INPUT* ei,
                          BEInput::BlockEntry* cf)
{

    int i, iph;


    ElectrodeBath* BG_ptr = ei->m_BG;

    /*
     *  Store the pointer to the BlockEntry structure that is used to parse the command file
     */
    ei->lastBlockEntryPtr_ = cf;

    /*
     * Load up the temperature and pressure
     */
    //electrodeA->T = ei->Temperature;
    //  electrodeA->Pres = ei->Pressure;

    electrodeA->setState_TP(ei->Temperature, ei->Pressure);
    /*
     *  Loop Over all phases in the PhaseList, adding these
     *  formally to the electrodeCell object.
     */
    for (iph = 0; iph < electrodeA->nPhases(); iph++) {

        ThermoPhase* tphase = &(electrodeA->thermo(iph));
        int nSpecies = tphase->nSpecies();

        // Find the name of the input block
        string phaseBath = "Bath Specification for Phase ";
        string phaseNm = tphase->name();
        phaseBath += phaseNm;

        /*
         * Set up the Bath Gas Conditions in each of
         * the ThermoPhase objects
         */


        BEInput::BlockEntry* pblock = cf->searchBlockEntry(phaseBath.c_str());
        /*
        bool molVecSpecified = false;
        if (pblock) {
            BEInput::BlockEntry* pbsmf = pblock->searchBlockEntry("Bath Species Mole Fraction");
            if (pbsmf) {
                if (pbsmf->get_NumTimesProcessed() > 0) {
                    molVecSpecified = true;
                }
            }
        }
        */

        double* molesSpecies = new double[nSpecies];

        bool molalVecSpecified = false;
        if (pblock) {
            BEInput::BlockEntry* pbsmm = pblock->searchBlockEntry("Bath Species Molalities");
            if (pbsmm) {
                if (pbsmm->get_NumTimesProcessed() > 0) {
                    molalVecSpecified = true;
                }
            }
        }
        if (molalVecSpecified) {
            MolalityVPSSTP* m_ptr = dynamic_cast<MolalityVPSSTP*>(tphase);
            if (m_ptr == 0) {
                printf("Dynamic cast failed for some reason\n");
                exit(-1);
            }
            m_ptr->setState_TPM(ei->Temperature, ei->Pressure,
                                BG_ptr->MolalitiesPLPhases[iph]);
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
        delete[] molesSpecies;

        //    tphase->getPartialMolarVolumes(&(electrodeA->VolPM[kstart]));
        //   tphase->getElectrochemPotentials(&(electrodeA->spElectroChemPot[kstart]));

        // transfer the voltages to the electrode struct
        double volts = tphase->electricPotential();
        electrodeA->setPhaseVoltage(iph, volts);
    }



    //            PRINT HEADER INFORMATION
    printf("\n");
    printf("            electrode_model_print()\n");
    printf("             Command file = %s\n", (ei->commandFile_).c_str());
    printf("\n");


    /*
     * Query whether an initial estimate has been made and then set iest.
     * Copy guess into vprobin
     */


    /*
     *          Printout the species information: PhaseID's and mole nums
     */
    printf("\n");
    print_char('-', 80);
    printf("\n");
    printf("             Phase IDs of species\n");
    printf("            species     phaseID        phaseName   ");
    printf(" Initial_Estimated_KMols\n");

    for (iph = 0; iph < electrodeA->nPhases(); iph++) {
        ThermoPhase* tphase = &(electrodeA->thermo(iph));
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
    printf("\n");
    print_char('-', 80);
    printf("\n");
    printf("             Information about phases\n");
    printf("  PhaseName    PhaseNum SingSpec GasPhase NumSpec");
    printf("  TMolesInert       TKmols\n");

    for (iph = 0; iph < electrodeA->nPhases(); iph++) {
        ThermoPhase* tphase = &(electrodeA->thermo(iph));
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
    printf("\n");
    print_char('-', 80);
    printf("\n");
    printf("             Information about Elements\n");
    printf("     ElementName  Abundance_Kmols\n");
    for (i = 0; i < electrodeA->nElements(); ++i) {
        string eName = electrodeA->elementName(i);
        printf("%12s ", eName.c_str());
        printf("  %11g \n", electrodeA->elementMoles(i));
    }


    printf("\n");
    print_char('=', 80);
    printf("\n");
    print_char('=', 20);
    printf(" mpequil: END OF PROBLEM STATEMENT ");
    print_char('=', 23);
    printf("\n");
    print_char('=', 80);
    printf("\n\n");


    return 0;
}


//=========================================================================================

