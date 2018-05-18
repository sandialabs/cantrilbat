/**
 *  @file Electrode_input.cpp
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode_input.h"
#include "Electrode_Factory.h"
#include "importPL.h"

#include "cantera/thermo/MolalityVPSSTP.h"

#include "BlockEntryGlobal.h"
#include "BE_MultiBlockVec.hpp"
#include "mdp_allo.h"

using std::cout; using std::endl;
using namespace BEInput;
using namespace ca_ab;
using namespace mdpUtil;
//==================================================================================================================================
//  Explicit instanteations of needed templated BE_MultiBlockVec classes

template class BEInput::BE_MultiBlockVec<Zuzax::EGRInput> ;
template class BEInput::BE_MultiBlockVec<Zuzax::ERSSpec> ;

//==================================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
EGRInput::EGRInput() :
    m_RSD_index(0),
    m_SS_KinSpeciesKindex(0),
    m_numElemReactions(0)
{
    m_ERSList.resize(1, nullptr);
    m_ERSList[0] = new ERSSpec();
}
//==================================================================================================================================
EGRInput::EGRInput(const EGRInput& right) :
    m_RSD_index(right.m_RSD_index),
    m_SS_KinSpeciesKindex( right.m_SS_KinSpeciesKindex ),
    m_numElemReactions(    right.m_numElemReactions )
{
    m_ERSList.resize(m_numElemReactions);
    for (size_t i = 0; i < static_cast<size_t>(m_numElemReactions); i++) {
        m_ERSList[i] = new ERSSpec(*(right.m_ERSList[i]));
    }
}
//==================================================================================================================================
EGRInput& EGRInput::operator=(const EGRInput& right)
{
    if (this == &right) {
        return *this;
    }
    for (ERSSpec* ptr : m_ERSList) {
        delete ptr;
    }
    m_RSD_index = right.m_RSD_index;
    m_SS_KinSpeciesKindex = right.m_SS_KinSpeciesKindex;
    m_numElemReactions = right.m_numElemReactions;
    m_ERSList.resize(m_numElemReactions, nullptr);
    for (size_t i = 0; i < static_cast<size_t>(m_numElemReactions); i++) {
        if (right.m_ERSList[i]) {
            m_ERSList[i] = new ERSSpec(*(right.m_ERSList[i]));
        }
    }
    return *this;
}
//==================================================================================================================================
EGRInput::~EGRInput()
{
    for (ERSSpec* ptr : m_ERSList) {
        delete ptr;
    }
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
ElectrodeBath::ElectrodeBath(PhaseList* pl) :
    m_pl(pl),
    XmolPLSpecVec(nullptr),
    MolalitiesPLSpecVec(nullptr),
    CapLeftCoeffPhases(nullptr),
    CapLeftCoeffSpecVec(nullptr),
    CapZeroDoDCoeffPhases(nullptr),
    CapZeroDoDCoeffSpecVec(nullptr)
{
}
//==================================================================================================================================
ElectrodeBath::ElectrodeBath(const ElectrodeBath &right) :
    ElectrodeBath(right.m_pl)
{
    operator=(right);
}
//==================================================================================================================================
 ElectrodeBath& ElectrodeBath::operator=(const ElectrodeBath& right)
 {
    if (this == &right) {
        return *this;
     }

    m_pl                         = right.m_pl;
    // Shallow pointer representation -> this is wrong and must be fixed up in parent routine.
    XmolPLSpecVec                = right.XmolPLSpecVec;
    MolalitiesPLSpecVec          = right.MolalitiesPLSpecVec;
    CapLeftCoeffPhases           = right.CapLeftCoeffPhases;
    CapLeftCoeffSpecVec          = right.CapLeftCoeffSpecVec;
    CapZeroDoDCoeffPhases        = right.CapZeroDoDCoeffPhases;
    CapZeroDoDCoeffSpecVec       = right.CapZeroDoDCoeffSpecVec;
 
    PhaseMoles                   = right.PhaseMoles;
    PhaseMass                    = right.PhaseMass;
 
    return *this;
 }
//==================================================================================================================================
ElectrodeBath::~ElectrodeBath() 
{
    mdpUtil::mdp_safe_free((void**) &XmolPLSpecVec);
    mdpUtil::mdp_safe_free((void**) &MolalitiesPLSpecVec);
    mdpUtil::mdp_safe_free((void**) &CapLeftCoeffPhases);
    mdpUtil::mdp_safe_free((void**) &CapLeftCoeffSpecVec);
    mdpUtil::mdp_safe_free((void**) &CapZeroDoDCoeffPhases);
    mdpUtil::mdp_safe_free((void**) &CapZeroDoDCoeffSpecVec);
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
OCV_Override_input::OCV_Override_input() :
    numTimes(0),
    surfacePhaseID(-1),
    OCVModel("Constant"),
    replacedSpeciesName(""),
    replacedGlobalSpeciesID(-1),
    replacedLocalSpeciesID(-1),
    replacedSpeciesPhaseID(-1),
    OCV_Format_(0),
    DoDSurrogateSpeciesName(""),
    MF_DoD_LocalSpeciesID(npos),
    rxnID(0),
    rxnID_deltaS(0),
    temperatureDerivType(0),
    temperatureBase(298.15),
    OCVTempDerivModel("Constant 0.0")
{
}
//==================================================================================================================================
OCV_Override_input::OCV_Override_input(const OCV_Override_input& right) :
    numTimes(right.numTimes),
    surfacePhaseID(right.surfacePhaseID),
    OCVModel(right.OCVModel),
    replacedSpeciesName(right.replacedSpeciesName),
    replacedGlobalSpeciesID(right.replacedGlobalSpeciesID),
    replacedLocalSpeciesID(right.replacedLocalSpeciesID),
    replacedSpeciesPhaseID(right.replacedSpeciesPhaseID),
    OCV_Format_(right.OCV_Format_),
    DoDSurrogateSpeciesName(right.DoDSurrogateSpeciesName),
    MF_DoD_LocalSpeciesID(right.MF_DoD_LocalSpeciesID),
    rxnID(right.rxnID),
    rxnID_deltaS(right.rxnID_deltaS),
    temperatureDerivType(right.temperatureDerivType),
    temperatureBase(right.temperatureBase),
    OCVTempDerivModel(right.OCVTempDerivModel)
{
}
//==================================================================================================================================
OCV_Override_input& OCV_Override_input::operator=(const OCV_Override_input& right)
{
    if (this == &right) {
       return *this;
    }
    numTimes                       = right.numTimes;
    surfacePhaseID                 = right.surfacePhaseID;
    OCVModel                       = right.OCVModel;
    replacedSpeciesName            = right.replacedSpeciesName;
    replacedGlobalSpeciesID        = right.replacedGlobalSpeciesID;
    replacedLocalSpeciesID         = right.replacedLocalSpeciesID;
    replacedSpeciesPhaseID         = right.replacedSpeciesPhaseID;
    OCV_Format_                    = right.OCV_Format_;
    DoDSurrogateSpeciesName        = right.DoDSurrogateSpeciesName;
    MF_DoD_LocalSpeciesID          = right.MF_DoD_LocalSpeciesID;
    rxnID                          = right.rxnID;
    rxnID_deltaS                   = right.rxnID_deltaS;
    temperatureDerivType           = right.temperatureDerivType;
    temperatureBase                = right.temperatureBase;
    OCVTempDerivModel              = right.OCVTempDerivModel;

    return *this;
}
//==================================================================================================================================
OCV_Override_input::~OCV_Override_input()
{
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
ELECTRODE_KEY_INPUT::ELECTRODE_KEY_INPUT(int printLvl) :
    m_pl(new PhaseList()),
    printLvl_(printLvl),
    commandFile_(""),
    NumberCanteraFiles(1),
    CanteraFileNames(0),
    Temperature(300.),
    Pressure(ZZCantera::OneAtm),
    m_BG(m_pl), 
    ProblemType(TP),
    SpeciesNames(0),
    PhaseNames(0),
    ElementNames(0),
    xmlStateInfoLevel(0),
    methodCapacityCalc(0), 
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
    nTotPhases(0),
    nTotSpecies(0),
    nTotElements(0),
    RelativeCapacityDischargedPerMole(-1.0),
    maxNumberSubGlobalTimeSteps(1000),
    relativeLocalToGlobalTimeStepMinimum(1.0E-3),
    extraRtolNonlinearSolver(1.0),
    doPolarizationAnalysis_(false)
{
    m_EGRList.resize(1, nullptr);
    m_EGRList[0] = new EGRInput();
}
//==================================================================================================================================
ELECTRODE_KEY_INPUT::ELECTRODE_KEY_INPUT(const ELECTRODE_KEY_INPUT &right) :
    m_pl(right.m_pl),
    printLvl_(right.printLvl_),
    commandFile_(""),
    NumberCanteraFiles(1),
    CanteraFileNames(0),
    Temperature(300.),
    Pressure(OneAtm),
    ProblemType(TP),
    SpeciesNames(0),
    PhaseNames(0),
    ElementNames(0),
    xmlStateInfoLevel(0),
    methodCapacityCalc(0),
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
    nTotPhases(0),
    nTotSpecies(0),
    nTotElements(0),
    RelativeCapacityDischargedPerMole(-1.0),
    maxNumberSubGlobalTimeSteps(1000),
    relativeLocalToGlobalTimeStepMinimum(1.0E-3),
    extraRtolNonlinearSolver(1.0),
    doPolarizationAnalysis_(false)
{   
    ELECTRODE_KEY_INPUT::operator=(right);
}
//==================================================================================================================================
ELECTRODE_KEY_INPUT&  ELECTRODE_KEY_INPUT::operator=(const ELECTRODE_KEY_INPUT& right)
{
     if (this == &right) {
        return *this;
     }

     if (m_pl && (m_pl != right.m_pl)) {
	 delete m_pl;
     }
     m_pl = new ZZCantera::PhaseList(*(right.m_pl));

     printLvl_                       = right.printLvl_;
     commandFile_                    = right.commandFile_;
     lastBlockEntryPtr_              = right.lastBlockEntryPtr_;
     NumberCanteraFiles              = right.NumberCanteraFiles;

     if (CanteraFileNames) {
        for (size_t i = 0; CanteraFileNames[i] != 0; i++) {
            free(CanteraFileNames[i]);
        }
        free(CanteraFileNames);
     }
     if (!right.CanteraFileNames) {
         CanteraFileNames = 0;
     } else {
         CanteraFileNames = (char**) mdpUtil::mdp_alloc_ptr_1(NumberCanteraFiles +1);
         for (int i = 0; i < NumberCanteraFiles; i++) {
             CanteraFileNames[i] = TKInput::copy_string(right.CanteraFileNames[i]);
         }
     }

     Temperature                      = right.Temperature;
     Pressure                         = right.Pressure;
     nTotSpecies                      = right.nTotSpecies;
     nTotPhases                       = right.nTotPhases;
     nTotElements                     = right.nTotElements;
    
     m_BG = right.m_BG;

     mdpUtil::mdp_realloc_dbl_1(&(m_BG.XmolPLSpecVec), nTotSpecies+2, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(m_BG.XmolPLSpecVec, (right.m_BG).XmolPLSpecVec, nTotSpecies);

     mdpUtil::mdp_realloc_dbl_1(&(m_BG.MolalitiesPLSpecVec), nTotSpecies+2, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(m_BG.MolalitiesPLSpecVec, (right.m_BG).MolalitiesPLSpecVec, nTotSpecies);

     mdpUtil::mdp_realloc_dbl_1(&(m_BG.CapLeftCoeffSpecVec), nTotSpecies+2, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(m_BG.CapLeftCoeffSpecVec, (right.m_BG).CapLeftCoeffSpecVec, nTotSpecies);

     mdpUtil::mdp_realloc_ptr_1((void ***) &(m_BG.CapLeftCoeffPhases), nTotPhases, 0);
     
     mdpUtil::mdp_realloc_dbl_1(&(m_BG.CapZeroDoDCoeffSpecVec), nTotSpecies+2, 0, 0.0);
     mdpUtil::mdp_copy_dbl_1(m_BG.CapZeroDoDCoeffSpecVec, (right.m_BG).CapZeroDoDCoeffSpecVec, nTotSpecies);
   
     mdpUtil::mdp_realloc_ptr_1((void ***) &(m_BG.CapZeroDoDCoeffPhases), nTotPhases, 0);

     size_t nVolPhases = right.m_pl->nVolPhases();
     for (size_t iph = 0; iph < static_cast<size_t>(nVolPhases); iph++) {
	 size_t kstart =  right.m_pl->globalSpeciesIndexVolPhaseIndex(iph);
	 m_BG.CapLeftCoeffPhases[iph]    = m_BG.CapLeftCoeffSpecVec + kstart;
	 m_BG.CapZeroDoDCoeffPhases[iph] = m_BG.CapZeroDoDCoeffSpecVec + kstart;
     }

     //---
   
     MoleNumber = right.MoleNumber;
     MoleFraction = right.MoleFraction;
     PotentialPLPhases = right.PotentialPLPhases;

     ProblemType = right.ProblemType;

     if (SpeciesNames) {
        for (size_t i = 0; i < nTotSpecies; i++) {
            free(SpeciesNames[i]);
        }
        free(SpeciesNames);
     }
     SpeciesNames = (char**) mdp_alloc_ptr_1(nTotSpecies);
     for (size_t i = 0; i < nTotSpecies; i++) {
	SpeciesNames[i] = TKInput::copy_string(right.SpeciesNames[i]);
     }

     if (PhaseNames) {
        for (size_t i = 0; i < nTotPhases; i++) {
            free(PhaseNames[i]);
        }
        free(PhaseNames);
     }
     PhaseNames = (char**) mdp_alloc_ptr_1(nTotPhases);
     for (size_t i = 0; i < nTotPhases; i++) {
	PhaseNames[i] = TKInput::copy_string(right.PhaseNames[i]);
     }

     if (ElementNames) {
        for (size_t i = 0; i < static_cast<size_t>(nTotElements); i++) {
            free(ElementNames[i]);
        }
        free(ElementNames);
     }
     ElementNames = (char**) mdp_alloc_ptr_1(nTotElements);
     for (size_t i = 0; i < nTotElements; i++) {
	ElementNames[i] = TKInput::copy_string(right.ElementNames[i]);
     }

     xmlStateInfoLevel                   = right.xmlStateInfoLevel;
     methodCapacityCalc                  = right.methodCapacityCalc;
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

     for (EGRInput* ep : m_EGRList) {
         delete ep;
     }
     m_EGRList.resize(numExtraGlobalRxns, nullptr);
     for (int i = 0; i < numExtraGlobalRxns  ; i++) {
	 m_EGRList[i] = new EGRInput(*(right.m_EGRList[i]));
     }

     RelativeCapacityDischargedPerMole   = right.RelativeCapacityDischargedPerMole;

     maxNumberSubGlobalTimeSteps         = right.maxNumberSubGlobalTimeSteps;
     relativeLocalToGlobalTimeStepMinimum = right.relativeLocalToGlobalTimeStepMinimum;
     extraRtolNonlinearSolver = right.extraRtolNonlinearSolver;
    doPolarizationAnalysis_ = right.doPolarizationAnalysis_;
   
     return *this;
}
//==================================================================================================================================
ELECTRODE_KEY_INPUT::~ELECTRODE_KEY_INPUT()
{
    if (CanteraFileNames) {
        for (int i = 0; CanteraFileNames[i] != 0; i++) {
            free(CanteraFileNames[i]);
        }
        free(CanteraFileNames);
    }

    free(SpeciesNames);
    free(PhaseNames);
    free(ElementNames);

    for (size_t j = 0; j < OCVoverride_ptrList.size(); j++) {
         delete OCVoverride_ptrList[j];
    }
    for (EGRInput* ep : m_EGRList) {
        delete ep;
    }
    delete m_pl;
    m_pl = 0;
}
//==================================================================================================================================
/*
 *  In this routine, we fill in the fields in the ELECTRODE_KEY_INPUT
 *  structure that will be needed to parse the keyword input.
 *  We also size the vectors in the structure appropriately, now we know the extent of the problem.
 *  The information to do this has already been gathered into the PhaseList structure
 *
 *  These are:
 *         nTotPhases
 *         nTotSpecies
 *         nTotElements
 *         SpeciesNames
 *         PhaseNames
 *         ElementNames
 *
 *  The sized fields include:
 *         MoleNumber
 *         MoleNumberIG
 */
void ELECTRODE_KEY_INPUT::InitForInput(const ZZCantera::PhaseList* const pl)
{
    nTotPhases  = pl->nPhases();
    nTotSpecies = pl->nGlobalSpecies();
    nTotElements = pl->nElements();

    MoleNumber.resize(nTotSpecies, 0.0);
    MoleFraction.resize(nTotSpecies, 0.0);
  
    PotentialPLPhases.resize(nTotPhases, 0.0);

    SpeciesNames = mdp_alloc_VecFixedStrings(nTotSpecies, MPEQUIL_MAX_NAME_LEN_P1);
    PhaseNames = mdp_alloc_VecFixedStrings(nTotPhases, MPEQUIL_MAX_NAME_LEN_P1);
    size_t kT = 0;
    for (size_t iphase = 0; iphase < static_cast<size_t>(nTotPhases); iphase++) {
        ThermoPhase* tPhase = &(pl->thermo(iphase));

        std::string id = tPhase->id();
        strncpy(PhaseNames[iphase], id.c_str(), MPEQUIL_MAX_NAME_LEN);
        size_t nspecies = tPhase->nSpecies();
        tPhase->getMoleFractions(MoleFraction.data() + kT);
        for (size_t k = 0; k < nspecies; k++) {
            std::string sname = tPhase->speciesName(k);
            strncpy(SpeciesNames[kT], sname.c_str(), MPEQUIL_MAX_NAME_LEN);
            kT++;
        }

    }

    ElementNames = mdp_alloc_VecFixedStrings(nTotElements, MPEQUIL_MAX_NAME_LEN_P1);
    const Elements* eObj = pl->globalElements();

    for (size_t e = 0; e < static_cast<size_t>(nTotElements); e++) {
        std::string eName = eObj->elementName(e);
        strncpy(ElementNames[e], eName.c_str(), MPEQUIL_MAX_NAME_LEN);
    }
}
//======================================================================================================================
//======================================================================================================================
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
    LE_OneStr* eName = new LE_OneStr("Electrode Name",  &(electrodeName), 15, 1, 0, "ElectrodeName");
    cf->addLineEntry(eName);

    /* --------------------------------------------------------------
     * Electrode Type = ["anode", "cathode"]
     *    (required)
     *    default = anode
     */
    const char* etype[2] = {"anode", "cathode"};
    LE_PickList* leptype = new LE_PickList("Electrode Type", &(electrodeCapacityType), etype, 2, 1, "ElectrodeType");
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
     *    (required)
     *    default = ""
     */
    LE_OneStr* eModelName = new LE_OneStr("Electrode Model Name",  &(electrodeModelName),
                                          1, 1, 1, "ElectrodeModelName");
    cf->addLineEntry(eModelName);

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
    LE_OneDblUnits* b5 = new LE_OneDblUnits("Pressure", &(Pressure), 0, "PO.Pressure", ucPres);
    b5->set_default(OneAtm);
    b5->set_limits(1.E20, 0.0);
    cf->addLineEntry(b5);

    /* ------------------------------------------------------------------
     * Maximum number of subGlobal time steps = [int]
     *  
     *     defaults to 1000.
     */
    LE_OneInt* m1 = new LE_OneInt("Maximum number of Subglobal time steps", &(maxNumberSubGlobalTimeSteps), 0,
				  "maxNumberSubGlobalTimeSteps");
    m1->set_default(1000);
    m1->set_limits(10000, 1);
    cf->addLineEntry(m1);

    BaseEntry::set_SkipUnknownEntries(3);

    /* --------------------------------------------------------------
     * Relative Local To Global Time Step Minimum -  double  [optional]  [default = 1.0E-3 ]
     *
     *   This modifies the lowest time step below which time step truncation errors aren't checked.
     *   Also the interval will start out at the minimum time step suggested by this value.
     *   The time step may decrease below this limit if there are issues with convergence.
     */
    LE_OneDbl* rl1 = new LE_OneDbl("Relative Local To Global Time Step Minimum", &(relativeLocalToGlobalTimeStepMinimum), 
                                  0, "relativeLocalToGlobalTimeStepMinimum");
    rl1->set_default(1.0E-3);
    rl1->set_limits(1.0, 0.0);
    cf->addLineEntry(rl1);

    /* --------------------------------------------------------------
     *   Extra Rtol Accuracy for Nonlinear Solve - double (default 1.0) (not-required)
     *
     *     Usually the nonlinear solve is only solved to the level of the time step truncation error.
     *     However, you can add additional rtol accuracy to the nonlinear solver (both on the
     *     residual tolerance and on the solution delta tolerance) by setting a multiplicative constant
     *     below 1.0 (1.0E-3 works well).
     */
    LE_OneDbl* rl2 = new LE_OneDbl("Extra Rtol Accuracy for Nonlinear Solve", &(extraRtolNonlinearSolver), 
                                  0, "extraRtolNonlinearSolver");
    rl2->set_default(1.0);
    rl2->set_limits(5.0, 1.0E-8);
    cf->addLineEntry(rl2);

}
//========================================================================================================================
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
    LE_MultiCStr* s1 = new LE_MultiCStr("Cantera File Name", &(CanteraFileNames), 1, 1,  0, "CanteraFileNames");
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
    BaseEntry::set_SkipUnknownEntries(3);
}
//========================================================================================================================
void  ELECTRODE_KEY_INPUT::setup_input_pass3(BlockEntry* cf)
{
    /*
     *  Store the pointer to the BlockEntry structure that is used to parse the command file
     */
    lastBlockEntryPtr_ = cf;

    PhaseList* pl = m_pl;

    /* ---------------------------------------------------------------------------
     * Particle diameter -
     *
     * Configure the characteristic particle size
     */
    BE_UnitConversion* ucLength = new BE_UnitConversionLength();
    LE_OneDblUnits* dpSize = new LE_OneDblUnits("Particle Diameter", &(particleDiameter), 1, "particleDiameter", ucLength);
    dpSize->set_default(1.0E-6);
    dpSize->set_limits(1.0, 0.0);
    cf->addLineEntry(dpSize);

    /* ---------------------------------------------------------------------------
     * Particle Number to Follow -
     *
     * Configure the number of particles to follow -> we are setting up an
     * extrinsic measure for the size of the system
     */
    LE_OneDbl* dpNum = new LE_OneDbl("Particle Number to Follow", &(particleNumberToFollow), 0, "particleNumbertoFollow");
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
    LE_OneDblUnits* eArea = new LE_OneDblUnits("Electrode Gross Area", &(electrodeGrossArea), 0, "electrodeGrossArea", ucLength2);
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
     * Electrode porosity = double
     *            (optional) (default = -1.0)
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

    /* ---------------------------------------------------------------------------
     *  Method for Electrode Capacity Calculation  = PickList [ Default,  Model, Coeffs , Coefficients]
     *                                                  (default = Default)
     *                                                  (optional)
     */
    const char* cMethCap[4] = {"Default", "Model", "Coeffs" , "Coefficients"};
    LE_PickList* lepMethCap = new LE_PickList("Method for Electrode Capacity Calculation", &(methodCapacityCalc),
                                              cMethCap, 4, 0, "methodCapacityCalc");
    lepMethCap->set_default(0);
    cf->addLineEntry(lepMethCap);

    /*
     *  Set up the Bath Gas BG object to receive input
     */
    int nVolPhases = pl->nVolPhases();


    m_BG.XmolPLSpecVec = mdpUtil::mdp_alloc_dbl_1(nTotSpecies + 2, 0.0);
    m_BG.MolalitiesPLSpecVec = mdpUtil::mdp_alloc_dbl_1(nTotSpecies + 2, 0.0);

    m_BG.CapLeftCoeffSpecVec = mdpUtil::mdp_alloc_dbl_1(nTotSpecies + 2, 0.0);
    m_BG.CapLeftCoeffPhases = (double**) mdpUtil::mdp_alloc_ptr_1(nVolPhases + pl->nSurPhases());
    m_BG.CapZeroDoDCoeffSpecVec = mdp_alloc_dbl_1(nTotSpecies + 2, 0.0);
    m_BG.CapZeroDoDCoeffPhases = (double**) mdp_alloc_ptr_1(nVolPhases + pl->nSurPhases());

    m_BG.PhaseMoles.resize(nVolPhases + pl->nSurPhases(), 0.0);
    m_BG.PhaseMass.resize(nVolPhases + pl->nSurPhases(), 0.0);

    for (size_t iph = 0; iph < (size_t) nVolPhases; iph++) {
        size_t kstart = pl->globalSpeciesIndexVolPhaseIndex(iph);
        m_BG.CapLeftCoeffPhases[iph] = m_BG.CapLeftCoeffSpecVec + kstart;
        m_BG.CapZeroDoDCoeffPhases[iph] = m_BG.CapZeroDoDCoeffSpecVec + kstart;
    }
    for (size_t iph = 0; iph < pl->nSurPhases(); iph++) {
        size_t tph = iph + pl->nVolPhases();
        size_t kstart = pl->globalSpeciesIndexSurPhaseIndex(iph);
        m_BG.CapLeftCoeffPhases[tph] = m_BG.CapLeftCoeffSpecVec + kstart;
        m_BG.CapZeroDoDCoeffPhases[tph] = m_BG.CapZeroDoDCoeffSpecVec + kstart;
    }

    // ---------------------------------------------------------------------------------
    /*
     *  Specify a block for each Phase to receive inputs on composition
     *  and voltage
     */
    for (size_t iph = 0; iph < pl->nVolPhases(); iph++) {
        std::string phaseBath = "Bath Specification for Phase ";
        ThermoPhase* tp = &(pl->volPhase(iph));
        std::string phaseNm = tp->name();
        int nSpecies = tp->nSpecies();
        phaseBath += phaseNm;
        /*
         *  create a section method description block and start writing
         *  line elements in it.
         */
        BlockEntry* bbathphase = new BlockEntry(phaseBath.c_str());
        cf->addSubBlock(bbathphase);
        size_t kstart = pl->globalSpeciesIndexVolPhaseIndex(iph);

        /* --------------------------------------------------------------
         * BG.PotentialPLPhases[iph]
         * 
         *        Voltage = double   [optional] [default = 0.0]
         *        Input the voltage for the phase
         */
        LE_OneDbl* iVolt = new LE_OneDbl("Voltage", &(PotentialPLPhases[iph]), 0,  "Voltage");
        iVolt->set_default(0.0);
        iVolt->set_limits(10., -10.);
        bbathphase->addLineEntry(iVolt);

        /* --------------------------------------------------------------
         * BG.PhaseMoles
         *  Input the number of moles for the phase
         */
        LE_OneDbl* iTMoles = new LE_OneDbl("Phase Moles", &(m_BG.PhaseMoles[iph]), 0, "PhaseMoles");
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
        LE_OneDbl* iTMass = new LE_OneDbl("Phase Mass", &(m_BG.PhaseMass[iph]), 0, "PhaseMass");
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
         * Create a PickList Line Element made out of the list of species
         */
        BE_MoleComp* bpmc = new BE_MoleComp("Bath Species Mole Fraction", &(m_BG.XmolPLSpecVec[kstart]), 0,
                                            SpeciesNames+kstart, nSpecies, 0, "XBathmol");
        bpmc->generateDefLE();
        bbathphase->addSubBlock(bpmc);

        /* --------------------------------------------------------------
         * BG.CapLeftCoeffPhases -
         *
         *  Create a StrDbl block vector made out of the list of species.
	 *  We store the capacity coefficients for the species here
         */
        BE_StrDbl* pclc= new BE_StrDbl("Capacity Left Coefficients", &(m_BG.CapLeftCoeffPhases[iph]), 0, 0,
                                       SpeciesNames+kstart, nSpecies, 0, "CapLeftCoeff");
        pclc->generateDefLE();
        bbathphase->addSubBlock(pclc);

        /* --------------------------------------------------------------
         * BG.CapZeroDoDCoeffPhases -
         *
         * Create a PickList Line Element made out of the list of  species
         */
        BE_StrDbl* pczc= new BE_StrDbl("Capacity At Zero DoD Coefficients",
                                       &(m_BG.CapZeroDoDCoeffPhases[iph]), 0, 0, SpeciesNames+kstart,
                                       nSpecies, 0, "CapZeroDoDCoeff");
        pczc->generateDefLE();
        bbathphase->addSubBlock(pczc);

        /* --------------------------------------------------------------
         * BG.BG.MolalitiesPLSpecVec -
         *
         *  Create a Molalities Block for entering molalities directly
         */
        if (tp->activityConvention() == cAC_CONVENTION_MOLALITY) {
            ZZCantera::MolalityVPSSTP* m_ptr = dynamic_cast<ZZCantera::MolalityVPSSTP*>(tp);
            if (m_ptr == 0) {
                printf("Dynamic cast failed for some reason\n");
                exit(-1);
            }
            int indS = m_ptr->solventIndex();
            double mwS = m_ptr->molecularWeight(indS);
            size_t kstart = pl->globalSpeciesIndexVolPhaseIndex(iph);
            BE_MolalityComp* bmolal = new BE_MolalityComp("Bath Species Molalities", &(m_BG.MolalitiesPLSpecVec[kstart]), 0,
                                                          SpeciesNames+kstart, nSpecies, indS, mwS, "MolalitiesBath");
            bbathphase->addSubBlock(bmolal);
        }
    }
    // ---------------------------------------------------------------------------------
    /*
     *  Specify a block for each surface phase to receive inputs on composition
     *  and voltage
     */
    for (size_t iphS = 0; iphS < pl->nSurPhases(); iphS++) {
        std::string phaseBath = "Bath Specification for Phase ";
        ThermoPhase* tp = &(pl->surPhase(iphS));
        std::string phaseNm = tp->name();
        size_t nSpecies = tp->nSpecies();
        phaseBath += phaseNm;
        /*
         *  create a section method description block and start writing
         *  line elements in it.
         */
        BlockEntry* bbathphase = new BlockEntry(phaseBath.c_str());
        cf->addSubBlock(bbathphase);
        size_t kstart = pl->globalSpeciesIndexSurPhaseIndex(iphS);

        /* --------------------------------------------------------------
         * BG.BathSpeciesMoleFractions -
         *
         * Create a PickList Line Element made out of the list of species
         */
        BE_MoleComp* bpmc = new BE_MoleComp("Bath Species Mole Fraction", &(m_BG.XmolPLSpecVec[kstart]), 0,
                                            SpeciesNames+kstart, nSpecies, 0, "XBathmol");
        bpmc->generateDefLE();
        bbathphase->addSubBlock(bpmc);
    }
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
    //  
    //  Specify a block for each surface phase to receive inputs on OCV Overrides
    //
    OCVoverride_ptrList.resize(pl->nSurPhases(), 0);
    for (size_t iphS = 0; iphS < (size_t)pl->nSurPhases(); iphS++) {
        std::string phaseBath = "Open Circuit Potential Override for interface ";
        ThermoPhase* tp = &(pl->surPhase(iphS));
        std::string phaseNm = tp->name();
        phaseBath += phaseNm;

        OCV_Override_input* ocv_input_ptr = OCVoverride_ptrList[iphS];
        if (ocv_input_ptr == 0) {
            OCVoverride_ptrList[iphS] = new OCV_Override_input();
            ocv_input_ptr = OCVoverride_ptrList[iphS]; 
            ocv_input_ptr->surfacePhaseID = iphS;
        }

        // 
        //   Create a section method description block and start writing
        //   line elements in it.
        // 
        BlockEntry* bOCVoverride = new BlockEntry(phaseBath.c_str());
        cf->addSubBlock(bOCVoverride);

        /* --------------------------------------------------------------
         *      ocv_input.OCVModel = string   (required)
         *
         *  Input a string name that will be input into a factory routine to pick the OCV model
         *  from a list of options. The first token will be interpreted as the string name.
         *  All other tokens will be interpreted as doubles and input as parameters for the model.
         *  
         *   default name = "Constant"
         */
        LE_OneStr* smodel = new LE_OneStr("Open Circuit Voltage Model", &(ocv_input_ptr->OCVModel), 10, 1, 1, "OCVModel");
        smodel->set_default("Constant");
        bOCVoverride->addLineEntry(smodel);

        /* ----------------------------------------------------------------------------------
         *    Replaced Species = string      (required)
         *
         *   Name of the replaced species
         */
        LE_OneStr* rspec = new LE_OneStr("Replaced Species", &(ocv_input_ptr->replacedSpeciesName), 1, 1, 1, "ReplacedSpeciesName");
        rspec->set_default("");
        bOCVoverride->addLineEntry(rspec);

        /* ----------------------------------------------------------------------------------
         *    Identify Reaction for OCV Model = int      (optional) (default = 0)
         *
         *   Integer name of the reaction. Since most models just have one reaction, we will start
         *   with assuming that it's the first reaction in the mechanism 
         */
        LE_OneInt* rid = new LE_OneInt("Identify Reaction for OCV Model", &(ocv_input_ptr->rxnID), 1, "OCVModel_rxnID");
        rid->set_default(0);
        bOCVoverride->addLineEntry(rid);

	/* ----------------------------------------------------------------------------------
         *    Identify Full Cell Reaction for OCV Model = int      (optional) (default = 0)
         *
         *   Integer name of the reaction. Since most models just have one reaction, we will start
         *   with assuming that it's the first reaction in the mechanism. 
	 *   If this card is not specified, there is an error exit. We recently added this
	 *   to specify   
         */
        LE_OneInt* ridFC = new LE_OneInt("Identify Full Cell Reaction for OCV Model", 
					 &(ocv_input_ptr->rxnID_deltaS), 1, "OCVModel_FCrxnID");
        ridFC->set_default(0);
        bOCVoverride->addLineEntry(ridFC);

        /* --------------------------------------------------------------
         *    Temperature Derivative = PickList = [Zero, Species, Model] (required)
         *                                  (default = Zero)
         *
         * Create a PickList Line Element made out of the list of
         * species
         *        zero    = The temperature derivative is set to zero.
         *        species = The temperature derivative is set to the value determined
         *                  by the thermo, even the thermo of the missing species.
         *        model   = The temperature derivative is set by a model, specified
         *                  further in the input deck.
         */
        const char* ccTD[3] = {"zero", "species", "model"};
        LE_PickList* plTD = new LE_PickList("Temperature Derivative", &(ocv_input_ptr->temperatureDerivType),
                                            ccTD, 3, 1, "OCVModel_tempDerivType");
        plTD->set_default(0);
        bOCVoverride->addLineEntry(plTD);

        /* --------------------------------------------------------------
         *    Temperature for OCV =  double [optional]  (default = 298.15 K)
         *
         *    Create an initial temperature for the OCV. If a temperature derivative is 
         *    entered, then this value will be used in a linear interpolation of the OCV.
         */
        LE_OneDbl* rtdval =  new LE_OneDbl("Temperature for OCV", &(ocv_input_ptr->temperatureBase),
                                           0, "OCVModel_tempBase");
        rtdval->set_default(298.15);
        bOCVoverride->addLineEntry(rtdval);

        /* --------------------------------------------------------------
         *      ocv_input.OCVTempDerivModel = string   (required)
         *
         *  Input a string name that will be inputted into a factor to pick the OCV model
         *  from a list of options.
         *  The first token will be interpreted as the string name for the model.
         *  All other tokens will be interpreted as doubles and input as parameters for the model.
         *  
         *   default name = "Constant"
         */
        LE_OneStr* rtdmodel = new LE_OneStr("Open Circuit Voltage Temperature Derivative Model", &(ocv_input_ptr->OCVTempDerivModel),
                                            10, 0, 1, "OCVTempDerivModel");
        rtdmodel->set_default("Constant 0.0");
        bOCVoverride->addLineEntry(rtdmodel);

	/* ----------------------------------------------------------------------------------
         *   Species identified as DoD Surrogate = string      (optional) (default = "")
         *
         *   Name of the  species whose mole fraction can be identifed with the DoD surrogate
         */
        LE_OneStr* dodspec = new LE_OneStr("Species identified as DoD Surrogate", &(ocv_input_ptr->DoDSurrogateSpeciesName), 
                                          1, 1, 0, "MF_DoD_LocalSpeciesID");
        dodspec->set_default("");
        bOCVoverride->addLineEntry(dodspec);

    }
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------

    BlockEntry* rb = new BlockEntry("Reaction Extent Limits", 0);
    cf->addSubBlock(rb);
    LE_OneDbl* iRETop = new LE_OneDbl("Reaction Extent Top Limit", &(RxnExtTopLimit), 1, "RxnExtTopLimit");
    iRETop->set_default(-1.0);
    rb->addLineEntry(iRETop);

    LE_OneDbl* iREBot = new LE_OneDbl("Reaction Extent Bottom Limit", &(RxnExtBotLimit), 1, "RxnExtBotLimit");
    iREBot->set_default(-1.0);
    rb->addLineEntry(iREBot);

    /*
     *   SPECIFY OTHER WAYS TO INPUT THE CONCENTRATION, without
     *   breaking it down via phases.
     */

    /* ------------------------------------------------------------------
     * Block Input For initial number of moles of species
     * -> optional block where not all entries must be specified in that block
     * -> Each entry must match an entry in SpeciesNames of which there are nTotSpecies.
     */
    BlockEntry* BESIM = new BE_StrVecDbl("Species Initial KMoles", &MoleNumber, 0, 0, SpeciesNames, nTotSpecies, 1, "MoleNumber");
    cf->addSubBlock(BESIM);

    /* ------------------------------------------------------------------
     * Specify the initial relative amount of capacity discharged
     */
    LE_OneDbl* iCapDis =  new LE_OneDbl("Capacity Discharged Per Initial Mole", &(RelativeCapacityDischargedPerMole), 0,
                                        "RelCapacityDischarged");
    iCapDis->set_default(-1.0);
    cf->addLineEntry(iCapDis);

    /* -----------------------------------------------------------------------------------------------------------
     *   Do Polarization Analysis = 1
     */
    LE_OneBool* eBool = new LE_OneBool("Do Polarization Analysis", &(doPolarizationAnalysis_), 0, "doPolarizationAnalysis");
    eBool->set_default(false);
    cf->addLineEntry(eBool);

    /* ------------------------------------------------------------------
     * Define a block that may occur multiples times
     *  The number of occurances of the block is .... PO.numExtraGlobalRxns
     *  All of the data defined in the block is
     *  in the pointer to the vector of structures .. PO.m_EGRList
     */
    BlockEntry* sbEGR = 
        new BEInput::BE_MultiBlockVec<Zuzax::EGRInput>("Extra Global Reaction", &(numExtraGlobalRxns), &(m_EGRList), 0);
 
    cf->addSubBlock(sbEGR);

    /*
     * OK now define items in the block
     */
    EGRInput* egr_ptr = m_EGRList[0];

    /* ------------------------------------------------------------------
     *  Define the reacting surface domain index
     *  -> optional, default is 0
     */
    LE_OneSizet* iRSDIndex = new LE_OneSizet("Reacting Surface Domain Index", &(egr_ptr->m_RSD_index), 0, "ReactingSurfaceDomainIndex");
    iRSDIndex->set_default(0);
    sbEGR->addLineEntry(iRSDIndex);

    /* --------------------------------------------------------------
     * EGR -> Special Species
     *
     * Create a PickList Line Element made out of the list of species
     */
    LE_PickList* egrSS = new LE_PickList("Special Species", &(egr_ptr->m_SS_KinSpeciesKindex), (const char**) SpeciesNames,
                                         nTotSpecies, 1, "Special_Species");
    egrSS->set_default(0);
    sbEGR->addLineEntry(egrSS);

    BlockEntry* sbERS = new BEInput::BE_MultiBlockVec<Zuzax::ERSSpec>("Elementary Reaction Specification",
                                                                       &(egr_ptr->m_numElemReactions), &(egr_ptr->m_ERSList), 0);
    sbEGR->addSubBlock(sbERS);
    /*
     * OK now define items in that block
     */
    ERSSpec* ers_ptr = egr_ptr->m_ERSList[0];

    /* ------------------------------------------------------------------
     *  Define a reaction index
     */
    LE_OneInt* iRxnIndex = new LE_OneInt("Reaction Index", &(ers_ptr->m_reactionIndex), 1, "ReactionIndex");
    iRxnIndex->set_default(-1);
    sbERS->addLineEntry(iRxnIndex);

    LE_OneDbl* iRxnMult = new LE_OneDbl("Reaction Multiplier", &(ers_ptr->m_reactionMultiplier), 1, "ReactionMultiplier");
    iRxnMult->set_default(0.0);
    sbERS->addLineEntry(iRxnMult);

    BaseEntry::set_SkipUnknownEntries(3);
}
//======================================================================================================================
/*
 *  printFlag 0 = no output
 *            1 = error message output
 *            2 = output -> processed lines
 *            3 = output -> original line and processed lines
 */
bool process_electrode_input(BlockEntry* cf, std::string fileName, int printFlag, int pass)
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
    const TKInput::TOKEN tok_in;
    TKInput::TOKEN tok_out;
    FILE* ifp = fopen(fileName.c_str(), "r");
    if (!ifp) {
        if (printFlag) {
            cout << "ERROR can't open file " << fileName<< endl;
        }
        return false;
    }
    if (printFlag > 1) {
        printf("==========================================================\n");
        printf(" STARTING PROCESSING COMMAND FILE %s, PASS # %d\n", fileName.c_str(), pass);
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
         * This catches the block input error messages
         */
        cout << bi.errorMessage() << endl;
        return false;
    }
    fclose(ifp);
    if (printFlag > 1) {
        printf("=========================================================\n");
        printf(" FINISHED PROCESSING COMMAND FILE %s, PASS # %d\n", fileName.c_str(), pass);
        printf("=========================================================\n");
    }
    return true;
}
//==================================================================================================================================
/*
 *  Set the bath conditions for the one ThermoPhase
 */
void setElectrodeBathSpeciesConditions(ZZCantera::thermo_t_double& tp, ELECTRODE_KEY_INPUT& EI, ElectrodeBath& BG, size_t iph, int printLvl)
{
    size_t nsp = tp.nSpecies();
    size_t kstart = BG.m_pl->globalSpeciesIndex(iph, 0);
    tp.setState_TPX(EI.Temperature, EI.Pressure, BG.XmolPLSpecVec + kstart);
    tp.setElectricPotential(EI.PotentialPLPhases[iph]);
    /*
     * Print out table summarizing bath gas conditions, if requested
     */
    if (printLvl) {
        double* C = new double [nsp];
        double* act = new double [nsp];
        tp.getConcentrations(C);
        tp.getActivities(act);
        print_char('=', 100);
        std::cout << "\n";
        dnt(0);
        std::cout << " SUMMARY OF SPECIES IN THE MECHANISM WITH A DESCRIPTION OF BATH COMPOSITION:\n\n"; 
        dnt(1);
        std::cout << "Total pressure = " << (EI.Pressure * 760. / OneAtm) << " torr \n";
        dnt(1);
        std::cout << "Temperature (where needed) = " << EI.Temperature << " Kelvin\n";
        dnt(1);
        std::cout << "Voltage (where needed) = " << EI.PotentialPLPhases[iph] << " Volts\n";
        dnt(4);
        std::cout << "Number       Name      Mole_fraction  Concentration (gmol/cm**3)   Activities\n";
        dnt(4);
        std::cout << "-----------------------------------------------------------------------------\n";
        std::string spN;
        for (size_t k = 0; k < nsp; k++) {
            dnt(4);
            pr_if(k+1, 5);
            spN = tp.speciesName(k);
            pr_sf(spN, 16);
            pr_df(BG.XmolPLSpecVec[kstart + k], 16, 4);
            pr_de(C[k] * 1.0E-3, 16, 3);
            pr_de(act[k], 16, 3);
            std::cout << "\n";
        }
        dnt(4);
        std::cout << "------------------------------------------------------------------------------";
        std::cout << std::endl;
        delete [] C;
        delete [] act;
    }
}
//==================================================================================================================================
int ELECTRODE_KEY_INPUT::electrode_input(std::string commandFile, BlockEntry* cf)
{
    int retn = 0;
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
    /*
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
        throw Electrode_Error("ELECTRODE_KEY_INPUT::electrode_input()", "Number of requested files differ");
    }

    /*
     * Read in all of the phase specifications from the cantera
     * input files into the PhaseList structure.
     */
    std::string fn;
    bool surNotFound = true;
    for (int i = 0; i < NumberCanteraFiles; i++) {
        fn = CanteraFileNames[i];
        importAllCTMLIntoPhaseList(m_pl, fn);
        if (surNotFound && (m_pl->nSurPhases() > 0)) {
            surNotFound = false;
        }
    }
    /*
     * Setup internally for next pass through the input file.
     */
    InitForInput(m_pl);

    /*
     * Setup and process the input deck for third time
     * -> Might have to print out the input and quit as well.
     */
    setup_input_pass3(cf);

    /*
     * Possibly change the print level for the last
     */
    Electrode_Types_Enum ieos = string_to_Electrode_Types_Enum(electrodeModelName);
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
//==================================================================================================================================
int ELECTRODE_KEY_INPUT::electrode_input_child(std::string commandFile, BlockEntry* cf)
{
    cf->clear();
    /*
     *  Redo the base input
     */
    electrode_input(commandFile, cf);
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
//==================================================================================================================================
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
    int nt;
    /*
     *  Loop Over all phases in the PhaseList, locating any default specifications
     *  of the mole numbers and mole fractions that exist in the input file.
     */
    for (size_t iph = 0; iph < m_pl->nPhases(); iph++) {
        size_t kstart = m_pl->globalSpeciesIndex(iph, 0);
        ThermoPhase* tphase = &(m_pl->thermo(iph));
        size_t nsp = tphase->nSpecies();
        // Find the name of the input block
        std::string phaseBath = "Bath Specification for Phase ";
        std::string phaseNm = tphase->name();
        phaseBath += phaseNm;
        //double* molF = m_BG->XmolPLPhases[iph];
        double* const molF = m_BG.XmolPLSpecVec + kstart;
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
            double kmol = m_BG.PhaseMoles[iph];

            BEInput::BlockEntry* pbsmf = pblock->searchBlockEntry("Bath Species Mole Fraction");
            if (pbsmf) {
                if (pbsmf->get_NumTimesProcessed() > 0) {
                    molVecSpecified = true;
                }
            }
            if (molVecSpecified) {
                for (size_t k = 0; k < nsp; k++) {
                    MoleFraction[kstart + k] = molF[k];
                }
                tphase->setMoleFractions(molF);
            } else {
                tphase->getMoleFractions(MoleFraction.data() + kstart);
            }
	    /*
	     *  Process the Phase Mass entry
	     */
            if (nt_pmoles == 0) {
                BEInput::LineEntry*  pbpm = pblock->searchLineEntry("Phase Mass");
                if (pbpm) {
                    nt = pbpm->get_NumTimesProcessed();
                    if (nt > 0) {
                        kmol = m_BG.PhaseMass[iph] / tphase->meanMolecularWeight();
                        m_BG.PhaseMoles[iph] = kmol;
                    }
                }
            }

            if (!molVecSpecified) {
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
                    m_ptr->setState_TPM(Temperature, Pressure, m_BG.MolalitiesPLSpecVec + kstart);
                    m_ptr->getMoleFractions(molF);
                    m_ptr->getMoleFractions(MoleFraction.data() + kstart);
                }
            }
            /*
             * From the bath specification, get the total number of moles of the phase
             */
            double totalMoles = m_BG.PhaseMoles[iph];
	    /*
	     *  Set the electric potential of the phase 
	     */
            tphase->setElectricPotential(PotentialPLPhases[iph]);
            /*
             *  Setup the global MoleNumber array in the Electrode object
             */
            for (size_t k = 0; k < nsp; k++) {
                MoleNumber[kstart + k] = totalMoles * molF[k];
            }
        } else {
            /*
             *  We have no information about the specification of the phase.
             *  All we can do is set the mole numbers of all of the species to zero
             *  We also read in the mole fraction vector from the ThermoPhase object
             */
            for (size_t k = 0; k < nsp; k++) {
                MoleNumber[kstart + k] = 0.0;
            }
            tphase->getMoleFractions(MoleFraction.data() + kstart);
        }

    }

    //
    // Post process the OCV override structure
    //
    for (size_t iphS = 0; iphS < m_pl->nSurPhases(); iphS++) {
	std::string phaseBath = "Open Circuit Potential Override for interface ";
        ThermoPhase* tp = &(m_pl->surPhase(iphS));
        std::string phaseNm = tp->name();
        phaseBath += phaseNm;

	BEInput::BlockEntry* pblock = cf->searchBlockEntry(phaseBath.c_str());
        int numTimes = 0;
        if (pblock) {
            numTimes = pblock->get_NumTimesProcessed();
        }
	if (numTimes > 1) {
	    throw Electrode_Error("Electrode_input::post_input_pass3", "block processed more than once");
	}
	if (numTimes == 1) {
	    OCV_Override_input* ocv_input_ptr = OCVoverride_ptrList[iphS];
            ocv_input_ptr->numTimes = 1;
            //
            // Discover the replacedSpeciesID
            //
            int kg = ocv_input_ptr->replacedGlobalSpeciesID = m_pl->globalSpeciesIndex(ocv_input_ptr->replacedSpeciesName);
            if (kg < 0) {
	        throw Electrode_Error("Electrode_input::post_input_pass3", "Species not found in phaselist : " 
	                            + ocv_input_ptr->replacedSpeciesName);
            }
            ocv_input_ptr->replacedGlobalSpeciesID = kg;
            
            size_t phaseID; 
            size_t localSpeciesIndex; 
            size_t kgs = kg;
            m_pl->getLocalIndecisesFromGlobalSpeciesIndex(kgs, phaseID, localSpeciesIndex);
            //
            // Store the phase index and local species index of the replaced species
            // 
            ocv_input_ptr->replacedSpeciesPhaseID = phaseID;
            ocv_input_ptr->replacedLocalSpeciesID = localSpeciesIndex;
            //
            // We can't check the validity of the kinetics reaction because we haven't yet set up the ReactingSurDomain,
            // And, we haven't yet read in the kinetics mechanism. We do this in the ReactingSurDomain setup.
            // 

	    double* CapZeroDoDCoeff = m_BG.CapZeroDoDCoeffPhases[phaseID];

	    double* CapLeftCoef = m_BG.CapLeftCoeffPhases[phaseID];

	    size_t nsp = m_pl->thermo(phaseID).nSpecies();
	    //
	    //  Find a species in the replacing phase whose mole fraction can be used as a simple DoD indicator.
	    //    (Confirmed to work for at least one case (MCMB -as an anode).
	    //
	    size_t ik = npos;
            bool foundNonzeroCapacity = false;
	    for (size_t k = 0; k <  nsp; k++) {
		if (CapZeroDoDCoeff[k] == 1.0) {
                    foundNonzeroCapacity = true;
		    if (CapLeftCoef[k] == 0.0) {
			if (ik != npos) {
			    // We don't have a simple model to employ be prepared to throw an error
			    ik = npos;
			    break;
			}
			ik = k;
		    }
		}
	    }
            if (!foundNonzeroCapacity) {
               throw Electrode_Error("post_input_pass3()",
                                     "CapZeroDodCoeff block is zero. Need to add the block to the input file");
            }
	    if (ocv_input_ptr->DoDSurrogateSpeciesName != "") {
		size_t k =  m_pl->thermo(phaseID).speciesIndex(ocv_input_ptr->DoDSurrogateSpeciesName);
		if (k != npos) {
		    ocv_input_ptr->MF_DoD_LocalSpeciesID = k;
		}
	    }
	    if (ocv_input_ptr->MF_DoD_LocalSpeciesID == npos) {
		if (ik != npos) {
		    ocv_input_ptr->MF_DoD_LocalSpeciesID = ik;
		} else {
		    //  Probably an error: but leaving this unused for the moment
		    throw Electrode_Error("post_input_pass3()",
		    			  "Couldn't determine easily the DoD variable. Specify it explicitly");
		}
	    }

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
 *  Virtual function that may get added into
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
} // end of namespace
//-----------------------------------------------------------------------------------------------------------------------------------
