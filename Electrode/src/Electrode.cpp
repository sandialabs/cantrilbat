/*
 * $Id: Electrode.cpp 593 2013-05-13 21:25:47Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"

#include "cantera/equilibrium.h"

#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/thermo/FixedChemPotSSTP.h"

#include "cantera/numerics/solveProb.h"
#include "cantera/numerics/BEulerInt.h"

#include "cantera/solvers.h"

#include "PhaseList.h"
#include "BlockEntry.h"
#include "LE_PickList.h"
#include "BE_MoleComp.h"
#include "BE_UnitConversionPressure.h"
#include "BE_MultiBlock.h"
#include "LE_OneDblUnits.h"
#include "LE_OneStr.h"
#include "LE_OneInt.h"
#include "LE_OneDbl.h"
#include "LE_OneBool.h"
#include "LE_MultiCStr.h"
#include "BE_MolalityComp.h"

#include "Electrode.h"
#include "ReactingSurDomain.h"
#include "importAllCTML.h"
#include "RxnMolChange.h"
#include "ExtraGlobalRxn.h"

#include "Electrode_input.h"
#include "Electrode_FuncCurrent.h"
#include "Electrode_Factory.h"

#include "ApplBase_print.h"

using namespace Cantera;
using namespace std;
using namespace BEInput;
using namespace TKInput;

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif
#ifndef SAFE_DELETE
#define SAFE_DELETE(x)  if ((x)) { delete (x) ; x = 0 ; }
#endif

// HKM debug code to turn it back to the bad way of doing things to see what the damage was
//#define OLD_FOLLOW 1

namespace Cantera {
//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode::Electrode() :
                PhaseList(),
                electrodeCapacityType_(CAPACITY_ANODE_ECT),
                pendingIntegratedStep_(0),

                prob_type(TP),
                numSurfaces_(0),
                followElectrolyteMoles_(0),
                electrolytePseudoMoles_(1.0),
                externalFieldInterpolationType_(T_FINAL_CONST_FIS),
                t_init_init_(0.0),
                t_final_final_(0.0),
                tinit_(0.0),
                tfinal_(0.0),
                deltaTsubcycleMax_(1.0E300),
                deltaTsubcycle_init_init_(1.0E300),
                deltaTsubcycleNext_(1.0E300),
                deltaTsubcycle_init_next_(1.0E300),
                choiceDeltaTsubcycle_init_(0),
                numIntegrationSubCycles_final_final_(1),
                temperature_(298.15),
                pressure_(1.0E5),
                ElectrodeSolidVolume_(0.0),
                phaseMolarVolumes_(0),
                sphaseMolarAreas_(0),
                VolPM_(0),
                spMoles_final_(0),
                spMoles_final_final_(0),
                spMoles_init_(0),
                spMoles_init_init_(0),
                spMoles_dot_(0),
                spMoles_predict_(0),
                spMf_init_(0),
                spMf_init_init_(0),
                spMf_final_(0),
                spElectroChemPot_(0),
                phaseVoltages_(0),
                RSD_List_(0),
                ActiveKineticsSurf_(0),
                phaseMoles_init_(0),
                phaseMoles_init_init_(0),
                phaseMoles_final_(0),
                phaseMoles_final_final_(0),
                phaseMoles_dot_(0),
                numExternalInterfacialSurfaces_(0),
                isExternalSurface_(0),
                surfaceAreaRS_init_(0),
                surfaceAreaRS_final_(0),
                surfaceAreaRS_init_init_(0),
                surfaceAreaRS_final_final_(0),
                spNetProdPerArea_List_(0, 0),
                spMoleIntegratedSourceTerm_(0),
                spMoleIntegratedSourceTermLast_(0),
                electrodeName_("Electrode"),
                numExtraGlobalRxns(0),
                m_EGRList(0),
                m_egr(0),
                m_rmcEGR(0),
                metalPhase_(-1),
                solnPhase_(-1),
                kElectron_(-1),
                kKinSpecElectron_sph_(0),
                deltaVoltage_(0.0),
                electronKmolDischargedToDate_(0.0),
                capacityLeftSpeciesCoeff_(0),
                capacityZeroDoDSpeciesCoeff_(0),
                Icurrent_(0.0),
                deltaG_(0),
                inputParticleDiameter_(1.0E-6),
                particleNumberToFollow_(-1.0),
                Radius_exterior_init_init_(5.0E-7),
                Radius_exterior_init_(5.0E-7),
                Radius_exterior_final_(5.0E-7),
                Radius_exterior_final_final_(5.0E-7),
                porosity_(0.0),
                molarAtol_(1.0E-16),
                xmlTimeIncrementData_(0),
                xmlTimeIncrementIntermediateData_(0),
                xmlExternalData_init_init_(0),
                xmlExternalData_init_(0),
                xmlExternalData_final_(0),
                xmlExternalData_final_final_(0),
                xmlStateData_init_init_(0),
                xmlStateData_init_(0),
                xmlStateData_final_(0),
                xmlStateData_final_final_(0),
                eState_final_(0),
                baseNameSoln_("soln"),
                electrodeModelType_(0),
                electrodeDomainNumber_(0),
                electrodeCellNumber_(0),
                counterNumberIntegrations_(0),
                counterNumberSubIntegrations_(0),
                printLvl_(4),
                printXMLLvl_(0),
                printCSVLvl_(0),
                detailedResidPrintFlag_(0),
                enableExtraPrinting_(false)
{
    //m_pl = new PhaseList();
    //m_rSurDomain = new ReactingSurDomain();
}

//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode::Electrode(const Electrode& right) :
    PhaseList(),
    electrodeCapacityType_(CAPACITY_ANODE_ECT),
    pendingIntegratedStep_(0),
                prob_type(TP),
                numSurfaces_(0),
                followElectrolyteMoles_(0),
                electrolytePseudoMoles_(1.0),
                externalFieldInterpolationType_(T_FINAL_CONST_FIS),
                t_init_init_(0.0),
                t_final_final_(0.0),
                tinit_(0.0),
                tfinal_(0.0),
                deltaTsubcycleMax_(1.0E300),
                deltaTsubcycle_init_init_(1.0E300),
                deltaTsubcycleNext_(1.0E300),
                deltaTsubcycle_init_next_(1.0E300),
                choiceDeltaTsubcycle_init_(0),
                numIntegrationSubCycles_final_final_(1),
                temperature_(298.15),
                pressure_(1.0E5),
                ElectrodeSolidVolume_(0.0),
                phaseMolarVolumes_(0),
                sphaseMolarAreas_(0),
                VolPM_(0),
                spMoles_final_(0),
                spMoles_final_final_(0),
                spMoles_init_(0),
                spMoles_init_init_(0),
                spMoles_dot_(0),
                spMoles_predict_(0),
                spMf_init_(0),
                spMf_init_init_(0),
                spMf_final_(0),
                spElectroChemPot_(0),
                phaseVoltages_(0),
                RSD_List_(0),
                ActiveKineticsSurf_(0),
                phaseMoles_init_(0),
                phaseMoles_init_init_(0),
                phaseMoles_final_(0),
                phaseMoles_final_final_(0),
                phaseMoles_dot_(0),
                numExternalInterfacialSurfaces_(0),
                isExternalSurface_(0),
                surfaceAreaRS_init_(0),
                surfaceAreaRS_final_(0),
                surfaceAreaRS_init_init_(0),
                surfaceAreaRS_final_final_(0),
                spNetProdPerArea_List_(0, 0),
                spMoleIntegratedSourceTerm_(0),
                spMoleIntegratedSourceTermLast_(0),
                electrodeName_("Electrode"),
                numExtraGlobalRxns(0),
                m_EGRList(0),
                m_egr(0),
                m_rmcEGR(0),
                metalPhase_(-1),
                solnPhase_(-1),
                kElectron_(-1),
                kKinSpecElectron_sph_(0),
                deltaVoltage_(0.0),
                electronKmolDischargedToDate_(0.0),
                capacityLeftSpeciesCoeff_(0),
                capacityZeroDoDSpeciesCoeff_(0),
                Icurrent_(0.0),
                deltaG_(0),
                inputParticleDiameter_(1.0E-6),
                particleNumberToFollow_(-1.0),
                Radius_exterior_init_init_(5.0E-7),
                Radius_exterior_init_(5.0E-7),
                Radius_exterior_final_(5.0E-7),
                Radius_exterior_final_final_(5.0E-7),
                porosity_(0.0),
                molarAtol_(1.0E-16),
                xmlTimeIncrementData_(0),
                xmlTimeIncrementIntermediateData_(0),
                xmlExternalData_init_init_(0),
                xmlExternalData_init_(0),
                xmlExternalData_final_(0),
                xmlExternalData_final_final_(0),
                xmlStateData_init_init_(0),
                xmlStateData_init_(0),
                xmlStateData_final_(0),
                xmlStateData_final_final_(0),
                eState_final_(0),
                baseNameSoln_("soln"),
                electrodeModelType_(0),
                electrodeDomainNumber_(0),
                electrodeCellNumber_(0),
                counterNumberIntegrations_(0),
                counterNumberSubIntegrations_(0),
                printLvl_(4),
                printXMLLvl_(0),
                printCSVLvl_(0),
                detailedResidPrintFlag_(0),
                enableExtraPrinting_(false)
{
    /*
     * Call the assignment operator.
     */
    operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode& Electrode::operator=(const Electrode& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    PhaseList::operator=(right);
    /*
     * We do a straight assignment operator on all of the
     * data. The vectors are copied.
     */

    /*
     * Copy over the ReactingSurDomain list
     */
    for (int i = 0; i < numSurfaces_; i++) {
        bool idHit = false;
        //if (right.m_rSurDomain == right.RSD_List_[i]) {
        //	idHit = true;
        //}
        if (RSD_List_[i]) {
            delete RSD_List_[i];
            RSD_List_[i] = 0;
        }
        if (right.RSD_List_[i]) {
            RSD_List_[i] = new ReactingSurDomain(*(right.RSD_List_[i]));
            if (idHit) {
                // m_rSurDomain = RSD_List_[i];
            }
        }
    }
    electrodeCapacityType_ = right.electrodeCapacityType_;
    pendingIntegratedStep_ = right.pendingIntegratedStep_;
    prob_type = right.prob_type;
    numSurfaces_ = right.numSurfaces_;
    followElectrolyteMoles_ = right.followElectrolyteMoles_;
    electrolytePseudoMoles_ = right.electrolytePseudoMoles_;
    externalFieldInterpolationType_ = right.externalFieldInterpolationType_;
    t_init_init_ = right.t_init_init_;
    t_final_final_ = right.t_final_final_;
    tinit_ = right.tinit_;
    tfinal_ = right.tfinal_;
    deltaTsubcycleMax_ = right.deltaTsubcycleMax_;
    deltaTsubcycle_init_init_ = right.deltaTsubcycle_init_init_;
    deltaTsubcycleNext_ = right.deltaTsubcycleNext_;
    deltaTsubcycle_init_next_ = right.deltaTsubcycle_init_next_;
    choiceDeltaTsubcycle_init_ = right.choiceDeltaTsubcycle_init_;
    numIntegrationSubCycles_final_final_ = right.numIntegrationSubCycles_final_final_;
    temperature_ = right.temperature_;
    pressure_ = right.pressure_;
    ElectrodeSolidVolume_ = right.ElectrodeSolidVolume_;
    phaseMolarVolumes_ = right.phaseMolarVolumes_;
    sphaseMolarAreas_ = right.sphaseMolarAreas_;
    VolPM_ = right.VolPM_;
    spMoles_final_ = right.spMoles_final_;
    spMoles_final_final_ = right.spMoles_final_final_;
    spMoles_init_ = right.spMoles_init_;
    spMoles_init_init_ = right.spMoles_init_init_;
    spMoles_dot_ = right.spMoles_dot_;
    spMoles_predict_ = right.spMoles_predict_;
    spMf_init_ = right.spMf_init_;
    spMf_init_init_ = right.spMf_init_init_;
    spMf_final_ = right.spMf_final_;
    spElectroChemPot_ = right.spElectroChemPot_;
    phaseVoltages_ = right.phaseVoltages_;
    ActiveKineticsSurf_ = right.ActiveKineticsSurf_;
    phaseMoles_init_ = right.phaseMoles_init_;
    phaseMoles_init_init_ = right.phaseMoles_init_init_;
    phaseMoles_final_ = right.phaseMoles_final_;
    phaseMoles_final_final_ = right.phaseMoles_final_final_;
    phaseMoles_dot_ = right.phaseMoles_dot_;
    justBornPhase_ = right.justBornPhase_;
    justDiedPhase_ = right.justDiedPhase_;
    numExternalInterfacialSurfaces_ = right.numExternalInterfacialSurfaces_;
    isExternalSurface_ = right.isExternalSurface_;
    surfaceAreaRS_init_ = right.surfaceAreaRS_init_;
    surfaceAreaRS_final_ = right.surfaceAreaRS_final_;
    surfaceAreaRS_init_init_ = right.surfaceAreaRS_init_init_;
    surfaceAreaRS_final_final_ = right.surfaceAreaRS_final_final_;
    spNetProdPerArea_List_ = right.spNetProdPerArea_List_;
    spMoleIntegratedSourceTerm_ = right.spMoleIntegratedSourceTerm_;
    spMoleIntegratedSourceTermLast_ = right.spMoleIntegratedSourceTermLast_;
    electrodeName_ = right.electrodeName_;
    numExtraGlobalRxns = right.numExtraGlobalRxns;
    Cantera::deepStdVectorPointerCopy<EGRInput>(right.m_EGRList, m_EGRList);
    Cantera::deepStdVectorPointerCopy<ExtraGlobalRxn>(right.m_egr, m_egr);
    Cantera::deepStdVectorPointerCopy<RxnMolChange>(right.m_rmcEGR, m_rmcEGR);
    metalPhase_ = right.metalPhase_;
    solnPhase_ = right.solnPhase_;
    kElectron_ = right.kElectron_;
    kKinSpecElectron_sph_ = right.kKinSpecElectron_sph_;
    deltaVoltage_ = right.deltaVoltage_;
    electronKmolDischargedToDate_ = right.electronKmolDischargedToDate_;
    capacityLeftSpeciesCoeff_ = right.capacityLeftSpeciesCoeff_;
    capacityZeroDoDSpeciesCoeff_ = right.capacityZeroDoDSpeciesCoeff_;
    capacityInitialZeroDod_ = right.capacityInitialZeroDod_;
    Icurrent_ = right.Icurrent_;
    deltaG_ = right.deltaG_;
    inputParticleDiameter_ = right.inputParticleDiameter_;
    particleNumberToFollow_ = right.particleNumberToFollow_;
    Radius_exterior_init_init_ = right.Radius_exterior_init_init_;
    Radius_exterior_init_ = right.Radius_exterior_init_;
    Radius_exterior_final_ = right.Radius_exterior_final_;
    Radius_exterior_final_final_ = right.Radius_exterior_final_final_;
    porosity_ = right.porosity_;
    molarAtol_ = right.molarAtol_;

    SAFE_DELETE(xmlTimeIncrementData_);
    xmlTimeIncrementData_ = new XML_Node(*right.xmlTimeIncrementData_);

    SAFE_DELETE(xmlTimeIncrementIntermediateData_);
    xmlTimeIncrementIntermediateData_ = new XML_Node(*right.xmlTimeIncrementIntermediateData_);

    SAFE_DELETE(xmlExternalData_init_init_);
    xmlExternalData_init_init_ = new XML_Node(*right.xmlExternalData_init_init_);

    SAFE_DELETE(xmlExternalData_init_);
    xmlExternalData_init_ = new XML_Node(*right.xmlExternalData_init_);

    SAFE_DELETE(xmlExternalData_final_);
    xmlExternalData_final_ = new XML_Node(*right.xmlExternalData_final_);

    SAFE_DELETE(xmlExternalData_final_final_);
    xmlExternalData_final_final_ = new XML_Node(*right.xmlExternalData_final_final_);

    SAFE_DELETE(xmlStateData_init_init_);
    xmlStateData_init_init_ = new XML_Node(*right.xmlStateData_init_init_);

    SAFE_DELETE(xmlStateData_init_);
    xmlStateData_init_ = new XML_Node(*right.xmlStateData_init_);

    SAFE_DELETE(xmlStateData_final_);
    xmlStateData_final_ = new XML_Node(*right.xmlStateData_final_);

    SAFE_DELETE(xmlStateData_final_final_);
    xmlStateData_final_final_ = new XML_Node(*right.xmlStateData_final_final_);

    SAFE_DELETE(eState_final_);
    eState_final_ = new EState();

    baseNameSoln_ = right.baseNameSoln_;

    electrodeModelType_ = right.electrodeModelType_;
    electrodeDomainNumber_ = right.electrodeDomainNumber_;
    electrodeCellNumber_ = right.electrodeCellNumber_;
    printLvl_ = right.printLvl_;
    printXMLLvl_ = right.printXMLLvl_;
    printCSVLvl_ = right.printCSVLvl_;
    detailedResidPrintFlag_ = right.detailedResidPrintFlag_;
    enableExtraPrinting_ = right.enableExtraPrinting_;
    counterNumberIntegrations_ = 0;
    counterNumberSubIntegrations_ = 0;

    eState_final_->initialize(this);

    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
/*
 *
 *  ELECTRODE_INPUT:destructor
 *
 * We need to manually free all of the arrays.
 */
Electrode::~Electrode()
{
    int is = m_EGRList.size();
    for (int i = 0; i < is; i++) {
        delete m_EGRList[i];
        m_EGRList[i] = 0;
    }

    for (int i = 0; i < numSurfaces_; i++) {
        if (RSD_List_[i]) {
            delete RSD_List_[i];
        }
    }
    //m_rSurDomain = 0;

    is = m_egr.size();
    for (int i = 0; i < is; i++) {
        if (m_egr[i]) {
            delete m_egr[i];
            m_egr[i] = 0;
        }
    }
    is = m_rmcEGR.size();
    for (int i = 0; i < is; i++) {
        if (m_rmcEGR[i]) {
            delete m_rmcEGR[i];
            m_rmcEGR[i] = 0;
        }
    }
    SAFE_DELETE(xmlTimeIncrementData_);
    SAFE_DELETE(xmlTimeIncrementIntermediateData_);
    SAFE_DELETE(xmlExternalData_init_init_);
    SAFE_DELETE(xmlExternalData_init_);
    SAFE_DELETE(xmlExternalData_final_);
    SAFE_DELETE(xmlExternalData_final_final_);

    SAFE_DELETE(xmlStateData_init_init_);
    SAFE_DELETE(xmlStateData_init_);
    SAFE_DELETE(xmlStateData_final_);
    SAFE_DELETE(xmlStateData_final_final_);

    SAFE_DELETE(eState_final_);

}
//======================================================================================================================
// Duplicator function
/*
 *  Duplicate the current electrode object, returning a base electrode pointer
 *
 *  (Virtual function from Electrode.h)
 */
Electrode* Electrode::duplMyselfAsElectrode() const
{
    return 0;
}
//======================================================================================================================
// Set the electrode ID information
void Electrode::setID(int domainNum, int cellNum)
{
    electrodeDomainNumber_ = domainNum;
    electrodeCellNumber_ = cellNum;
}
//======================================================================================================================
Electrode_Types_Enum Electrode::electrodeType() const
{
    return BASE_TYPE_ET;
}
//======================================================================================================================
// Add additional Keylines for child electrode objects, and then read them in
/*
 *   (virtual function from Electrode)
 *   (overload virtual function - replaces parent members)
 *
 *   This function will replace the ELECTRODE_KEY_INPUT structure with an expanded
 *   child member structure containing the extra information.
 *
 *   If the command file has been read before, it will then reparse the command file
 *   storring the new information in the ELECTRODE_KEY_INPUT structure.
 *
 *    @param ei  Handle to the ELECTRODE_KEY_INPUT base pointer. This handle may change
 *               as the child class of  ELECTRODE_KEY_INPUT gets malloced.
 *
 *    @return  0 successful but no change in ei
 *             1 Successful and ei has changed
 *            -1 unsuccessful fatal error of some kind.
 */
int Electrode::electrode_input_child(ELECTRODE_KEY_INPUT** ei)
{
    return 0;
}
//======================================================================================================================
void ErrorModelType(int pos, std::string actual, std::string expected)
{
    throw Cantera::CanteraError("Electrode::electrode_model_create() model id",
                                "At pos " + int2str(pos) + ", expected phase " + expected + " but got phase " + actual);
}
void pmatch(std::vector<std::string>& pn, int pos, std::string expected)
{
    if (pos <= ((int) pn.size() - 1)) {
        if (pn[pos] != expected) {
            ErrorModelType(pos, pn[pos], expected);
        }
    }
}
//======================================================================================================================
//  Setup the electrode
/*
 * @param ei    ELECTRODE_KEY_INPUT pointer object
 */
int Electrode::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{

    int i, iph;

    Electrode_Types_Enum ieos = electrodeType();
    Electrode_Types_Enum seos = string_to_Electrode_Types_Enum(ei->electrodeModelName);

    if (ieos != seos) {
        throw CanteraError(
                " Electrode::electrode_model_create()",
                "Electrode Object Type, " + Electrode_Types_Enum_to_string(ieos) + ", is different than requested type, "
                        + ei->electrodeModelName);
    }

    // use the assignment operator to transfer for now
    // May get more sophisticated
    PhaseList::operator=(*(ei->m_pl));

    // resize m_NumTotPhases = total number of phase vectors
    phaseMoles_init_.resize(m_NumTotPhases, 0.0);
    phaseMoles_init_init_.resize(m_NumTotPhases, 0.0);
    phaseMoles_final_.resize(m_NumTotPhases, 0.0);
    phaseMoles_final_final_.resize(m_NumTotPhases, 0.0);
    phaseMoles_dot_.resize(m_NumTotPhases, 0.0);
    phaseVoltages_.resize(m_NumTotPhases, 0.0);
    phaseMolarVolumes_.resize(m_NumTotPhases, 0.0);
    justDiedPhase_.resize(m_NumTotPhases, 0);

    // resize volume phase vectors

    // Identify the number of surfaces with the number of surface phases for now
    numSurfaces_ = m_NumSurPhases;
    // resize surface phase vectors
    isExternalSurface_.resize(numSurfaces_, true);
    surfaceAreaRS_init_.resize(numSurfaces_, 0.0);
    surfaceAreaRS_final_.resize(numSurfaces_, 0.0);
    surfaceAreaRS_init_init_.resize(numSurfaces_, 0.0);
    surfaceAreaRS_final_final_.resize(numSurfaces_, 0.0);
    RSD_List_.resize(numSurfaces_, 0);
    numRxns_.resize(numSurfaces_, 0);
    ActiveKineticsSurf_.resize(numSurfaces_, 0);
    sphaseMolarAreas_.resize(numSurfaces_, 0.0);

    // Resize species vectors
    spMoles_init_.resize(m_NumTotSpecies, 0.0);
    spMoles_final_.resize(m_NumTotSpecies, 0.0);
    spMoles_final_final_.resize(m_NumTotSpecies, 0.0);
    spMoles_init_init_.resize(m_NumTotSpecies, 0.0);
    spMoles_dot_.resize(m_NumTotSpecies, 0.0);
    spMoles_predict_.resize(m_NumTotSpecies, 0.0);
    spMf_final_.resize(m_NumTotSpecies, 0.0);
    spMf_init_.resize(m_NumTotSpecies, 0.0);
    spMf_init_init_.resize(m_NumTotSpecies, 0.0);
    spMf_final_final_.resize(m_NumTotSpecies, 0.0);
    VolPM_.resize(m_NumTotSpecies, 0.0);
    spElectroChemPot_.resize(m_NumTotSpecies, 0.0);

    capacityLeftSpeciesCoeff_.resize(m_NumTotSpecies, 0.0);
    capacityZeroDoDSpeciesCoeff_.resize(m_NumTotSpecies, 0.0);

    capacityLeftSpeciesCoeffPlat_.resize(5, m_NumTotSpecies);
    capacityZeroDoDSpeciesCoeffPlat_.resize(5, m_NumTotSpecies);

    /*
     * OK, Find the first kinetics object
     */
    int isphfound = -1;
    for (int isph = 0; isph < m_NumSurPhases; isph++) {
        if (SurPhaseHasKinetics[isph]) {
            isphfound = isph;
            break;
        }
    }
    if (isphfound == -1) {
        for (i = 0; i < NumVolPhases_; i++) {
            if (VolPhaseHasKinetics[i]) {
                break;
            }
        }
    }

    /*
     * Assign a reacting surface to each surface in the PhaseList object
     */
    for (int isurf = 0; isurf < numSurfaces_; isurf++) {
        // Right now there is a correspondence between isph and isurf for the initial read in.
        // However this may break in the future. Therefore, I've put in this logic for future expansion
        int isph = isurf;
        if (SurPhaseHasKinetics[isph]) {
            ReactingSurDomain* rsd = new ReactingSurDomain();
            int ok = rsd->importFromPL(this, -1, isph);
            if (!ok) {
                throw CanteraError("cttables main:", "rSurDomain returned an error");
            }
            // We carry a list of pointers to ReactingSurDomain
            RSD_List_[isurf] = rsd;
            numRxns_[isurf] = rsd->nReactions();
            ActiveKineticsSurf_[isurf] = 1;
        }
    }

    /*
     * Calculate the number of external interfacial surfaces
     */
    numExternalInterfacialSurfaces_ = 1;

    /*
     * Resize the species production vector for the electrode
     *  -> per area and total net
     */
    spMoleIntegratedSourceTerm_.resize(m_NumTotSpecies, 0.0);
    spMoleIntegratedSourceTermLast_.resize(m_NumTotSpecies, 0.0);
    spNetProdPerArea_List_.resize(m_NumTotSpecies, numSurfaces_, 0.0);

    /*
     * Load up the temperature and pressure
     */
    temperature_ = ei->Temperature;
    pressure_ = ei->Pressure;

    /*
     *  Loop Over all phases in the PhaseList, adding these
     *  formally to the Electrode object.
     */
    int nspecies = 0;
    for (iph = 0; iph < m_NumTotPhases; iph++) {
	// Get the pointer to the ThermoPhase object
        ThermoPhase* tphase = &(thermo(iph));
        int nSpecies = tphase->nSpecies();

        // Find the name of the input block
        string phaseBath = "Bath Specification for Phase ";
        string phaseNm = tphase->name();
        phaseBath += phaseNm;

        /*
         * Search the ReactingSurDomain to see if the current phase is in
         * the object
         */
        string pname = tphase->id();

        int kstart = nspecies;
        nspecies += nSpecies;

        /*
         *  We copy the mole numbers and mole fraction information from the bath gas section of the input file
         *  here. Note, the total quantitites of electrode material and electrolyte material may be overwritten
         *  later.
	 *  We only really need to populate the _final_ state. We will set all states equal to each other at the end 
	 *  of this process.
         */
        double sum = 0.0;
        for (int k = 0; k < nSpecies; k++) {
            spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] = ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
            spMoles_init_[m_PhaseSpeciesStartIndex[iph] + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k];
            sum += ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
        }

        if (sum > 1.0E-300) {
            for (int k = 0; k < nSpecies; k++) {
                spMf_final_[kstart + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] / sum;
                spMf_init_[kstart + k] = spMf_final_[kstart + k];
                spMf_init_init_[kstart + k] = spMf_final_[kstart + k];
                spMf_final_final_[kstart + k] = spMf_final_[kstart + k];
            }
	    /*
	     *  If the phase formally exists, we set the mole fraction vector from the usual calculation
	     *  and post it down to the ThermoPhase object
	     */
            tphase->setMoleFractions(&(spMf_final_[kstart]));
        } else {
	    /*
	     *  If the phase doesn't formally exist, we get the mole fraction vector from the
	     *  ThermoPhase object itself. This is important as the ThermoPhase object has an
	     *  idea of what a reasonable initial state would be, while we do not have any idea.
	     */
            tphase->getMoleFractions(&(spMf_final_[kstart]));
            tphase->getMoleFractions(&(spMf_init_[kstart]));
        }
	/*
	 *  set the voltage from the thermophase object 
	 */
        phaseVoltages_[iph] = tphase->electricPotential();
    }

    /* set the absolute tolerance to 1e-5 the total number of moles */
    // HEWSON -- 12/1/10 -- this particular tolerance value
    // and all magic tolerance values throughout are subject to revision
    double tMoles = 0.0;
    for (int k = 0; k < m_NumTotSpecies; k++) {
        tMoles += spMoles_final_[k];
    }
    molarAtol_ = tMoles * 1.0E-5;

    /*
     * Now that spMoles_final_ is sized appropriately, we can call updatePhaseNumbers()
     */
    for (iph = 0; iph < m_NumTotPhases; iph++) {
        updatePhaseNumbers(iph);
    }

    electrodeName_ = ei->electrodeName;

    if (ei->electrodeCapacityType == 0) {
        electrodeCapacityType_ = CAPACITY_ANODE_ECT;
    } else {
        electrodeCapacityType_ = CAPACITY_CATHODE_ECT;
    }

    numExtraGlobalRxns = ei->numExtraGlobalRxns;
    if (numExtraGlobalRxns > 0) {
        m_EGRList.resize(numExtraGlobalRxns, 0);

        for (int i = 0; i < numExtraGlobalRxns; i++) {
            AssertTrace(!m_EGRList[i]);
            m_EGRList[i] = new EGRInput(*(ei->m_EGRList[i]));
        }
    }

    for (iph = 0; iph < NumVolPhases_; iph++) {
        phaseMoles_init_[iph] = phaseMoles_final_[iph];
        phaseMoles_init_init_[iph] = phaseMoles_final_[iph];
    }

    /*
     * Identify the phases which will ID the voltage drops
     *  metalPhase_ = Phase where the electron exists.
     */
    metalPhase_ = -1;
    solnPhase_ = -1;
    kElectron_ = -1;
    //   InterfaceKinetics *iK = m_rSurDomain;
    ReactingSurDomain* rsd = RSD_List_[0];
    std::vector<RxnMolChange*>& rmcV = rsd->rmcVector;
    int nr = rmcV.size();
    for (iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tp = &(thermo(iph));
        int nSpecies = tp->nSpecies();
        int nElements = tp->nElements();
        int eElectron = tp->elementIndex("E");
        if (eElectron >= 0) {
            for (int k = 0; k < nSpecies; k++) {
                if (tp->nAtoms(k, eElectron) == 1) {
                    int ifound = 1;
                    for (int e = 0; e < nElements; e++) {
                        if (tp->nAtoms(k, e) != 0.0) {
                            if (e != eElectron) {
                                ifound = 0;
                            }
                        }
                    }
                    if (ifound == 1) {
                        metalPhase_ = iph;
                        kElectron_ = m_PhaseSpeciesStartIndex[iph] + k;
                    }
                }
            }
        }
        if (iph != metalPhase_) {
            int jph = rsd->PLtoKinPhaseIndex_[iph];
            if (jph >= 0) {
                for (int i = 0; i < nr; i++) {
                    RxnMolChange* rmc = rmcV[i];
                    if (rmc->m_phaseChargeChange[jph] != 0) {
                        if (rmc->m_phaseDims[jph] == 3) {
                            solnPhase_ = iph;
                            break;
                        }
                    }
                }
            }
        }
    }
    if (solnPhase_ == -1) {
        throw CanteraError("Electrode::electrode_model_create()", "Couldn't find soln phase");
    }
    if (metalPhase_ == -1) {
        throw CanteraError("Electrode::electrode_model_create()", "Couldn't find metal phase");
    }

    kKinSpecElectron_sph_.resize(numSurfaces_, -1);
    for (int isurf = 0; isurf < numSurfaces_; isurf++) {
        kKinSpecElectron_sph_[isurf] = -1;
        ReactingSurDomain* rsd = RSD_List_[isurf];
        if (rsd) {
            std::vector<RxnMolChange*>& rmcV = rsd->rmcVector;
            if (rsd) {
                InterfaceKinetics* iK = rsd;
                for (iph = 0; iph < m_NumTotPhases; iph++) {
                    ThermoPhase* tp = &(thermo(iph));
                    int nSpecies = tp->nSpecies();
                    int nElements = tp->nElements();
                    int eElectron = tp->elementIndex("E");
                    if (eElectron >= 0) {
                        for (int k = 0; k < nSpecies; k++) {
                            if (tp->nAtoms(k, eElectron) == 1) {
                                int ifound = 1;
                                for (int e = 0; e < nElements; e++) {
                                    if (tp->nAtoms(k, e) != 0.0) {
                                        if (e != eElectron) {
                                            ifound = 0;
                                        }
                                    }
                                }
                                if (ifound == 1) {
                                    metalPhase_ = iph;
                                    kElectron_ = m_PhaseSpeciesStartIndex[iph] + k;
                                    kKinSpecElectron_sph_[isurf] = iK->kineticsSpeciesIndex(k, iph);
                                }
                            }
                        }
                    }
                    if (iph != metalPhase_) {
                        int jph = rsd->PLtoKinPhaseIndex_[iph];
                        if (jph >= 0) {
                            for (int i = 0; i < nr; i++) {
                                RxnMolChange* rmc = rmcV[i];
                                if (rmc->m_phaseChargeChange[jph] != 0) {
                                    if (rmc->m_phaseDims[jph] == 3) {
                                        solnPhase_ = iph;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];

    /*
     *  Find the maximum number of reactions in any reacting domain
     */
    int mR = spMoles_final_.size();
    for (int isk = 0; isk < (int) RSD_List_.size(); isk++) {
        ReactingSurDomain* rsd = RSD_List_[isk];
        if (rsd) {
            int nr = rsd->nReactions();
            if (mR < nr) {
                mR = nr;
            }
        }
    }
    deltaG_.resize(mR);

    /*
     *  Set up the particleDiameter_ and particleNumberToFollow_ fields.
     *   If any geometry information is supplied in the input file, then the geometry is
     *   used to rescale the input mole numbers of species. If no geometry information
     *   is input, the the partcleNumbertoFollow_ variable is calculated to be consistent
     *   with the input mole numbers.
     */
    // particle diameter is a required parameter
    if (ei->particleDiameter <= 0.0) {
        throw CanteraError("Electrode::electrode_model_create()", "Particle diameter is <= 0");
    }
    inputParticleDiameter_ = ei->particleDiameter;
    Radius_exterior_final_ = inputParticleDiameter_ / 2.0;
    Radius_exterior_init_ = Radius_exterior_final_;
    Radius_exterior_init_init_ = Radius_exterior_final_;
    Radius_exterior_final_final_ = Radius_exterior_final_;
    double partVol = Pi * inputParticleDiameter_ * inputParticleDiameter_ * inputParticleDiameter_ / 6.0;
    bool rewriteMoleNumbers = false;
    /*
     *  Gross volume is the total superficial volume of the electrode. It is the sum of the
     *  solid volume and the liquid electrylte volume.
     *    - Compute the gross volume if this information has been supplied through the input file
     */
    double grossVolume = 0.0;
    if (ei->electrodeGrossArea > 0.0 && ei->electrodeGrossThickness > 0.0) {
        grossVolume = ei->electrodeGrossArea * ei->electrodeGrossThickness;
    } else if (ei->electrodeGrossDiameter > 0.0 && ei->electrodeGrossThickness > 0.0) {
        grossVolume = Pi * 0.25 * ei->electrodeGrossDiameter * ei->electrodeGrossDiameter * ei->electrodeGrossThickness;
    }

    /* Consider four cases
     * 1) specify the particle number and porosity --> calculate total superficial volume, renormalize solid mole numbers
     * 2) Specify the Total superficial volume of electrode and porosity --> calculate solid volume, renormalize solid mole numbers
     * 3) specify the mass/moles electrodes and Total superficial volume --> calculate porosity
     * 4) specify the mass/moles electrodes and porosity --> calculate total superficial volume
     */
    if (ei->particleNumberToFollow > 0.0) {
        particleNumberToFollow_ = ei->particleNumberToFollow;
        rewriteMoleNumbers = true;
        ElectrodeSolidVolume_ = particleNumberToFollow_ * partVol;
        if (ei->porosity > 0.0) {
            porosity_ = ei->porosity;
            grossVolume = ElectrodeSolidVolume_ / (1.0 - porosity_);
        } else {
            if (grossVolume > 0.0) {
                porosity_ = (grossVolume - ElectrodeSolidVolume_) / grossVolume;
            } else {
                throw CanteraError(" Electrode::electrode_model_create() ",
                                   " Need to either specify the total superficial volume or the porosity");
            }
            porosity_ = 0.3;
        }
        printf("Electrode::electrode_model_create \n\tWARNING: reset moles/masses "
               "based on \'Particle Number to Follow\' keyword.\n");
    } else if ((ei->porosity > 0.0) && (grossVolume > 0.0)) {
        porosity_ = ei->porosity;
        ElectrodeSolidVolume_ = grossVolume * (1.0 - porosity_);
        particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        rewriteMoleNumbers = true;
    } else if (Electrode::SolidVol() > 0.0) {
        if (grossVolume > 0.0) {
            ElectrodeSolidVolume_ = Electrode::SolidVol();
            porosity_ = (grossVolume - ElectrodeSolidVolume_) / grossVolume;
            if (porosity_ < 0.0) {
                throw CanteraError(" Electrode::electrode_model_create() ",
                                   " porosity_ is negative " + fp2str(porosity_));
            }
            particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        } else if (ei->porosity > 0.0) {
            ElectrodeSolidVolume_ = Electrode::SolidVol();
            particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
            grossVolume = ElectrodeSolidVolume_ / (1.0 - ei->porosity);
        } else {
            throw CanteraError(
                    " Electrode::electrode_model_create() ",
                    " no specification of cell volume through \'Electrode Porosity\' "
                            "or \'Electrode Thickness\', etc., keywords. grossVol gives " + fp2str(grossVolume));
        }
    } else {
        throw CanteraError(
                " Electrode::electrode_model_create() ",
                " no specification of solid volume through \'Particle Number to Follow\', "
                        " Porosity and Superficial Volume, "
                        "\'Phase Moles\', or \'Phase Mass\' keyword.  SolidVol() gives "
                        + fp2str(Electrode::SolidVol()));
    }

    /*
     *  Calculate the total number of moles of species in the electrode that are consistent
     *  with the input parameters, when the input parameters have specified the total solid volume of the electrode.
     */
    if (rewriteMoleNumbers) {
        Electrode::resizeMoleNumbersToGeometry();
    }

    phaseMoles_init_ = phaseMoles_final_;
    phaseMoles_init_init_ = phaseMoles_final_;
    phaseMoles_final_final_ = phaseMoles_final_;

    /*
     *  Resize the value of the surface areas to ones appropriate for the input geometry.
     */
    Electrode::resizeSurfaceAreasToGeometry();

    Electrode::updateState_OnionOut();

    spMoles_init_init_ = spMoles_final_;
    spMoles_final_final_ = spMoles_final_;

    ElectrodeBath* BG = ei->m_BG;

    for (iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();

        double* vv = BG->CapZeroDoDCoeffPhases[iph];
        double* vl = BG->CapLeftCoeffPhases[iph];
        for (int k = 0; k < nspPhase; k++) {
            int iGlobSpeciesIndex = kStart + k;

            capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = vl[k];
            capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = vv[k];

        }

        //ei->SpeciesNames+kstart,
        //nSpecies
    }

    for (iph = 0; iph < m_NumTotPhases; iph++) {

        if (PhaseNames_[iph] == "LiFeS_X_combo") {
            electrodeModelType_ = 4;
            break;
        }
        if (PhaseNames_[iph] == "FeS2(S)") {
            electrodeModelType_ = 2;
        }
        if (PhaseNames_[iph] == "Li13Si4(S)") {
            if (electrodeModelType_ != 3) {
                electrodeModelType_ = 1;
            }
        }
        if (PhaseNames_[iph] == "Li7Si3_Interstitial") {
            electrodeModelType_ = 3;
        }
        if (PhaseNames_[iph] == "MCMB_Interstitials_anode") {
            electrodeModelType_ = 5;
        }
    }

    /*
     *  Now that we have identified the model type, make sure our list of phases is what
     *  we expect from the model type.
     */
    int nRSD = RSD_List_.size();
    if (electrodeModelType_ == 1) {
        pmatch(PhaseNames_, 2, "Li13Si4(S)");
        pmatch(PhaseNames_, 3, "Li7Si3(S)");
        if (nRSD > 1) {
            pmatch(PhaseNames_, 4, "Li12Si7(S)");
        }
        if (nRSD > 2) {
            pmatch(PhaseNames_, 5, "Si(S)");
        }
    }
    if (electrodeModelType_ == 2) {
        pmatch(PhaseNames_, 2, "FeS2(S)");
        pmatch(PhaseNames_, 3, "Li3Fe2S4(S)");
        if (nRSD > 1) {
            pmatch(PhaseNames_, 4, "Li[2+x]Fe[1-x]S2(S)");
            pmatch(PhaseNames_, 5, "Fe[1-x]S(S)");
        }
        if (nRSD > 2) {
            pmatch(PhaseNames_, 6, "Li2S(S)");
            pmatch(PhaseNames_, 7, "Fe(S)");
        }
    }
    if (electrodeModelType_ == 4) {
        pmatch(PhaseNames_, 2, "FeS2(S)");
        pmatch(PhaseNames_, 3, "Li3Fe2S4(S)");
        if (nRSD > 1) {
            pmatch(PhaseNames_, 4, "LiFeS_X_combo");
        }
        if (nRSD > 2) {
            pmatch(PhaseNames_, 5, "Li2S(S)");
            pmatch(PhaseNames_, 6, "Fe(S)");
        }
    }

    if (electrodeModelType_ == 1) {
        setCapacityCoeff_LiSi();
    } else if (electrodeModelType_ == 2) {
        setCapacityCoeff_FeS2();
    } else if (electrodeModelType_ == 3) {
        setCapacityCoeff_LiSi_Li();
    } else if (electrodeModelType_ == 4) {
        setCapacityCoeff_FeS2_Combo();
    } else if (electrodeModelType_ == 5) {
        setCapacityCoeff_MCMB();
    }

    if (ei->RelativeCapacityDischargedPerMole != -1) {

        // setRelativeCapacityDischarged(ei->RelativeCapacityDischarged);
    }

    capacityInitialZeroDod_ = capacity();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];

    return 0;
}
//====================================================================================================================
int Electrode::setInitialConditions(ELECTRODE_KEY_INPUT* ei)
{
    int nspecies = 0;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {

        ThermoPhase* tphase = &(thermo(iph));
        int nSpecies = tphase->nSpecies();
        int kstart = nspecies;
        nspecies += nSpecies;

        /*
         *  We copy the mole numbers and mole fraction information from the bath gas section of the input file
         *  here. Note, the total quantitites of electrode material and electrolyte material may be overwritten
         *  later.
         */
        double sum = 0.0;
        for (int k = 0; k < nSpecies; k++) {
            spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] = ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
            sum += ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
        }

        if (sum > 0.0) {
            for (int k = 0; k < nSpecies; k++) {
                spMf_final_[kstart + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] / sum;
            }
            tphase->setMoleFractions(&(spMf_final_[kstart]));
        } else {
            tphase->getMoleFractions(&(spMf_final_[kstart]));
        }
        phaseVoltages_[iph] = tphase->electricPotential();
    }

    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];

    /* set the absolute tolerance to 1e-5 the total number of moles */
    double tMoles = 0.0;
    for (int k = 0; k < m_NumTotSpecies; k++) {
        tMoles += spMoles_final_[k];
    }
    molarAtol_ = tMoles * 1.0E-5;

    /*
     * Now that spMoles_final_ is sized appropriately, we can call updatePhaseNumbers()
     */
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        Electrode::updatePhaseNumbers(iph);
    }

    /*
     *  Set up the particleDiameter_ and particleNumberToFollow_ fields.
     *   If any geometry information is supplied in the input file, then the geometry is
     *   used to rescale the input mole numbers of species. If no geometry information
     *   is input, the the partcleNumbertoFollow_ variable is calculated to be consistent
     *   with the input mole numbers.
     */
    // particle diameter is a required parameter
    if (ei->particleDiameter <= 0.0) {
        throw CanteraError("Electrode::electrode_model_create()", "Particle diameter is <= 0");
    }
    inputParticleDiameter_ = ei->particleDiameter;
    Radius_exterior_final_ = inputParticleDiameter_ / 2.0;

    double partVol = Pi * inputParticleDiameter_ * inputParticleDiameter_ * inputParticleDiameter_ / 6.0;
    bool rewriteMoleNumbers = false;
    /*
     *  Gross volume is the total superficial volume of the electrode. It is the sum of the
     *  solid volume and the liquid electrylte volume.
     *    - Compute the gross volume if this information has been supplied through the input file
     */
    double grossVolume = 0.0;
    if (ei->electrodeGrossArea > 0.0 && ei->electrodeGrossThickness > 0.0) {
        grossVolume = ei->electrodeGrossArea * ei->electrodeGrossThickness;
    } else if (ei->electrodeGrossDiameter > 0.0 && ei->electrodeGrossThickness > 0.0) {
        grossVolume = Pi * 0.25 * ei->electrodeGrossDiameter * ei->electrodeGrossDiameter * ei->electrodeGrossThickness;
    }

    /* Consider four cases
     * 1) Specify the particle number and porosity --> calculate total superficial volume, renormalize solid mole numbers
     * 2) Specify the Total superficial volume of electrode and porosity --> calculate solid volume, renormalize solid mole numbers
     * 3) Specify the mass/moles electrodes and Total superficial volume --> calculate porosity
     * 4) Specify the mass/moles electrodes and porosity --> calculate total superficial volume
     */
    if (ei->particleNumberToFollow > 0.0) {
        particleNumberToFollow_ = ei->particleNumberToFollow;
        rewriteMoleNumbers = true;
        ElectrodeSolidVolume_ = particleNumberToFollow_ * partVol;
        if (ei->porosity > 0.0) {
            porosity_ = ei->porosity;
            grossVolume = ElectrodeSolidVolume_ / (1.0 - porosity_);
        } else {
            if (grossVolume > 0.0) {
                porosity_ = (grossVolume - ElectrodeSolidVolume_) / grossVolume;
            } else {
                throw CanteraError(" Electrode::electrode_model_create() ",
                                   " Need to either specify the total superficial volume or the porosity");
            }
            porosity_ = 0.3;
        }
        printf("Electrode::electrode_model_create \n\tWARNING: reset moles/masses "
               "based on \'Particle Number to Follow\' keyword.\n");
    } else if ((ei->porosity > 0.0) && (grossVolume > 0.0)) {
        porosity_ = ei->porosity;
        ElectrodeSolidVolume_ = grossVolume * (1.0 - porosity_);
        particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        rewriteMoleNumbers = true;
    } else if (SolidVol() > 0.0) {
        if (grossVolume > 0.0) {
            ElectrodeSolidVolume_ = SolidVol();
            porosity_ = (grossVolume - ElectrodeSolidVolume_) / grossVolume;
            if (porosity_ < 0.0) {
                throw CanteraError(" Electrode::electrode_model_create() ",
                                   " porosity_ is negative " + fp2str(porosity_));
            }
            particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        } else if (ei->porosity > 0.0) {
            ElectrodeSolidVolume_ = SolidVol();
            particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
            grossVolume = ElectrodeSolidVolume_ / (1.0 - ei->porosity);
        } else {
            throw CanteraError(
                    " Electrode::setInitialConditions() ",
                    " no specification of cell volume through \'Electrode Porosity\' "
                            "or \'Electrode Thickness\', etc., keywords. grossVol gives " + fp2str(grossVolume));
        }
    } else {
        throw CanteraError(
                " Electrode::setInitialConditions() ",
                " no specification of solid volume through \'Particle Number to Follow\', "
                        " Porosity and Superficial Volume, "
                        "\'Phase Moles\', or \'Phase Mass\' keyword.  SolidVol() gives " + fp2str(SolidVol()));
    }

    /*
     *  Calculate the total number of moles of species in the electrode that are consistent
     *  with the input parameters, when the input parameters have specified the total solid volume of the electrode.
     */
    if (rewriteMoleNumbers) {
        resizeMoleNumbersToGeometry();
    }

    /*
     *  Resize the value of the surface areas to ones appropriate for the input geometry.
     */
    resizeSurfaceAreasToGeometry();

    capacityInitialZeroDod_ = capacity();

    Electrode::updateState();
    Electrode::updateSurfaceAreas();

    /*
     *  New addition -> we should be able to put this here, but let's run test suite.
     */
    if (ei->RelativeCapacityDischargedPerMole != -1) {
        setRelativeCapacityDischargedPerMole(ei->RelativeCapacityDischargedPerMole);
    }

    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    /*
     *  Now transfer the final state to the init and init init states.
     */
    Electrode::setInitStateFromFinal_Oin(true);

    return 0;
}
//====================================================================================================================
int Electrode::electrode_stateSave_create(ELECTRODE_KEY_INPUT* ei)
{
    return 0;
}
//====================================================================================================================
// Set the sizes of the electrode from the input parameters
/*
 *  We resize all of the information within the electrode from the input parameters
 *
 * @param electrodeArea   Area of the electrode
 * @param electrodeThickness  Width of the electrode
 * @param porosity        Volume of the electrolyte phase
 */
void Electrode::setElectrodeSizeParams(doublereal electrodeArea, doublereal electrodeThickness, doublereal porosity)
{

    porosity_ = porosity;
    double grossVolume = electrodeArea * electrodeThickness;
    double partVol = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * Radius_exterior_final_ / 3.0;
    if (grossVolume > 0.0) {
        ElectrodeSolidVolume_ = (1.0 - porosity_) * grossVolume;
        particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        resizeMoleNumbersToGeometry();
        resizeSurfaceAreasToGeometry();
    } else {
        throw CanteraError("Electrode::setElectrodeSizeParams()", "current volume is zero or less");
    }
}
//====================================================================================================================
// Resize the solid phase and electrolyte mole numbers within the object
/*
 *
 *  (Initialization safe, because this function doesn't call child routines)
 *
 *  This routine uses particleDiameter_ , particleNumberToFollow_, and porosity_ to recalculate
 *  all the mole numbers in the electrode. This is done by rescaling all of the numbers.
 *  At the end of the process, the total volume of the electrode object is
 *
 *    grossVol = SolidVol() / ( 1.0 - porosity_)
 *
 *  where the SolidVol() is equal to
 *
 *   SolidVol() =  particleNumberToFollow_  Pi *  particleDiameter_**3 / 6.0;
 *
 */
void Electrode::resizeMoleNumbersToGeometry()
{

    double partVol = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * Radius_exterior_final_ / 3.0;
    double targetSolidVol = partVol * particleNumberToFollow_;
    double currentSolidVol = Electrode::SolidVol();
    if (currentSolidVol <= 0.0) {
        throw CanteraError("Electrode::resizeMoleNumbersToGeometry()", "current solid volume is zero or less");
    }
    double ratio = targetSolidVol / currentSolidVol;
    double tMoles = 0.0;
    for (int k = 0; k < m_NumTotSpecies; k++) {
        spMoles_final_[k] *= ratio;
        spMoles_init_[k] = spMoles_final_[k];
        spMoles_init_init_[k] = spMoles_final_[k];
        tMoles += spMoles_final_[k];
    }
    molarAtol_ = tMoles * 1.0E-5;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        Electrode::updatePhaseNumbers(iph);
    }
    phaseMoles_init_ = phaseMoles_final_;
    phaseMoles_init_init_ = phaseMoles_final_;
    phaseMoles_final_final_ = phaseMoles_final_;

    /*
     * Resize the electrolyte so that the total volume of the electrolyte is consistent with the given
     * porosity, porosity_
     */
    // Electrode::resizeSolutionNumbersToPorosity();
    currentSolidVol = Electrode::SolidVol();
    double grossVol = currentSolidVol / (1.0 - porosity_);
    double targetSolnVol = grossVol - currentSolidVol;
    double totalVol = Electrode::TotalVol();
    double currentSolnVol = totalVol - currentSolidVol;
    ratio = targetSolnVol / currentSolnVol;

    int istart = m_PhaseSpeciesStartIndex[solnPhase_];
    int nspSoln = m_PhaseSpeciesStartIndex[solnPhase_ + 1] - m_PhaseSpeciesStartIndex[solnPhase_];
    for (int kk = 0; kk < nspSoln; kk++) {
        int k = istart + kk;
        spMoles_final_[k] *= ratio;
        spMoles_init_[k] = spMoles_final_[k];
        spMoles_init_init_[k] = spMoles_final_[k];
        spMoles_final_final_[k] = spMoles_final_[k];
    }
    Electrode::updatePhaseNumbers(solnPhase_);
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];

    currentSolidVol = Electrode::SolidVol();
    totalVol = Electrode::TotalVol();
    double calcPor = (totalVol - currentSolidVol) / totalVol;
    if (fabs(calcPor - porosity_) > 1.0E-6) {
        throw CanteraError("Electrode::resizeMoleNumbersToGeometry() Error",
                           "Couldn't set the porosity correctly: " + fp2str(calcPor) + " vs " + fp2str(porosity_));
    }

    /*
     *  Update all of the states again
     */
    //Electrode::updateState();
}
//====================================================================================================================
//   Resize the electrolyte mole numbers so that it consistent with the porosity
/*
 *  (serial virtual function -> calls children)
 *
 *  This routine resizes the electrolyte mole numbers within the object to be consistent
 *  with the porosity. This has the effect of changing the total volume of the
 *  object.
 *
 *  @param porosityReset   If positive, this resets the porosity within the electrode
 *                         to the given value.
 */
void Electrode::resizeSolutionNumbersToPorosity(doublereal porosityReset)
{
    if (porosityReset > 0.0) {
        porosity_ = porosityReset;
    }
    double currentSolidVol = SolidVol();
    double grossVol = currentSolidVol / (1.0 - porosity_);
    double targetSolnVol = grossVol - currentSolidVol;
    double currentSolnVol = TotalVol() - currentSolidVol;
    double ratio = targetSolnVol / currentSolnVol;

    int istart = m_PhaseSpeciesStartIndex[solnPhase_];
    int nspSoln = m_PhaseSpeciesStartIndex[solnPhase_ + 1] - m_PhaseSpeciesStartIndex[solnPhase_];
    for (int kk = 0; kk < nspSoln; kk++) {
        int k = istart + kk;
        spMoles_final_[k] *= ratio;
        spMoles_init_[k] = spMoles_final_[k];
        spMoles_init_init_[k] = spMoles_final_[k];
    }
    updatePhaseNumbers(solnPhase_);

}
//====================================================================================================================
//    Resize the surface areas according to the input geometry.
/*
 *  We resize the surface areas of the Reacting Surfaces to a value which is
 *  equal to the particle exterior surface area multiplied by the number of particles.
 *
 *  Note we can't do anything better in this routine because we do not know the specifics
 *  of the particles.
 */
void Electrode::resizeSurfaceAreasToGeometry()
{
    double totalSurfArea = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * particleNumberToFollow_;

    double total = 0.0;
    double maxVal = -1.0;
    for (int i = 0; i < numSurfaces_; i++) {
        if (maxVal < surfaceAreaRS_final_[i]) {
            maxVal = surfaceAreaRS_final_[i];
        }
        total += surfaceAreaRS_final_[i];
    }

    if (total == 0.0) {
        for (int i = 0; i < numSurfaces_; i++) {
            surfaceAreaRS_final_[i] = totalSurfArea;
            surfaceAreaRS_init_[i] = totalSurfArea;
            surfaceAreaRS_init_init_[i] = totalSurfArea;
            surfaceAreaRS_final_final_[i] = totalSurfArea;

            if (RSD_List_[i]) {
                ActiveKineticsSurf_[i] = 1;
            } else {
                ActiveKineticsSurf_[i] = 0;
            }
        }
    } else {
        double ratio = totalSurfArea / maxVal;
        for (int i = 0; i < numSurfaces_; i++) {
            surfaceAreaRS_final_[i] *= ratio;
            surfaceAreaRS_init_[i] = surfaceAreaRS_final_[i];
            surfaceAreaRS_init_init_[i] = surfaceAreaRS_final_[i];
            surfaceAreaRS_final_final_[i] = surfaceAreaRS_final_[i];
            ActiveKineticsSurf_[i] = 0;
            if (surfaceAreaRS_final_[i] > 0.0) {
                if (RSD_List_[i]) {
                    ActiveKineticsSurf_[i] = 1;
                }
            }
        }
    }

}
//====================================================================================================================
// We change the particleNumbertoFollow_ field to comply with the number of moles and the particle diameter
void Electrode::resizeParticleNumbersToMoleNumbers()
{
    double partVol = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * Radius_exterior_final_ / 3.0;
    double currentSolidVol = SolidVol();
    particleNumberToFollow_ = currentSolidVol / partVol;
    resizeSurfaceAreasToGeometry();
}
//====================================================================================================================
Cantera::ReactingSurDomain* Electrode::currOuterReactingSurface()
{
    for (int isk = 0; isk < numSurfaces_; isk++) {
        if (ActiveKineticsSurf_[isk]) {
            return RSD_List_[isk];
        }
    }
    return 0;
}
//====================================================================================================================
Cantera::ReactingSurDomain* Electrode::reactingSurface(int iSurf)
{
    return RSD_List_[iSurf];
}
//================================================================================================
// Set the mole numbers in the electrolyte phase
/*
 *  We set the mole numbers of the electrolyte phase separately from
 *  the rest of the phases.
 *
 *  We always make sure that mole numbers are positive by clipping. We
 *  always make sure that mole fractions sum to one.
 *
 *  If we are not following mole numbers in the electrode, we set the
 *  total moles to the internal constant, electrolytePseudoMoles_, while
 *  using this vector to set the mole fractions.
 *
 * @param electrolyteMoleNum vector of mole numbers of the species in the
 *                     electrolyte phase.
 *                     units = kmol
 *                     size = number of species in the electrolyte phase
 *
 * @param setInitial   Boolean indicating that we should set the initial values of the
 *                     electrolyte mole numbers as well. 
 */
void Electrode::setElectrolyteMoleNumbers(const double* const electrolyteMoleNum, bool setInitial)
{
    int istart = m_PhaseSpeciesStartIndex[solnPhase_];
    ThermoPhase& tp = thermo(solnPhase_);
    int nsp = tp.nSpecies();
    AssertTrace(nsp == (m_PhaseSpeciesStartIndex[solnPhase_+1] - m_PhaseSpeciesStartIndex[solnPhase_]));
    double tmp = 0.0;
    for (int k = 0; k < nsp; k++) {
        spMoles_final_[istart + k] = MAX(electrolyteMoleNum[k], 0.0);
        tmp += spMoles_final_[istart + k];
    }
    phaseMoles_final_[solnPhase_] = tmp;

    if (tmp > 1.0E-200) {
        for (int k = 0; k < nsp; k++) {
            spMf_final_[istart + k] = spMoles_final_[istart + k] / tmp;
        }
        /*
         *  @TODO: take this section out. If we are setting the moles at the interface,
         *         then we should set the moles. See the explanation in the Electrode.h file.
         */
        if (!followElectrolyteMoles_) {
            for (int k = 0; k < nsp; k++) {
                spMoles_final_[istart + k] *= electrolytePseudoMoles_ / tmp;
            }
            phaseMoles_final_[solnPhase_] = electrolytePseudoMoles_;
        }
    }
    /*
     *  Update the internal state 
     */
    updateState();

    /*
     *  Set other states retroactively
     */
    if (setInitial) {
	if (pendingIntegratedStep_) {
	    throw CanteraError("Electrode::setElectrolyteMoleNumbers ERROR",
			       "Trying to set initial electroltye mole numbers");
	}
        for (int k = 0; k < nsp; k++) {
            spMoles_init_[istart + k] = spMoles_final_[istart + k];
            spMoles_init_init_[istart + k] = spMoles_final_[istart + k];
        }
        phaseMoles_init_[solnPhase_] = phaseMoles_final_[solnPhase_];
        phaseMoles_init_init_[solnPhase_] = phaseMoles_final_[solnPhase_];
        phaseMoles_final_final_[solnPhase_] = phaseMoles_final_[solnPhase_];
    

	/*
	 *  Set the init and the init_init state and the final_final state from the final state
	 */
	setInitStateFromFinal();
	setFinalFinalStateFromFinal();
    }
}
//================================================================================================
// Set the mole numbers in a single phase
/*
 *  We set the mole numbers of a single phase separately from
 *  the rest of the phases.
 *
 *  We always make sure that mole numbers are positive by clipping. We
 *  always make sure that mole fractions sum to one.
 *
 *  If we are not following mole numbers in the electrode, we set the
 *  total moles to the internal constant, electrolytePseudoMoles_, while
 *  using this vector to set the mole fractions.
 *
 * @param iph     Phase id.
 * @param moleNum vector of mole numbers of the species in the
 *                     electrolyte phase.
 *                     units = kmol
 *                     size = number of species in the electrolyte phase
 */
void Electrode::setPhaseMoleNumbers(int iph, const double* const moleNum)
{
    if (pendingIntegratedStep_) {
        if (iph != solnPhase_) {
            throw CanteraError("Electrode::setPhaseMoleNumbers", "called when there is a pending step");
        }
    }
    int istart = m_PhaseSpeciesStartIndex[iph];
    int nsp = m_PhaseSpeciesStartIndex[iph + 1] - istart;
    for (int k = 0; k < nsp; k++) {
        spMoles_final_[istart + k] = MAX(moleNum[k], 0.0);
        spMoles_init_[istart + k] = spMoles_final_[istart + k];
        spMoles_init_init_[istart + k] = spMoles_final_[istart + k];
    }
    updatePhaseNumbers(iph);

    spMoles_init_init_ = spMoles_final_;
    spMoles_init_ = spMoles_final_;
    spMoles_final_final_ = spMoles_final_;
    phaseMoles_init_init_ = phaseMoles_final_;
    phaseMoles_init_ = phaseMoles_final_;
    phaseMoles_final_final_ = phaseMoles_final_;

    /*
     *  Change the particleNumbersToFollow_ field
     */
    resizeParticleNumbersToMoleNumbers();

    /*
     *  Call the virtual function updateState. This may fill in further state information that
     *  is not known to this function. 
     */
    updateState();

    /*
     *  Set the init and the init_init state and the final_final state from the final state
     */
    setInitStateFromFinal(true);
    setFinalFinalStateFromFinal();
}

// -----------------------------------------------------------------------------------------------------------------
// ------------------------------------ CARRY OUT INTEGRATION OF EQUATIONS -----------------------------------------
// -----------------------------------------------------------------------------------------------------------------

//================================================================================================
void Electrode::setTime(double time)
{
    if (pendingIntegratedStep_) {
        throw CanteraError("Electrode::setTime()", "called when there is a pending step");
    }
    t_init_init_ = time;
    tfinal_ = t_init_init_;
    tinit_ = t_init_init_;
    t_final_final_ = t_init_init_;
}
//================================================================================================
double Electrode::getFinalTime() const
{
    return tfinal_;
}
//================================================================================================
// Update all mole numbers in the object from the mole numbers in the spMoles_final_[] vector
/*
 *  We use the field spMoles_final_[] to set the field phaseMoles_final_[].
 *  We set the mole numbers of a single phase separately from the rest of the phases.
 *  We do not clip the mole numbers to be positive. We allow negative mole numbers.
 *  We make sure that the mole fractions sum to one.
 *  The following fields are state variables
 *
 *  Summary: State Variables
 *              spMoles_final_[kGlobal]  : kGlobal in solid phase species
 *
 *              phaseMoles_final_[iph]  : When the phase moles are negative
 *              spMf_final_[kGlobal]
 *
 *
 *           Independent Variables
 *              spMoles_final_[kGlobal]  :    kGlobal in electrolyte phase
 *              phaseVoltages_[iph] :         iph in solid phase electrodes
 *              Temperature
 *              Pressure
 *
 *            spElectroChemPot_[]
 *
 *   The following fields in this object are set:
 *            phaseMoles_final_[iph]          iph in solid phase electrode
 *            spMf_final_[kGlobal]            kGlobal in solid phase species
 *            VolPM_[kGlobal]                 kGlobal in solid phase species
 *            spElectroChemPot_[kGlobal]      kGlobal in solid phase species
 *            ThermoPhase[iph]                iph in solid phase electrode
 *           	phaseMolarVolumes_[iph]         iph in solid phase electrode
 *
 *  If we are not following  the mole numbers in the electrode, we set the
 *  total moles to the internal constant, electrolytePseudoMoles_, while
 *  using this vector to set the mole fractions, using the ThermoPhase object
 *  to get the mole fractions.
 *
 * A Note about negative Phase Numbers.
 * -----------------------------------------------
 *     The algorithms shouldn't rely on this routine to handle negative mole numbers or negative phase numbers.
 *  The algorithms should already have arranged this issue before calling this routine.
 *  However, this routine does have a method for handling this. For positive mole numbers, the spMoles vector
 *  is treated as the state vector. For negative mole numbers,
 *
 *
 * @param iph     Phase id.
 */
void Electrode::updatePhaseNumbers(int iph)
{
    int istart = m_PhaseSpeciesStartIndex[iph];
    ThermoPhase& tp = thermo(iph);
    int nsp = m_PhaseSpeciesStartIndex[iph + 1] - istart;
    double tmp = 0.0;
    for (int k = 0; k < nsp; k++) {
        tmp += spMoles_final_[istart + k];
    }
    phaseMoles_final_[iph] = tmp;
    if (tmp > 1.0E-200) {
        for (int k = 0; k < nsp; k++) {
            spMf_final_[istart + k] = spMoles_final_[istart + k] / tmp;
        }
#ifdef OLD_FOLLOW
        if (iph == solnPhase_) {
            if (!followElectrolyteMoles_) {
                ThermoPhase& tp = thermo(solnPhase_);
                tp.getMoleFractions(&spMf_final_[istart]);
                for (int k = 0; k < nsp; k++) {
                    spMoles_final_[istart + k] = electrolytePseudoMoles_ * spMf_final_[istart + k];
                }
                phaseMoles_final_[iph] = electrolytePseudoMoles_;
            }
        }
#endif
    } else if (tmp < -1.0E-200) {
        /*
         *  We are here when the phase moles are less than zero. This is a permissible condition when we
         *  are calculating a phase death, as we solve for the time when the phase moles is exactly zero.
         *  letting the phaseMoles go negative during the solve.
         */
        //if (phaseMoles_final_[iph] < 0.0) {
        //    throw CanteraError("Electrode::updatePhaseNumbers()", "unexpected condition");
        //}
        /*
         * For the current method the spMf_final_[] contains the valid information.
         */
        for (int k = 0; k < nsp; k++) {
            tmp = spMf_final_[istart + k];
            if (tmp < 0.0 || tmp > 1.0) {
                throw CanteraError("Electrode::updatePhaseNumbers()",
                                   "Mole fractions out of bounds:" + int2str(k) + " " + fp2str(tmp));
            }
            // TODO: this may be wrong. Experiment with leaving this alone
            spMoles_final_[istart + k] = 0.0;
        }
    } else {
        // We are here when the mole numbers of the phase are zero. In this case, we respect the
	// mole fraction vector in spMf_final_[]. There are cases where we are solving for the
	// mole fraction vector when the phase is zero, and we are just trying to get a 1 on the
	// diagonal. In this case we need to respect  spMF_final_[].
	// However we check for out of bounds values.
        for (int k = 0; k < nsp; k++) {
            tmp = spMf_final_[istart + k];
            if (tmp < 0.0 || tmp > 1.0) {
                throw CanteraError("Electrode::updatePhaseNumbers()",
                                   "Mole fractions out of bounds:" + int2str(k) + " " + fp2str(tmp));
            }
            spMoles_final_[istart + k] = 0.0;
            phaseMoles_final_[iph] = 0.0;
        }

    }

    // Here we set the state within the phase object
    tp.setState_TPX(temperature_, pressure_, &spMf_final_[istart]);
    tp.setElectricPotential(phaseVoltages_[iph]);
    tp.getPartialMolarVolumes(&(VolPM_[istart]));
    tp.getElectrochemPotentials(&(spElectroChemPot_[istart]));
    if (iph < NumVolPhases_) {
        phaseMolarVolumes_[iph] = tp.molarVolume();
    } else {
        phaseMolarVolumes_[iph] = 0.0;
        int isurf = iph - NumVolPhases_;
        sphaseMolarAreas_[isurf] = tp.molarVolume();
    }
}
//================================================================================================
void Electrode::updatePhaseNumbersTmp(vector<doublereal>& spMoles_tmp, vector<doublereal>& phaseMoles_tmp,
        vector<doublereal>& spMf_tmp)
{
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        int istart = m_PhaseSpeciesStartIndex[iph];
        // ThermoPhase &tp = thermo(solnPhase_);
        int nsp = m_PhaseSpeciesStartIndex[iph + 1] - istart;
        double tmp = 0.0;
        for (int k = 0; k < nsp; k++) {
            tmp += spMoles_tmp[istart + k];
        }
        phaseMoles_tmp[iph] = tmp;
        if (tmp > 1.0E-200) {
            for (int k = 0; k < nsp; k++) {
                spMf_tmp[istart + k] = spMoles_tmp[istart + k] / tmp;
            }
            if (iph == solnPhase_) {
#ifdef OLD_FOLLOW
                if (!followElectrolyteMoles_) {
                    ThermoPhase& tp = thermo(solnPhase_);
                    tp.getMoleFractions(&spMf_tmp[istart]);
                    for (int k = 0; k < nsp; k++) {
                        spMoles_tmp[istart + k] = electrolytePseudoMoles_ * spMf_tmp[istart + k];
                    }
                    phaseMoles_tmp[iph] = electrolytePseudoMoles_;
                }
#endif
            }
        } else {
            // We are here when the mole numbers of the phase are zero. In this case, we still need
            // a valid mole fraction vector. This is kept in the ThermoPhase object.
            ThermoPhase& tp = thermo(iph);
            tp.getMoleFractions(&(spMf_tmp[istart]));
        }
    }
}
//====================================================================================================================
// Calculates the change in the surface area of all external and internal interfaces within the electrode
/*
 *  (virtual)
 *  variables to be potentially altered
 *   surfaceAreaRS_[];
 *   isExternalSurface[]
 *   numExternalInterfacialSurfaces_;
 */
double Electrode::calcSurfaceAreaChange(double deltaT)
{
    double sa_final = surfaceAreaRS_init_[0];
    return sa_final;
}

//=====================================================================================================================
double Electrode::updateElectrolytePseudoMoles()
{
    /*
     * Call updateState to update spMf_final_ and phaseMoles_final_
     *   We also save the value of   electrolytePseudoMoles_  here.
     */
    if (pendingIntegratedStep_) {
        throw CanteraError(" Electrode::updateElectrolytePseudoMoles()",
                           " call is illegal when there is a pending step");
    }
    int origFEM = followElectrolyteMoles_;
    followElectrolyteMoles_ = 1;
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    double currentSolidVol = SolidVol();
    double totalVol = TotalVol();
    double currentSolnVol = totalVol - currentSolidVol;
    phaseMoles_final_[solnPhase_] = currentSolnVol / phaseMolarVolumes_[solnPhase_];
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    followElectrolyteMoles_ = 0;
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];

    followElectrolyteMoles_ = 0;
    followElectrolyteMoles_ = origFEM;
    return electrolytePseudoMoles_;
}
//====================================================================================================================
void Electrode::turnOffFollowElectrolyteMoles()
{
    //if (followElectrolyteMoles_ == 1) {
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    //}

    followElectrolyteMoles_ = 0;
    if (electrolytePseudoMoles_ <= 0.0) {
        throw CanteraError("Electrode::turnOffFollowElectrolyteMoles()", "electrolyte moles set negative or zero");
    }
    updatePhaseNumbers(solnPhase_);
}
//====================================================================================================================
void Electrode::turnOnFollowElectrolyteMoles()
{
    if (followElectrolyteMoles_ == 0) {
        phaseMoles_final_[solnPhase_] = electrolytePseudoMoles_;
    }

    followElectrolyteMoles_ = 1;
    updatePhaseNumbers(solnPhase_);

    if (phaseMoles_final_[solnPhase_] <= 0.0) {
        if (electrolytePseudoMoles_ > 0.0) {
            followElectrolyteMoles_ = 0;
            updatePhaseNumbers(solnPhase_);
            followElectrolyteMoles_ = 1;
        }
        if (phaseMoles_final_[solnPhase_] <= 0.0) {
            throw CanteraError("Electrode::turnOnFollowElectrolyteMoles()", "electrolyte moles set negative or zero");
        }
    }

}
//====================================================================================================================
// Get the ID of the external surface, if there is one
const std::vector<bool>& Electrode::getExternalSurfaceBooleans() const
{
    return isExternalSurface_;
}
//====================================================================================================================
doublereal Electrode::phaseMoles(int iph) const
{
    return phaseMoles_final_[iph];
}
//================================================================================================
doublereal Electrode::elementMoles(int ieG) const
{
    const Elements* elemList = getGlobalElements();
    std::string eN = elemList->elementName(ieG);
    return elementMoles(eN);
}
//================================================================================================
doublereal Electrode::elementMoles(std::string eN) const
{
    double sum = 0.0;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        const ThermoPhase* tp_ptr = &thermo(iph);
        size_t ie = tp_ptr->elementIndex(eN);
        if (ie != npos) {
            int nsp = thermo(iph).nSpecies();
            int kStart = m_PhaseSpeciesStartIndex[iph];
            for (int ik = 0; ik < nsp; ik++) {
                int na = thermo(iph).nAtoms(ik, ie);
                sum += na * spMoles_final_[kStart + ik];
            }
        }
    }
    return sum;
}
//================================================================================================
doublereal Electrode::elementSolidMoles(int ieG) const
{
    const Elements* elemList = getGlobalElements();
    std::string eN = elemList->elementName(ieG);
    return elementMoles(eN);
}
//================================================================================================
doublereal Electrode::elementSolidMoles(std::string eN) const
{
    double sum = 0.0;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph != solnPhase_) {
            const ThermoPhase* tp_ptr = &thermo(iph);
            size_t ie = tp_ptr->elementIndex(eN);
            if (ie != npos) {
                int nsp = thermo(iph).nSpecies();
                int kStart = m_PhaseSpeciesStartIndex[iph];
                for (int ik = 0; ik < nsp; ik++) {
                    double na = thermo(iph).nAtoms(ik, ie);
                    sum += na * spMoles_final_[kStart + ik];
                }
            }
        }
    }
    return sum;
}
//================================================================================================
double Electrode::speciesElectrochemPotential(int iGlobalSpIndex) const
{
    // DANGER -> Check to see this is updated

    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tphase = &(thermo(iph));
        int kStart = m_PhaseSpeciesStartIndex[iph];
        tphase->getElectrochemPotentials(&(spElectroChemPot_[kStart]));
    }
    return spElectroChemPot_[iGlobalSpIndex];
}
//================================================================================================
double Electrode::speciesChemPotential(int iGlobalSpIndex) const
{
    // DANGER -> Check to see this is updated

    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tphase = &(thermo(iph));
        int kStart = m_PhaseSpeciesStartIndex[iph];
        tphase->getChemPotentials(&(spElectroChemPot_[kStart]));
    }
    double tmp = spElectroChemPot_[iGlobalSpIndex];
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tphase = &(thermo(iph));
        int kStart = m_PhaseSpeciesStartIndex[iph];
        tphase->getElectrochemPotentials(&(spElectroChemPot_[kStart]));
    }
    return tmp;
}
//================================================================================================
void Electrode::getMoleFractions(doublereal* const x) const
{
    std::copy( spMf_final_.begin(), spMf_final_.end(), x );
}
//================================================================================================
void Electrode::getMoleNumSpecies(doublereal* const n) const
{
    std::copy( spMoles_final_.begin(), spMoles_final_.end(), n );
}
//================================================================================================
const std::vector<doublereal> & Electrode::getMoleNumSpecies() const
{
  return spMoles_final_;
}
//================================================================================================
void Electrode::getMoleNumPhases(doublereal* const np) const
{
    std::copy( phaseMoles_final_.begin(), phaseMoles_final_.end(), np);
}
//================================================================================================
double Electrode::moleFraction(int globalSpeciesIndex) const
{
    return spMf_final_[globalSpeciesIndex];
}
//================================================================================================
double Electrode::moleNumSpecies(int globalSpeciesIndex) const
{
    return spMoles_final_[globalSpeciesIndex];
}
//================================================================================================
// Returns the index of a phase in the ReactionSurfaceDomain object
// given the index of that phase in the PhaseList object
/*
 * @param PLph index of the phase in the PhaseList object, which is also the
 *             Electrode_Model object.
 *
 *  @return  Returns the index of the phase in the current ReactingSurDomain
 *           object. A value of -1 in this slot means that the phase doesn't
 *           participate in the  current ReactingSurDomain object
 */
int Electrode::ReactingSurfacePhaseIndex(int isk, int PLph) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        throw CanteraError("ReactingSurfacePhaseIndex", "Reacting surface not found");
    }
    return rsd->PLtoKinPhaseIndex_[PLph];
}
//================================================================================================
void Electrode::setState_TP(doublereal temperature, doublereal pressure)
{
    temperature_ = temperature;
    pressure_ = pressure;
    updateState();
}
//================================================================================================
double Electrode::temperature() const
{
    return temperature_;
}
//================================================================================================
double Electrode::pressure() const
{
    return pressure_;
}
//====================================================================================================================
double Electrode::SolidVol() const
{
    double vol = 0.0;
    for (int iph = 0; iph < NumVolPhases_; iph++) {
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        if (iph != solnPhase_) {
            for (int k = 0; k < nspPhase; k++) {
                vol += spMoles_final_[kStart + k] * VolPM_[kStart + k];
            }
        }
    }
    return vol;
}
//====================================================================================================================
double Electrode::SolidTotalMoles() const
{
    double tot = 0.0;
    for (int iph = 0; iph < NumVolPhases_; iph++) {
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        if (iph != solnPhase_ && (iph != metalPhase_)) {
            for (int k = 0; k < nspPhase; k++) {
                tot += spMoles_final_[kStart + k];
            }
        }
    }
    return tot;
}
//====================================================================================================================
double Electrode::TotalVol(bool ignoreErrors) const
{
    double vol = 0.0;

    for (int iph = 0; iph < NumVolPhases_; iph++) {
        double psum = 0.0;
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        for (int k = 0; k < nspPhase; k++) {
            vol += spMoles_final_[kStart + k] * VolPM_[kStart + k];
            psum += spMoles_final_[kStart + k] * VolPM_[kStart + k];
        }
        double mv = tp.molarVolume();
        double palt = mv * phaseMoles_final_[iph];
        if (!ignoreErrors) {
            if (palt < -1.0E-15) {
                throw CanteraError(" Electrode::TotalVol() ", " phase volume is negative " + fp2str(palt));
            }
            if (psum < -1.0E-15) {
                throw CanteraError(" Electrode::TotalVol() ", " phase volume is negative " + fp2str(psum));
            }
        }
        double denom = palt + psum + 1.0E-9;
        if (tp.eosType() != cLattice) {
            if (!ignoreErrors && 0) {
                if (fabs((palt - psum) / denom) > 1.0E-4) {
                    throw CanteraError(" Electrode::TotalVol() ",
                                       " internal inconsistency " + fp2str(palt) + " " + fp2str(psum));
                }
            }
        }
    }
    return vol;
}
//====================================================================================================================
void Electrode::getPhaseVol(double* const phaseVols) const
{
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        phaseVols[iph] = 0.0;
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        for (int k = 0; k < nspPhase; k++) {
            phaseVols[iph] += spMoles_final_[kStart + k] * VolPM_[kStart + k];
        }
        if (iph >= NumVolPhases_) {
            phaseVols[iph] = 0.0;
        }
    }
}
//====================================================================================================================
int Electrode::nSurfaces() const
{
    return numSurfaces_;
}
//====================================================================================================================
int Electrode::nReactions(int isk) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        return 0;
    } else {
        return rsd->rmcVector.size();
    }
}
//====================================================================================================================
void Electrode::getSurfaceAreas(double* const surfArea) const
{
    for (int i = 0; i < numSurfaces_; i++) {
        surfArea[i] = surfaceAreaRS_final_[i];
    }
}
//====================================================================================================================
// ----------------------------------- QUERY AND SET VOLTAGES --------------------------------------------------------
//====================================================================================================================
// This sets the metal and solution voltage
/*
 *  This sets the metal and the electrolyte voltages
 *  This assumes that there are only two voltages in the system.
 *   The voltage of the interface is defined as VOLTS = phiMetal - phiElectrolyte
 *
 *  @param[in]  phiMetal     Potential of the metal
 *  @param[in]  phElectrolyte Potential of the electrolyte
 */
void Electrode::setVoltages(const double phiMetal, const double phiElectrolyte)
{
    phaseVoltages_[metalPhase_] = phiMetal;
    phaseVoltages_[solnPhase_] = phiElectrolyte;
    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];
    updateState();
}
//====================================================================================================================
double Electrode::voltage() const
{
  return deltaVoltage_;
}
//====================================================================================================================
// Return the voltage of a phase
/*
 * @param iph  Phase id
 *
 * @return Returns the voltage in volts
 */
double Electrode::phaseVoltage(int iph) const
{
    return phaseVoltages_[iph];
}
//====================================================================================================================
// Set the voltage of a phase
/*
 * @param iph  Phase id
 * @param volts volts
 */
void Electrode::setPhaseVoltage(int iph, double volts)
{
    phaseVoltages_[iph] = volts;
    updateState();
}
//====================================================================================================================
//   Take the state (i.e., the final state) within the Electrode object and push it down
//   to the ThermoPhase objects and other variables that are part of the Electrode object
/*
 *  (virtual function from Electrode)
 *  This virtual function should be written so that child routines are not called by parent routines.
 *  It may be the case that child routines will surplant parent routines, in order to create efficiency.
 *
 *  We take the values of spMoles_final_[], and the number of particles, and their size,
 *  which are the default specification of the state variables within the Electrode object,
 *  and propagate them down to the ThermoPhase objects in the electrode.
 *  We also calculate the volumetric properties of the Electrode, the phase moles,
 *  and the mole fractions, and the external radius of the particle
 *
 *  All of these properties are defined for the _final_ state.
 *
 *  Thus, this is the main routine that reconciles all of the state information within the object.
 *  At the end of this routine, all aspects of the final state are consistent with each other.
 *
 *  prerequisites: The object must have been already been created.
 *
 *
 *  Summary: State Variables
 *              spMoles_final_[kGlobal]  : kGlobal in solid phase species
 *              particleNumberToFollow_
 *
 *           Independent Variables
 *              spMoles_final_[kGlobal]  :   kGlobal in electrolyte phase
 *              phaseVoltages_[iph] :        iph in solid phase electrodes
 *              Temperature
 *              Pressure
 *
 *   Dependent StateVariables
 *              phaseMoles_final_[iph]       iph in solid phase electrode
 *              spMf_final_[kGlobal]         kGlobal in solid phase species
 *              VolPM_[kGlobal]              kGlobal in solid phase species
 *              spElectroChemPot_[kGlobal]   kGlobal in solid phase species
 *              ThermoPhase[iph]             iph in solid phase electrode
 *        	    phaseMolarVolumes_[iph]      iph in solid phase electrode
 *              ElectrodeSolidVolume_        (this is called from SolidVolume();
 *              Radius_exterior_final_       (this is calculated from ElectrodeSolidVolume_)
 *
 */
void Electrode::updateState()
{
    int iph;

    /*
     * Loop over all phases in the object, volume and surface phases
     */
    for (iph = 0; iph < m_NumTotPhases; iph++) {
        updatePhaseNumbers(iph);
    }
    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];
    /*
     * Calculate the volume of the electrode phase. This is the main routine to do this.
     */
    ElectrodeSolidVolume_ = SolidVol();

    double vol = ElectrodeSolidVolume_ / particleNumberToFollow_;

    /*
     *  Calculate the external radius of the particle (in meters) assuming that all particles are the same
     *  size and are spherical
     */
    Radius_exterior_final_ = pow(vol * 3.0 / (4.0 * Pi), 0.3333333333333333);
}
//====================================================================================================================
//
/*
 *  Note this version doesn't call child member functions
 */
void Electrode::updateState_OnionOut()
{
    int iph;

    /*
     * However, I want to make sure mole fractions are
     * consistent with final moles.
     */
    for (int i = 0; i < m_NumTotPhases; i++) {
        Electrode::updatePhaseNumbers(i);
    }

    /*
     * Loop over all phases in the object
     */
    for (iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tphase = &(thermo(iph));
        std::string pName = tphase->id();
        int kStart = m_PhaseSpeciesStartIndex[iph];
        /*
         * A direct call to tphase to set the state is unnecessary,
         * because m_mp will do it anyway.
         */
        tphase->setState_TPX(temperature_, pressure_, &spMf_final_[kStart]);
        tphase->setElectricPotential(phaseVoltages_[iph]);

        /*
         * Ok, we have set the state. Now upload the Vol and mu's.
         */
        tphase->getPartialMolarVolumes(&(VolPM_[kStart]));
        tphase->getElectrochemPotentials(&(spElectroChemPot_[kStart]));

        if (iph < NumVolPhases_) {
            phaseMolarVolumes_[iph] = tphase->molarVolume();
        } else {
            phaseMolarVolumes_[iph] = 0.0;
        }

    }
    /*
     * Calculate the volume of the electrode phase. This is the main routine to do this.
     */
    ElectrodeSolidVolume_ = Electrode::SolidVol();

    double vol = ElectrodeSolidVolume_ / particleNumberToFollow_;

    /*
     *  Calculate the external radius of the particle (in meters) assuming that all particles are the same
     *  size and are spherical
     */
    Radius_exterior_final_ = pow(vol * 3.0 / (4.0 * Pi), 0.3333333333333333);
}
//====================================================================================================================
//  Recalcualte the surface areas of the surfaces for the final state
/*
 * (virtual function from Electrode)
 *
 *    Here, we don't know anything about the morphology of the particle. Models that do know will override this
 *    function.
 *    Hwere we assume that the surface area is equal to the exterior surface area multiplied by the number
 *    of particles.
 *    We also respect keeping the surface area as zero, if it is initially set to zero.
 */
void Electrode::updateSurfaceAreas()
{
    double totalSurfArea = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * particleNumberToFollow_;
    for (int i = 0; i < numSurfaces_; i++) {
        if (surfaceAreaRS_final_[i] != 0.0) {
            surfaceAreaRS_final_[i] = totalSurfArea;
        }
    }
    /*
     * Here we assume that the surface just creates species enough to fill the available surface sites
     */
    int i = 0;
    for (int iph = NumVolPhases_; iph < m_NumTotPhases; iph++, i++) {
        ThermoPhase* tphase = &(thermo(iph));
        sphaseMolarAreas_[i] = tphase->molarVolume();
        int kstart = m_PhaseSpeciesStartIndex[iph];
        int nsp = m_PhaseSpeciesStartIndex[iph + 1] - kstart;
        phaseMoles_final_[iph] = surfaceAreaRS_final_[i] / sphaseMolarAreas_[i];
        for (int k = 0; k < nsp; k++) {
            spMoles_final_[k + kstart] = spMf_final_[k] * phaseMoles_final_[iph];
        }
    }
}
//====================================================================================================================
//   This is used to set the phase information that is implicit but not set by a restart or an initialization
/*
 *  (virtual function from Electrode)
 *
 *  This is called immediately after the restart file's contents are loaded into the electrode object.
 *  We then call this function to calculate the internal flags. Then, we call updateState() to make sure
 *  all information about the state is self-consistent.
 *
 *  @param flagErrors If true any changes in the current flags caused by a mismatch between the state
 *                    and the values of the flags will cause an error exit.
 */
bool Electrode::stateToPhaseFlagsReconciliation(bool flagErrors)
{
    /*
     *  We probably need to set the phase existance for Reacting Surfaces.
     */
    setPhaseExistenceForReactingSurfaces(false);
    return false;
}
//====================================================================================================================
//   Reactant stoichiometric coefficient
/*
 * Get the reactant stoichiometric coefficient for the kth global species
 * in the ith reaction of the reacting surface domain with index isk.
 */
double Electrode::reactantStoichCoeff(const int isk, int kGlobal, int i) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    int krsd = rsd->PLtoKinSpeciesIndex_[kGlobal];
    if (krsd == -1) {
        return 0.0;
    }
    double rst = rsd->reactantStoichCoeff(krsd, i);
    return rst;
}
//====================================================================================================================
//   Reactant stoichiometric coefficient
/*
 * Get the reactant stoichiometric coefficient for the kth global species
 * in the ith reaction of the reacting surface domain with index isk.
 */
double Electrode::productStoichCoeff(const int isk, int kGlobal, int i) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    int krsd = rsd->PLtoKinSpeciesIndex_[kGlobal];
    if (krsd == -1) {
        return 0.0;
    }
    double rst = rsd->productStoichCoeff(krsd, i);
    return rst;
}

//====================================================================================================================
// Get the net production rates of all species in the electrode object
// at the current conditions
/*
 *
 *  This routine assumes that the underlying objects have been updated
 */
void Electrode::getNetProductionRates(const int isk, doublereal* const net) const
{
    std::fill_n(net, m_NumTotSpecies, 0.);

    /*
     *  This routine basically translates between species lists for the reacting surface
     *  domain and the Electrode.
     *  Later, when we have more than one reacting surface domain in the electrode object,
     *  this will do a lot more
     */

    /*
     *  For each Reacting surface
     */
    if (ActiveKineticsSurf_[isk]) {
        /*
         *  Get the species production rates for the reacting surface
         */
        const vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();

        /*
         *  loop over the phases in the reacting surface
         */
        int nphRS = RSD_List_[isk]->nPhases();
        int jph, kph;
        int kIndexKin = 0;
        for (kph = 0; kph < nphRS; kph++) {
            jph = RSD_List_[isk]->kinOrder[kph];
            int istart = m_PhaseSpeciesStartIndex[jph];
            int nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
            for (int k = 0; k < nsp; k++) {
                net[istart + k] += rsSpeciesProductionRates[kIndexKin];
                kIndexKin++;
            }
        }
    }
}
//====================================================================================================================
void Electrode::getIntegratedProductionRates(doublereal* const net) const
{
    if (pendingIntegratedStep_ != 1) {
        throw CanteraError(" Electrode:integratedProductionRatesCurrent()",
                           " No integration step has been carried out");
    }
    double invDelT = 1.0;
    if (t_final_final_ > t_init_init_) {
        invDelT = 1.0 / (t_final_final_ - t_init_init_);
    }
    for (int k = 0; k < m_NumTotSpecies; k++) {
        net[k] = invDelT * spMoleIntegratedSourceTerm_[k];
    }
}
//====================================================================================================================
//  Returns the current and the net production rates per unit surface area of surface isk
/*
 *    columb sec-1 m-2
 */
double Electrode::getNetProductionRatesCurrent(const int isk, doublereal* const net) const
{
    getNetProductionRates(isk, net);
    // kmol sec-1 m-2
    double Eprod = net[kElectron_];
    // coulomb / kmol
    return Eprod * Faraday;
}
//====================================================================================================================
//  Return the global species index of the electron in the Electrode object
int Electrode::kSpecElectron() const
{
    return kElectron_;
}
//====================================================================================================================
double Electrode::getIntegratedProductionRatesCurrent(doublereal* const net) const
{
    if (pendingIntegratedStep_ != 1) {
        double* netOne = &deltaG_[0];
        std::fill_n(net, m_NumTotSpecies, 0.);
        for (int isk = 0; isk < numSurfaces_; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                getNetProductionRates(isk, netOne);
                for (int k = 0; k < m_NumTotSpecies; k++) {
                    net[k] += netOne[k] * surfaceAreaRS_final_[isk];
                }
            }
        }
        double Eprod = net[kElectron_];
        return Eprod * Faraday;
    }
    getIntegratedProductionRates(net);
    // kmol sec-1 m-2
    double Eprod = net[kElectron_];
    // coulomb / kmol
    return Eprod * Faraday;
}
//====================================================================================================================
double Electrode::integratedCurrent() const
{
    if (pendingIntegratedStep_ != 1) {
        std::vector<double> net(m_NumTotSpecies, 0.0);
        double* netOne = &deltaG_[0];
        for (int isk = 0; isk < numSurfaces_; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                getNetProductionRates(isk, netOne);
                for (int k = 0; k < m_NumTotSpecies; k++) {
                    net[k] += netOne[k] * surfaceAreaRS_final_[isk];
                }
            }
        }
        double Eprod = net[kElectron_];
        // coulomb / kmol
        return Eprod * Faraday;
    }
    double invDelT = 1.0;
    if (t_final_final_ > t_init_init_) {
        invDelT = 1.0 / (t_final_final_ - t_init_init_);
    }
    // kmol sec-1 m-2
    double Eprod = invDelT * spMoleIntegratedSourceTerm_[kElectron_];
    // coulomb / kmol
    return Eprod * Faraday;
}
//====================================================================================================================
double Electrode::integratedLocalCurrent() const
{
    if (pendingIntegratedStep_ != 1) {
        std::vector<double> net(m_NumTotSpecies, 0.0);
        double* netOne = &deltaG_[0];
        for (int isk = 0; isk < numSurfaces_; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                getNetProductionRates(isk, netOne);
                for (int k = 0; k < m_NumTotSpecies; k++) {
                    net[k] += netOne[k] * surfaceAreaRS_final_[isk];
                }
            }
        }
        double Eprod = net[kElectron_];
        // coulomb / kmol
        return Eprod * Faraday;
    }
    double invDelT = 1.0;
    if (t_final_final_ > t_init_init_) {
        invDelT = 1.0 / (tfinal_ - tinit_);
    }
    // kmol sec-1
    double Eprod = invDelT * spMoleIntegratedSourceTermLast_[kElectron_];
    // returns coulomb / sec
    return Eprod * Faraday;
}
//====================================================================================================================
//  Returns the current and the net production rates of the phases in kg/m2/s from a single surface
/*
 *  Returns the net production rates of all phases from reactions on a single surface
 *
 *  @param isk Surface ID to get the fluxes from.
 *  @param phaseMassFlux  Returns the mass fluxes of the phases
 */
void Electrode::getPhaseMassFlux(const int isk, doublereal* const phaseMassFlux)
{
    std::fill_n(phaseMassFlux, m_NumTotPhases, 0.);

    if (ActiveKineticsSurf_[isk]) {
        const vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();

        int nphRS = RSD_List_[isk]->nPhases();
        int kIndexKin = 0;
        for (int kph = 0; kph < nphRS; kph++) {
            int jph = RSD_List_[isk]->kinOrder[kph];
            ThermoPhase& tp = thermo(jph);
            int istart = m_PhaseSpeciesStartIndex[jph];
            int nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
            for (int k = 0; k < nsp; k++) {
                double net = rsSpeciesProductionRates[kIndexKin];
                double mw = tp.molecularWeight(k);
                phaseMassFlux[jph] += net * mw;
                kIndexKin++;
            }
        }
    }
}
//================================================================================================
//  Returns the current and the net production rates of the phases in kmol/m2/s from a single surface
/*
 *  Returns the net production rates of all phases from reactions on a single surface
 *
 *  @param isk Surface ID to get the fluxes from.
 *  @param phaseMassFlux  Returns the mass fluxes of the phases
 */
void Electrode::getPhaseMoleFlux(const int isk, doublereal* const phaseMoleFlux)
{
    std::fill_n(phaseMoleFlux, m_NumTotPhases, 0.);
    if (ActiveKineticsSurf_[isk]) {
        const vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
        int nphRS = RSD_List_[isk]->nPhases();
        int PLph, RSph;
        int kIndexKin = 0;
        for (RSph = 0; RSph < nphRS; RSph++) {
            PLph = RSD_List_[isk]->kinOrder[RSph];
            int istart = m_PhaseSpeciesStartIndex[PLph];
            int nsp = m_PhaseSpeciesStartIndex[PLph + 1] - istart;
            for (int k = 0; k < nsp; k++) {
                double net = rsSpeciesProductionRates[kIndexKin];
                phaseMoleFlux[PLph] += net;
                kIndexKin++;
            }
        }
    }
}
//================================================================================================

void Electrode::getIntegratedPhaseMoleTransfer(doublereal* const phaseMolesTransfered)
{

    if (!pendingIntegratedStep_) {
        throw CanteraError(" Electrode::getIntegratedPhaseMoleTransfer", "no pending integration step");
    }
    double sum = 0.0;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        phaseMolesTransfered[iph] = 0.0;
        ThermoPhase& tp = thermo(iph);
        string pname = tp.id();
        int istart = m_PhaseSpeciesStartIndex[iph];
        int nsp = tp.nSpecies();
        for (int ik = 0; ik < nsp; ik++) {
            int k = istart + ik;
            phaseMolesTransfered[iph] += spMoleIntegratedSourceTerm_[k];
        }
        sum += fabs(phaseMolesTransfered[iph]);
    }
    /*
     *  TO DO : We can put in a check here on phaseMoles_Final vs phaseMoles_init_
     */
    /* THIS DOESNT WORK FOR ELECTRON PHASE and perhaps toher phases */
    /*
     double deltaPhaseMolesIntegrated = 0.0;

     sum += 1.0E-14;
     for (int iph = 0; iph < m_NumTotPhases; iph++) {
     deltaPhaseMolesIntegrated = phaseMoles_final_[iph] - phaseMoles_init_init_[iph];
     if (fabs(deltaPhaseMolesIntegrated - phaseMolesTransfered[iph])/ sum > 1.0E-7) {
     throw CanteraError(" Electrode::getIntegratedPhaseMoleTransfer",
     " mismatched phase deltas");
     }
     }
     */
}
//=================================================================================================
// Returns the current potential drop across the electrode
/*
 *   This is equal to phiMetal - phiSoln
 */
double Electrode::potentialDrop() const
{
    // check
    double phiM = phaseVoltages_[metalPhase_];
    double phiS = phaseVoltages_[solnPhase_];
    double dd = phiM - phiS;
    if (fabs(dd - deltaVoltage_) > 1.0E-7) {
        throw CanteraError(" Electrode::potentialDrop()", "inconsistency");
    }
    return deltaVoltage_;
}

//=================================================================================================
// Returns the standard state OCV for the selected ReactingSurfaceDomain.
/*
 *  When there is more than a single reaction,
 *  pick open circuit potential for reaction that is
 *  closest to equilibrium given the cell voltage since this one
 *  is the one for which open circuit is most relevant.
 *  Note that it will be possible for the standard state OCV
 *  to be computed for a different reaction relative to the
 *  method openCircuitVoltage(isk) that computes OCV from
 *  the current concentrations because it looks for the reaction
 *  closest to equilibrium given the current cell voltage.
 */
double Electrode::openCircuitVoltageSSRxn(int isk, int iReaction) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        return 0.0;
    }
    int rxnIndex = 0;
    if (iReaction >= 0) {
        rxnIndex = iReaction;
    }
    rsd->getDeltaSSGibbs(DATA_PTR(deltaG_));
    double length = rsd->rmcVector.size();
    double PhiRxn = 0.0;

    int metalPhaseRS = rsd->PLtoKinPhaseIndex_[metalPhase_];
    if (metalPhaseRS >= 0) {
        RxnMolChange* rmc = rsd->rmcVector[rxnIndex];
        /*
         *  Compute open circuit potential from first reaction
         */
        double nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
        if (nStoichElectrons != 0.0) {
            PhiRxn = deltaG_[rxnIndex] / Faraday / nStoichElectrons;
        }
        if (iReaction >= 0) {
            return PhiRxn;
        }

        /*
         *  Compute open circuit potential from other reactions.
         *  Pick open circuit potential for reaction that is
         *  closest to equilibrium given cell voltage since this one
         *  is the one for which open circuit is most relevant.
         */
        if (length > 1) {
            for (int i = 1; i < length; i++) {
                rmc = rsd->rmcVector[i];
                nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
                if (nStoichElectrons != 0.0) {
                    if (fabs(deltaVoltage_ - PhiRxn) > fabs(deltaVoltage_ - deltaG_[i] / Faraday / nStoichElectrons)) {
                        rxnIndex = i;
                        PhiRxn = deltaG_[i] / Faraday / nStoichElectrons;
                    }
                }
            }
        }
        return PhiRxn;
    }
    return 0.0;
}
//=================================================================================================
// Returns the equilibrium OCV for the selected ReactingSurfaceDomain and current conditions.
/*
 *  When there is more than a single reaction,
 *  pick open circuit potential for reaction that is
 *  closest to equilibrium given the cell voltage since this one
 *  is the one for which open circuit is most relevant.
 */
double Electrode::openCircuitVoltageRxn(int isk, int iReaction) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        return 0.0;
    }
    int rxnIndex = 0;
    if (iReaction >= 0) {
        rxnIndex = iReaction;
    }

    rsd->getDeltaGibbs(DATA_PTR(deltaG_));

    int nR = rsd->nReactions();
    double PhiRxn = 0.0; // store open circuit voltage here

    int metalPhaseRS = rsd->PLtoKinPhaseIndex_[metalPhase_];
    if (metalPhaseRS >= 0) {

        RxnMolChange* rmc = rsd->rmcVector[rxnIndex];
        /*
         *  Compute open circuit potential from first reaction
         */
        double nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
        if (nStoichElectrons != 0.0) {
            PhiRxn = deltaG_[rxnIndex] / Faraday / nStoichElectrons;
        }
        if (iReaction >= 0) {
            return PhiRxn;
        }
        /*
         *  Compute open circuit potential from other reactions.
         *  Pick open circuit potential for reaction that is
         *  closest to equilibrium given cell voltage since this one
         *  is the one for which open circuit is most relevant.
         */
        if (nR > 1) {
            for (int i = 1; i < nR; i++) {
                rmc = rsd->rmcVector[i];
                nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
                if (nStoichElectrons != 0.0) {
                    if (fabs(deltaVoltage_ - PhiRxn) > fabs(deltaVoltage_ - deltaG_[i] / Faraday / nStoichElectrons)) {
                        rxnIndex = i;
                        PhiRxn = deltaG_[i] / Faraday / nStoichElectrons;
                    }
                }
            }
        }
        return PhiRxn;
    }
    return 0.0;
}

//======================================================================================================
// A calculation of the open circuit potential of a surface
/*!
 *  This routine uses a root finder to find the voltage at which there
 *  is zero net electron production.  It leaves the object unchanged. However, it
 *  does change the voltage of the phases during the calculation, so this is a nonconst function.
 */
double Electrode::openCircuitVoltage(int isk)
{
    int printDebug = 30;
    static int oIts = 0;
    double nStoichElectrons = 0.0;
    double deltaT = 0.005;
    int rxnIndex = 0;
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        return 0.0;
    }

  startOver: ;
    /*
     *  Get the current value of deltaG_[] for reactions defined on the current reacting surface
     */
    rsd->getDeltaGibbs(DATA_PTR(deltaG_));
    if (printDebug < 0) {
        printf("we are here, oIts = %d\n", oIts);
    }
    RxnMolChange* rmc = 0;
    int nR = rsd->nReactions();
    if (nR == 0) {
        return 0.0;
    }
    int nP = rsd->nPhases();
    double phiRxn = 0.0; // store open circuit voltage here

    vector<int> phaseExistsInit(nP, 1);
    vector<int> phaseStabInit(nP, 1);

    for (int iph = 0; iph < nP; iph++) {
        phaseExistsInit[iph] = rsd->phaseStability(iph);
        phaseStabInit[iph] = rsd->phaseStability(iph);
        if (printLvl_ > printDebug) {
            ThermoPhase& tp = rsd->thermo(iph);
            string pname = tp.name();
            printf(" %20s %d  %d \n", pname.c_str(), phaseExistsInit[iph], phaseStabInit[iph]);
        }
    }

    int metalPhaseRS = rsd->PLtoKinPhaseIndex_[metalPhase_];
    if (metalPhaseRS < 0) {
        return 0.0;
    }
    double phiMin = +100.;
    double phiMax = -100.;
    int nT = 0;
    vector<int> reactionTurnedOn(nR, 1);
    bool thereIsOneTurnedOn = false;
    for (int i = 0; i < nR; i++) {
        int rxnIndex = i;
        rmc = rsd->rmcVector[rxnIndex];
        /*
         *  Compute open circuit potential for the current reaction reaction
         */
        nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
        if (nStoichElectrons != 0.0) {
            nT++;
            phiRxn = deltaG_[rxnIndex] / Faraday / nStoichElectrons;
            phiMax = std::max(phiMax, phiRxn);
            phiMin = std::min(phiMin, phiRxn);

            for (int jph = 0; jph < nP; jph++) {
                if (rmc->m_phaseReactantMoles[jph] > 0.0 || rmc->m_phaseProductMoles[jph] > 0.0) {
                    if (phaseExistsInit[jph] == 0 || phaseStabInit[jph] == 0) {
                        reactionTurnedOn[i] = 0;
                    }
                }
            }

        } else {
            phiRxn = 0.0;
            reactionTurnedOn[i] = 0;
        }

        if (printLvl_ > printDebug) {
            string rstring = rsd->reactionString(i);
            printf("%d %100.100s %g  \n", i, rstring.c_str(), phiRxn);
        }
        if (reactionTurnedOn[i]) {
            thereIsOneTurnedOn = true;
        }
    }

    if (printLvl_ > printDebug) {
        printf("    deltaVoltage_ = %g\n", deltaVoltage_);
    }
    /*
     *  If there is only one participating electron reaction then the
     *  open circuit potential is given by the Nernst equation for that reaction.
     *  We just calculated it above. Therefore, we can do a quick return.
     */
    if (nT <= 1) {
        return phiRxn;
    }
    double phiMetalInit = phaseVoltages_[metalPhase_];

    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];

    vector<int> phaseExists(phaseExistsInit);
    vector<int> phaseStab(phaseStabInit);

    double phiRxnBest = 1.0E300;
    for (int i = 0; i < nR; i++) {
        bool electReact = false;
        bool electProd = false;
        /*
         *  get the extra information for the reaction
         */
        rmc = rsd->rmcVector[i];
        /*
         *  First figure out if electron is a reactant or a product
         */
        nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
        if (nStoichElectrons > 0) {
            electReact = true;
        }
        if (nStoichElectrons < 0) {
            electProd = true;
        }
        if (nStoichElectrons != 0.0) {

            phiRxn = deltaG_[i] / Faraday / nStoichElectrons;
            if (fabs(deltaVoltage_ - phiRxnBest) > fabs(deltaVoltage_ - phiRxn)) {
                phiRxnBest = phiRxn;
            }
            if (!thereIsOneTurnedOn) {
                if (electProd) {
                    if (phiRxn > deltaVoltage_) {
                        for (int jph = 0; jph < nP; jph++) {
                            if (rmc->m_phaseReactantMoles[jph] > 0.0 || rmc->m_phaseProductMoles[jph] > 0.0) {
                                phaseExists[jph] = 1;
                                phaseStab[jph] = 1;
                            }
                        }
                    }
                }
                if (electReact) {
                    if (phiRxn <= deltaVoltage_) {
                        for (int jph = 0; jph < nP; jph++) {
                            if (rmc->m_phaseReactantMoles[jph] > 0.0 || rmc->m_phaseProductMoles[jph] > 0.0) {
                                phaseExists[jph] = 1;
                                phaseStab[jph] = 1;
                            }
                        }
                    }
                }
            }
        }
    }
    if (!thereIsOneTurnedOn) {
        rmc = rsd->rmcVector[rxnIndex];
        for (int jph = 0; jph < nP; jph++) {
            if (rmc->m_phaseReactantMoles[jph] > 0.0 || rmc->m_phaseProductMoles[jph] > 0.0) {
                phaseExists[jph] = 1;
                phaseStab[jph] = 1;
            }
        }
    }

    for (int iph = 0; iph < nP; iph++) {
        rsd->setPhaseExistence(iph, phaseExists[iph]);
        rsd->setPhaseStability(iph, phaseStab[iph]);
        if (printLvl_ > printDebug) {
            ThermoPhase& tp = rsd->thermo(iph);
            string pname = tp.name();
            printf(" %20s %d  %d \n", pname.c_str(), phaseExists[iph], phaseStab[iph]);
        }
    }

    /*
     *  bump the values
     */
    phiMin = phiMin - fabs(0.001 + fabs((phiMax - phiMin) * 0.01));
    phiMax = phiMax + fabs(0.001 + fabs((phiMax - phiMin) * 0.01));

    RSD_ElectronProduction ep(rsd, metalPhase_, deltaT);

    if (printLvl_ > printDebug) {
        ep.printLvl_ = 15;
    }

    RootFind rf(&ep);

    phiRxn = openCircuitVoltageSSRxn(isk);
    phiRxn = phiRxnBest;
    double phiMetalRxn = phiRxn + phaseVoltages_[solnPhase_];
    if (printLvl_ > printDebug) {
        rf.setPrintLvl(15);
    }
    rf.setTol(1.0E-5, 1.0E-10);
    rf.setFuncIsGenerallyIncreasing(true);
    rf.setDeltaX(0.01);

    double funcTargetValue = 0.0;
    double phiMetalMin = phiMin + phaseVoltages_[solnPhase_];
    double phiMetalMax = phiMax + phaseVoltages_[solnPhase_];

    int retn = rf.solve(phiMetalMin, phiMetalMax, 100, funcTargetValue, &phiMetalRxn);

    if (retn != 0) {
        printf(" Electrode::openCircuitVoltage: rootSolver returned an error = %d\n", retn);
        printf("            phiMetalMin = %g, phiMetalMax = %g, phiMetalRxn = %g, funcTargetValue = %g\n", phiMetalMin,
               phiMetalMax, phiMetalRxn, funcTargetValue);
        printf("            Redoing calculating with printing turned on, then terminating, oIts = %d\n", oIts);
        if (printDebug >= 0) {
            printDebug = -1;
            goto startOver;
        }
        exit(-1);
    }
    phaseVoltages_[metalPhase_] = phiMetalInit;
    phiRxn = phiMetalRxn - phaseVoltages_[solnPhase_];

    for (int iph = 0; iph < rmc->m_nPhases; iph++) {
        rsd->setPhaseExistence(iph, phaseExistsInit[iph]);
        rsd->setPhaseStability(iph, phaseStabInit[iph]);
        if (printLvl_ > printDebug) {
            ThermoPhase& tp = rsd->thermo(iph);
            string pname = tp.name();
            printf(" %20s %d  %d \n", pname.c_str(), phaseExistsInit[iph], phaseStabInit[iph]);
        }
    }

    updateState();

    oIts++;

    return phiRxn;
}

//=================================================================================================
// Returns the equilibrium OCV for the selected ReactingSurfaceDomain and current conditions.
/*
 *  When there is more than a single reaction,
 *  pick open circuit potential for reaction that is
 *  closest to equilibrium given the cell voltage since this one
 *  is the one for which open circuit is most relevant.
 */
void Electrode::getOpenCircuitVoltages(int isk, double* ocv) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        ocv[0] = 0.0;
        return;
    }
    /*
     * reaction Gibbs free energy divided by number of electrons and Faraday constant
     */
    rsd->getDeltaGibbs(DATA_PTR(deltaG_));

    //loop over reactions
    double length = rsd->rmcVector.size();

    int metalPhaseRS = rsd->PLtoKinPhaseIndex_[metalPhase_];
    if (metalPhaseRS >= 0) {
        RxnMolChange* rmc;
        double nStoichElectrons;
        for (int i = 0; i < length; i++) {
            rmc = rsd->rmcVector[i];
            nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
            if (nStoichElectrons != 0.0) {
                ocv[i] = deltaG_[i] / Faraday / nStoichElectrons;
            } else {
                ocv[i] = 0.0;
            }
        }

    } else {
        // no metalPhase_
        ocv[0] = 0.0;
    }
    return;
}
//====================================================================================================================
// Returns the overpotential for the current conditions
/*
 *  The overpotential is the current voltage minus the open circuit voltage.
 */
double Electrode::overpotential(int isk)
{
    double Erxn = openCircuitVoltage(isk);
    return (deltaVoltage_ - Erxn);
}
//====================================================================================================================
// Returns the overpotential for the current conditions based on a particular reaction
/*
 *  The overpotential is the current voltage minus the open circuit voltage.
 *  This routine calculates the open circuit voltage for the surface by calling a
 *  openCircuitVoltageRxn() routine.
 *
 *   @param isk   reacting surface number
 *   @param irxn  Reaction number
 */
double Electrode::overpotentialRxn(int isk, int irxn)
{
    double Erxn = openCircuitVoltageRxn(isk, irxn);
    return (deltaVoltage_ - Erxn);
}

//====================================================================================================================
//=================================================================================================
// Returns the vector of exchange current densities for given surface
/*
 * For an exchange current to be relevant, we need the reaction to involve an electron.
 * However, the InterfaceKinetics version of this method was implemented to let
 * i0 = kf / Kc^0.5 for non-electrochemical reactions.
 */
double Electrode::getExchangeCurrentDensity(int isk, int irxn) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        return 0.0;
    } else {
        double nstoich;
        double ocv;
        double io;
        double nu;
        rsd->getExchangeCurrentFormulation(irxn, &nstoich, &ocv, &io, &nu);
        return io;
    }
}
//====================================================================================================================
int Electrode::kKinSpecElectron(int isurf) const
{
    return kKinSpecElectron_sph_[isurf];
}

//====================================================================================================================
int Electrode::metalPhaseIndex() const
{
    return metalPhase_;
}
//====================================================================================================================
int Electrode::solnPhaseIndex() const
{
    return solnPhase_;
}
//====================================================================================================================
int Electrode::numSolnPhaseSpecies() const
{
    return getPhase(phaseName(solnPhase_).c_str())->nSpecies();
}

//====================================================================================================================
// Return the number of extra print tables
int Electrode::getNumPrintTables() const
{
    return 0;
}
//====================================================================================================================
//! Get the values that are printed in tables for the 1D code.
void Electrode::getPrintTable(int itable, std::vector<std::string>& colNames, std::vector<double>& colValues) const
{

}
//====================================================================================================================
// Set the current depth of discharge
/*
 * This is roughly equal to the total number of electrons that has been discharged from a fully charged state
 * divided by the total moles
 *
 * @return  returns the depth of discharge in amps sec
 * @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 */
void Electrode::setRelativeCapacityDischargedPerMole(double relDischargedPerMole, int platNum)
{
    if (platNum != -1) {
        throw CanteraError(" Electrode::setRelativeCapacityDischargedPerMole()", "platNum not -1");
    }
    if (pendingIntegratedStep_) {
        throw CanteraError(" Electrode::setRelativeCapacityDischargedPerMole()", "pending step");
    }
    throw CanteraError(" Electrode::setRelativeCapacityDischargedPerMole()", "not implemented");
}

//====================================================================================================================
RxnMolChange* Electrode::rxnMolChangesEGR(int iegr)
{
    return m_rmcEGR[iegr];
}
//====================================================================================================================
void Electrode::addExtraGlobalRxn(EGRInput* egr_ptr)
{
    InterfaceKinetics* iKA = RSD_List_[0];
    int nReactionsA = iKA->nReactions();

    Cantera::ExtraGlobalRxn* egr = new ExtraGlobalRxn(iKA);
    double* RxnVector = new double[nReactionsA];
    for (int i = 0; i < nReactionsA; i++) {
        RxnVector[i] = 0.0;
    }

    for (int iErxn = 0; iErxn < egr_ptr->m_numElemReactions; iErxn++) {
        struct ERSSpec* ers_ptr = egr_ptr->m_ERSList[iErxn];
        RxnVector[ers_ptr->m_reactionIndex] = ers_ptr->m_reactionMultiplier;
    }
    egr->setupElemRxnVector(RxnVector, egr_ptr->m_SS_KinSpeciesKindex);

    RxnMolChange* rmcEGR = new RxnMolChange(iKA, egr);

    m_egr.push_back(egr);
    m_rmcEGR.push_back(rmcEGR);

    delete[] RxnVector;
}
//======================================================================================================================
int Electrode::processExtraGlobalRxnPathways()
{
    for (int i = 0; i < numExtraGlobalRxns; i++) {
        addExtraGlobalRxn(m_EGRList[i]);
    }
    return numExtraGlobalRxns;
}
//====================================================================================================================

Cantera::ExtraGlobalRxn* Electrode::extraGlobalRxnPathway(int iegr)
{
    return m_egr[iegr];
}
//===================================================================================================================
int Electrode::numExtraGlobalRxnPathways() const
{
    return numExtraGlobalRxns;
}
//===================================================================================================================
Electrode::phasePop_Resid::phasePop_Resid(Electrode* ee, int iphaseTarget, double* const Xmf_stable,
        double deltaTsubcycle) :
                ResidEval(),
                ee_(ee),
                iphaseTarget_(iphaseTarget),
                Xmf_stable_(Xmf_stable),
                deltaTsubcycle_(deltaTsubcycle)
{
}
int Electrode::phasePop_Resid::evalSS(const doublereal t, const doublereal* const y, doublereal* const r)
{
    int retn = ee_->phasePopResid(iphaseTarget_, y, deltaTsubcycle_, r);
    return retn;
}
int Electrode::phasePop_Resid::getInitialConditions(const doublereal t0, doublereal* const y, doublereal* const ydot)
{
    int ne = nEquations();
    for (int k = 0; k < ne; k++) {
        y[k] = Xmf_stable_[k];
    }
    return 1;
}
int Electrode::phasePop_Resid::nEquations() const
{
    ThermoPhase* tp = &(ee_->thermo(iphaseTarget_));
    int nsp = tp->nSpecies();
    return nsp;
}
//===================================================================================================================
int Electrode::phasePopResid(int iphaseTarget, const double* const Xf_phase, double deltaTsubcycle, double* const resid)
{
    int k;
    vector<doublereal> spMoles_tmp(m_NumTotSpecies, 0.0);

    ThermoPhase* tptarget = PhaseList_[iphaseTarget];
    int nspPhase = tptarget->nSpecies();
    if (nspPhase == 1) {
        resid[0] = Xf_phase[0] - 1.0;
        return 0;
    }

    int kstartTarget = m_PhaseSpeciesStartIndex[iphaseTarget];
    /*
     *   Set the internal objects to the correct conditions
     *    -> This will be the final conditions.
     *    -> Because we do this within an iterator loop, we are essentially trying to emulate
     *       an implicit algorithm
     */
    updateState();
    tptarget->setState_PX(pressure_, (double*) Xf_phase);

    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    for (int isk = 0; isk < numSurfaces_; isk++) {
        /*
         *  Loop over phases, figuring out which phases have zero moles.
         *  Volume phases exist if the initial or final mole numbers are greater than zero
         *  Surface phases exist if the initial or final surface areas are greater than zero.
         */
        if (ActiveKineticsSurf_[isk]) {
            ReactingSurDomain* rsd = RSD_List_[isk];
            int nph = rsd->nPhases();
            for (int jph = 0; jph < nph; jph++) {
                int iph = rsd->kinOrder[jph];
                if (iph == metalPhase_) {
                    continue;
                }
                double mm = phaseMoles_init_[iph];
                double mmf = phaseMoles_final_[iph];
                ThermoPhase& tp = thermo(iph);
                int nsp = tp.nSpecies();
                if (iph >= NumVolPhases_) {
                    // we are in a surface phase
                    int isur = iph - NumVolPhases_;
                    double sa_init = surfaceAreaRS_init_[isur];
                    double sa_final = surfaceAreaRS_final_[isur];
                    if (sa_init > 0.0 || sa_final > 0.0) {
                        rsd->setPhaseExistence(jph, true);
                    } else {
                        rsd->setPhaseExistence(jph, false);
                    }
                } else {
                    if (mm <= 0.0 && mmf <= 0.0) {
                        rsd->setPhaseExistence(jph, false);
                        if (nsp == 1 || nsp == 2) {
                            rsd->setPhaseStability(jph, true);
                        }
                    } else {
                        rsd->setPhaseExistence(jph, true);
                    }
                }
                if (iph == iphaseTarget) {
                    rsd->setPhaseExistence(jph, true);
                }
            }
        }
    }
    /*
     *  This routine basically translates between species lists for the reacting surface
     *  domain and the Electrode.
     *  Later, when we have more than one reacting surface domain in the electrode object,
     *  this will do a lot more
     */

    for (int isk = 0; isk < numSurfaces_; isk++) {
        // Loop over phases, figuring out which phases have zero moles.

        if (ActiveKineticsSurf_[isk]) {

            /*
             *  For each Reacting surface
             *
             *  Get the species production rates for the reacting surface
             */
            //    m_rSurDomain->getNetProductionRates(&RSSpeciesProductionRates_[0]);
            const vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();

            double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
            /*
             *  loop over the phases in the reacting surface
             *  Get the net production vector
             */
            std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
            int nphRS = RSD_List_[isk]->nPhases();
            int jph, kph;
            int kIndexKin = 0;
            for (kph = 0; kph < nphRS; kph++) {
                jph = RSD_List_[isk]->kinOrder[kph];
                int istart = m_PhaseSpeciesStartIndex[jph];
                int nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
                for (int k = 0; k < nsp; k++) {
                    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                    if (rsSpeciesProductionRates[kIndexKin] > 0.0) {
                        if (phaseMoles_init_[jph] <= 0.0) {
                            if (nsp > 1) {
                                int bornMultiSpecies = jph;
                                if (iphaseTarget != bornMultiSpecies) {
                                    throw CanteraError("", "two multispecies phases");
                                }
                            }
                        }
                    }
                    kIndexKin++;
                }
            }
        }
    }
    /*
     *  Find the initial surface area to use
     */
    double sa_init = surfaceAreaRS_init_[0];

    /*
     *  Find the final surface area to use
     */
    double sa_final = sa_init;

    /*
     *  Calculate the change in the moles of all of the species
     */

    for (k = 0; k < m_NumTotSpecies; k++) {
        spMoles_tmp[k] = spMoles_init_[k];
        for (int isk = 0; isk < numSurfaces_; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                sa_init = surfaceAreaRS_init_[isk];
                sa_final = surfaceAreaRS_final_[isk];
                double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
                spMoles_tmp[k] += 0.5 * deltaTsubcycle * (sa_init + sa_final) * spNetProdPerArea[k];
            }
        }
    }

    /*
     * Calculate the total moles
     */
    double ptotal = 0.0;
    for (int kp = 0; kp < nspPhase; kp++) {
        k = kp + kstartTarget;
        ptotal += spMoles_tmp[k];
    }
    if (ptotal <= 0.0) {
        return -1;
    }
    int kmax = 0;
    double xmax = Xf_phase[0];
    double xsum = 0.0;
    for (int kp = 0; kp < nspPhase; kp++) {
        k = kp + kstartTarget;
        resid[kp] = Xf_phase[kp] - spMoles_tmp[k] / ptotal;
        if (Xf_phase[kp] > xmax) {
            kmax = kp;
            xmax = Xf_phase[kp];
        }
        xsum += Xf_phase[kp];
    }
    resid[kmax] = 1.0 - xsum;

    return 0;
}
//===================================================================================================================
int Electrode::phasePop(int iphaseTarget, double* const Xmf_stable, double deltaTsubcycle)
{
    int k;
    vector<doublereal> Xf_phase(m_NumTotSpecies, 0.0);
    int retn = 1;
    ThermoPhase* tptarget = PhaseList_[iphaseTarget];
    int nspPhase = tptarget->nSpecies();
    if (nspPhase == 1) {
        Xmf_stable[0] = 1.0;
        return 0;
    }
    for (k = 0; k < nspPhase; k++) {
        Xf_phase[k] = Xmf_stable[k];
    }
    tptarget->setMoleFractions(DATA_PTR(Xf_phase));

    int kstartTarget = m_PhaseSpeciesStartIndex[iphaseTarget];

    /*
     * Check starting conditions
     */
    for (int kp = 0; kp < nspPhase; kp++) {
        k = kp + kstartTarget;
        if (spMoles_init_[k] != 0.0) {
            throw CanteraError(" Electrode::phasePop", "spMoles_init_[k] != 0.0");
        }
        if (spMoles_final_[k] != 0.0) {
            throw CanteraError(" Electrode::phasePop", "spMoles_final_[k] != 0.0");
        }
    }

    /*
     * Load with an initial guess
     */
    Electrode::phasePop_Resid* pSolve_Res = new phasePop_Resid(this, iphaseTarget, Xmf_stable, deltaTsubcycle);

    Cantera::solveProb* pSolve = new solveProb(pSolve_Res);
    pSolve->m_ioflag = 10;
    pSolve->setAtolConst(1.0E-11);
    retn = pSolve->solve(SOLVEPROB_RESIDUAL, 1.0E-6, 1.0E-3);
    if (retn == 0) {
        pSolve->reportState(Xmf_stable);
        vector_fp resid(nspPhase);
        phasePopResid(iphaseTarget, Xmf_stable, deltaTsubcycle, DATA_PTR(resid));
    }

    delete pSolve;
    delete pSolve_Res;
    return retn;
}

//====================================================================================================================

double Electrode::reportStateVariableIntegrationError(int& numSV, double* const errorVector) const
{
    numSV = 0;
    throw CanteraError("Electrode::reportStateVariableIntegrationError()", "Base Class Called");
    return 0.0;
}
//====================================================================================================================
// Return the number of equations in the Nonlinear equation system used to solve the system
// at every time step
/*
 *  This is also equal to the number of state variables in the problem
 */
int Electrode::nEquations() const
{
    throw CanteraError("Electrode::nEquations()", "Base Class Called");
}
//====================================================================================================================
Electrode::integrate_ResidJacEval::integrate_ResidJacEval(Electrode* ee) :
                ResidJacEval(),
                ee_(ee)
{
}
//  -----------------------------------------------------------------------------------------------------------------
int Electrode::integrate_ResidJacEval::evalResidNJ(doublereal t, const doublereal deltaT, const doublereal* const y,
        const doublereal* const ydot, doublereal* const resid, const ResidEval_Type_Enum evalType, int id_x,
        doublereal delta_x)
{
    int retn = ee_->integrateResid(t, deltaT, y, ydot, resid, evalType, id_x, delta_x);
    return retn;
}
//  -----------------------------------------------------------------------------------------------------------------
int Electrode::integrate_ResidJacEval::getInitialConditions(const doublereal t0, doublereal* const y,
        doublereal* const ydot)
{
    int ne = nEquations();
    for (int k = 0; k < ne; k++) {
        y[k] = 0.0;
    }
    return 1;
}
//  -----------------------------------------------------------------------------------------------------------------
int Electrode::integrate_ResidJacEval::nEquations() const
{
    return neq_;
}
//===================================================================================================================
//  Residual calculation for the solution of the Nonlinear integration problem
/*
 *
 * @param t             Time                    (input)
 * @param delta_t       The current value of the time step (input)
 * @param y             Solution vector (input, do not modify)
 * @param ydot          Rate of change of solution vector. (input, do not modify)
 * @param resid         Value of the residual that is computed (output)
 * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
 * @param id_x          Index of the variable that is being numerically differenced to find
 *                      the jacobian (defaults to -1, which indicates that no variable is being
 *                      differenced or that the residual doesn't take this issue into account)
 * @param delta_x       Value of the delta used in the numerical differencing
 */
int Electrode::integrateResid(const doublereal tfinal, const doublereal deltaTsubcycle, const doublereal* const y,
        const doublereal* const ydot, doublereal* const resid, const ResidEval_Type_Enum evalType, const int id_x,
        const doublereal delta_x)
{

    throw CanteraError(" Electrode::integrateResid()", "Base class called");
    return 0;
}
//====================================================================================================================
// The internal state of the electrode must be kept for the initial and
// final times of an integration step.
/*
 *  This function advances the initial state to the final state that was calculated
 *  in the last integration step. If the initial time is input, then the code doesn't advance
 *  or change anything.
 *
 * @param Tinitial   This is the New initial time. This time is compared against the "old"
 *                   final time, to see if there is any problem.
 */
void Electrode::resetStartingCondition(double Tinitial, bool doTestsAlways)
{
    int i;

    /*
     *  If we haven't done a time step then we don't have anything to update.
     */
    if (pendingIntegratedStep_ != 1) {
#ifdef DEBUG_ELECTRODE
        // printf(" Electrode::resetStartingCondition WARNING: resetStartingCondition called with no pending integration step\n");
#endif
        return;
    }

    /*
     * If this routine is called with Tinitial = t_init_init_, then we should return without doing anything
     * We have already advanced the time step to the new time.
     */
    double tbase = MAX(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-9 * tbase) && !doTestsAlways) {
        return;
    }

    /*
     *  The final_final time must be equal to the new Tinitial time
     */ 
    tbase = MAX(Tinitial, tbase);
    tbase = MAX(tbase, t_final_final_);
    if (fabs(Tinitial - t_final_final_) > (1.0E-9 * tbase)) {
        throw CanteraError("Electrode::resetStartingCondition()",
                           "Tinitial " + fp2str(Tinitial) + " is not compatible with t_final_final_ " + fp2str(t_final_final_));
    }

    /*
     *  Here is where we store the electrons discharged
     */
    if (pendingIntegratedStep_ == 1) {
        electronKmolDischargedToDate_ += spMoleIntegratedSourceTerm_[kElectron_];
    }

    t_init_init_ = Tinitial;

    // reset surface quantities
    for (i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_init_[i] = surfaceAreaRS_final_[i];
        surfaceAreaRS_init_init_[i] = surfaceAreaRS_final_[i];
    }

    // Reset total species quantities
    for (int k = 0; k < m_NumTotSpecies; k++) {
        spMoles_init_[k] = spMoles_final_[k];
        spMoles_init_init_[k] = spMoles_final_[k];
    }

    std::fill(spMoleIntegratedSourceTerm_.begin(), spMoleIntegratedSourceTerm_.end(), 0.);
    std::fill(spMoleIntegratedSourceTermLast_.begin(), spMoleIntegratedSourceTermLast_.end(), 0.);

    // Reset the total phase moles quantities
    for (i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_init_[i] = phaseMoles_final_[i];
        phaseMoles_init_init_[i] = phaseMoles_final_[i];
    }

    // Reset the particle size
    Radius_exterior_init_init_ = Radius_exterior_final_final_;

    /*
     *  Change the initial subcycle time delta here. Note, we should not change it during the integration steps
     *  because we want jacobian calculations to mainly use the same time step history, so that the answers are
     *  comparible irrespective of the time step truncation error.
     */;
    if (deltaTsubcycle_init_next_ < 1.0E299) {
        deltaTsubcycle_init_init_ = deltaTsubcycle_init_next_;
    }
    deltaTsubcycle_init_next_ = 1.0E300;

    pendingIntegratedStep_ = 0;
}
//====================================================================================================================
// Set the internal initial intermediate and initial global state from the internal final state
/*
 *  (non-virtual function)  -> function should onionize in-first.
 *
 *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 * @param setInitInit   Boolean indicating whether you should set the init_init state as well
 */
void Electrode::setInitStateFromFinal_Oin(bool setInitInit)
{
    int i;
    // reset surface quantities
    for (i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_init_[i] = surfaceAreaRS_final_[i];
    }
    if (setInitInit) {
        for (i = 0; i < numSurfaces_; i++) {
            surfaceAreaRS_init_init_[i] = surfaceAreaRS_final_[i];
        }
    }

    // Reset total species quantities
    for (int k = 0; k < m_NumTotSpecies; k++) {
        spMoles_init_[k] = spMoles_final_[k];
        spMf_init_[k] = spMf_final_[k];
        spMoles_final_final_[k] = spMoles_final_[k];
        if (setInitInit) {
            spMf_init_init_[k] = spMf_final_[k];
            spMoles_init_init_[k] = spMoles_final_[k];
        }
    }

    // Reset the total phase moles quantities
    for (i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_init_[i] = phaseMoles_final_[i];
        if (setInitInit) {
            phaseMoles_init_init_[i] = phaseMoles_final_[i];
        }
    }

    // Reset the particle size
    if (setInitInit) {
        if (pendingIntegratedStep_) {
            throw CanteraError("Electrode::setInitStateFromFinal(true)",
                               "Function called to overwrite init_init during a pending step");
        }
        Radius_exterior_init_init_ = Radius_exterior_final_;
    }
    Radius_exterior_init_ = Radius_exterior_final_;
    Radius_exterior_final_final_ = Radius_exterior_final_;

    if (xmlStateData_final_) {
        SAFE_DELETE(xmlStateData_init_);
        xmlStateData_init_ = new XML_Node(*xmlStateData_final_);
        SAFE_DELETE(xmlStateData_final_final_);
        xmlStateData_final_final_ = new XML_Node(*xmlStateData_final_);
        if (setInitInit) {
            SAFE_DELETE(xmlStateData_init_init_);
            xmlStateData_init_init_ = new XML_Node(*xmlStateData_final_);
        }
    }

    if (xmlExternalData_final_) {
        SAFE_DELETE(xmlExternalData_init_);
        xmlStateData_init_ = new XML_Node(*xmlExternalData_final_);
        SAFE_DELETE(xmlExternalData_final_final_);
        xmlExternalData_final_final_ = new XML_Node(*xmlExternalData_final_);
        if (setInitInit) {
            SAFE_DELETE(xmlExternalData_init_init_);
            xmlExternalData_init_init_ = new XML_Node(*xmlExternalData_final_);
        }
    }

    tinit_ = tfinal_;
    //t_final_final_ = tfinal_;
    if (setInitInit) {
        t_init_init_ = tfinal_;
    }
}
//====================================================================================================================
// Set the internal initial intermediate and initial global state from the internal final state
/*
 *  (virtual function)
 *
 *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 * @param setInitInit   Boolean indicating whether you should set the init_init state as well
 */
void Electrode::setInitStateFromFinal(bool setInitInit)
{
    setInitStateFromFinal_Oin(setInitInit);
}
//====================================================================================================================
//! Set the internal initial intermediate and initial global state from the internal final_final state
/*!
 *  (virtual function)
 *
 *  Set the intial  and init_int state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 */
void Electrode::setInitInitStateFromFinalFinal()
{
    int i;
    // reset surface quantities
    for (i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_init_[i] = surfaceAreaRS_final_final_[i];
        surfaceAreaRS_init_init_[i] = surfaceAreaRS_final_final_[i];
    }

    // Reset total species quantities
    for (int k = 0; k < m_NumTotSpecies; k++) {
        spMoles_init_[k] = spMoles_final_final_[k];
        spMf_init_[k] = spMf_final_final_[k];
        spMf_init_init_[k] = spMf_final_final_[k];
        spMoles_init_init_[k] = spMoles_final_final_[k];
    }

    // Reset the total phase moles quantities
    for (i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_init_[i] = phaseMoles_final_final_[i];
        phaseMoles_init_init_[i] = phaseMoles_final_final_[i];
    }

    if (pendingIntegratedStep_) {
        throw CanteraError("Electrode::setInitInitStateFromFinalFinal(true)",
                           "Function called to overwrite init_init during a pending step");
    }

    Radius_exterior_init_init_ = Radius_exterior_final_final_;
    Radius_exterior_init_ = Radius_exterior_final_final_;

    if (xmlStateData_final_final_) {
        SAFE_DELETE(xmlStateData_init_);
        xmlStateData_init_ = new XML_Node(*xmlStateData_final_final_);
        SAFE_DELETE(xmlStateData_init_init_);
        xmlStateData_init_init_ = new XML_Node(*xmlStateData_final_final_);
    }

    if (xmlExternalData_final_final_) {
        SAFE_DELETE(xmlExternalData_init_);
        xmlStateData_init_ = new XML_Node(*xmlExternalData_final_final_);
        SAFE_DELETE(xmlExternalData_init_init_);
        xmlExternalData_init_init_ = new XML_Node(*xmlExternalData_final_final_);
    }

    tinit_ = t_final_final_;
    //t_final_final_ = tfinal_;

    t_init_init_ = t_final_final_;

}
//====================================================================================================================
// Revert the object's conditions to the initial conditions
/*!
 *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init.
 *  We get rid of the pending flag here as well.
 *
 * @param revertToInitInit revert to the t_init_init solution
 *                         Defaults to true.
 */
void Electrode::revertToInitialTime(bool revertToInitInit)
{
    if (revertToInitInit) {
       setInitStateFromInitInit(true);
       pendingIntegratedStep_ = 0;
    } else {
       setFinalStateFromInit();
    }
}
//====================================================================================================================
// Set the internal final intermediate from the internal init state
/*
 *  (nonvirtual function onionized from outer to inner)
 *
 *  Set the final state from the init state.
 *
 */
void Electrode::setFinalStateFromInit_Oin()
{
    int i;
    // reset surface quantities
    for (i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_final_[i] = surfaceAreaRS_init_[i];
    }
    // Reset total species quantities
    for (i = 0; i < m_NumTotSpecies; i++) {
        spMoles_final_[i] = spMoles_init_[i];
        spMf_final_[i] = spMf_init_[i];
    }
    // Reset the total phase moles quantities
    for (i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_final_[i] = phaseMoles_init_[i];
    }
    // Reset the particle size
    Radius_exterior_final_ = Radius_exterior_init_;

    if (xmlStateData_init_) {
        SAFE_DELETE(xmlStateData_final_);
        xmlStateData_final_ = new XML_Node(*xmlStateData_init_);
    }
    if (xmlExternalData_init_) {
        SAFE_DELETE(xmlExternalData_final_);
        xmlExternalData_final_ = new XML_Node(*xmlExternalData_init_);
    }
    tfinal_ = tinit_;
}
//====================================================================================================================
void Electrode::setFinalStateFromInit()
{
    setFinalStateFromInit_Oin();
//    updateState();
}
//====================================================================================================================
//    Set the internal initial intermediatefrom the internal initial global state
/*
 *  Set the intial state from the init init state. We also can set the final state from this
 *  routine as well.
 *
 *  The final_final is not touched.
 *
 * @param setFinal   Boolean indicating whether you should set the final as well
 */
void Electrode::setInitStateFromInitInit(bool setFinal)
{
    int i, k;

    // reset surface quantities
    for (i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_init_[i] = surfaceAreaRS_init_init_[i];
    }
    // Reset total species quantities
    if (setFinal) {
        for (k = 0; k < m_NumTotSpecies; k++) {
            spMoles_init_[k] = spMoles_init_init_[k];
            spMoles_final_[k] = spMoles_init_init_[k];
            spMf_init_[k] = spMf_init_init_[k];
            spMf_final_[k] = spMf_init_init_[k];
        }
    } else {
        for (k = 0; k < m_NumTotSpecies; k++) {
            spMoles_init_[k] = spMoles_init_init_[k];
            spMf_init_[k] = spMf_init_init_[k];
        }

    }
    // Reset the total phase moles quantities
    for (i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_init_[i] = phaseMoles_init_init_[i];
    }
    if (setFinal) {
        for (i = 0; i < m_NumTotPhases; i++) {
            phaseMoles_init_init_[i] = phaseMoles_init_init_[i];
        }
    }
    // Reset the particle size
    Radius_exterior_init_ = Radius_exterior_init_init_;
    if (setFinal) {
        Radius_exterior_final_ = Radius_exterior_init_init_;
    }
    // Copy XML data structures
    if (xmlStateData_init_init_) {
        SAFE_DELETE(xmlStateData_init_);
        xmlStateData_init_ = new XML_Node(*xmlStateData_init_init_);
        if (setFinal) {
            SAFE_DELETE(xmlStateData_final_);
            xmlStateData_final_ = new XML_Node(*xmlStateData_init_init_);
        }
    }
    if (xmlExternalData_init_init_) {
        SAFE_DELETE(xmlExternalData_init_);
        xmlStateData_init_ = new XML_Node(*xmlExternalData_init_init_);
        if (setFinal) {
            SAFE_DELETE(xmlExternalData_final_);
            xmlExternalData_final_ = new XML_Node(*xmlExternalData_init_init_);
        }
    }
    // reset the time
    tinit_ = t_init_init_;
    if (setFinal) {
        tfinal_ = t_init_init_;
    }
}
//====================================================================================================================
// Set the internal final global state from the internal final intermediate state
/*
 *  (virtual function from Electrode)
 *
 *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
 */
void Electrode::setFinalFinalStateFromFinal()
{
    Electrode::setFinalFinalStateFromFinal_Oin();
}
//====================================================================================================================
//   Set the internal final global state from the internal final intermediate state
/*
 *  (non-virtual function from Electrode)
 *
 *  Set the final_final state from the final state. This is commonly called at the end of successful base integration.
 *  This is an onionize in function
 */
void Electrode::setFinalFinalStateFromFinal_Oin()
{
    for (int i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_final_final_[i] = surfaceAreaRS_final_[i];
    }

    for (int i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_final_final_[i] = phaseMoles_final_[i];
    }

    for (int k = 0; k < m_NumTotSpecies; k++) {
        spMoles_final_final_[k] = spMoles_final_[k];
        spMf_final_final_[k] = spMf_final_[k];
    }

    Radius_exterior_final_final_ = Radius_exterior_final_;

    if (xmlStateData_final_) {
        SAFE_DELETE(xmlStateData_final_final_);
        xmlStateData_final_final_ = new XML_Node(*xmlStateData_final_);
    }

    if (xmlExternalData_final_) {
        SAFE_DELETE(xmlExternalData_final_final_);
        xmlExternalData_final_final_ = new XML_Node(*xmlExternalData_final_);
    }
    // We may not have reached the final time
    // 
    t_final_final_ = tfinal_;
}
//====================================================================================================================
//  Calculate the change in the state of the system when integrating from Tinitial to Tfinal
//  at constant current
/*
 *
 *  @param current    Current in amps. On return it returns the current actually obtained
 *  @param deltaT     DeltaT for the integration step.
 *  @param phiMax     Maximum value of the voltage. Defaults to  100.
 *  @param phiMin     Minimum value of the voltage. Defaults to -100.
 *  @param maxIntegrationSteps  defaults to 5000
 *
 *  @return Return the voltage used to obtain the current
 */
double Electrode::integrateConstantCurrent(doublereal& current, double& deltaT, double phiMax, double phiMin,
					   int maxIntegrationSteps)
{
    int status;
    double currentNeeded = current;
    Electrode_ECurr ec(this, deltaT);
    //printf("current needed = %g\n", currentNeeded);

    /*
     *  We take our cues about the best voltages from the electrode object itself
     */
    double phiM = phaseVoltages_[metalPhase_];
    double phimax = phiMax;
    double phimin = phiMin;

    if (phiMax == 100. && phiMin == -100.) {
        phimax = phiM + 2.0;
        phimin = phiM - 2.0;
    }
    if (phiMax < phiM) {
	throw CanteraError("Electrode::integrateConstantCurrent()",
			   "phiMax , " + fp2str(phiMax) + ", is less than starting phi, " + fp2str(phiM));
    }
    if (phiMin > phiM) {
	throw CanteraError("Electrode::integrateConstantCurrent()",
			   "phiMin , " + fp2str(phiMin) + ", is greater than starting phi, " + fp2str(phiM));
    }

    RootFind rf(&ec);
    rf.setTol(1.0E-5, 1.0E-10);
    rf.setFuncIsGenerallyIncreasing(true);
    rf.setDeltaX(0.01);
    double xbest = phiM;
    int oldP = printLvl_;
    double deltaT_curr = deltaT;
    // Turn down printing for lower levels
    printLvl_ = MAX(0, printLvl_ - 3);
    int printLvlLocal = oldP;
    if (printLvlLocal > 0) {
        rf.setPrintLvl(printLvlLocal);
        ec.printLvl_ = printLvlLocal;
    } else {
        rf.setPrintLvl(0);
        ec.printLvl_ = 0;
    }
    Electrode_Integrator* eei = dynamic_cast<Electrode_Integrator*>(this);
  
    int numSteps = integrate(deltaT, 1.0E-3,  T_FINAL_CONST_FIS, BASE_TIMEINTEGRATION_SIR);
    if (numSteps > maxIntegrationSteps) {
	if (eei) {
	    SubIntegrationHistory& sih = eei->timeHistory();
	    if (printLvlLocal > 2) {
		sih.print(5);
	    } 
	    int nn = 3 * maxIntegrationSteps / 4 + 1;
	    TimeStepHistory& tsh = sih.TimeStepList_[nn];
	    deltaT_curr = tsh.t_final_calc_;
	} else {
	    deltaT_curr = 0.5 * deltaT *  maxIntegrationSteps /  numSteps;
	}
    } else {
	if (printLvlLocal > 1) {
	    if (eei) {
		SubIntegrationHistory& sih = eei->timeHistory(); 
		sih.print(5);
	    }
	}
    }


    /*
     *  We put the rootfinder in a loop here. Sometimes the rootfinder will fail. The algorithm is to cut
     *  the time in half and retry the problem again.
     */
    for (int iTrial = 0; iTrial < 5; iTrial++) {
        ec.set_deltaT(deltaT_curr);
        double currentObtained = currentNeeded;
        xbest = phiM;
        status = rf.solve(phimin, phimax, 100, currentObtained, &xbest);

        /*
         *  Status = 0 means success. Return with a success indicator.
         */
        if (status == 0) {
            if (printLvl_ > 1) {
                printf("Electrode::integrateConstantCurrent(): Volts (%g amps) = %g at deltaT = %g\n", currentObtained,
                       xbest, deltaT_curr);
            }
	    numSteps = ec.nIntegrationSteps();
	    if (eei) {
		SubIntegrationHistory& sih = eei->timeHistory(); 
		if (printLvl_ > 2) {
		    sih.print(5);
		}
	    
		if (numSteps > maxIntegrationSteps) {
		    int nn = 3 * maxIntegrationSteps / 4 + 1;
		    TimeStepHistory& tsh = sih.TimeStepList_[nn];
		    deltaT_curr = tsh.t_final_calc_ - t_init_init_;
		    deltaT = deltaT_curr;
		    current = currentObtained;
		    continue;
		}
	    }
            deltaT = deltaT_curr;
            current = currentObtained;
            printLvl_ = oldP;
            return xbest;
        }
        /*
         *  Any other status indicator, we print out a warning message
         */
        if (printLvlLocal > 1) {
            printf("Electrode::integrateConstantCurrent(): bad status = %d Volts (%g amps) = %g at deltaT = %g\n",
                   status, currentObtained, xbest, deltaT_curr);
        }
        /*
         *  Deal with the returned voltage being at the top limit voltage
         */
        if (xbest >= phiMax) {
            /*
             *  If we have tried the rootfinder 3 times with decreasing time steps and we are at the max voltage,
             *  then give up early. We are probably at an end-of-electron condition. With this condition, we need to
             *  exit gracefully and handle the situation at a higher level.
             */
            if (iTrial >= 2) {
                if (printLvlLocal > 1) {
                    printf("                            volts at top = %g : returning with voltMax solution\n", phiMax);
                }
                current = currentObtained;
                deltaT = deltaT_curr;
                printLvl_ = oldP;
                return xbest;
            } else {
                if (printLvlLocal > 1) {
		    if (fabs(currentObtained ) < 1.0E-10) {
			printf("                            volts at top = %g probably because no electrons left, curr = %g\n", phiMax, currentObtained);
		    } else {
			printf("                            volts at top = %g\n", phiMax);
		    }
                }
		if (fabs(currentObtained ) < 1.0E-10) {
		    current = currentObtained;
		    deltaT = deltaT_curr;
		    printLvl_ = oldP;
		    return xbest;
		}
            }
        }
        /*
         *  Deal with the returned voltage being at the bottom limit voltage
         */
        if (xbest <= phiMin) {
            /*
             *  If we have tried the rootfinder 3 times with decreasing time steps and we are at the min voltage,
             *  then give up early. We are probably at an end-of-electron condition. With this condition, we need to
             *  exit gracefully and handle the situation at a higher level.
             */
            if (iTrial >= 2) {
                if (printLvlLocal > 1) {
                    printf("                            volts at bottom = %g: returning with voltMin solution\n",
                           phiMin);
                }
                current = currentObtained;
                deltaT = deltaT_curr;
                printLvl_ = oldP;
                return xbest;
            } else {
		if (printLvlLocal > 1) {
		    if (fabs(currentObtained ) < 1.0E-6) {
			printf("                            volts at bottom = %g probably because no electrons left, curr = %g\n", phiMin, currentObtained);
		    } else {
			printf("                            volts at bottom = %g, current = %g\n", phiMin,  currentObtained);
		    }
                }
		if (fabs(currentObtained ) < 1.0E-6) {
		    current = currentObtained;
		    deltaT = deltaT_curr;
		    printLvl_ = oldP;
		    return xbest;
		}
            }
        }
        /*
         *  We are here for ealy cases of max or min voltages, where we are running out of electrons. We are
         *  also here when the rootfinder fails for indeterminate reasons. For this latter case, we will try
         *  reducing the time step to see if this helps solve the system. In some cases, this has actually worked,
         *  while in other cases, it has not worked.
         */
        /*
         *  Fill in the return variables, current and deltaT, with the last quantities actually obtained from
         *  the rootfinder. Then reduce the time step by a factor of two and resolve the system.
         *  Maybe we can find the requested current at a smaller time step increment.  This often
         *  happens when we are running out of electrons in the electrode.
         */
        current = currentObtained;
        deltaT = deltaT_curr;
        deltaT_curr *= 0.5;
    }
    /*
     *  Return with failed indicator values. We return with the last voltage, xbest, the last current, current,
     *  and the last deltaT, deltaT, actually returned from the rootfinder.
     */
    if (printLvl_ > 1) {
        printf("                            Returning with failed rootfinder values\n");
    }
    printLvl_ = oldP;
    return xbest;
}
//====================================================================================================================
// Report the integrated source term for the electrode over an interval in time.
/*
 *  This is the net change in the moles of species defined in the electrode over that
 *  interval of time. The conditions at the end of the interval are used to carry
 *  out the integrations.
 *
 *  @param spMoleDelta The end result in terms of the change in moles of species in the
 *                     electrode.
 *
 *  @return Tfinal    Final time to integrate to.
 */
double Electrode::integratedSourceTerm(doublereal* const spMoleDelta)
{
    if (tfinal_ == tinit_) {
        throw CanteraError(" Electrode::integratedSourceTerm()", "tfinal == tinit");
    }
    /*
     *  We may do more here to ensure that the last integration is implicit
     */
    std::copy(spMoleIntegratedSourceTerm_.begin(), spMoleIntegratedSourceTerm_.end(), spMoleDelta);
    return tfinal_;
}
//====================================================================================================================
//! Report the energy source term for the electrode over an interval in time
/*!
 *  Sum over phases ( enthalpy phase * (phaseMoles_final_ - phaseMoles_init_init_) )
 *  This should only be called after integrate() has finished running.
 */
double Electrode::energySourceTerm()
{
    doublereal energySource = 0.0;
    for (int iph = 0; iph < NumVolPhases_; iph++) {
        ThermoPhase& tp = thermo(iph);
        //int nspPhase = tp.nSpecies();
        energySource += tp.enthalpy_mole() * (phaseMoles_final_[iph] - phaseMoles_init_init_[iph]);
    }
    return energySource;
}
//====================================================================================================================
double Electrode::getSourceTerm(SOURCES sourceType)
{
  double result = 0.0;

  int species_index = 0;
  int solnSpeciesStart = m_PhaseSpeciesStartIndex[solnPhase_];
  int nsp = thermo(solnPhase_).nSpecies();
  if( sourceType > SPECIES_SOURCE )
  {
    species_index = sourceType - SPECIES_SOURCE;
    AssertThrow(species_index < nsp, "Electrode::getSourceTerm");
    sourceType = SPECIES_SOURCE;
  }

  switch (sourceType) {
  case CURRENT_SOURCE:
    result = spMoleIntegratedSourceTerm_[kElectron_];
    break;
  case ELECTROLYTE_PHASE_SOURCE:
    for (int ik = 0; ik < nsp; ik++) {
      int k = solnSpeciesStart + ik;
      result += spMoleIntegratedSourceTerm_[k];
    }
    break;
  case ENTHALPY_SOURCE:
    result = energySourceTerm();
    break;
  case SPECIES_SOURCE:
    result = spMoleIntegratedSourceTerm_[solnSpeciesStart + species_index];
    break;
  case MAX_SOURCE:
    result = 0.0;
    break;
  }
  return result;
}
//====================================================================================================================
// Calculate the integrated source term for the electrode over an interval in time.
/*
 *  This is the net change in the moles of species defined in the electrode over that
 *  interval of time. The conditions at the beginning of the interval are used to carry
 *  out the integrations. This may be used at the start of the time step, since no
 *  unknown conditions are required.
 *
 *  @param deltaT     time to integrate
 *  @param spMoleDelta The end result in terms of the change in moles of species in the
 *                     electrode.
 *
 *  @return Tfinal    Final time to integrate to.
 */
double Electrode::integrateAndPredictSourceTerm(doublereal deltaT, doublereal* const spMoleDelta)
{
    integrate(deltaT);
    std::copy(spMoleIntegratedSourceTerm_.begin(), spMoleIntegratedSourceTerm_.end(), spMoleDelta);
    return t_final_final_;
}
//====================================================================================================================
// Set the phase existence flag in the electrode kinetics object so that kinetics
// are calculated correctly
/*
 *    Flags are set in the kinetics object to tell the kinetics object which phases
 *    have zero moles.  The zero mole indicator is taken from phaseMoles_final_[]. Therefore,
 *    the final state is queried.
 *    There is a special case. Phases that are in justBornPhase_[] vector are allowed to
 *    be set to exist even if their phase moles are zero.
 *
 * @param assumeStableSingleSpeciesPhases Assume that single phases are stable. This
 *                         allows their production rates to be calculated
 */
void Electrode::setPhaseExistenceForReactingSurfaces(bool assumeStableSingleSpeciesPhases)
{
    for (int isk = 0; isk < numSurfaces_; isk++) {
        /*
         *  Loop over phases, figuring out which phases have zero moles.
         *  Volume phases exist if the initial or final mole numbers are greater than zero
         *  Surface phases exist if the initial or final surface areas are greater than zero.
         */
        if (RSD_List_[isk]) {
            ReactingSurDomain* rsd = RSD_List_[isk];
            int nph = rsd->nPhases();
            for (int jph = 0; jph < nph; jph++) {
                int iph = rsd->kinOrder[jph];
                if (iph == metalPhase_) {
                    continue;
                }
                double mm = phaseMoles_init_[iph];
                double mmf = phaseMoles_final_[iph];
                ThermoPhase& tp = thermo(iph);
                int nsp = tp.nSpecies();
                if (iph >= NumVolPhases_) {
                    // we are in a surface phase
                    int isur = iph - NumVolPhases_;
                    double sa_init = surfaceAreaRS_init_[isur];
                    double sa_final = surfaceAreaRS_final_[isur];
                    if (sa_init > 0.0 || sa_final > 0.0) {
                        rsd->setPhaseExistence(jph, true);
                    } else {
                        rsd->setPhaseExistence(jph, false);
                    }
                } else {
                    if (mm <= 0.0 && mmf <= 0.0) {
                        rsd->setPhaseExistence(jph, false);
                        if (nsp == 1 && assumeStableSingleSpeciesPhases) {
                            rsd->setPhaseStability(jph, true);
                        }
                    } else {
                        rsd->setPhaseExistence(jph, true);
                    }
                }

                for (int iiph = 0; iiph < (int) justBornPhase_.size(); iiph++) {
                    if (iph == justBornPhase_[iiph]) {
                        rsd->setPhaseExistence(jph, true);
                    }
                }
            }
        }
    }
}
//====================================================================================================================
void Electrode::printElectrodePhaseList(int pSrc, bool subTimeStep)
{
    printf("\n");
    printf("         PName                   MoleNum      molarVol     Volume       FractVol     Voltage   \n");
    printf("     ============================================================================================\n");
    double egv = TotalVol();
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        std::string pname = PhaseNames_[iph];
        ThermoPhase& tp = thermo(iph);
        double mv = tp.molarVolume();
        if (iph >= NumVolPhases_) {
            mv = 0.0;
        }
        double pv = mv * phaseMoles_final_[iph];
        printf("     ");
        ca_ab::pr_sf_lj(pname, 24, 1);
        printf(" %12.3E", phaseMoles_final_[iph]);
        printf(" %12.3E", mv);
        printf(" %12.3E", pv);
        printf(" %12.3E", pv / egv);
        if (iph == metalPhase_ || iph == solnPhase_) {
            printf(" %12.5E", tp.electricPotential());
        } else {
            printf("             ");
        }
        printf("\n");
    }
    printf("     ============================================================================================\n");
}
//====================================================================================================================
// Print conditions of the electrode for the current integration step to stdout
/*
 *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
 *                       to the final_final time.
 *                       The default is to print out the source terms
 *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
 *                       time step. The default is to print out the global values
 */
void Electrode::printElectrode(int pSrc, bool subTimeStep)
{
    int iph;
    double egv = TotalVol();
    double tsm = SolidTotalMoles();
    printf("   ==============================================================================================\n");
    if (subTimeStep) {
        printf("      Electrode at intermediate-step time final = %g\n", tfinal_);
        printf("                   intermediate-step time init  = %g                deltaT = %g\n", tinit_,
               deltaTsubcycle_);
        printf("                    Model Type = %3d , DomainNumber = %2d , CellNumber = %2d , SubIntegrationCounter = %d\n",
               electrodeModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberSubIntegrations_);
    } else {
        printf("      Electrode at time final = %g\n", t_final_final_);
        printf("                   time init  = %g                         deltaTglobal = %g\n", t_init_init_,
               t_final_final_ - t_init_init_);
        printf("                    Model Type = %3d , DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n",
               electrodeModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberIntegrations_);
    }
    printf("   ==============================================================================================\n");
    printf("\n");
    printf("          Number of external surfaces = %d\n", numExternalInterfacialSurfaces_);
    printf("          Solid Volume = %11.4E m**3\n", ElectrodeSolidVolume_);
    printf("          Total Volume = %11.4E m**3\n", egv);
    if (egv > 0.0) {
        printf("          Porosity     = %10.3E\n", (egv - ElectrodeSolidVolume_) / egv);
    } else {
        printf("          Porosity     =       NA\n");
    }
    printf("          Temperature = %g K\n", temperature_);
    printf("          Pressure = %g Pa\n", pressure_);
    printf("          Total Solid Moles = %11.4E kmol\n", tsm);
    if (tsm > 0.0) {
        printf("          Molar Volume of Solid = %11.4E cm3 gmol-1\n", ElectrodeSolidVolume_ / tsm * 1.0E3);
    } else {
    }
    printf("          Particle Number to Follow = %11.4E\n", particleNumberToFollow_);
    if (subTimeStep) {
        double curr = integratedLocalCurrent();
        printf("          Current = %12.5E amps\n", curr);
    } else {
        double curr = integratedCurrent();
        printf("          Current = %12.5E amps\n", curr);
    }
    printf("          Voltage (phiMetal - phiElectrolyte) = %12.5E volts\n", deltaVoltage_);

    printf("          followElectrolyteMoles = %d\n", followElectrolyteMoles_);
    printf("          ElectrolytePseudoMoles = %g\n", electrolytePseudoMoles_);
    printf("\n");

    printElectrodeCapacityInfo(pSrc, subTimeStep);
    printf("\n");
    printElectrodePhaseList(pSrc, subTimeStep);

    int m = m_NumTotPhases;
    if (numSurfaces_ > m_NumSurPhases) {
        m = numSurfaces_ + NumVolPhases_;
    }
    for (iph = 0; iph < m; iph++) {
        printElectrodePhase(iph, pSrc, subTimeStep);
    }

}
//===================================================================================================================
void Electrode::printElectrodeCapacityInfo(int pSrc, bool subTimeStep)
{
    double capacd = capacityDischarged();
    printf("          Capacity Discharged Since Start = %12.6g coulombs = %12.6g Ah\n", capacd, capacd / 3600.);
    double dod = depthOfDischarge();
    printf("          Depth of Discharge (Current)    = %12.6g coulombs\n", dod);
    double capLeft = capacityLeft();
    double capZero = capacity();
    printf("          Capacity Left                   = %12.6g coulombs         Capacity at Zero DOD = %12.6g coulombs\n",
           capLeft, capZero);
}
//===================================================================================================================
void Electrode::printElectrodePhase(int iph, int pSrc, bool subTimeStep)
{
    int isurf;
    double* netROP = new double[m_NumTotSpecies];
    ThermoPhase& tp = thermo(iph);
    string pname = tp.id();
    int istart = m_PhaseSpeciesStartIndex[iph];
    int nsp = tp.nSpecies();
    if (printLvl_ <= 1) {
        return;
    }
    printf("     ============================================================================================\n");
    printf("          PHASE %d %s \n", iph, pname.c_str());
    printf("                Total Moles  = %11.5E kmol\n", phaseMoles_final_[iph]);
    double mv = tp.molarVolume();
    if (iph >= NumVolPhases_) {
        printf("                Molar Volume = %11.5E cm3 gmol-1\n", 0.0);
        printf("                Molar Area   = %11.5E cm2 gmol-1\n", mv * 10.);
    } else {
        printf("                Molar Volume = %11.5E cm3 gmol-1\n", mv * 1.0E3);
    }
    if (iph == metalPhase_) {
        double deltaT = t_final_final_ - t_init_init_;
        if (subTimeStep) {
            deltaT = tfinal_ - tinit_;
        }
        if (deltaT > 1.0E-200) {
            double amps = spMoleIntegratedSourceTerm_[istart] / deltaT * Faraday;
            if (subTimeStep) {
                amps = spMoleIntegratedSourceTermLast_[istart] / deltaT * Faraday;
            }
            printf("                Current = %g amps \n", amps);
        } else {
            printf("                Current = NA amps \n");
        }
    }
    if (iph == metalPhase_ || iph == solnPhase_) {
        printf("                Voltage = %g volts\n", tp.electricPotential());
    }
    /*
     * Do specific surface phase printouts
     */
    double radius;
    if (iph >= NumVolPhases_) {
        isurf = iph - NumVolPhases_;
        radius = sqrt(surfaceAreaRS_final_[isurf] / (4.0 * Pi * particleNumberToFollow_));
        radius *= 1.0E6;
        printf("                surface area (final) = %13.5E m**2      Radius(final) = %13.5E um\n",
               surfaceAreaRS_final_[isurf], radius);
        if (subTimeStep) {
            radius = sqrt(surfaceAreaRS_init_[isurf] / (4.0 * Pi * particleNumberToFollow_));
            radius *= 1.0E6;
            printf("                surface area (init)  = %13.5E m**2      Radius(init)  = %13.5E um\n",
                   surfaceAreaRS_init_[isurf], radius);
        } else {
            radius = sqrt(surfaceAreaRS_init_init_[isurf] / (4.0 * Pi * particleNumberToFollow_));
            radius *= 1.0E6;
            printf("                surface area (init)  = %13.5E m**2      Radius(init)  = %13.5E um\n",
                   surfaceAreaRS_init_init_[isurf], radius);
        }
        int ttt = isExternalSurface_[isurf];
        printf("                IsExternalSurface = %d\n", ttt);
        ttt = ActiveKineticsSurf_[isurf];
        printf("                ActiveKineticsSurf = %d\n", ttt);
        //double oc = openCircuitVoltage(isurf);
        double oc = openCircuitVoltage(isurf);
        if (oc != 0.0) {
            printf("                Open Circuit Voltage = %g Volts\n", oc);
        }
    }
    if (printLvl_ >= 3) {
        printf("\n");
        printf("                Name              MoleFrac_final kMoles_final kMoles_init SrcTermIntegrated(kmol)\n");
        for (int k = 0; k < nsp; k++) {
            string sname = tp.speciesName(k);
            if (pSrc) {
                if (subTimeStep) {
                    printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(),
                           spMf_final_[istart + k], spMoles_final_[istart + k], spMoles_init_[istart + k],
                           spMoleIntegratedSourceTermLast_[istart + k]);
                } else {
                    printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(),
                           spMf_final_[istart + k], spMoles_final_[istart + k], spMoles_init_init_[istart + k],
                           spMoleIntegratedSourceTerm_[istart + k]);
                }
            } else {
                if (subTimeStep) {
                    printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                           spMoles_final_[istart + k], spMoles_init_[istart + k]);
                } else {
                    printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                           spMoles_final_[istart + k], spMoles_init_init_[istart + k]);
                }
            }
        }
    }
    if (printLvl_ >= 4) {
        if (iph >= NumVolPhases_) {
            const vector<double>& rsSpeciesProductionRates = RSD_List_[isurf]->calcNetProductionRates();
            RSD_List_[isurf]->getNetRatesOfProgress(netROP);

            doublereal* spNetProdPerArea = (doublereal*) spNetProdPerArea_List_.ptrColumn(isurf);
            std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
            int nphRS = RSD_List_[isurf]->nPhases();
            int kIndexKin = 0;
            for (int kph = 0; kph < nphRS; kph++) {
                int jph = RSD_List_[isurf]->kinOrder[kph];
                int istart = m_PhaseSpeciesStartIndex[jph];
                int nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
                for (int k = 0; k < nsp; k++) {
                    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                    kIndexKin++;
                }
            }
            printf("\n");
            printf("                           spName                  SourceRateLastStep (kmol/m2/s) \n");
            for (int k = 0; k < m_NumTotSpecies; k++) {
                string ss = speciesName(k);
                printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
            }
        }
    }
    printf("     ============================================================================================\n");
    delete[] netROP;

}
//====================================================================================================================
// Determines the level of printing for each step.
/*
 *   0 -> absolutely nothing is printed for a single call to integrate.
 *   1 -> One line summary per integrate call
 *   2 -> short description, points of interest: Table of nonlinear solve - one line per iteration
 *   3 -> Table is included -> More printing per nonlinear iteration (default) that occurs during the table
 *   4 -> Summaries of the nonlinear solve iteration as they are occurring -> table no longer printed
 *   5 -> Algorithm information on the nonlinear iterates are printed out
 *   6 -> Additional info on the nonlinear iterates are printed out
 *   7 -> Additional info on the linear solve is printed out.
 *   8 -> Info on a per iterate of the linear solve is printed out.
 */
void Electrode::setPrintLevel(int printLvl)
{
    printLvl_ = printLvl;
}
//====================================================================================================================
//   Set the deltaT used for the subcycle step
/*
 *  The default setting for this is infinite. If it's infinite then deltaTSubcycle is set to the
 *  deltaT of the integration step. However, there are many times where a smaller natural delta T is
 *  required to solve the equation system
 */
void Electrode::setDeltaTSubcycle(doublereal deltaTsubcycle)
{
    deltaTsubcycleNext_ = deltaTsubcycle;
    deltaTsubcycle_init_next_ = deltaTsubcycle;
    deltaTsubcycle_init_init_ = deltaTsubcycle;
}
//====================================================================================================================
//   Set the maximum deltaT used for the subcycle step
void Electrode::setDeltaTSubcycleMax(doublereal deltaTsubcycle)
{
    deltaTsubcycleMax_ = deltaTsubcycle;
}
//====================================================================================================================
//  Calculate the time derivatives of the mole numbers at the current subcycle step
/*
 *   This may be used as a virtual function
 */
void Electrode::calculateTimeDerivatives(doublereal deltaTsubcycle)
{
    if (deltaTsubcycle <= 0.0) {
        deltaTsubcycle = 1.0;
    }
    double invT = 1.0 / deltaTsubcycle;
    for (int k = 0; k < m_NumTotSpecies; k++) {
        spMoles_dot_[k] = invT * (spMoles_final_[k] - spMoles_init_[k]);
    }
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        phaseMoles_dot_[iph] = invT * (phaseMoles_final_[iph] - phaseMoles_init_[iph]);
    }
}
//====================================================================================================================
// Write out CSV tabular data on the integrations
/*
 *  The idea is to print out tabular data about each intermediate run and about each
 *  global run
 *
 *  @param itype Type of the data
 *            - 0      Initialization information - overwrite file
 *            - 1      Normal intermediate information
 *            - 2      Normal global information at the end of a global time step
 *            - -1     Failed intermediate calculation, a failure from a nonlinear solver step
 *            - -2     Failed calculation from the predictor step - not necessarily significant.
 */
void Electrode::writeCSVData(int itype)
{
    //if (printLvl_ < 2) return;
    int k;
    int kstart;
    static std::string globOutputName = "";
    static std::string intOutputName = "";
    static FILE* fpG = 0;
    static FILE* fpI = 0;
    static int runNumber = -1;
    bool ignoreErrors = false;
    if (itype < 0) {
        ignoreErrors = true;
    }

    if (itype == 0 || fpI == 0) {
        runNumber++;
        if (runNumber == 0) {
            intOutputName = (electrodeName_ + "_intResults_" + int2str(electrodeDomainNumber_) + "_"
                    + int2str(electrodeCellNumber_) + ".csv");
        } else {
            intOutputName = (electrodeName_ + "_intResults_" + int2str(electrodeDomainNumber_) + "_"
                    + int2str(electrodeCellNumber_) + "_" + int2str(runNumber) + ".csv");
        }
        fpI = fopen(intOutputName.c_str(), "w");

        fprintf(fpI, "         Tinit ,        Tfinal ,");

        fprintf(fpI, "      Volts_Soln ,  Volts_Electrode ,");

        fprintf(fpI, "      Current ,");

        fprintf(fpI, "  CapDischarged ,");
        fprintf(fpI, "  CapacityLeft ,");
        fprintf(fpI, "  Capacity ,");

        for (k = 0; k < m_NumTotSpecies; k++) {
            string sss = speciesName(k);
            fprintf(fpI, " MN_%-20.20s,", sss.c_str());
        }

        for (k = 0; k < m_NumTotSpecies; k++) {
            string sss = speciesName(k);
            fprintf(fpI, " SRC_%-20.20s,", sss.c_str());
        }

        fprintf(fpI, " iType,");

        fprintf(fpI, " SolidVol,");
        fprintf(fpI, " ElectrolyteVol,");
        fprintf(fpI, " TotalVol,");
        fprintf(fpI, " Porosity");

        fprintf(fpI, "\n");
        fclose(fpI);
    }

    if (itype == 0 || fpG == 0) {
        if (runNumber == 0) {
            globOutputName = (electrodeName_ + "_globalResults_" + int2str(electrodeDomainNumber_) + "_"
                    + int2str(electrodeCellNumber_) + ".csv");
        } else {
            intOutputName = (electrodeName_ + "_globalResults_" + int2str(electrodeDomainNumber_) + "_"
                    + int2str(electrodeCellNumber_) + "_" + int2str(runNumber) + ".csv");
        }
        fpG = fopen(globOutputName.c_str(), "w");
        fprintf(fpG, "         Tinit ,        Tfinal ,");

        fprintf(fpG, "      Volts_Soln ,  Volts_Electrode ,");

        fprintf(fpG, "      Current ,");

        fprintf(fpG, "  CapDischarged ,");
        fprintf(fpG, "  CapacityLeft ,");
        fprintf(fpG, "  Capacity ,");

        for (k = 0; k < m_NumTotSpecies; k++) {
            string sss = speciesName(k);
            fprintf(fpG, " MN_%-20.20s,", sss.c_str());
        }

        for (k = 0; k < m_NumTotSpecies; k++) {
            string sss = speciesName(k);
            fprintf(fpG, " SRC_%-20.20s,", sss.c_str());
        }

        fprintf(fpG, " iType,");

        fprintf(fpG, " SolidVol,");
        fprintf(fpG, " ElectrolyteVol,");
        fprintf(fpG, " TotalVol,");
        fprintf(fpG, " Porosity");
        fprintf(fpG, "\n");
        fclose(fpG);
    }

    // Provide a hook here for the debugger
    if (itype == -1) {
        kstart = m_PhaseSpeciesStartIndex[metalPhase_];
    }

    if (itype == 1 || itype == 2) {
        fpI = fopen(intOutputName.c_str(), "a");
        fprintf(fpI, "  %12.5E ,  %12.5E ,", tinit_, tfinal_);

        double phiS = phaseVoltages_[solnPhase_];
        double phiM = phaseVoltages_[metalPhase_];
        fprintf(fpI, "    %12.5E ,     %12.5E ,", phiS, phiM);

        kstart = m_PhaseSpeciesStartIndex[metalPhase_];
        double deltaT = tfinal_ - tinit_;
        double amps = 0.0;
        if (deltaT > 1.0E-200) {
            amps = spMoleIntegratedSourceTermLast_[kstart] / deltaT * Faraday;
        }
        fprintf(fpI, " %12.5E ,", amps);

        double cap = capacityDischarged();
        fprintf(fpI, "   %12.5E ,", cap);

        cap = capacityLeft();
        fprintf(fpI, "   %12.5E ,", cap);

        cap = capacity();
        fprintf(fpI, "   %12.5E ,", cap);

        for (k = 0; k < m_NumTotSpecies; k++) {
            fprintf(fpI, " %12.5E           ,", spMoles_final_[k]);
        }

        for (k = 0; k < m_NumTotSpecies; k++) {
            fprintf(fpI, " %12.5E            ,", spMoleIntegratedSourceTermLast_[k]);
        }

        fprintf(fpI, "   %3d ,", itype);

        double solidVol = SolidVol();
        double grossVol = TotalVol(ignoreErrors);
        double lyteVol = grossVol - solidVol;
        double porosity = lyteVol / grossVol;
        //  This will fail until we figure out how we want to handle this
        // if (fabs(porosity - porosity_) > 1.0E-6) {
        //	throw CanteraError("writeCSVData() ERROR", "porosity isn't consistent");
        //}
        fprintf(fpI, " %12.5E ,", solidVol);
        fprintf(fpI, " %12.5E ,", lyteVol);
        fprintf(fpI, " %12.5E ,", grossVol);
        fprintf(fpI, " %12.5E ", porosity);
        fprintf(fpI, " \n");

        fclose(fpI);
    }

    if (itype == 2) {
        fpG = fopen(globOutputName.c_str(), "a");
        fprintf(fpG, "  %12.5E ,  %12.5E ,", t_init_init_, t_final_final_);

        double phiS = phaseVoltages_[solnPhase_];
        double phiM = phaseVoltages_[metalPhase_];
        fprintf(fpG, "    %12.5E ,     %12.5E ,", phiS, phiM);

        kstart = m_PhaseSpeciesStartIndex[metalPhase_];
        double deltaT = t_final_final_ - t_init_init_;
        double amps = 0.0;
        if (deltaT > 1.0E-200) {
            amps = spMoleIntegratedSourceTerm_[kstart] / deltaT * Faraday;
        }
        fprintf(fpG, " %12.5E ,", amps);

        double cap = capacityDischarged();
        fprintf(fpG, "   %12.5E ,", cap);

        cap = capacityLeft();
        fprintf(fpG, "   %12.5E ,", cap);

        cap = capacity();
        fprintf(fpG, "   %12.5E ,", cap);

        for (k = 0; k < m_NumTotSpecies; k++) {
            fprintf(fpG, " %12.5E           ,", spMoles_final_[k]);
        }

        for (k = 0; k < m_NumTotSpecies; k++) {
            fprintf(fpG, " %12.5E            ,", spMoleIntegratedSourceTerm_[k]);
        }

        fprintf(fpG, "   %3d,", itype);
        double solidVol = SolidVol();
        double grossVol = TotalVol();
        double lyteVol = grossVol - solidVol;
        double porosity = lyteVol / grossVol;

        fprintf(fpG, " %12.5E ,", solidVol);
        fprintf(fpG, " %12.5E ,", lyteVol);
        fprintf(fpG, " %12.5E ,", grossVol);
        fprintf(fpG, " %12.5E ", porosity);
        fprintf(fpG, " \n");
        fclose(fpG);
    }

}

//====================================================================================================================
}// End of namespace Cantera
//======================================================================================================================
