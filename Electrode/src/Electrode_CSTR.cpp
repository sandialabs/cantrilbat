/**
 *  @file Electrode_CSTR.cpp
 *     Definitions for the  Electrode_CSTR class, used to model 
 *     Electrode processes in particles with no transport limits
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode_CSTR Electrode_CSTR\endlink).
 */

#include "Electrode_CSTR.h"

#include "LE_OneDbl.h"
#include "BE_BlockEntry.h"

#include "cantera/numerics/NonlinearSolver.h"
#include "cantera/base/vec_functions.h"

using namespace std;
using namespace BEInput;

#ifndef MAX
//! define a fast max operator for doubles
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif
#ifndef MIN
//! define a fast min operator for doubles
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

#ifndef SAFE_DELETE
//! Delete and zero a pointer
#define SAFE_DELETE(x)  if ((x)) { delete (x) ; x = nullptr ; }
#endif

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
ELECTRODE_CSTR_KEY_INPUT::ELECTRODE_CSTR_KEY_INPUT(int printLvl) :
    ELECTRODE_KEY_INPUT(printLvl),
    boundaryResistance_(0.0)
{
}
//==================================================================================================================================
ELECTRODE_CSTR_KEY_INPUT::~ELECTRODE_CSTR_KEY_INPUT()
{

}
//==================================================================================================================================
void ELECTRODE_CSTR_KEY_INPUT::setup_input_child1(BEInput::BlockEntry* cf)
{
    /*
     * Obtain the number of cantera files to be read
     */
    LE_OneDbl* d1 = new LE_OneDbl("Boundary Resistance", &(boundaryResistance_), 0, "boundaryResistance");
    d1->set_default(0.0);
    cf->addLineEntry(d1);
    BaseEntry::set_SkipUnknownEntries(3);
}
//====================================================================================================================
void ELECTRODE_CSTR_KEY_INPUT::setup_input_child2(BEInput::BlockEntry* cf)
{
    BaseEntry::set_SkipUnknownEntries(3);
}
//======================================================================================================================
//======================================================================================================================
//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_CSTR::Electrode_CSTR() :
    Electrode_Integrator(),

    electrodeType_(ELECTRODETYPE_ANODE),
    RelativeExtentRxn_init_(0.0),
    RelativeExtentRxn_init_init_(0.0),
    RelativeExtentRxn_final_(0.0),
    RelativeExtentRxn_final_final_(0.0),
    RelativeExtentRxn_NormalizationFactor_(1.0),
    RelativeExtentRxn_RegionBoundaries_(0),
    onRegionBoundary_init_(-1),
    onRegionBoundary_final_(-1),
    onRegionBoundary_init_init_(-1),
    onRegionBoundary_final_final_(-1),
    xRegion_init_(-1),
    xRegion_init_init_(-1),
    xRegion_final_(-1),
    xRegion_final_final_(-1),
    goNowhere_(0),
    deltaT_RegionBoundaryCollision_(1.0E300),
    atolBaseResid_(1.0E-12),
    DspMoles_final_(0),
    SrcDot_RxnExtent_final_(0.0),
    deltaSpMoles_(0),
    minPH_(npos)
{
}
//======================================================================================================================
Electrode_CSTR::Electrode_CSTR(const Electrode_CSTR& right) :
    Electrode_Integrator(),
    electrodeType_(right.electrodeType_),
    RelativeExtentRxn_init_(0.0),
    RelativeExtentRxn_init_init_(0.0),
    RelativeExtentRxn_final_(0.0),
    RelativeExtentRxn_final_final_(0.0),
    RelativeExtentRxn_NormalizationFactor_(1.0),
    RelativeExtentRxn_RegionBoundaries_(0),
    onRegionBoundary_init_(-1),
    onRegionBoundary_final_(-1),
    onRegionBoundary_init_init_(-1),
    onRegionBoundary_final_final_(-1),
    xRegion_init_(-1),
    xRegion_init_init_(-1),
    xRegion_final_(-1),
    xRegion_final_final_(-1),
    goNowhere_(0),
    deltaT_RegionBoundaryCollision_(1.0E300),
    atolBaseResid_(1.0E-12),
    DspMoles_final_(0),
    SrcDot_RxnExtent_final_(0.0),

    deltaSpMoles_(0),
    minPH_(npos)
{
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
}
//======================================================================================================================
Electrode_CSTR& Electrode_CSTR::operator=(const Electrode_CSTR& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode_Integrator::operator=(right);

    electrodeType_                       = right.electrodeType_;
    RelativeExtentRxn_init_              = right.RelativeExtentRxn_init_;
    RelativeExtentRxn_init_init_         = right.RelativeExtentRxn_init_init_;
    RelativeExtentRxn_final_             = right.RelativeExtentRxn_final_;
    RelativeExtentRxn_final_final_       = right.RelativeExtentRxn_final_final_;
    RelativeExtentRxn_NormalizationFactor_=right.RelativeExtentRxn_NormalizationFactor_;
    RelativeExtentRxn_RegionBoundaries_  = right.RelativeExtentRxn_RegionBoundaries_;
    onRegionBoundary_init_               = right.onRegionBoundary_init_;
    onRegionBoundary_final_              = right.onRegionBoundary_final_;
    onRegionBoundary_init_init_          = right.onRegionBoundary_init_init_;
    onRegionBoundary_final_final_        = right.onRegionBoundary_final_final_;
    xRegion_init_                        = right.xRegion_init_;
    xRegion_init_init_                   = right.xRegion_init_init_;
    xRegion_final_                       = right.xRegion_final_;
    xRegion_final_final_                 = right.xRegion_final_final_;
    goNowhere_                           = right.goNowhere_;

    deltaT_RegionBoundaryCollision_      = right.deltaT_RegionBoundaryCollision_;

    atolBaseResid_                       = right.atolBaseResid_;

    ROP_                                 = right.ROP_;
    DspMoles_final_                      = right.DspMoles_final_;

    SrcDot_RxnExtent_final_              = right.SrcDot_RxnExtent_final_;

    phaseIndexSolidPhases_               = right.phaseIndexSolidPhases_;
    numSpecInSolidPhases_                = right.numSpecInSolidPhases_;
    deltaSpMoles_                        = right.deltaSpMoles_;
    minPH_                               = right.minPH_;

    phaseMFBig_                          = right.phaseMFBig_;
    justDied_                            = right.justDied_;
    phaseMoles_final_lagged_             = right.phaseMoles_final_lagged_;
    DphMoles_final_                      = right. DphMoles_final_;

    return *this;
}
//======================================================================================================================
/*
 *
 *  ELECTRODE_INPUT:destructor
 *
 * We need to manually free all of the arrays.
 */
Electrode_CSTR::~Electrode_CSTR()
{

}
//======================================================================================================================
Electrode_Types_Enum Electrode_CSTR::electrodeType() const
{
    return CSTR_ET;
}

//====================================================================================================================
int Electrode_CSTR::electrode_input_child(ELECTRODE_KEY_INPUT** ei_ptr)
{
    /*
     *  malloc an expanded child input
     */
    ELECTRODE_CSTR_KEY_INPUT* ei_mp = new ELECTRODE_CSTR_KEY_INPUT();
    /*
     *  Find the command file
     */
    ELECTRODE_KEY_INPUT* ei = *(ei_ptr);
    string commandFile = ei->commandFile_;
    BlockEntry* cf = ei->lastBlockEntryPtr_;
    ei_mp->printLvl_ = ei->printLvl_;
    /*
     *  Parse the complete child input file
     */
    int ok = ei_mp->electrode_input_child(commandFile, cf);
    if (ok == -1) {
        return -1;
    }
    /*
     * Switch the pointers around so that the child input file is returned.
     * Delete the original pointer.
     */
    delete ei;
    *ei_ptr = ei_mp;
    return 0;
}
//======================================================================================================================
int Electrode_CSTR::electrode_model_create(ELECTRODE_KEY_INPUT* eibase)
{

    int flag = Electrode_Integrator::electrode_model_create(eibase);
    if (flag != 0) {
        return flag;
    }

    /*
     *  Downcast the Key input to make sure we are being fed the correct child object
     */
    ELECTRODE_CSTR_KEY_INPUT* ei = dynamic_cast<ELECTRODE_CSTR_KEY_INPUT*>(eibase);
    if (!ei) {
        throw ZuzaxError(" Electrode_CSTR::electrode_model_create()",
                         " Expecting a child ELECTRODE_CSTR_KEY_INPUT object and didn't get it");
    }

    /*
     *  We determine two arrays here that are useful.
     *
     *       phaseIndexXSolidPhases_[]
     *       numSpecInSolidPhases_[]
     *
     *  phaseIndexSolidPhases_ lists the volume phases in the PhaseList object
     *  which are actually solid phases that are part of the electrode. Right now, everything that is not
     *  the metal phase or the electrolyte phase is a solid phase.
     *  numSpecInSolidPhases_[] is the number of species in each of those solid phases.
     */
    phaseIndexSolidPhases_.clear();
    numSpecInSolidPhases_.clear();
    for (size_t ph = 0; ph < m_NumVolPhases; ph++) {
        ThermoPhase* tp = VolPhaseList_[ph];
        size_t iph = globalPhaseIndex(tp);
        if (iph == metalPhaseIndex() || iph == solnPhaseIndex()) {
            //do nothing
        } else {
            phaseIndexSolidPhases_.push_back(iph);
            numSpecInSolidPhases_.push_back(tp->nSpecies());
        }
    }

    /*
     * Initialize the arrays in this object now that we know the number of equations
     */
    init_sizes();

    /*
     *  Calculate the Normalization factor. It is the factor that relates the change in the relative extent of
     *  reaction (which is dimensionless) wrt solid phase mole number. The default treatment is that
     *  there is a strict proportionality, i.e., that the value is equal to the total solid phase moles.
     *  This value can be overriden in child objects, but so far, this is not necessary.
     */
    double solidMoles = 0.0;
    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        solidMoles += phaseMoles_final_[iph];
    }
    RelativeExtentRxn_NormalizationFactor_ = solidMoles;

    if (ei->RxnExtTopLimit >= 0.0) {
        RelativeExtentRxn_RegionBoundaries_.resize(2);
        RelativeExtentRxn_RegionBoundaries_[1] = ei->RxnExtTopLimit;
        RelativeExtentRxn_RegionBoundaries_[0] = ei->RxnExtBotLimit;
    }

    /*
     *  We assume that the current state of the electrode is that it is located in the 0 region for the
     *  extent of reaction variable.
     */
    xRegion_init_ = 0;
    xRegion_init_init_ = 0;
    xRegion_final_ = 0;
    xRegion_final_final_ = 0;

    return 0;
}
//======================================================================================================================
//   Set the electrode initial conditions from the input file.
/*
 *   (virtual from Electrode)
 *   (This is a serial virtual function or an overload function)
 *
 *    This is one of the most important routines. It sets up the initial conditions of the electrode
 *    from the input file. The electrode itself has been set up from a call to electrode_model_create().
 *    After the call to this routine, the electrode should be internally ready to be integrated and reacted.
 *    It takes its input from an ELECTRODE_KEY_INPUT object which specifies the setup of the electrode
 *    object and the initial state of that object.
 *
 *    The routine works like an onion initialization. The parent object is initialized before the
 *    child. This means the child object first calls the parent, before it does its own initializations.
 *
 */
int Electrode_CSTR::setInitialConditions(ELECTRODE_KEY_INPUT* eibase)
{

    ELECTRODE_CSTR_KEY_INPUT* ei = dynamic_cast<ELECTRODE_CSTR_KEY_INPUT*>(eibase);
    if (!ei) {
        throw Electrode_Error(" Electrode_CSTR::electrode_model_create()",
                           " Expecting a child ELECTRODE_KEY_INPUT object and didn't get it");
    }



    int flag = Electrode::setInitialConditions(ei);
    if (flag != 0) {
        return flag;
    }



    /*
     *  Calculate the Normalization factor. It is the factor that relates the change in the relative extent of
     *  reaction (which is dimensionless) wrt solid phase mole number. The default treatment is that
     *  there is a strict proportionality, i.e., that the value is equal to the total solid phase moles.
     *  This value can be overriden in child objects, but so far, this is not necessary.
     */
    double solidMoles = 0.0;
    for (int ph = 0; ph < (int) phaseIndexSolidPhases_.size(); ph++) {
        int iph = phaseIndexSolidPhases_[ph];
        solidMoles += phaseMoles_final_[iph];
    }
    RelativeExtentRxn_NormalizationFactor_ = solidMoles;



    if (ei->RxnExtTopLimit >= 0.0) {
        RelativeExtentRxn_RegionBoundaries_.resize(2);
        RelativeExtentRxn_RegionBoundaries_[1] = ei->RxnExtTopLimit;
        RelativeExtentRxn_RegionBoundaries_[0] = ei->RxnExtBotLimit;
    }

    /*
     *  We assume that the current state of the electrode is that it is located in the 0 region for the
     *  extent of reaction variable.
     */
    xRegion_init_ = 0;
    xRegion_init_init_ = 0;
    xRegion_final_ = 0;
    xRegion_final_final_ = 0;

    /*
     *   We keep a record of the extent of reaction
     */
    RelativeExtentRxn_final_ = calcRelativeExtentRxn_final();
    //
    //  Calculate  onRegionBoundary_final_ 
    // 
    stateToPhaseFlagsReconciliation(false);


    setInitStateFromFinal(true);

    return 0;
}
//====================================================================================================================
int
Electrode_CSTR::electrode_stateSave_create()
{
    eState_final_ = new EState();
    int rr = eState_final_->initialize(this);
    if (rr >= 0) {
        rr = 0;
    }
    return rr;
}
//======================================================================================================================
//   local routine to resize arrays that this object is responsible for
void Electrode_CSTR::init_sizes()
{
    neq_ = nEquations_calc();
    DspMoles_final_.resize(m_NumTotSpecies, 0.0);

    int maxNumRxns = RSD_List_[0]->nReactions();
    ROP_.resize(maxNumRxns, 0.0);

    deltaSpMoles_.resize(m_NumTotSpecies, 0.0);

    phaseMFBig_.resize(m_NumVolPhases, 0);

    phaseMoles_final_lagged_.resize(m_NumTotPhases, 0.0);
    DphMoles_final_.resize(m_NumTotPhases, 0.0);
    phaseMFBig_.resize(m_NumTotPhases, 0);
    justDied_.resize(m_NumTotPhases, 0);
    justBornPhase_.resize(m_NumTotPhases, 0);
    //soln_predict_.resize(neq_, 0.0);
    //spMoles_predict_.resize(m_NumTotSpecies, 0.0);
}

//======================================================================================================================
int Electrode_CSTR::create_solvers()
{
    int neqNLS = nEquations_calc();
    Electrode_Integrator::create_solvers();
    atolNLS_.resize(neqNLS, 1.0E-12);
    atolResidNLS_.resize(neqNLS, 1.0E-12);
    return neqNLS;
}
//======================================================================================================================
// Set the sizes of the electrode from the input parameters
/*
 *  We resize all of the information within the electrode from the input parameters
 *
 */
void Electrode_CSTR::setElectrodeSizeParams(double electrodeArea, double
        electrodeThickness, double porosity)
{
    Electrode_Integrator::setElectrodeSizeParams(electrodeArea, electrodeThickness, porosity);
    /*
     *  Relcalculate the normalization factor
     */
    double solidMoles = 0.0;
    for (int ph = 0; ph < (int) phaseIndexSolidPhases_.size(); ph++) {
        int iph = phaseIndexSolidPhases_[ph];
        solidMoles += phaseMoles_final_[iph];
    }
    RelativeExtentRxn_NormalizationFactor_ = solidMoles;
}
//======================================================================================================================
// Resize the solid phase and electrolyte mole numbers within the object
/*
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
void Electrode_CSTR::resizeMoleNumbersToGeometry()
{
    Electrode::resizeMoleNumbersToGeometry();
    /*
     *  Recalculate the normalization factor
     */
    double solidMoles = 0.0;
    for (int ph = 0; ph < (int) phaseIndexSolidPhases_.size(); ph++) {
        int iph = phaseIndexSolidPhases_[ph];
        solidMoles += phaseMoles_final_[iph];
    }
    RelativeExtentRxn_NormalizationFactor_ = solidMoles;
}
//====================================================================================================================
void Electrode_CSTR::updateState()
{
    /*
     *  We rely on the base implementation to update most of the fields.
     *  This includes calculation of Radius_exterior_final_ 
     */
    Electrode::updateState();
    /*
     *   We keep a record of the extent of reaction
     */
    RelativeExtentRxn_final_ = calcRelativeExtentRxn_final();

    /*
     *  Update the surface areas
     */
    updateSurfaceAreas();
}
//====================================================================================================================
//  Recalculate the surface areas of the surfaces for the final state
/*
 *    (virtual function from Electrode)
 *
 *    We used the internal variable locationOfReactingSurface_ to determine the behavior.
 *    A value of zero indicates that the surface 0 follows the reaction front as it goes from outer to inner as
 *    a function of the % though the plateau.
 *    A value of locationOfReactingSurface_ = 1 indicates that the surface 0 follows the exterior surface of the particle
 *
 *    We also assume that the surface area is equal to the particle surface area multiplied by the numbers of particles.
 *
 *    Dependent StateVariables Used
 *         Radius_exterior_final_;
 *         particleNumberToFollow_
 *
 *    Dependent StateVariables Calculated
 *          surfaceAreaRS_final_[]
 *          Radius_internal_final_
 */
void Electrode_CSTR::updateSurfaceAreas()
{
    double totalSA = 4. * Pi * Radius_exterior_final_ * Radius_exterior_final_ * particleNumberToFollow_;
    surfaceAreaRS_final_[0] = totalSA;
}
//==================================================================================================================================
bool Electrode_CSTR::stateToPhaseFlagsReconciliation(bool flagErrors)
{
    double tol = 1.0E-6;
    int onRegionBoundary_final = -1;

    bool retn = Electrode::stateToPhaseFlagsReconciliation(flagErrors);
    int sz = (int) RelativeExtentRxn_RegionBoundaries_.size();
    int xRegion_final = sz;
    for (int i = 0; i < sz; i++) {
        if (RelativeExtentRxn_final_ < RelativeExtentRxn_RegionBoundaries_[i]) {
            xRegion_final = i - 1;
        }
        if (fabs(RelativeExtentRxn_final_ - RelativeExtentRxn_RegionBoundaries_[i]) < tol) {
            onRegionBoundary_final = RelativeExtentRxn_RegionBoundaries_[i];
        }
    }
    // If we are before the initial region boundary, we set it to be on the initial boundary
    if (xRegion_final < 0) {
        xRegion_final = 0;
        onRegionBoundary_final = 0;
    }
    // If we are beyond the final boundary, we set it to be on the final boundary.
    if (xRegion_final >= sz) {
        xRegion_final = sz - 1;
        onRegionBoundary_final = sz - 1;
    }

    if (flagErrors) {
        if (xRegion_final != xRegion_final_) {
            printf("Electrode_CSTR::stateToPhaseFlagsReconciliation() WARNING\n"
                   "                xRegin_final %d is different than storred value %d\n", xRegion_final,  xRegion_final_);
            retn = true;
        }
        if (onRegionBoundary_final != onRegionBoundary_final_) {
            printf("Electrode_CSTR::stateToPhaseFlagsReconciliation() WARNING\n"
                   "                onRegionBoundary_final %d is different than storred value %d\n",
                   onRegionBoundary_final, onRegionBoundary_final_);
            retn = true;
        }
    }
    onRegionBoundary_final_ = onRegionBoundary_final;
    xRegion_final_ = xRegion_final;
    return retn;
}
//==================================================================================================================================
//  Value of the standard state open circuit voltage for the standard state conditions for region xRegion
//  at the final_ conditions.
/*
 *  This is the standard state open circuit potential at the current _final_ conditions.
 *  The value is dependent on the region input. The region input, currently is identified
 *  with the surface reacting phase object.
 *  Therefore, this call is basically a wrapper around openCircuitVoltageSS(isk) which
 *  calculates the open circuit standard state voltage for the isk reacting surface.
 *
 *  Additions include extening the region values to include false values for DoD = 0 and 1
 *  conditions.
 *
 */
double Electrode_CSTR::openCircuitVoltageSS_Region(int xRegion) const
{
    int doRegion = xRegion;
    double vOffset = 0.0;
    if (xRegion < 0) {
        doRegion = 0;
        if (electrodeType_ == 0) {
            // For anodes setting voltage negative of open circuit at DoD=0 will cause no reactions to occur
            vOffset = -5.0;
        } else {
            // For cathodes setting voltage positive of open circuit at DoD=0 will cause no reactions to occur
            vOffset =  5.0;
        }
    }
    if (xRegion >= (int) RelativeExtentRxn_RegionBoundaries_.size() - 1) {
        doRegion = RelativeExtentRxn_RegionBoundaries_.size() - 2;
        if (electrodeType_ == 0) {
            // For anodes setting voltage positve of open circuit at DoD=1 will cause no reactions to occur
            vOffset = 5.0;
        } else {
            // For cathodes setting voltage negative of open circuit at DoD=1 will cause no reactions to occur
            vOffset =  -5.0;
        }
    }
    /*
     *  go get deltaG_ss for the current surface
     */
    int isk = doRegion;
    double ocss = openCircuitVoltageSSRxn(isk);
    return (vOffset + ocss);
}
//====================================================================================================================
// Value of the open circuit voltage for region xRegion at the final_ conditions.
/*
 *  This is the open circuit potential at the current _final_ conditions.
 *  The value is dependent on the region input. The region input, currently is identified
 *  with the surface reacting phase object.
 *  Therefore, this call is basically a wrapper around openCircuitVoltage(isk) which
 *  calculates the open circuit voltage for the isk reacting surface.
 *
 *  Additions include extending the region values to include false values for DoD = 0 and 1
 *  conditions.
 *
 */
double Electrode_CSTR::openCircuitVoltage_Region(int xRegion, bool comparedToReferenceElectrode) const
{
    int doRegion = xRegion;
    double vOffset = 0.0;
    if (xRegion < 0) {
        doRegion = 0;
        if (electrodeType_ == 0) {
            // For anodes setting voltage negative of open circuit at DoD=0 will cause no reactions to occur
            vOffset = -5.0;
        } else {
            // For cathodes setting voltage positive of open circuit at DoD=0 will cause no reactions to occur
            vOffset =  5.0;
        }
    }
    if (xRegion >= (int) RelativeExtentRxn_RegionBoundaries_.size() - 1) {
        doRegion = RelativeExtentRxn_RegionBoundaries_.size() - 2;
        if (electrodeType_ == 0) {
            // For anodes setting voltage positve of open circuit at DoD=1 will cause no reactions to occur
            vOffset = 5.0;
        } else {
            // For cathodes setting voltage negative of open circuit at DoD=1 will cause no reactions to occur
            vOffset =  -5.0;
        }
    }
    /*
     *  go get deltaG_ss and openCircuit for the current surface
     */
    int isk = doRegion;
    double ocss = openCircuitVoltageRxn(isk, comparedToReferenceElectrode);
    return (vOffset + ocss);
}
//====================================================================================================================
double Electrode_CSTR::capacity(int platNum) const
{
    if (electrodeChemistryModelType_ == 2) {
        setCapacityCoeff_FeS2();
    }
    double capRaw = capacityRaw(platNum);
    double fac = capRaw / (RelativeExtentRxn_NormalizationFactor_ * Faraday);
    double unFillable = 0.0;
    if (RelativeExtentRxn_RegionBoundaries_.size() >= 2) {
        if (RelativeExtentRxn_RegionBoundaries_[0] > 0.0) {
            unFillable = RelativeExtentRxn_RegionBoundaries_[0];
            unFillable  *= RelativeExtentRxn_NormalizationFactor_;
        }
    }
    double unExtractable = 0.0;
    if (RelativeExtentRxn_RegionBoundaries_.size() >= 2) {
        if (RelativeExtentRxn_RegionBoundaries_[1] < fac) {
            unExtractable = fac - RelativeExtentRxn_RegionBoundaries_[1];
            unExtractable  *= RelativeExtentRxn_NormalizationFactor_;
        }
    }
    double capZeroDoD = capRaw - Faraday * (unFillable + unExtractable);
    return capZeroDoD;
}
//==================================================================================================================================
double Electrode_CSTR::capacityRaw(int platNum) const
{
    if (electrodeChemistryModelType_ == 2) {
        setCapacityCoeff_FeS2();
    }
    double capZeroDoD = 0.0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        for (size_t k = 0; k < NumSpeciesList_[iph]; k++) {
            capZeroDoD += spMoles_final_[kStart + k] * capacityZeroDoDSpeciesCoeff_[kStart + k];
        }
    }
    double tmp = capZeroDoD * Faraday;
    return tmp;
}
//==================================================================================================================================
// Amount of charge that the electrode that has available to be discharged
/*
 *   We report the number in terms of Amp seconds = coulombs
 */
double Electrode_CSTR::capacityLeft(int platNum, double voltsMax, double voltsMin) const
{
    if (electrodeChemistryModelType_ == 2) {
        setCapacityCoeff_FeS2();
    }
    double capLeftRaw = capacityLeftRaw(platNum);
    double capRaw = capacityRaw(platNum);
    double fac = capRaw / (RelativeExtentRxn_NormalizationFactor_ * Faraday);

    int lastRegionBoundary = RelativeExtentRxn_RegionBoundaries_.size() - 1;
    double unExtractable = 0.0;
    if (RelativeExtentRxn_RegionBoundaries_.size() >= 2) {
        if (RelativeExtentRxn_RegionBoundaries_[lastRegionBoundary] < fac) {
            unExtractable = fac - RelativeExtentRxn_RegionBoundaries_[1];
            unExtractable  *= RelativeExtentRxn_NormalizationFactor_;
        }
    }
    double capLeft = capLeftRaw / Faraday;
    capLeft -= unExtractable;
    /*
     *  For print out purposes and for stability we switch to integer logic at the end of discharge.
     *  We do it here and we should do it everywhere
     */
    if (onRegionBoundary_final_ == lastRegionBoundary) {
        if (fabs(capLeft) > 1.0E-6 * RelativeExtentRxn_NormalizationFactor_) {
            printf("Electrode_CSTR::capacityLeft() - Logic error maybe\n");
            exit(-1);
        }
        capLeft = 0.0;
    }

    if (capLeft < 0.0) {
        capLeft = 0.0;
    }
    return capLeft * Faraday;
}
//====================================================================================================================
double Electrode_CSTR::capacityLeftRaw(int platNum) const
{
    if (electrodeChemistryModelType_ == 2) {
        setCapacityCoeff_FeS2();
    }
    double capLeft = 0.0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph == solnPhase_ || iph ==  metalPhase_) {
            continue;
        }
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        for (size_t k = 0; k < NumSpeciesList_[iph]; k++) {
            capLeft += spMoles_final_[kStart + k] * capacityLeftSpeciesCoeff_[kStart + k];
        }
    }
    return capLeft * Faraday;
}
//====================================================================================================================
// Set the relative current capacity discharged per Mole

void Electrode_CSTR::setRelativeCapacityDischargedPerMole(double relDischargedPerMole, int platNum)
{
    /*
     *  Needs to be converted to the relative extent of reaction, which is in terms of mole fractions
     *  and has an upper and a lower limit.
     */
    if (RelativeExtentRxn_RegionBoundaries_.size() != 2) {
        throw Electrode_Error("Electrode_CSTR::setRelativeCapacityDischargedPerMole",
                           "confused");
    }

    double relExtentRxn =  RelativeExtentRxn_RegionBoundaries_[0] + relDischargedPerMole;

    if (platNum == 0 || platNum == -1) {
        setState_relativeExtentRxn(relExtentRxn);
    } else {
        throw Electrode_Error("Electrode_CSTR::setRelativeCapacityDischargedPerMole",
                           " wrong platNum");
    }
}
//====================================================================================================================
//  Extract the ROP of the multiple reaction fronts from Cantera within this routine
/*
 *  In this routine we calculate the rates of progress of reactions and species on all active reacting surfaces.
 *
 *        ROP_[jRxn]
 *        spNetProdPerArea_List_[isk][kIndexKin]
 */
void Electrode_CSTR::extractInfo()
{

    int maxNumRxns = RSD_List_[0]->nReactions();
    std::vector<double> netROP(maxNumRxns, 0.0);
    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    for (size_t isk = 0; isk < numSurfaces_; isk++) {
        double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
        std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
        /*
         *   Only loop over surfaces that have kinetics objects associated with them
         */
        if (ActiveKineticsSurf_[isk]) {
            /*
             *  Sometimes the ROP for goNowhere is nonzero, when it should be zero.
             */
            if (goNowhere_) {
                for (int i = 0; i < maxNumRxns; i++) {
                    ROP_[i] = 0.0;
                }
                continue;
            } else {
                /*
                 *  Get the species production rates for the reacting surface.
                 *  This is used to integrate cell model equations.
                 */
                getNetSurfaceProductionRates(isk, spNetProdPerArea);
                /*
                 *  Get reaction net rates of progress for the reacting surface.
                 *  This does not appear to be used (cmtenne 2012.06.06).
                 */
                RSD_List_[isk]->getNetRatesOfProgress(DATA_PTR(netROP));
                for (int i = 0; i < maxNumRxns; i++) {
                    ROP_[i] = netROP[i];
                }
            }
        }
    }
}
//====================================================================================================================
// Calculate the production rate of species in the electrode at the final time of the time step
/*
 *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
 *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
 *   all species in the electrode.
 */
void Electrode_CSTR::updateSpeciesMoleChangeFinal()
{
    double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(0);
    std::fill(DspMoles_final_.begin(), DspMoles_final_.end(), 0.0);
    double mult = (surfaceAreaRS_init_[0] + surfaceAreaRS_final_[0]) / 2.0;
    for (size_t i = 0; i < m_NumVolSpecies; i++) {
        DspMoles_final_[i] += mult * spNetProdPerArea[i];
    }

    // Also need to update DphMoles_final_ since it is used in the residual calculation.
    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        double DphaseMoles = 0.0;
        for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; ++sp) {
          size_t isp = globalSpeciesIndex(iph, sp);
          DphaseMoles += DspMoles_final_[isp];
        }
        DphMoles_final_[iph] = DphaseMoles;
    }

    /*
     * We define the srcdot for the Extent of reaction as the net electron loss out of the electrode. Thus,
     * we put a negative sign here.
     */
    if (electrodeType_ ==  ELECTRODETYPE_ANODE) {
        SrcDot_RxnExtent_final_  =   DspMoles_final_[kElectron_];
    } else if (electrodeType_ ==  ELECTRODETYPE_CATHODE) {
        SrcDot_RxnExtent_final_  = - DspMoles_final_[kElectron_];
    }
    double fac = capacityRaw() / (RelativeExtentRxn_NormalizationFactor_ * Faraday);
    fac = 1.0;
    double tmp = - capacityLeftDot() / Faraday / fac;
    SrcDot_RxnExtent_final_ = tmp;
   
    if (fabs(tmp -  SrcDot_RxnExtent_final_) / RelativeExtentRxn_NormalizationFactor_ > 1.0E-6) {
	printf("logic error %g %g\n", tmp ,  SrcDot_RxnExtent_final_);
        //exit(-1);
    }
}
//================================================================================================================
double Electrode_CSTR::capacityLeftDot(int platNum) const
{
    double capLeftDot = 0.0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        for (size_t k = 0; k < NumSpeciesList_[iph]; k++) {
            double ll =  DspMoles_final_[kStart + k];
            capLeftDot += ll * capacityLeftSpeciesCoeff_[kStart + k];
        }
    }
    return capLeftDot * Faraday;
}
//================================================================================================================
double Electrode_CSTR::capacityDot(int platNum) const
{
    double capDot = 0.0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph ==  solnPhase_ || iph == metalPhase_) {
            continue;
        }
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        for (size_t k = 0; k < NumSpeciesList_[iph]; k++) {
            double ll = DspMoles_final_[kStart + k];
            capDot += ll * capacityZeroDoDSpeciesCoeff_[kStart + k];
        }
    }
    return capDot * Faraday;
}
//================================================================================================================
double Electrode_CSTR::RxnExtentDot(int platNum) const
{
    double capLeft = capacityLeft(platNum, 50., -50.);
    double capZero = capacity(platNum);
    double capLeftDot = capacityLeftDot(platNum);
    double capDot = capacityDot(platNum);
    double reDot = - capLeftDot / capZero + capLeft * capDot / ( capZero * capZero);
    return reDot;
}
//================================================================================================================
void Electrode_CSTR::speciesProductionRates(double* const spMoleDot)
{     
    zeroD(m_NumTotSpecies, spMoleDot);
    //
    //  For non-pending we calculate the instantaneous value
    // 
    if (pendingIntegratedStep_ == 1) {
	printf("WARNING: speciesProductionRate called with pendingIntegratedStep_ = 1\n");
    }
  
    //
    // Look over active kinetics surfaces
    //
    for (size_t isk = 0; isk < numSurfaces_; isk++) {
	if (ActiveKineticsSurf_[isk]) {
	    /*
	     *  For each Reacting surface
	     *      (  m_rSurDomain->getNetProductionRates(&RSSpeciesProductionRates_[0]);
	     *  Get the species production rates for the reacting surface
	     */
	    // TODO: Check this logic for end of region conditions and goNowhere issues
	    const vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetSurfaceProductionRateDensities();
 
	    /*
	     *  loop over the phases in the reacting surface
	     *  Get the net production vector
	     */
	    size_t nphRS = RSD_List_[isk]->nPhases();
	    size_t jph, kph;
	    size_t kIndexKin = 0;
	    for (kph = 0; kph < nphRS; kph++) {
		jph = RSD_List_[isk]->kinOrder[kph];
		size_t istart = m_PhaseSpeciesStartIndex[jph];
		size_t nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
		for (size_t k = 0; k < nsp; k++) {
		    spMoleDot[istart + k] += rsSpeciesProductionRates[kIndexKin] * surfaceAreaRS_final_[isk];
		    kIndexKin++;
		}
	    }
	}
    }
}
//====================================================================================================================
void Electrode_CSTR::setState_relativeExtentRxn(double relExtentRxn)
{
    throw Electrode_Error("Electrode_CSTR::setState_relativeExtentRxn()",
                       "Base class called");
}
//====================================================================================================================
// Predict the solution
/*
 * Ok at this point we have a time step deltalimiTsubcycle_
 * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
 * We now calculate predicted solution components from these conditions.
 *
 */
int  Electrode_CSTR::predictSoln()
{

    double vleft, vright;
    double vtop, vbot;

    // predict that the calculated deltaT is equal to the input deltaT
    deltaTsubcycleCalc_ = deltaTsubcycle_;

    int redoSteps = 0;
    int reDo = 0;
    /* ---------------------------------------------------------------------------------------
     *                                    Loop over the predictor process
     * --------------------------------------------------------------------------------------- */
    do {
        redoSteps++;
        /*
         * Copy the initial state to the final state
         */
        std::fill(justBornPhase_.begin(), justBornPhase_.end(), 0);
        RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
        std::copy(spMf_init_.begin(), spMf_init_.end(), spMf_final_.begin());
        std::copy(spMoles_init_.begin(), spMoles_init_.end(), spMoles_final_.begin());
        std::copy(phaseMoles_init_.begin(), phaseMoles_init_.end(), phaseMoles_final_.begin());
        std::copy(surfaceAreaRS_init_.begin(), surfaceAreaRS_init_.end(), surfaceAreaRS_final_.begin());
        deltaTsubcycleCalc_ = deltaTsubcycle_;
	solidMoles_init_ = SolidTotalMoles();

        std::vector<double> deltaSpMoles(m_NumTotSpecies, 0.0);
        onRegionBoundary_final_ = onRegionBoundary_init_;
        /*
         * Get the top and bottom voltages
         */
        int nRegions = RelativeExtentRxn_RegionBoundaries_.size() - 1;
        onRegionBoundary_final_ = onRegionBoundary_init_;


	for (int itsP = 0; itsP < 2; itsP++) {
        /* ---------------------------------------------------------------------------------------
         *                                Update the Internal State of ThermoPhase Objects
         * --------------------------------------------------------------------------------------- */
        updateState();

        /* ---------------------------------------------------------------------------------------
         *                                 Handle the case where we start on a Region Boundary
         * --------------------------------------------------------------------------------------- */
        goNowhere_ = 0;
        if (onRegionBoundary_init_ >= 0) {

          // FIXME: this if block that sets vleft, vright, vtop, vbot will always have them overwritten
          // by the following block on lines 1026-1034. Should those lines be in an else clause?
            if (onRegionBoundary_init_ == 0) {
                vleft = openCircuitVoltage_Region(-1);
                vright = openCircuitVoltage_Region(0);
                if (electrodeType_ == ELECTRODETYPE_CATHODE) {
                    vtop = vleft;
                    vbot = vright;
                } else if (electrodeType_ == ELECTRODETYPE_ANODE) {
                    vtop = vright;
                    vbot = vleft;
                }
            }
            vleft  = openCircuitVoltage_Region(onRegionBoundary_init_ - 1);
            vright = openCircuitVoltage_Region(onRegionBoundary_init_);
            if (electrodeType_ == ELECTRODETYPE_CATHODE) {
                vtop = vleft;
                vbot = vright;
            } else { // electrodeType_ == ELECTRODETYPE_ANODE
                vtop = vright;
                vbot = vleft;
            }


            if (deltaVoltage_ > vtop) {
                if (electrodeType_ == ELECTRODETYPE_ANODE) {
                    if (onRegionBoundary_init_ < nRegions) {
                        xRegion_init_ = onRegionBoundary_init_;
                        xRegion_final_ = onRegionBoundary_init_;
                        onRegionBoundary_final_ = -1;
                        goNowhere_ = 0;
                    } else {
                        xRegion_init_ = nRegions - 1;
                        xRegion_final_ = nRegions - 1;
                        onRegionBoundary_final_ = nRegions;
                        goNowhere_ = 1;
                    }
                } else if (electrodeType_ == ELECTRODETYPE_CATHODE) {
                    if (onRegionBoundary_init_ > 0) {
                        xRegion_init_ = onRegionBoundary_init_ -1;
                        xRegion_final_ = onRegionBoundary_init_ -1;
                        onRegionBoundary_final_ = -1;
                        goNowhere_ = 0;
                    } else {
                        xRegion_init_ = 0;
                        xRegion_final_ = 0;
                        onRegionBoundary_final_ = 0;
                        goNowhere_ = 1;
                    }
                }
            } else if (deltaVoltage_ < vbot) {
                if (electrodeType_ == ELECTRODETYPE_ANODE) {
                    if (onRegionBoundary_init_ > 0) {
                        xRegion_init_ = onRegionBoundary_init_ - 1;
                        xRegion_final_ = onRegionBoundary_init_ - 1;
                        onRegionBoundary_final_ = -1;
                        goNowhere_ = 0;
                    } else {
                        xRegion_init_ = 0;
                        xRegion_final_ = 0;
                        onRegionBoundary_final_ = 0;
                        goNowhere_ = 1;
                    }
                } else if (electrodeType_ == ELECTRODETYPE_CATHODE) {
                    if (onRegionBoundary_init_ < nRegions) {
                        xRegion_init_ = onRegionBoundary_init_;
                        xRegion_final_ = onRegionBoundary_init_;
                        onRegionBoundary_final_ = -1;
                        goNowhere_ = 0;
                    } else {
                        xRegion_init_ = nRegions - 1;
                        xRegion_final_ = nRegions - 1;
                        onRegionBoundary_final_ = nRegions;
                        goNowhere_ = 1;
                    }
                }
            } else {
                xRegion_final_ = xRegion_init_;
                onRegionBoundary_final_ = onRegionBoundary_init_;
                goNowhere_ = 1;
            }
        }

        /* ---------------------------------------------------------------------------------------
         *                             Update the Phase Existance Flags
         * --------------------------------------------------------------------------------------- */
        /*
        * At the start of this step the _final_ values are equal to the _init_ values
        * - we may want to put a check here that this is true. However, for now we will assume it
        *
        * Set the phase existence -> this depends on SpMoles_init_ and SpMoles_final_
        */
        if (redoSteps == 1) {
            setPhaseExistenceForReactingSurfaces(false);
        } else {
            setPhaseExistenceForReactingSurfaces(false);
        }
        /* ---------------------------------------------------------------------------------------
         *                           Update the Rates of Progress and Derivative properties
         * --------------------------------------------------------------------------------------- */
        /*
         * Calculate ROP_inner and ROP_outer, and justBornPhase_[];
         */
        extractInfo();
        /*
         *  Calculate DspMoles_final_[]
         */
        updateSpeciesMoleChangeFinal();

        /* ---------------------------------------------------------------------------------------
         *                         Estimate the final total moles
         * --------------------------------------------------------------------------------------- */


        /* ---------------------------------------------------------------------------------------
         *                         Figure out the minimum death time for a phase to zero
         * --------------------------------------------------------------------------------------- */
        //    First assume that no phase dies.
        minPH_ = npos;
        double minDT = 1.0E300;

	//
	//  Then we assign minPH_ as the phase that dies first
	//  and minDT as the time at which that phase dies.
	//
        double DT, DTmin, tmp;
        int hasAPos, hasANeg;
        for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
            hasAPos = 0;
            hasANeg = 0;
            int iph = phaseIndexSolidPhases_[ph];
            justDiedPhase_[iph] = 0;
            DTmin = 1.0E+300;
            for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                size_t isp = globalSpeciesIndex(iph, sp);
                deltaSpMoles_[isp] = DspMoles_final_[isp] * deltaTsubcycle_;
                tmp = spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycle_;
                if (spMoles_init_[isp] > 0.0) {
                    if (tmp <= 0.0) {
                        hasANeg = 1;
                        DT = - spMoles_init_[isp] / DspMoles_final_[isp];
                        DTmin = MIN(DT, DTmin);
                    } else if (tmp > 0.0) {
                        hasAPos = 1;
                    }
                }
            }
            if (hasANeg && (!hasAPos)) {
                DT = DTmin;
                if (DT < minDT) {
                    minPH_ = iph;
                    minDT = DT;
                }
            }
        }
        /* ---------------------------------------------------------------------------------------
         *                         Figure out if we are going to run into a plateau boundary
         * --------------------------------------------------------------------------------------- */
        /*
         *   Here we figure out if we are going to run into a plateau boundary. If so, we reduce
         *   deltaTsubcycleCalc_ value to the estimated collision time.
         */
        deltaT_RegionBoundaryCollision_ = 1.0E300;
        double deltax = 0.0;
        if (goNowhere_) {
            RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
        } else {
            RelativeExtentRxn_final_ = RelativeExtentRxn_init_
                                       + SrcDot_RxnExtent_final_ * deltaTsubcycle_ / RelativeExtentRxn_NormalizationFactor_;

            deltax =  RelativeExtentRxn_final_ - RelativeExtentRxn_init_;

            double xBd = RelativeExtentRxn_RegionBoundaries_[xRegion_init_ + 1];
            /*
             *   Look for collisions with boundary regions
             */
            if ((SrcDot_RxnExtent_final_ > 1.0E-200) &&
                    (RelativeExtentRxn_final_ > RelativeExtentRxn_RegionBoundaries_[xRegion_init_ + 1])) {
                // Forward boundary collision -> decrease expected deltaT time.
                deltax = xBd - RelativeExtentRxn_init_;
                deltaT_RegionBoundaryCollision_ = deltax * RelativeExtentRxn_NormalizationFactor_ / SrcDot_RxnExtent_final_;
                onRegionBoundary_final_ = xRegion_init_ + 1;
                deltaTsubcycleCalc_ = deltaT_RegionBoundaryCollision_;
                RelativeExtentRxn_final_  = RelativeExtentRxn_RegionBoundaries_[xRegion_init_ + 1];

            } else  if ((SrcDot_RxnExtent_final_ < -1.0E-200) &&
                        (RelativeExtentRxn_final_ < RelativeExtentRxn_RegionBoundaries_[xRegion_init_])) {
                // backwards boundary collision -> decrease expected deltaT time.
                xBd = RelativeExtentRxn_RegionBoundaries_[xRegion_init_];
                deltax = xBd - RelativeExtentRxn_init_;
                deltaT_RegionBoundaryCollision_ = deltax * RelativeExtentRxn_NormalizationFactor_ / SrcDot_RxnExtent_final_;
                onRegionBoundary_final_ = xRegion_init_;
                deltaTsubcycleCalc_ = deltaT_RegionBoundaryCollision_;
                RelativeExtentRxn_final_  = RelativeExtentRxn_RegionBoundaries_[xRegion_init_];
            }
        }

        /* ---------------------------------------------------------------------------------------
         *                         Figure out what phases are going to die, now that we have a deltaTsubcycleCalc_
         * --------------------------------------------------------------------------------------- */
	//
	//  Here we calculate justDiedPhase_[] vector
	//
        if (minDT < 1.15 * deltaTsubcycleCalc_) {
            for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
                hasAPos = 0;
                hasANeg = 0;
                DTmin = 1.0E+300;
                int iph = phaseIndexSolidPhases_[ph];
                for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                    size_t isp = globalSpeciesIndex(iph, sp);
                    deltaSpMoles_[isp] = DspMoles_final_[isp] * deltaTsubcycleCalc_;
                    tmp = spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_;
                    if (spMoles_init_[isp] > 0.0) {
                        if (tmp <= 1.0E-40) {
                            hasANeg = 1;
                            DT = - spMoles_init_[isp] / DspMoles_final_[isp];
                            DTmin = MIN(DT, DTmin);
                        } else if (tmp > 0.0) {
                            hasAPos = 1;
                        }
                    }
                }
                if (hasANeg && (!hasAPos)) {
                    DT = DTmin;
                    if (DT < minDT + 0.1 *minDT) {
                        justDiedPhase_[iph] = 1;
                    }
                }
            }
	    //
	    // We redo the calculation of deltaTsubcycleCalc_
	    //
            if (minPH_ != npos) {
                deltaTsubcycleCalc_ = minDT;
                //   If we had previously found that a collision with a region boundary occurs,
                //   This is now negated if the time scales indicate that the phase disappearance isn't equal
                //   to the boundary collision. (not been exp)
                if (onRegionBoundary_final_ >= 0) {
                    if (fabs(deltaT_RegionBoundaryCollision_ - deltaTsubcycleCalc_) >  0.01 * deltaTsubcycleCalc_) {
                        onRegionBoundary_final_  = -1;
                    }
                }
            }
        }

        /* ---------------------------------------------------------------------------------------
         *                        Predict Soln
         * --------------------------------------------------------------------------------------- */
        for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
            size_t iph = phaseIndexSolidPhases_[ph];
            /*
             *  If we predict the phase will die, then we set the mole fractions equal to the
             *  previous step, we set the mole numbers to zero, and we predict deltaTsubcycleCalc_
             *  based on the ratio of the initial moles to the rate of mole destruction.
             */
            if (justDiedPhase_[iph]) {
                phaseMoles_final_[iph] = 0.0;
                for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                    size_t isp = globalSpeciesIndex(iph,sp);
                    spMf_final_[isp] = spMf_init_[isp];
		    spMoles_final_[isp] = 0.0;
                }
            } else {
                double sum = 0.0;
                for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                    size_t isp = globalSpeciesIndex(iph,sp);
                    spMoles_final_[isp] = spMoles_init_[isp] +  DspMoles_final_[isp] * deltaTsubcycleCalc_;
                    if (spMoles_final_[isp] < 0.0) {
                        spMoles_final_[isp] = 0.1 * spMoles_init_[isp];
                    }
                    sum +=  spMoles_final_[isp];
                }
                phaseMoles_final_[iph]= sum;
		if ( sum > 1.0E-300) {
		    for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
			size_t isp = globalSpeciesIndex(iph,sp);
			spMf_final_[isp] =  spMoles_final_[isp] /  phaseMoles_final_[iph];
		    }
		} else {
		    for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
			size_t isp = globalSpeciesIndex(iph, sp);
			spMf_final_[isp] = spMf_init_[isp];
		    }
		}
            }
        }

        /* ---------------------------------------------------------------------------------------
         *                        Predict the relative extent of reaction variable
         * --------------------------------------------------------------------------------------- */

        RelativeExtentRxn_final_ = RelativeExtentRxn_init_
                                   + SrcDot_RxnExtent_final_ * deltaTsubcycleCalc_ / RelativeExtentRxn_NormalizationFactor_;
	}

    } while (reDo);
    /* ---------------------------------------------------------------------------------------
     *                            Keep a copy of the estimated final total moles
     * --------------------------------------------------------------------------------------- */

    packNonlinSolnVector(DATA_PTR(soln_predict_));

    return 1;
}
//====================================================================================================================
// Unpack the soln vector
/*
 *  This function unpacks the solution vector into deltaTsubcycleCalc_ and RelativeExtentRxn_final_
 */
void  Electrode_CSTR::unpackNonlinSolnVector(const double* const y)
{
    size_t index = 0;
    deltaTsubcycleCalc_ = y[0];
    tfinal_ = tinit_ + deltaTsubcycleCalc_;
    index++;

    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        phaseMoles_final_[iph] = y[index];
        index++;
    }
    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        size_t isp = globalSpeciesIndex(iph, 0);
        size_t iStart = isp;
        if (numSpecInSolidPhases_[ph] == 1) {
            spMf_final_[isp] = 1.0;
            spMoles_final_[isp] = phaseMoles_final_[iph];
        }  else {
            double BigMF = 1.0;
            for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                if (sp != phaseMFBig_[iph]) {
                    spMf_final_[isp] = y[index];
                    spMoles_final_[isp] = spMf_final_[isp] * phaseMoles_final_[iph];
                    BigMF -= spMf_final_[isp];
                    index++;
                }
                isp++;
            }
            isp = iStart + phaseMFBig_[iph];
            spMf_final_[isp] = BigMF;
            spMoles_final_[isp] = spMf_final_[isp] * phaseMoles_final_[iph];
        }
    }
    /*
     * Calculate the relative extents of reaction
     */
    double tmp  = calcRelativeExtentRxn_final();
    RelativeExtentRxn_final_ = tmp;

    checkStillOnRegionBoundary();
}
//==================================================================================================================================
void Electrode_CSTR::checkStillOnRegionBoundary()
{
  const double upperBound = RelativeExtentRxn_RegionBoundaries_[xRegion_init_ + 1];
  const double lowerBound = RelativeExtentRxn_RegionBoundaries_[xRegion_init_];
  if (onRegionBoundary_final_ >= 0) {
     if (lowerBound < RelativeExtentRxn_final_ && RelativeExtentRxn_final_ < upperBound) {
        //printf("Caution: in valid region, even though onRegionBoundary_final_ = %d\n", onRegionBoundary_final_);
     }

  } 
  //if (lowerBound <= RelativeExtentRxn_final_ && RelativeExtentRxn_final_ <= upperBound ) {
    // onRegionBoundary_final_ = -1;
  //}
}
//==================================================================================================================================
void Electrode_CSTR::setOnRegionBoundary()
{
  const double upperBound = RelativeExtentRxn_RegionBoundaries_[xRegion_init_ + 1];
  const double lowerBound = RelativeExtentRxn_RegionBoundaries_[xRegion_init_];
  if (RelativeExtentRxn_final_ >= upperBound) {
    onRegionBoundary_final_ = xRegion_init_ + 1;
    setState_relativeExtentRxn(upperBound);
  } else if (RelativeExtentRxn_final_ <= lowerBound) {
    onRegionBoundary_final_ = xRegion_init_;
    setState_relativeExtentRxn(lowerBound);
  } else {
    onRegionBoundary_final_ = -1;
  }
}
//==================================================================================================================================
void  Electrode_CSTR::packNonlinSolnVector(double* const y) const
{

    size_t index = 0;
    y[0] = deltaTsubcycleCalc_;
    // tfinal_ = tinit_ + deltaTsubcycleCalc_;
    index++;

    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        y[index] = phaseMoles_final_[iph];
        index++;
    }
    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        size_t isp = globalSpeciesIndex(iph, 0);
        if (numSpecInSolidPhases_[ph] == 1) {
            //spMf_final_[isp] = 1.0;
            //spMoles_final_[isp] = phaseMoles_final_[iph];
        }  else {
            for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                if (sp != phaseMFBig_[iph]) {
                    y[index] = spMf_final_[isp];
                    index++;
                }
                isp++;
            }
        }
    }
    y[index] = RelativeExtentRxn_final_;
    y[index+1] = onRegionBoundary_final_;
}
//==================================================================================================================
//! Pack the nonlinear solver proplem
/*!
 *  formulate the nonlinear solver problem to be solved.
 *     Fields to be filled in
 *             yvalNLS_
 *             ylow
 *             yhigh
 *             ydeltaBoundsMagnitudes
 *             ysType
 *
 *
 *  Let's formulate the solution vector. The solution vector will consist of the following form
 *
 *     y =   deltaTsubcycleCalc_
 *           phaseMoles_final[iph]    for iph = phaseIndexSolidPhases[0]
 *           phaseMoles_final[iph]    for iph = phaseIndexSolidPhases[1]
 *           . . .
 *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[0]
 *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[0]
 *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[0]
 *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[0]
 *           ...
 *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[1]
 *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[1]
 *           ...
 */
void  Electrode_CSTR::initialPackSolver_nonlinFunction()
{
    // int nTop = (int) RegionBoundaries_ExtentRxn_.size() - 1;
    //double delBounds = 0.01 * fabs(RegionBoundaries_ExtentRxn_[nTop]);
    // yvalNLS_.clear();
    size_t index = 0;
    /*
     *  Set up the detaTsubcycle variable
     */
    yvalNLS_[0] = deltaTsubcycleCalc_;
    ylowNLS_[0] = 1.0E-40;
    yhighNLS_[0] = 1.0E300;
    deltaBoundsMagnitudesNLS_[0] = 1.0E300;
    index++;

    /*
     *  Calculate the atolVal that will be used in the calculation of phase moles.
     *  Note, from experience we cannot follow within the equil solver the phases with mole number that
     *  are additively insignificant compared to the total number of moles. This is a basic
     *  limitation. However, they must be followed kinetically up to the level that they contribute
     *  electrons. So, we will set atolBaseResid_ to 1.0E-12 with possible additions later.
     */
    //
    // Pick a solidMoles value that is representative of the average number of moles in the nonlinear solution
    // This number shouldn't ever go small, even if the moles in the electrode goes to zero
    //
    double solidMoles = SolidTotalMoles();
    solidMoles = std::max(solidMoles_init_, solidMoles);
    solidMoles = std::max(solidMoles,RelativeExtentRxn_NormalizationFactor_); 
    solidMoles = std::max(solidMoles, 1.0E-10);

    double atolVal = solidMoles * atolBaseResid_;

    /*
     * Add the unknowns for each phase moles
     */
    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        yvalNLS_[index] = phaseMoles_final_[iph];
        /*
         *  HKM Experimental bounds check on total number of moles in a phase
         *    -> This is looking like it is a bad idea to bounds check phaseMoles_final_[] values,
         *       a very bad idea.  However, it needs to be pursued. The idea is to move all births
         *       and deaths to the predictor. PURSUE!!!
         */
#ifdef STRICT_BOUNDS_NP
        if (justDiedPhase_[iph]) {
            ylowNLS_[index] = -0.1 * solidMoles;
        } else {
            ylowNLS_[index] = 0.0;
        }
#else
        /*
         * If we are birthing a phase, we decide to do this in the predictor. We search for
         * a solution where this is the case. If this doesn't happen, we produce a hard error.
         */
        if (justBornPhase_[iph]) {
            ylowNLS_[index] = atolVal;
        } else {
            ylowNLS_[index] = (-0.0001 * solidMoles);
        }
#endif
        yhighNLS_[index] = (10 * solidMoles);

        deltaBoundsMagnitudesNLS_[index] = solidMoles * 1.0E-6;
        index++;
    }

    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
            size_t iph = phaseIndexSolidPhases_[ph];
            size_t isp = globalSpeciesIndex(iph,sp);
            if (sp != phaseMFBig_[iph]) {
                yvalNLS_[index]  = spMf_final_[isp];
                ylowNLS_[index]  = 0.0;
                yhighNLS_[index] = 1.0;
                deltaBoundsMagnitudesNLS_[index] = 1.0E-16;
                index++;
            }
        }
    }
    if (static_cast<int>(index) != nEquations()) {
        printf("we have a prob\n");
        exit(-1);
    }
}
//====================================================================================================================
size_t Electrode_CSTR::nEquations_calc() const
{
    size_t neq = phaseIndexSolidPhases_.size();
    neq++;
    for (size_t i = 0; i < phaseIndexSolidPhases_.size(); i++) {
        if (numSpecInSolidPhases_[i] > 1) {
            neq += numSpecInSolidPhases_[i] - 1;
        }
    }
    return neq;
}
//====================================================================================================================
//  Return a vector of delta y's for calculation of the numerical Jacobian
/*
 *   There is a default algorithm provided.
 *
 *        delta_y[i] = atol[i] + 1.0E-6 ysoln[i]
 *        delta_y[i] = atol[i] + MAX(1.0E-6 ysoln[i] * 0.01 * solnWeights[i])
 *
 * @param t             Time                    (input)
 * @param y             Solution vector (input, do not modify)
 * @param ydot          Rate of change of solution vector. (input, do not modify)
 * @param delta_y       Value of the delta to be used in calculating the numerical jacobian
 * @param solnWeights   Value of the solution weights that are used in determining convergence (default = 0)
 *
 * @return Returns a flag to indicate that operation is successful.
 *            1  Means a successful operation
 *            0  Means an unsuccessful operation
 */
int Electrode_CSTR::calcDeltaSolnVariables(const double t, const double* const ySoln,
        const double* const ySolnDot, double* const deltaYSoln,
        const double* const solnWeights)
{
    if (!solnWeights) {
        for (size_t i = 0; i < static_cast<size_t>(neq_); i++) {
            deltaYSoln[i] = 1.0E-12 + fabs(1.0E-6 * ySoln[i]);
        }
    } else {
        for (size_t i = 0; i < static_cast<size_t>(neq_); i++) {
            deltaYSoln[i] = fmax(1.0E-2 * solnWeights[i], 1.0E-6 * fabs(ySoln[i]));
        }
    }
    return 1;
}
//====================================================================================================================
//   Determine the big mole fraction in the phase
void  Electrode_CSTR::determineBigMoleFractions()
{
    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        phaseMFBig_[iph] = 0;
        size_t iStart = globalSpeciesIndex(iph, 0);
        double xBig = spMf_final_[iStart];
        for (size_t sp = 1; sp < numSpecInSolidPhases_[ph]; sp++) {
            size_t isp = globalSpeciesIndex(iph, sp);
            if (spMf_final_[isp] > xBig) {
                phaseMFBig_[iph] = sp;
                xBig = spMf_final_[isp];
            }
        }
    }
}
//==================================================================================================================================
//  Residual calculation for the solution of the Nonlinear integration problem
/*
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
int Electrode_CSTR::integrateResid(const double t, const double delta_t, const double* const y, const double* const ySolnDot,
                                   double* const resid, const ResidEval_Type evalType, const int id_x, const double delta_x)
{
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printf("\t\t===============================================================================================================================\n");
        printf("\t\t  EXTRA PRINTING FROM NONLINEAR RESIDUAL: ");
        if (evalType ==  ResidEval_Type::Base_ResidEval) {
            printf(" BASE RESIDUAL");
        } else if (evalType == ResidEval_Type::JacBase_ResidEval) {
            printf(" BASE JAC RESIDUAL");
        } else  if (evalType == ResidEval_Type::JacDelta_ResidEval) {
            printf(" DELTA JAC RESIDUAL");
            printf(" var = %d delta_x = %12.4e Y_del = %12.4e Y_base = %12.4e", id_x, delta_x, y[id_x], y[id_x] - delta_x);
        } else  if (evalType == ResidEval_Type::Base_ShowSolution) {
            printf(" BASE RESIDUAL - SHOW SOLUTION");
        }
        printf(" DomainNumber = %2d , CellNumber = %2d , SubIntegrationCounter = %d\n",
               electrodeDomainNumber_, electrodeCellNumber_, counterNumberSubIntegrations_);
    }
    if ((evalType != ResidEval_Type::JacDelta_ResidEval)) {
        std::fill(justDied_.begin(), justDied_.end(), 0);
    }
    /*
     *  UNPACK THE SOLUTION VECTOR
     */
    unpackNonlinSolnVector(y);

    if (evalType != ResidEval_Type::JacDelta_ResidEval && (evalType != ResidEval_Type::Base_LaggedSolutionComponents)) {
        std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_final_lagged_.begin());
    }
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        if (minPH_ != npos) {
            if (deltaTsubcycleCalc_  <= deltaTsubcycle_) {
                printf("\t\t\tTfinal = %12.4e Tinit = %12.4e DeltaTsubcycle = %12.4e  DeltaTsubcycleCalc = %12.4e SPECIAL CASE DECREASED deltaT\n",
                       tfinal_, tinit_,deltaTsubcycle_, deltaTsubcycleCalc_);
            }
            if (deltaTsubcycleCalc_  > deltaTsubcycle_) {
                printf("\t\t\tQWARNING Tfinal = %12.4e Tinit = %12.4e DeltaTsubcycle = %12.4e  DeltaTsubcycleCalc = %12.4e SPECIAL CASE INCREASED deltaT\n",
                       tfinal_, tinit_,deltaTsubcycle_, deltaTsubcycleCalc_);
            }
        } else if (onRegionBoundary_final_ >= 0 && (onRegionBoundary_init_ != onRegionBoundary_final_)) {
            if (deltaTsubcycleCalc_  <= deltaTsubcycle_) {
                printf("\t\t\tTfinal = %12.4e Tinit = %12.4e DeltaTsubcycle = %12.4e  DeltaTsubcycleCalc = %12.4e SPECIAL CASE DECREASED deltaT\n",
                       tfinal_, tinit_,deltaTsubcycle_, deltaTsubcycleCalc_);
            }
            if (deltaTsubcycleCalc_  > deltaTsubcycle_) {
                printf("\t\t\tQWARNING Tfinal = %12.4e Tinit = %12.4e DeltaTsubcycle = %12.4e  DeltaTsubcycleCalc = %12.4e SPECIAL CASE INCREASED deltaT\n",
                       tfinal_, tinit_,deltaTsubcycle_, deltaTsubcycleCalc_);
            }
            printf("\t\t\t          Solving for Solution on Region Boundary %d\n", onRegionBoundary_final_);
        } else if (onRegionBoundary_final_ >= 0 && onRegionBoundary_init_ >= 0) {

            printf("\t\t\tTfinal = %12.4e Tinit = %12.4e DeltaTsubcycle = %12.4e  DeltaTsubcycleCalc = %12.4e \n",
                   tfinal_, tinit_,deltaTsubcycle_, deltaTsubcycleCalc_);
            printf("\t\t\t          Solving for GoNoWhere, boundary %d\n", onRegionBoundary_final_);

        } else {
            printf("\t\t\tTfinal = %12.4e Tinit = %12.4e DeltaTsubcycle = %12.4e\n", tfinal_, tinit_, deltaTsubcycle_);
        }
        printf("\t\t\tDeltaVoltage = %12.4e\n", deltaVoltage_);
    }

REDO:
    updateState();
    extractInfo();
    /*
     *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *   all species in the electrode.
     */
    updateSpeciesMoleChangeFinal();
    /*
     * Calculate the residual
     */
    calcResid(resid, evalType);
    /*
     * Change the problem specification for the nonlinear solve in certain circumstances when the solver
     *  is calculating the base residual of a jacobian
     */
    if (evalType == ResidEval_Type::JacBase_ResidEval) {
#ifdef DEBUG_HKM
        if (deltaTsubcycleCalc_ < 0.0) {
            printf("we are here deltTsybcycleCalc_ = %g\n",deltaTsubcycleCalc_);
            enableExtraPrinting_ = 1;
            detailedResidPrintFlag_  = 14;
            setPrintLevel(14);
            pSolve_->setPrintLvl(14);
        }
#endif
        if (deltaTsubcycleCalc_ > 1.2 * deltaTsubcycle_) {
            if (onRegionBoundary_final_ >= 0 && minPH_ == npos) {
                if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
                    printf("\t\t PROBLEM CHANGE: Stretching too far to find a boundary. Getting rid of boundary specification\n");
                    printf("\t\t 	               Current Calc deltaT = %g, nominal deltaT = %g\n", deltaTsubcycleCalc_, deltaTsubcycle_);
                    printf("\t\t                 RER_final = %g, RER_init = %g, RER_bd = %g",  RelativeExtentRxn_final_,
                           RelativeExtentRxn_init_,  RelativeExtentRxn_RegionBoundaries_[onRegionBoundary_final_]);
                }
                onRegionBoundary_final_ = -1;
                goto REDO ;
            }
        }
        if (deltaTsubcycleCalc_ > 1.5 * deltaTsubcycle_) {
            if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
                printf("\t\tDeltaTCalc is too big. %g %g BAILING FROM RESIDUAL\n", deltaTsubcycleCalc_, deltaTsubcycle_);
            }
            return -2;
        }
        for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
            size_t iph = phaseIndexSolidPhases_[ph];
            if (justDied_[iph] != justDiedPhase_[iph]) {
                size_t index = 1 + ph;
                std::vector<double>& ylow = pSolve_->lowBoundsConstraintVector();
                if (justDied_[iph] == 1) {
#ifdef DEBUG_ELECTRODE
                    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
                        printf("\t\tPhase %d has just died according to the residual\n", iph);
                    }
#endif
                    justDiedPhase_[iph] = 1;
                    if (ylow[index] >= 0.0) {
                        ylow[index] = -1.0E30;
                    }
                } else {
                    // HKM  -> We have turned off this code, because it leads to multiple errors.
                    //         If a phase has bee slated to be killed off, you need to let the nonlinear solver try to kill
                    //         it off.
                    // if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
                    // printf("\t\tPhase %d is alive even though it is assigned to die\n", iph);
                    // }
                    // justDiedPhase_[iph] = 0;
#ifdef STRICT_BOUNDS_NP
                    if (ylow[index] < 0.0) {
                        ylow[index] = 0.0;
                    }
#endif
                    if (DphMoles_final_[iph] > 0.0) {
                        if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
                            printf("\t\tPhase %d is growing even though it is assigned to die. BAILING FROM RESIDUAL\n", 
                                   static_cast<int>(iph));
                        }
                        return -3;

                    }
                }
            }
        }
    }

    int index = 1;
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {

        if (minPH_ != npos) {
            printf("\t\t minPH = %3d  delTcalc = %16.7E  Pinit/-DelP = %16.7E  Res = %16.7E   \n",
                   static_cast<int>(minPH_),  deltaTsubcycleCalc_, -phaseMoles_init_[minPH_] / DphMoles_final_[minPH_],  resid[0]);
        }
        printf("\t\t        PhaseName        Moles_Init    Moles_final     KExists  |   Src_Moles  Pred_Moles_Final  |    Resid     |\n");
        for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
            double src =   DphMoles_final_[iph] * deltaTsubcycleCalc_;
            printf("\t\t %20.20s  %12.4e  %12.4e            %2d | %12.4e %12.4e", PhaseNames_[iph].c_str(), phaseMoles_init_[iph],
                   phaseMoles_final_[iph], 1, src,  phaseMoles_init_[iph] + src);
            bool found = false;
            for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
                size_t jph = phaseIndexSolidPhases_[ph];
                if (jph == iph) {
                    found = true;
                    index++;
                    double res =  phaseMoles_final_[iph] - (phaseMoles_init_[iph] + src);
                    printf("      |%12.4e  | ",  res);
                }
            }
            if (!found) {
                printf("      |              | ");
            }
            if (justDiedPhase_[iph]) {
                if (justDiedPhase_[2] && justDiedPhase_[3]) {
                    printf("we are here 2 & 3 - how?\n");
                }
#ifdef DEBUG_ELECTRODE
                if (justDiedPhase_[3]) {
                    printf("we are here - 3 death\n");
                }
#endif
                printf("   --- PHASE DEATH ---");
                if (!justDied_[iph]) {
                    printf(" (with resurrection)");
                }
            } else {
                if (justDied_[iph]) {
                    printf("   --- JUST DIED ---");
                }
            }
            printf("\n");
        }
        printf("\t\t      SpeciesName        Moles_Init    Moles_final   SpMF       |   Src_Moles  Pred_Moles_Final\n");
        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            string ss = speciesName(k);
            double src =  DspMoles_final_[k] * deltaTsubcycleCalc_;
            printf("\t\t %20s  %12.4e  %12.4e  %12.4e | %12.4e %12.4e", ss.c_str(), spMoles_init_[k], spMoles_final_[k],
                   spMf_final_[k], src, spMoles_init_[k] + src);
            bool found = false;
            for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
                size_t iph = phaseIndexSolidPhases_[ph];
                if (numSpecInSolidPhases_[ph] > 1) {
                    size_t iStart = globalSpeciesIndex(iph, 0);
                    for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                        if (sp != phaseMFBig_[iph]) {
                            size_t i = iStart + sp;
                            if (i == k) {
                                found = true;
                                double res =  spMoles_final_[k] - (spMoles_init_[k] + src);
                                printf("      |%12.4e  | ", res);
                            }
                        }
                    }
                }
            }
            if (!found) {
                printf("      |              | ");
            }
            printf("\n");
        }
    }

    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printf("\t\t===============================================================================================================================\n");
    }
    return 1;
}
//====================================================================================================================
// Main internal routine to calculate the rate constant
/*
 *  This routine calculates the functional at the current stepsize, deltaTsubcycle_.
 *  A new stepsize, deltaTsubcycleCalc_, is calculated within this routine for changes
 *  in topology.
 *
 *  This routine calcules yval_retn, which is the calculated value of the residual for the
 *  nonlinear function
 *
 *   resid[i] = y[i] - yval_retn[i]
 *
 *  The formulation of the solution vector is as follows. The solution vector will consist of the following form
 *
 *     y =   phaseMoles_final[iph]    for iph = phaseIndexSolidPhases[0]
 *           phaseMoles_final[iph]    for iph = phaseIndexSolidPhases[1]
 *           . . .
 *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[0]
 *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[0]
 *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[0]
 *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[0]
 *           ...
 *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[1]
 *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[1]
 *
 */
int Electrode_CSTR::calcResid(double* const resid, const ResidEval_Type evalType)
{

    int index = 0;

    size_t sp;
    std::vector<double> phaseMoles_pos_(m_NumTotPhases, 0.0);
    // The first residual is the size of the time step
    resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
    double minDT = deltaTsubcycle_;
    double minDTpossible = deltaTsubcycle_;
    double tmp;
    double deltaTdeath = deltaTsubcycle_;
    if (minPH_ != npos) {
        if (deltaTsubcycleCalc_ >  deltaTsubcycle_ * 0.8333) {
            deltaTdeath =  deltaTsubcycleCalc_  * 1.2;
        }
    }


    /*
     *  Find a solid phase that may disappear. Mark the phase that disappears the fastest with the
     *  variable minPH_. Calculate the deltaT that zeroes the phase by the variable minDT.
     *
     * HKM -> CONSIDER WHACKING THIS SECTION - chose minPH_ at the predictor level and don't let it change?
     */
    size_t minPHtmp = minPH_;
    if (evalType != ResidEval_Type::JacDelta_ResidEval) {

        for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
            size_t iph = phaseIndexSolidPhases_[ph];

            if (phaseMoles_init_[iph] > 0.0) {
                if ((DphMoles_final_[iph] < 0.0) &&
                        (phaseMoles_init_[iph] / (-DphMoles_final_[iph]) < deltaTdeath)) {
                    bool allNeg = true;
                    if (numSpecInSolidPhases_[ph] > 1) {
                        for (sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                            size_t isp = globalSpeciesIndex(iph, sp);
                            tmp = spMoles_init_[isp] + DspMoles_final_[isp] * deltaTdeath;
                            if (spMoles_init_[isp] > 0.0 ||  DspMoles_final_[isp] > 0.0) {
                                if (tmp > 0.0) {
                                    allNeg = false;
                                    break;
                                }
                            }

                        }
                    }
                    // We are here if all predicted mole numbers in the phase are less than zero Or
                    // if we are solving modified equations for a phase death. In that case, if the predicted
                    // moles are greater than zero then that's ok.
                    if (allNeg || justDiedPhase_[iph]) {
                        // Find the predicted time that the phase goes away
                        minDTpossible = phaseMoles_init_[iph]/(-DphMoles_final_[iph]);
                        // We don't want to signal that the phase is dead, until the actual solution uknown is three orders
                        // of magnitude lower that the initial value and the deltaT for its death is close to the deltaT
                        // of the step -> previously without these limitations we had false positives.
                        if (minDTpossible < 1.2 * deltaTsubcycleCalc_ &&  minDTpossible < 1.2 * minDT) {
                            if (phaseMoles_final_[iph] < 1.0E-6 * phaseMoles_init_[iph] || phaseMoles_init_[iph] < 1.0E-22) {
                                justDied_[iph] = 1;
                            }
                        }

                        // solve for the smallest deltaT
                        if (minDTpossible < minDT) {
                            if (justDiedPhase_[iph] || justDied_[iph]) {
                                minDT = minDTpossible;
                                minPHtmp = iph;
                            }
                        }
                    }
                }
            }
        }

        if (evalType == ResidEval_Type::JacBase_ResidEval) {
            if (minPH_ == npos && minPH_ != minPHtmp) {
                minPH_ = minPHtmp;
            }
        }

    }

    double xBd = RelativeExtentRxn_RegionBoundaries_[xRegion_init_ + 1];
    double deltax = 1.0E300;
    if (onRegionBoundary_final_ >= 0) {
        if (SrcDot_RxnExtent_final_ > 1.0E-200) {
            deltax = xBd - RelativeExtentRxn_init_;
            deltaT_RegionBoundaryCollision_ = deltax * RelativeExtentRxn_NormalizationFactor_ / SrcDot_RxnExtent_final_;
        } else  if (SrcDot_RxnExtent_final_ < -1.0E-200) {
            xBd = RelativeExtentRxn_RegionBoundaries_[xRegion_init_ ];
            deltax = xBd - RelativeExtentRxn_init_;
            deltaT_RegionBoundaryCollision_ = deltax * RelativeExtentRxn_NormalizationFactor_ / SrcDot_RxnExtent_final_;
        } else {
            deltax = xBd - RelativeExtentRxn_init_;
            deltaT_RegionBoundaryCollision_ = 1000. * deltaTsubcycle_;
        }
    }

    /*
     *  We are no longer guarranteed to have deltaTsubcycleCalc_ less than or equal to deltaTsubcycle_
     *  That is OK.
     */
    if (goNowhere_) {
        resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
    } else {
        if (onRegionBoundary_final_ >= 0) {
            resid[0] = deltaTsubcycleCalc_ - deltax * RelativeExtentRxn_NormalizationFactor_ / SrcDot_RxnExtent_final_;
        } else if (minPH_ != npos) {
            resid[0] = deltaTsubcycleCalc_ + phaseMoles_init_[minPH_] / DphMoles_final_[minPH_];
        } else {
            resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
        }
    }
    index++;

    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        if (numSpecInSolidPhases_[ph] > 1) {
            phaseMoles_pos_[iph] = 0.0;
            for (sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                size_t isp = globalSpeciesIndex(iph, sp);
                tmp = spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_;
                if (tmp > 0.0) {
                    phaseMoles_pos_[iph] += tmp;
                } else {
                    phaseMoles_pos_[iph] += phaseMoles_final_[iph] * spMf_final_[isp];
                }
            }
        }
    }

    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        double tmp3 = phaseMoles_init_[iph] + DphMoles_final_[iph] * deltaTsubcycleCalc_;
        resid[index] = phaseMoles_final_[iph] - tmp3;
        index++;
    }


    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        if (numSpecInSolidPhases_[ph] > 1) {
            for (sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                size_t isp = globalSpeciesIndex(iph,sp);
                double frac;
                if (DphMoles_final_[iph] == 0.0 && phaseMoles_init_[iph] == 0.0) {
                    frac = spMf_init_[isp];
                } else if ((minPH_ == iph) || (justDiedPhase_[iph])) {
                    //  This approach to killing off multispecies phases is robust and stable.
                    frac = DspMoles_final_[isp] / DphMoles_final_[iph];
                } else {
                    if (phaseMoles_final_lagged_[iph] > 0.0) {
                        // HKM -> Note all of these approaches were tried and failed. Only the lagged approach seemed to be robust.
                        // frac = (spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_) / phaseMoles_pos_[iph];
                        // frac = (spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_) / phaseMoles_final_[iph];
                        // frac = (spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_) / phaseMoles_tmp_[iph];
                        frac = (spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_) / phaseMoles_final_lagged_[iph];
                    } else if (phaseMoles_pos_[iph] > 0.0) {
                        frac = (spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_) / phaseMoles_pos_[iph];
                    } else {
                        double ptmp = 0.0;
                        for (size_t ssp = 0; ssp < numSpecInSolidPhases_[ph]; ssp++) {
                            size_t iisp = globalSpeciesIndex(iph, ssp);
                            tmp = spMoles_init_[iisp] + DspMoles_final_[iisp] * deltaTsubcycleCalc_;
                            if (tmp > 0.0) {
                                ptmp += tmp;
                            }
                        }
                        if (ptmp > 0.0) {
                            frac = (spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_) / ptmp;
                        } else {
                            if (phaseMoles_init_[iph] > 0.0) {
                                frac = (spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_) / phaseMoles_init_[iph];
                            } else {
                                if (spMoles_init_[isp] > 0.0) {
                                    frac = (spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_) / (4.0 * spMoles_init_[isp]);
                                } else {
                                    frac = spMf_init_[isp];
                                }
                            }
                        }
                    }
                }
                if (sp != phaseMFBig_[iph]) {
                    resid[index] = spMf_final_[isp] - frac;
#ifdef DEBUG_MODE
                    checkFinite(resid[index]);
#endif
                    index++;
                }
            }
        }
    }

    return 1;
}
//==================================================================================================================================
int Electrode_CSTR::GFCEO_evalResidNJ(const double tdummy, const double delta_t_dummy,
                                const double* const y, const double* const ySolnDot,
                                double* const resid, const ResidEval_Type evalType,
                                const int id_x, const double delta_x)
{

   return 0;
}
//==================================================================================================================
// Main internal routine to calculate the rate constant
/*
 *  This routine calculates the functional at the current stepsize, deltaTsubcycle_.
 *  A new stepsize, deltaTsubcycleCalc_, is calculated within this routine for changes
 *  in topology.
 *
 *  This routine calcules yval_retn, which is the calculated value of the residual for the
 *  nonlinear function
 *
 *   resid[i] = y[i] - yval_retn[i]
 *
 *  The formulation of the solution vector is as follows. The solution vector will consist of the following form
 *
 *     y =   phaseMoles_final[iph]    for iph = phaseIndexSolidPhases[0]
 *           phaseMoles_final[iph]    for iph = phaseIndexSolidPhases[1]
 *           . . .
 *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[0]
 *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[0]
 *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[0]
 *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[0]
 *           ...
 *           Xmol[k = 0]              for iph = phaseIndexSolidPhases[1]
 *           Xmol[k = nSpecies()-1]   for iph = phaseIndexSolidPhases[1]
 *
 *    @param yval_retn calculated return vector whose form is described above
 *
 *
 *  @return  1 Means a good calculation that produces a valid result
 *           0 Bad calculation that means that the current nonlinear iteration should be terminated
 */
int Electrode_CSTR::GFCEO_calcResid(double* const resid, const ResidEval_Type evalType)
{
    size_t index = 0;
    // This has to be scaled
    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
	resid[index] = phaseMoles_dot_[iph] - DphMoles_final_[iph];
        index++;
    }
    return 1;
}
//==================================================================================================================
//   Set the Residual absolute error tolerances
/*
 *  (virtual from Electrode_Integrator)
 *
 *   Set the absolute error tolerances for the residuals for the nonlinear solvers. This is called at the top
 *   of the integrator() routine.
 *
 *   Calculates residAtolNLS_[]
 *   Calculates atolNLS_[]
 */
void  Electrode_CSTR::setResidAtolNLS()
{
    double deltaT = t_final_final_ - t_init_init_;
    deltaT = MAX(deltaT, 1.0E-3);
    double solidMoles = SolidTotalMoles();
    /*
     *  We don't care about differences that are 1E-6 less than the global time constant.
     *  residual has the same units as soln.
     */
    atolResidNLS_[0] = deltaT * 1.0E-6;
    /*
     *  Calculate the atolVal that will be used in the calculation of phase moles.
     *  Note, from experience we cannot follow within the equil solver the phases with mole number that
     *  are additively insignificant compared to the total number of moles. This is a basic
     *  limitation. However, they must be followed kinetically up to the level that they contribute
     *  electrons. So, we will set atolBaseResid_ to 1.0E-25 with possible additions later.
     */
    double atolVal = solidMoles * atolBaseResid_;

    atolNLS_.clear();
    atolNLS_.push_back(1.0E-50);
    int index = 0;
    index++;

    for (int ph = 0; ph < (int) phaseIndexSolidPhases_.size(); ph++) {
        // int iph = phaseIndexSolidPhases_[ph];
        atolNLS_.push_back(atolVal);
        atolResidNLS_[index] = atolNLS_[index];
        index++;
    }

    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {

        for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
            size_t iph = phaseIndexSolidPhases_[ph];
            //int isp = globalSpeciesIndex(iph,sp);
            if (sp != phaseMFBig_[iph]) {
                atolNLS_.push_back(1.0E-14);
                atolResidNLS_[index] = atolNLS_[index];
                index++;
            }
        }

    }
    if (atolNLS_.size() !=  neq_) {
        printf("ERROR\n");
        exit(-1);
    }

    determineBigMoleFractions();

}
//====================================================================================================================

void  Electrode_CSTR::printElectrodeCapacityInfo(int pSrc, bool subTimeStep)
{

    double capacd = capacityDischarged();
    printf("          Capacity Discharged Since Start = %12.6g coulombs = %12.6g Ah\n", capacd, capacd / 3600.);
    double dod = depthOfDischarge();
    printf("          Depth of Discharge (Current)    = %12.6g coulombs\n", dod);
    double capLeft = capacityLeft();
    double capZero = capacity();
    printf("          Capacity Left                   = %12.6g coulombs         Capacity at Zero DOD = %12.6g coulombs\n",
           capLeft, capZero);
    if (subTimeStep) {
        printf("          RelativeExtentRxn_final              = %g",      RelativeExtentRxn_final_);
        if (onRegionBoundary_final_ >= 0) {
            printf("        onRegionBoundaryFinal = %d\n",onRegionBoundary_final_);
        }
        printf("\n");
        printf("          RelativeExtentRxn_init               = %g", RelativeExtentRxn_init_);
        if (onRegionBoundary_init_ >= 0) {
            printf("        onRegionBoundaryInit = %d\n", onRegionBoundary_init_);
        }
        printf("\n");
        printf("                                 xRegion_init_ = %d xRegion_final = %d\n", xRegion_init_, xRegion_final_);
        if (goNowhere_) {
            printf("          goNowhere                            = true\n");
        }
    } else {
        printf("          RelativeExtentRxn_final_final         = %g",      RelativeExtentRxn_final_final_);
        if (onRegionBoundary_final_final_ >= 0) {
            printf("        onRegionBoundary_final_final = %d\n", onRegionBoundary_final_final_);
        }
        printf("\n");
        printf("          RelativeExtentRxn_init_init           = %g",      RelativeExtentRxn_init_init_);
        if (onRegionBoundary_init_init_ >= 0) {
            printf("        onRegionBoundary_init_init = %d\n", onRegionBoundary_init_init_);
        }
        printf("\n");
        printf("                             xRegion_init_init_ = %d xRegion_final_final = %d\n", xRegion_init_init_,
               xRegion_final_final_);
        if (goNowhere_) {
            printf("          goNowhere (at End)                    = true\n");
        }
    }
    printf("          RelativeExtentRxn_NormalizationFactor = %g kmol\n", RelativeExtentRxn_NormalizationFactor_);
    double tmp = RelativeExtentRxn_NormalizationFactor_ * Faraday;
    printf("                                                = %g coulombs\n", tmp);
    printf("          RelativeExtentRxn_bd[0]               = %g\n",     RelativeExtentRxn_RegionBoundaries_[0]);
    printf("          RelativeExtentRxn_bd[1]               = %g\n",    RelativeExtentRxn_RegionBoundaries_[1]);



}
//===================================================================================================================
void Electrode_CSTR::printElectrodePhase(size_t iph, int pSrc, bool subTimeStep)
{
    size_t isph = npos;
    double* netROP = new double[m_NumTotSpecies];
    ThermoPhase& tp = thermo(iph);
    std::string pname = tp.id();
    size_t istart = m_PhaseSpeciesStartIndex[iph];
    size_t nsp = tp.nSpecies();
    if (printLvl_ <= 1) {
        return;
    }
    printf("     ============================================================================================\n");
    printf("          PHASE %d %s \n", static_cast<int>(iph), pname.c_str());
    printf("                Total moles = %g\n", phaseMoles_final_[iph]);
    double mv = tp.molarVolume();
    if (iph < m_NumVolPhases) {
        printf("                Molar Volume = %11.5E cm3 gmol-1\n", mv * 1.0E3);
    } else {
        double ma = tp.molarArea();
        printf("                Molar Volume = %11.5E cm3 gmol-1    Molar Area = %11.5E cm2 gmol-1\n", 
               mv * 1.0E3, ma * 10.);
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
        printf("                Electric Potential = %g\n", tp.electricPotential());
    }
    /*
     * Do specific surface phase printouts
     */
    if (iph >= m_NumVolPhases) {
        isph = iph - m_NumVolPhases;
        printf("                surface area (final) = %11.5E m2\n",  surfaceAreaRS_final_[isph]);
        printf("                surface area (init)  = %11.5E m2\n",  surfaceAreaRS_init_[isph]);
        int ttt = isExternalSurface_[isph];
        printf("                IsExternalSurface = %d\n", ttt);
        double oc = openCircuitVoltage(isph);
        if (oc != 0.0) {
            printf("                 Open Circuit Voltage = %g\n", oc);
        }
    }
    if (printLvl_ >= 3) {
        printf("\n");
        printf("                Name              MoleFrac_final kMoles_final kMoles_init SrcTermIntegrated(kmol)\n");
        for (size_t k = 0; k < nsp; k++) {
            std::string sname = tp.speciesName(k);
            if (pSrc) {
                if (subTimeStep) {
                    printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                           spMoles_final_[istart + k], spMoles_init_[istart + k],
                           spMoleIntegratedSourceTermLast_[istart + k]);
                } else {
                    printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                           spMoles_final_[istart + k], spMoles_init_init_[istart + k],
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
        if (iph >= m_NumVolPhases) {
            const std::vector<double>& rsSpeciesProductionRates = RSD_List_[isph]->calcNetSurfaceProductionRateDensities();
            RSD_List_[isph]->getNetRatesOfProgress(netROP);

            double* spNetProdPerArea = (double*) spNetProdPerArea_List_.ptrColumn(isph);
            std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
            if (!goNowhere_) {
                size_t nphRS = RSD_List_[isph]->nPhases();
                size_t kIndexKin = 0;
                for (size_t kph = 0; kph < nphRS; kph++) {
                    size_t jph = RSD_List_[isph]->kinOrder[kph];
                    size_t istart = m_PhaseSpeciesStartIndex[jph];
                    size_t nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
                    for (size_t k = 0; k < nsp; k++) {
                        spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                        kIndexKin++;
                    }
                }
            }
            printf("\n");
            printf("                           spName                  SourceRateLastStep (kmol/m2/s) \n");
            for (size_t k = 0; k < m_NumTotSpecies; k++) {
                string ss = speciesName(k);
                printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
            }
        }
    }
    printf("     ============================================================================================\n");
    delete [] netROP;
}
//====================================================================================================================
// Evaluate the residual function
/*
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
int Electrode_CSTR::evalResidNJ(const double tdummy, const double delta_t_dummy,
                                const double* const y,
                                const double* const ySolnDot,
                                double* const resid,
                                const ResidEval_Type evalType,
                                const int id_x,
                                const double delta_x)
{
    int retn =  integrateResid(tdummy, delta_t_dummy, y, ySolnDot, resid, evalType, id_x, delta_x);
    return retn;
}
//====================================================================================================================
int Electrode_CSTR::getInitialConditions(const double t0, double* const ySoln, double* const ySolnDot)
{
    zeroD(neq_, ySoln);
    return 1;
}
//==============================================================================================================
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
void Electrode_CSTR::resetStartingCondition(double Tinitial, bool doTestsAlways)
{
    bool resetToInitInit = false;
    /*
     * If the initial time input from the parameter list, Tinitial, is the same as the current initial time,
     * Then, we don't advance the time step.
     */
    double tbase = MAX(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-13 * tbase) && !doTestsAlways) {
        resetToInitInit = true;
    }
    /*
     *  Call the base class resetStarting condition
     */
    Electrode_Integrator::resetStartingCondition(Tinitial);

    // Copy The final Extent of reaction to the beginning extent
    if (!resetToInitInit) {
    RelativeExtentRxn_init_init_ =   RelativeExtentRxn_final_final_;
    RelativeExtentRxn_init_      =   RelativeExtentRxn_final_final_;

    xRegion_init_init_ = xRegion_final_final_;
    xRegion_init_      = xRegion_final_final_;

    onRegionBoundary_init_init_  = onRegionBoundary_final_final_;
    onRegionBoundary_init_       = onRegionBoundary_final_final_;

    if (eState_final_) {
        SAFE_DELETE(xmlStateData_init_init_);
        xmlStateData_init_init_ =   xmlStateData_final_;
        SAFE_DELETE(xmlStateData_init_);
        xmlStateData_init_ = new XML_Node(*xmlStateData_final_);
        xmlStateData_final_ = 0;
        SAFE_DELETE(xmlStateData_final_final_);
    }
    }
}
//====================================================================================================================
// Check to see that the preceding step is a successful one
/*
 *   We check to see if the preceding step is a successful one.
 *
 *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
 */
bool  Electrode_CSTR::checkSubIntegrationStepAcceptable() const
{
    //
    // If we have solved the alternate problem where we find the time for phase death, then we expect to be at the
    // phase death condition. This is a valid throw I think.
    //
    if (onRegionBoundary_final_ >= 0) {
        if (fabs(RelativeExtentRxn_final_  - RelativeExtentRxn_RegionBoundaries_[onRegionBoundary_final_]) > 1.0E-5) {
            throw Electrode_Error("Electrode_CSTR::integrate() ERROR: cell " + int2str(electrodeCellNumber_) +
                               ", integrationCounter " + int2str(counterNumberIntegrations_),
                               " RelativeExtentRxn_final_ not on boundary and onRegionBoundary_final_ is. RelativeExtentRxn_final_ = "
                               + fp2str(RelativeExtentRxn_final_));
        }
    }

    if (deltaTsubcycle_  <= 0.0) {
        throw Electrode_Error("Electrode_CSTR::integrate() ERROR: cell " + int2str(electrodeCellNumber_) +
                           ", integrationCounter " + int2str(counterNumberIntegrations_),
                           "negative deltaTsubcycle_ = " + fp2str(deltaTsubcycle_));
    }
    if (deltaTsubcycleCalc_  <= 0.0) {
        throw Electrode_Error("Electrode_CSTR::integrate() ERROR: cell " + int2str(electrodeCellNumber_) +
                           ", integrationCounter " + int2str(counterNumberIntegrations_),
                           "negative deltaTsubcycleCalc_ = " + fp2str(deltaTsubcycleCalc_));
    }
    return true;
}
//====================================================================================================================
// Possibly change the solution due to phase births and deaths.
/*
 *   (virtual from Electrode_Integrator)
 *
 *  @return  Returns true if the solution step is good. It returns false if there is a problem.
 */
bool Electrode_CSTR::changeSolnForBirthDeaths()
{
    if (onRegionBoundary_final_ >= 0) {
        if (fabs(RelativeExtentRxn_final_  - RelativeExtentRxn_RegionBoundaries_[onRegionBoundary_final_]) > 1.0E-5) {
            throw Electrode_Error("Electrode_CSTR::integrate() ERROR: cell " + int2str(electrodeCellNumber_) +
                               ", integrationCounter " + int2str(counterNumberIntegrations_),
                               " RelativeExtentRxn_final_ not on boundary and onRegionBoundary_final_ is");
        } else {
            RelativeExtentRxn_final_  = RelativeExtentRxn_RegionBoundaries_[onRegionBoundary_final_];
            setState_relativeExtentRxn(RelativeExtentRxn_final_);
        }
    }

    /*
     *  Calculate the atolVal that will be used in the calculation of phase moles.
     *  Note, from experience we cannot follow within the equil solver the phases with mole number that
     *  are additively insignificant compared to the total number of moles. This is a basic
     *  limitation. However, they must be followed kinetically up to the level that they contribute
     *  electrons. So, we will set atolBaseResid_ to 1.0E-25 with possible additions later.
     */
    double solidMoles = SolidTotalMoles();
    double atolVal = solidMoles * atolBaseResid_;

    bool stepAcceptable = true;


    for (size_t ph = 0; ph < phaseIndexSolidPhases_.size(); ph++) {
        size_t iph = phaseIndexSolidPhases_[ph];
        std::string ppp = phase_name(iph);
        if (justDiedPhase_[iph]) {
            if (!justDied_[iph]) {
                if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
                    printf(" XXX integrate(): Phase didn't die like it should %d %s\n", static_cast<int>(iph), ppp.c_str());
                    printf("                  Final moles = %g, initial moles = %g\n", phaseMoles_final_[iph], phaseMoles_init_[iph]);
                    printf("                  Tfinal = %g Tinitial = %g, deltaTsubcycle_ = %g, deltaTsubcycleCalc_ = %g\n",
                           tfinal_, tinit_, deltaTsubcycle_, deltaTsubcycleCalc_);
                }
                stepAcceptable = false;
                break;
            }
            for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                size_t isp = globalSpeciesIndex(iph,sp);
                if (spMoles_final_[isp] > 0.0) {
                    if (spMoles_final_[isp] > spMoles_init_[isp] * 1.0E-6) {
                        if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
                            printf(" XXX integrate(): species in zeroed phase nonzero %d = %g\n", static_cast<int>(isp), spMoles_final_[isp]);
                            printf(" XXX integrate(): Phase didn't die like it should have %d %s\n", static_cast<int>(iph), ppp.c_str());
                        }
                        spMoles_final_[isp] = 0.0;
                        stepAcceptable = false;
                        break;
                    } else {
                        spMoles_final_[isp] = 0.0;
                    }
                } else if (spMoles_final_[isp] < 0.0) {
                    if (fabs(spMoles_final_[isp]) > spMoles_init_[isp] * 1.0E-6) {
                        if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
                            printf(" XXX integrate(): species in zeroed phase nonzero %lu = %g\n", 
                                   isp, spMoles_final_[isp]);
                        }
                        spMoles_final_[isp] = 0.0;
                        stepAcceptable = false;
                        break;
                    } else {
                        spMoles_final_[isp] = 0.0;
                    }
                }
            }
        } else {
            for (size_t sp = 0; sp < numSpecInSolidPhases_[ph]; sp++) {
                size_t isp = globalSpeciesIndex(iph,sp);
                if (spMoles_final_[isp] < 0.0) {
                    if (spMoles_final_[isp] > -atolVal) {
                        spMoles_final_[isp] = 0.0;
                    } else {
                        printf(" XXX integrate(): species in nonzeroed phase below zero %lu = %g\n", 
                               isp, spMoles_final_[isp]);
                        stepAcceptable = false;
                        break;
                    }
                }
            }
        }
    }
    return stepAcceptable;
}
//====================================================================================================================
void  Electrode_CSTR::manageBirthDeathSuccessfulStep()
{
}
//==================================================================================================================
//   Calculate the integrated source terms and do other items now that we have a completed time step
/*
 *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
 *  for the current step, and we are calculating post-processed quantities like source terms.
 */
void Electrode_CSTR::calcSrcTermsOnCompletedStep()
{
    bool doOneWay = false;
    if (doOneWay) {
        /*
         *  Calculate the integrated source term
         *       An alternative would be to redo the residual calculation. However, here we assume that
         *       the residual calculation has been done and the results are in _final_
         */
        for (size_t i = 0; i < m_NumTotSpecies; i++) {
            spMoleIntegratedSourceTermLast_[i] = spMoles_final_[i] - spMoles_init_[i];
        }
    } else {
        extractInfo();
        updateSpeciesMoleChangeFinal();
        for (size_t isp = 0; isp < m_NumTotSpecies; isp++) {
            spMoleIntegratedSourceTermLast_[isp] = DspMoles_final_[isp] * deltaTsubcycleCalc_;
        }
    }
    if (doThermalPropertyCalculations_) {
        integratedThermalEnergySourceTermLast_ = thermalEnergySourceTerm_EnthalpyFormulation_SingleStep();
        /*
         * these last two are needed for informative output only
         */
        integratedThermalEnergySourceTerm_overpotential_Last_ = thermalEnergySourceTerm_Overpotential_SingleStep();
        integratedThermalEnergySourceTerm_reversibleEntropy_Last_ = thermalEnergySourceTerm_ReversibleEntropy_SingleStep();
    }
}
//====================================================================================================================
void Electrode_CSTR::gatherIntegratedSrcPrediction()
{
    extractInfo();
    updateSpeciesMoleChangeFinal();
    IntegratedSrc_Predicted.resize(m_NumTotSpecies);
    for (size_t isp = 0; isp < m_NumTotSpecies; isp++) {
        IntegratedSrc_Predicted[isp] = DspMoles_final_[isp] * deltaTsubcycleCalc_;
    }
}

//====================================================================================================================
//  Calculate the norm of the difference between the predicted answer and the final converged answer
//  for the current time step
/*
 *  (virtual from Electrode_Integrator)
 *
 *   The norm calculated by this routine is used to determine whether the time step is accurate enough.
 *
 *  @return    Returns the norm of the difference. Normally this is the L2 norm of the difference
 */
double Electrode_CSTR::predictorCorrectorWeightedSolnNorm(const std::vector<double>& yval)
{
    double pnorm = l0norm_PC_NLS(soln_predict_, yval, neq_, atolNLS_, rtolNLS_);
    return pnorm;
}
//====================================================================================================================
void Electrode_CSTR::predictorCorrectorGlobalSrcTermErrorVector()
{

}
//====================================================================================================================
void Electrode_CSTR::predictorCorrectorPrint(const std::vector<double>& yval, double pnormSrc, double pnormSoln) const
{
    double atolVal =  1.0E-8;
    double denom;
    double tmp;
    double RelativeExtentRxnPredict = soln_predict_[neq_];
    int onRegionPredict = soln_predict_[neq_ + 1];
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf(" PREDICTOR_CORRECTOR  SubIntegrationCounter = %7d       t_init = %12.5E,       t_final = %12.5E\n",
           counterNumberSubIntegrations_, tinit_, tfinal_);
    printf("                         IntegrationCounter = %7d  t_init_init = %12.5E, t_final_final = %12.5E\n",
           counterNumberIntegrations_, t_init_init_, t_final_final_);
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("                               Initial       Prediction      Actual          Difference         Tol   Contrib      |\n");
    printf("onRegionPredict          |         %3d            %3d          %3d      |                                          |\n",
           onRegionBoundary_init_, onRegionPredict, onRegionBoundary_final_);
    denom = MAX(fabs(yval[0]), fabs(soln_predict_[0]));
    denom = MAX(denom, atolVal);
    tmp = fabs((yval[0] - soln_predict_[0])/ denom);
    printf("DeltaT                   | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E |\n",
           deltaTsubcycle_, soln_predict_[0],  yval[0], yval[0] - soln_predict_[0], atolVal, tmp);
    denom = MAX(fabs(RelativeExtentRxn_final_), fabs(RelativeExtentRxnPredict));
    denom = MAX(denom, atolVal);
    tmp = fabs((RelativeExtentRxn_final_ - RelativeExtentRxnPredict)/ denom);
    printf("RelativeExtentRxn        | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E | \n",
           RelativeExtentRxn_init_, RelativeExtentRxnPredict,  RelativeExtentRxn_final_, RelativeExtentRxn_final_ - RelativeExtentRxnPredict, atolVal, tmp);
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("                                                                                                        %10.3E\n",
           pnormSoln);

}
//====================================================================================================================
double Electrode_CSTR::predictorCorrectorGlobalSrcTermErrorNorm()
{
    return 0.0;
}

//====================================================================================================================
void Electrode_CSTR::check_final_state()
{
    /*
     * Reconcile RelativeExtentRxn_final with mole numbers
     */
    double relEcalc = calcRelativeExtentRxn_final();
    if (fabs(relEcalc - RelativeExtentRxn_final_) > 1.0E-6) {
        printf("problems\n");
        exit(-1);
    }
    setOnRegionBoundary();
}
//====================================================================================================================
void Electrode_CSTR::setInitStateFromFinal_Oin(bool setInitInit)
{
    Electrode_Integrator::setInitStateFromFinal(setInitInit);

    RelativeExtentRxn_init_        = RelativeExtentRxn_final_;
    RelativeExtentRxn_final_final_ = RelativeExtentRxn_final_;
    onRegionBoundary_init_         = onRegionBoundary_final_;
    onRegionBoundary_final_final_  = onRegionBoundary_final_;
    xRegion_init_                  = xRegion_final_;
    xRegion_final_final_           = xRegion_final_;

    if (setInitInit) {
        RelativeExtentRxn_init_init_   = RelativeExtentRxn_final_;
        onRegionBoundary_init_init_    = onRegionBoundary_final_;
        xRegion_init_init_             = xRegion_final_;
    }
}
//====================================================================================================================
// Set the internal initial intermediate from the internal initial global state
/*
 *  Set the intial state from the init init state. We also can set the final state from this
 *  routine as well.
 *
 *  The final_final is not touched.
 *
 * @param setFinal   Boolean indicating whether you should set the final as well
 */
void Electrode_CSTR::setInitStateFromInitInit(bool setFinal)
{  
    Electrode_Integrator::setInitStateFromInitInit(setFinal);

    RelativeExtentRxn_init_        = RelativeExtentRxn_init_init_;
    onRegionBoundary_init_         = onRegionBoundary_init_init_;
    xRegion_init_                  = xRegion_init_init_;
    if (setFinal) {
	RelativeExtentRxn_final_       = RelativeExtentRxn_init_init_;
	onRegionBoundary_final_        = onRegionBoundary_init_init_;
	xRegion_final_                 = xRegion_init_init_;
    }
}
//====================================================================================================================
void Electrode_CSTR::setInitInitStateFromFinalFinal()
{
    Electrode_Integrator::setInitInitStateFromFinalFinal();

    RelativeExtentRxn_init_init_   = RelativeExtentRxn_final_final_;
    RelativeExtentRxn_init_        = RelativeExtentRxn_final_final_;
    onRegionBoundary_init_init_    = onRegionBoundary_final_final_;
    onRegionBoundary_init_         = onRegionBoundary_final_final_;
    xRegion_init_                  = xRegion_final_final_;
    xRegion_init_init_             = xRegion_final_final_;
}
//====================================================================================================================
// Set the internal final intermediate and from the internal init state
/*
 *  (virtual function from Electrode)
 *
 *  Set the final state from the init state. This is commonly called during a failed time step
 *
 */
void  Electrode_CSTR::setFinalStateFromInit()
{
    Electrode_Integrator::setFinalStateFromInit_Oin();
    /*
     * Do stuff not done in base class
     */
    RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
    xRegion_final_ = xRegion_init_;
    onRegionBoundary_final_ = onRegionBoundary_init_;


    /*
     *  If this object becomes a parent, then we will have to create a setFinalStateFromInit_Oin() to avoid
     *  calling this routine.
     */
    updateState();
}
//====================================================================================================================
//! Set the internal final global state from the internal final intermediate state
/*!
 *  (virtual function from Electrode)
 *
 *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
 */
void Electrode_CSTR::setFinalFinalStateFromFinal()
{
    Electrode_Integrator::setFinalFinalStateFromFinal_Oin();

    RelativeExtentRxn_final_final_ = RelativeExtentRxn_final_;
    xRegion_final_final_           = xRegion_final_;
    onRegionBoundary_final_final_  = onRegionBoundary_final_;
}
//====================================================================================================================
//   Set the internal initial intermediate and initial global state from the internal final state
/*
 *  (virtual function)
 *
 *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 * @param setInitInit   Boolean indicating whether you should set the init_init state as well
 */
void Electrode_CSTR::setInitStateFromFinal(bool setInitInit)
{
    setInitStateFromFinal_Oin(setInitInit);
}
//==================================================================================================================
void Electrode_CSTR::setNLSGlobalSrcTermTolerances(double rtolResid)
{
    double sum = SolidTotalMoles();
    double val = 1.0E-14 * sum;

    for (size_t i = 0; i < numIntegratedSrc_; i++) {
        atol_IntegratedSrc_global_[i] = val;
    }
    rtol_IntegratedSrc_global_ = rtolResid;
}
//======================================================================================================================
double Electrode_CSTR::calcRelativeExtentRxn_final() const
{
    double capInit = capacityRaw();
    double capLeft = capacityLeftRaw();
    double tmp = (capInit - capLeft) / (Faraday * RelativeExtentRxn_NormalizationFactor_);
    return tmp;
}


//====================================================================================================================
} // End of ZZCantera namespace
//----------------------------------------------------------------------------------------------------------------------------------

