/*
 * $Id: Electrode_MP_RxnExtent.cpp 604 2013-05-24 16:27:35Z hkmoffa $
 */

#include "Electrode_MP_RxnExtent.h"
#include "Electrode_input.h"

#include "BlockEntryGlobal.h"

#include "cantera/thermo.h"

using namespace std;
using namespace BEInput;
using namespace TKInput;

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
#ifdef DEBUG_MODE_PREDICTION
double predictSave[30];
#endif
//====================================================================================================================
ELECTRODE_MP_RxnExtent_KEY_INPUT::ELECTRODE_MP_RxnExtent_KEY_INPUT(int printLvl) :
    ELECTRODE_KEY_INPUT(printLvl),
    numRegions_(1),
    solidDiffusionModel_(0),
    locationOfReactingSurface_(0),
    rxnPerturbRegions_(0)
{

}
//====================================================================================================================
ELECTRODE_MP_RxnExtent_KEY_INPUT::~ELECTRODE_MP_RxnExtent_KEY_INPUT()
{

}
//====================================================================================================================
ELECTRODE_MP_RxnExtent_KEY_INPUT::ELECTRODE_MP_RxnExtent_KEY_INPUT(const ELECTRODE_MP_RxnExtent_KEY_INPUT& right) :
    ELECTRODE_KEY_INPUT(right),
    numRegions_(right.numRegions_),
    solidDiffusionModel_(right.solidDiffusionModel_),
    locationOfReactingSurface_(right.locationOfReactingSurface_)
{
    diffusionCoeffRegions_          = right.diffusionCoeffRegions_;
    molarVolumeRegionBoundaries_    = right.molarVolumeRegionBoundaries_;
    rxnPerturbRegions_              = right.rxnPerturbRegions_;
}
//====================================================================================================================
ELECTRODE_MP_RxnExtent_KEY_INPUT&
ELECTRODE_MP_RxnExtent_KEY_INPUT::operator=(const ELECTRODE_MP_RxnExtent_KEY_INPUT& right)
{
    if (this == &right) {
        return *this;
    }

    ELECTRODE_KEY_INPUT::operator=(right);

    numRegions_                     = right.numRegions_;
    solidDiffusionModel_            = right.solidDiffusionModel_;
    diffusionCoeffRegions_          = right.diffusionCoeffRegions_;
    molarVolumeRegionBoundaries_    = right.molarVolumeRegionBoundaries_;
    locationOfReactingSurface_      = right.locationOfReactingSurface_;
    rxnPerturbRegions_              = right.rxnPerturbRegions_;

    return *this;
}
//====================================================================================================================
void ELECTRODE_MP_RxnExtent_KEY_INPUT::setup_input_child1(BEInput::BlockEntry* cf)
{
    /*
     * Obtain the number of regions
     */
    LE_OneInt* s1 = new LE_OneInt("Number of Regions", &(numRegions_), 0, "numRegions");
    s1->set_default(1);
    cf->addLineEntry(s1);
    BaseEntry::set_SkipUnknownEntries(3);

    /*
     * Obtain the number of regions
     */
    LE_OneInt* sdm = new LE_OneInt("Solid Diffusion Model", &(solidDiffusionModel_), 0, "solidDiffusionModel");
    sdm->set_default(0);
    cf->addLineEntry(sdm);
    BaseEntry::set_SkipUnknownEntries(3);


}
//====================================================================================================================
void ELECTRODE_MP_RxnExtent_KEY_INPUT::setup_input_child2(BEInput::BlockEntry* cf)
{
    LE_StdVecDblVarLength* v1 = new LE_StdVecDblVarLength("Molar Volumes Region Boundaries",
            &(molarVolumeRegionBoundaries_),
            numRegions_ + 1, numRegions_ + 1,
            "molarVolumesRegionBoundaries");
    v1->set_default(50.0);
    v1->set_limits(1000., 0.001);
    cf->addLineEntry(v1);

    /*
     * Obtain the number of regions
     */
    int reqd = 0;
    if (solidDiffusionModel_) {
        reqd = 1;
    }
    LE_OneInt* s3 = new LE_OneInt("Location of Rate Limiting Surface", &(locationOfReactingSurface_), reqd,
                                  "locationOfReactingSurface");
    s3->set_default(0);

    cf->addLineEntry(s3);

    /*
     * Make the diffusion coefficients optional unless the solid diffusion model is turned on
     */
    reqd = 0;
    if (solidDiffusionModel_) {
        reqd = numRegions_  + 1;
    }

    LE_StdVecDblVarLength* d1 = new LE_StdVecDblVarLength("Diffusion Coefficients for Regions",
            &(diffusionCoeffRegions_),
            numRegions_ + 1, reqd, "diffusionCoeffRegions_");
    d1->set_default(0.0);
    d1->set_limits(1000., 1.0E-30);
    cf->addLineEntry(d1);

    reqd = 0;

    LE_StdVecDblVarLength* rr1 = new LE_StdVecDblVarLength("Reaction Rate Constant Perturbations for Regions",
            &(rxnPerturbRegions_),
            numRegions_, reqd, "rxnPerturbRegions_");
    rr1->set_default(1.0);
    rr1->set_limits(1.0E30, 1.0E-30);
    cf->addLineEntry(rr1);

    BaseEntry::set_SkipUnknownEntries(0);
}
//====================================================================================================================
Electrode_MP_RxnExtent::Electrode_MP_RxnExtent() :
    Electrode_Integrator(),
    numRegions_(0),
    RelativeExtentRxn_init_(0.0),
    RelativeExtentRxn_init_init_(0.0),
    RelativeExtentRxn_final_(0.0),
    RelativeExtentRxn_final_final_(0.0),
    RelativeExtentRxn_tmp_(0.0),
    xRegion_init_(0),
    xRegion_init_init_(0),
    xRegion_final_(0),
    xRegion_final_final_(0),
    // deltaTsubcycleCalc_(0.0),
    RegionBoundaries_ExtentRxn_(0),
    SrcDot_ExtentRxn_final_(0.0),
    deltaTdeath_(0.0),
    ROP_(0),
    DspMoles_final_(0),
    phaseIndexSolidPhases_(0),
    numSpecInSolidPhases_(0),
    Li_liq_(0),
    Hf_B_base_(0.0),
    Hf_B_current_(0.0),
    volts_OCV_SS_final_(0.0),
    volts_OCV_final_(0.0),
    spMoles_FeS2_Normalization_(0.0),
    onRegionBoundary_init_(0),
    onRegionBoundary_final_(0),
    goNowhere_(0),
    Radius_internal_final_(0.0),
    Radius_internal_final_final_(0.0),
    Radius_internal_init_(0.0),
    Radius_internal_init_init_(0.0),
    locationOfReactingSurface_(0),
    ROPModificationType_(0),
    indexOfReactingSurface_(0),
    solidDiffusionModel_(0),
    diffusionCoeffRegions_(0),
    molarVolumeRegions_(0),
    rxnPerturbRegions_(0),
    actEquilInterstitialsRegions_(0),
    kf_id_(0),
    kf_dir_(0),
    kfExt_(0.0),
    krExt_(0.0),
    kfInner_(0.0),
    krInner_(0.0),
    DaOuter_(0.0),
    DaOuter_Bar_(0.0),
    DaInner_(0.0),
    Lin_(0.0),
    Lout_(0.0),
    kfAB_(0.0),
    krAB_(0.0),
    ca_Lip_(0.0),
    limitingEquationBehavior_(0),
    innerDaTreatmentType_(2)
{
    neq_ = 2;
    ROP_.resize(10);
    Li_liq_ = newPhase("Li_Liquid.xml","Li(L)");
    if (!Li_liq_) {
        throw CanteraError("", "");
    }
    setVoltageVsExtent_FeS2();
    soln_predict_.resize(3);
}
//====================================================================================================================
//====================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_MP_RxnExtent::Electrode_MP_RxnExtent(const Electrode_MP_RxnExtent& right) :
    Electrode_Integrator(),
    numRegions_(0),
    RelativeExtentRxn_init_(0.0),
    RelativeExtentRxn_init_init_(0.0),
    RelativeExtentRxn_final_(0.0),
    RelativeExtentRxn_final_final_(0.0),
    RelativeExtentRxn_tmp_(0.0),
    xRegion_init_(0),
    xRegion_init_init_(0),
    xRegion_final_(0),
    xRegion_final_final_(0),
    // deltaTsubcycleCalc_(0.0),
    RegionBoundaries_ExtentRxn_(0),
    SrcDot_ExtentRxn_final_(0.0),
    deltaTdeath_(0.0),
    ROP_(0),
    DspMoles_final_(0),
    phaseIndexSolidPhases_(0),
    numSpecInSolidPhases_(0),
    Li_liq_(0),
    Hf_B_base_(0.0),
    Hf_B_current_(0.0),
    volts_OCV_SS_final_(0.0),
    volts_OCV_final_(0.0),
    spMoles_FeS2_Normalization_(0.0),
    //    soln_predict_(0),
    onRegionBoundary_init_(0),
    onRegionBoundary_final_(0),
    goNowhere_(0),
    molarVolume_final_(0.0),
    molarVolume_final_final_(0.0),
    molarVolume_init_(0.0),
    molarVolume_init_init_(0.0),
    Radius_internal_final_(0.0),
    Radius_internal_final_final_(0.0),
    Radius_internal_init_(0.0),
    Radius_internal_init_init_(0.0),
    locationOfReactingSurface_(0),
    ROPModificationType_(0),
    indexOfReactingSurface_(0),
    solidDiffusionModel_(0),
    diffusionCoeffRegions_(0),
    molarVolumeRegions_(0),
    rxnPerturbRegions_(0),
    actEquilInterstitialsRegions_(0),
    kf_id_(0),
    kf_dir_(0),
    kfExt_(0.0),
    krExt_(0.0),
    kfInner_(0.0),
    krInner_(0.0),
    DaOuter_(0.0),
    DaOuter_Bar_(0.0),
    DaInner_(0.0),
    Lin_(0.0),
    Lout_(0.0),
    kfAB_(0.0),
    krAB_(0.0),
    ca_Lip_(0.0),
    limitingEquationBehavior_(0),
    innerDaTreatmentType_(2)
{
    neq_ = 2;
    ROP_.resize(10);
    Li_liq_ = newPhase("Li_Liquid.xml","Li(L)");
    if (!Li_liq_) {
        throw CanteraError("", "");
    }
    *this = operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_MP_RxnExtent& Electrode_MP_RxnExtent::operator=(const Electrode_MP_RxnExtent& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode_Integrator::operator=(right);

    numRegions_                    = right.numRegions_;
    RelativeExtentRxn_init_        = right.RelativeExtentRxn_init_;
    RelativeExtentRxn_init_init_   = right.RelativeExtentRxn_init_init_;
    RelativeExtentRxn_final_       = right.RelativeExtentRxn_final_;
    RelativeExtentRxn_final_final_ = right.RelativeExtentRxn_final_final_;
    RelativeExtentRxn_tmp_         = right.RelativeExtentRxn_tmp_;

    xRegion_init_                  = right.xRegion_init_;
    xRegion_init_init_             = right.xRegion_init_init_;
    xRegion_final_                 = right.xRegion_final_;
    xRegion_final_final_           = right.xRegion_final_final_;
    deltaTsubcycleCalc_            = right.deltaTsubcycleCalc_;


    RegionBoundaries_ExtentRxn_    = right.RegionBoundaries_ExtentRxn_;
    SrcDot_ExtentRxn_final_        = right.SrcDot_ExtentRxn_final_;
    deltaTdeath_                   = right.deltaTdeath_;

    detailedResidPrintFlag_        = right.detailedResidPrintFlag_;
    enableExtraPrinting_           = right.enableExtraPrinting_;

    ROP_                           = right.ROP_;

    DspMoles_final_                = right.DspMoles_final_;
    phaseIndexSolidPhases_         = right.phaseIndexSolidPhases_;
    numSpecInSolidPhases_          = right.numSpecInSolidPhases_;

    delete Li_liq_;
    Li_liq_ = right.Li_liq_->duplMyselfAsThermoPhase();
    Hf_B_base_                     = right.Hf_B_base_;
    Hf_B_current_                  = right.Hf_B_current_;
    volts_OCV_SS_final_            = right.volts_OCV_SS_final_;
    volts_OCV_final_               = right.volts_OCV_final_;
    spMoles_FeS2_Normalization_    = right.spMoles_FeS2_Normalization_;

    onRegionBoundary_init_         = right.onRegionBoundary_init_;
    onRegionBoundary_final_        = right.onRegionBoundary_final_;
    goNowhere_                     = right.goNowhere_;

    molarVolume_final_             = right.molarVolume_final_;
    molarVolume_final_final_       = right.molarVolume_final_final_;
    molarVolume_init_              = right.molarVolume_init_;
    molarVolume_init_init_         = right.molarVolume_init_init_;

    Radius_internal_final_         = right.Radius_internal_final_;
    Radius_internal_final_final_   = right.Radius_internal_final_final_;
    Radius_internal_init_          = right.Radius_internal_init_;
    Radius_internal_init_init_     = right.Radius_internal_init_init_;
    locationOfReactingSurface_     = right.locationOfReactingSurface_;
    ROPModificationType_           = right.ROPModificationType_;
    indexOfReactingSurface_        = right.indexOfReactingSurface_;

    solidDiffusionModel_           = right.solidDiffusionModel_;
    diffusionCoeffRegions_         = right.diffusionCoeffRegions_;
    molarVolumeRegions_            = right.molarVolumeRegions_;
    rxnPerturbRegions_             = right.rxnPerturbRegions_;
    actEquilInterstitialsRegions_  = right.actEquilInterstitialsRegions_;
    kf_id_                         = right.kf_id_;
    kf_dir_                        = right.kf_dir_;

    kfExt_                         = right.kfExt_;
    krExt_                         = right.krExt_;
    kfInner_                       = right.kfInner_;
    krInner_                       = right.krInner_;
    DaOuter_                       = right.DaOuter_;
    DaOuter_Bar_                   = right.DaOuter_Bar_;
    DaInner_                       = right.DaInner_;
    Lin_                           = right.Lin_;
    Lout_                          = right.Lout_;
    kfAB_                          = right.kfAB_;
    krAB_                          = right.krAB_;
    ca_Lip_                        = right.ca_Lip_;
    limitingEquationBehavior_      = right.limitingEquationBehavior_;
    innerDaTreatmentType_          = right.innerDaTreatmentType_;


    return *this;
}
//=======================================================================================================
Electrode_MP_RxnExtent::~Electrode_MP_RxnExtent()
{

    delete Li_liq_;
    Li_liq_ = 0;


}
//=======================================================================================================
// Return the type of electrode
/*
 *  Returns the enum type of the electrode. This is used in the factory routine.
 *
 *  @return Returns an enum type, called   Electrode_Types_Enum
 */
Electrode_Types_Enum  Electrode_MP_RxnExtent::electrodeType() const
{
    return MP_RXNEXTENT_ET;
}
//====================================================================================================================
int Electrode_MP_RxnExtent::electrode_input_child(ELECTRODE_KEY_INPUT** ei_ptr)
{
    /*
     *  malloc an expanded child input
     */
    ELECTRODE_MP_RxnExtent_KEY_INPUT* ei_mp = new ELECTRODE_MP_RxnExtent_KEY_INPUT();
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
//====================================================================================================================
//  Setup the electrode
/*
 * @param ei    ELECTRODE_KEY_INPUT pointer object
 */
int
Electrode_MP_RxnExtent::electrode_model_create(ELECTRODE_KEY_INPUT* eibase)
{

    int flag;

    /*
     *  Downcast the Key input to make sure we are being fed the correct child object
     */
    ELECTRODE_MP_RxnExtent_KEY_INPUT* ei = dynamic_cast<ELECTRODE_MP_RxnExtent_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_MP_RxnExtent::electrode_model_create()",
                           " Expecting a child ELECTRODE_KEY_INPUT object and didn't get it");
    }
    radiusSmallBound_  = 1.0E-12;
    /*
     *  Call the base routine.
     *   During this routine typcally it reads a field called "Capacity Discharged per Mole = 1.7"
     *   This sets the RelativeExtentRxn_final_ value to an initial value, RelativeExtentRxn_final_ is the main
     *   state of this object.
     */
    flag = Electrode_Integrator::electrode_model_create(ei);

    if (flag != 0) {
        return flag;
    }

    DspMoles_final_.resize(m_NumTotSpecies, 0.0);

    numRegions_ = ei->numRegions_;
    molarVolumeRegions_.resize(numRegions_ + 1, 0.0);
    molarVolumeRegions_ = ei->molarVolumeRegionBoundaries_;
    rxnPerturbRegions_ = ei->rxnPerturbRegions_;
    if ((int) rxnPerturbRegions_.size() != numRegions_) {
        rxnPerturbRegions_.resize(numRegions_, 1.0);
    }
    molarVolume_final_ = molarVolumeRegions_[0];
    if ((int) molarVolumeRegions_.size() != numRegions_ + 1) {
        throw CanteraError(" Electrode_MP_RxnExtent::electrode_model_create()",
                           " Wrong size of molarVolumeRegions_");
    }

    diffusionCoeffRegions_ = ei->diffusionCoeffRegions_;
    if ((int) diffusionCoeffRegions_.size() != numRegions_ + 1) {
        throw CanteraError(" Electrode_MP_RxnExtent::electrode_model_create()",
                           " Wrong size of  diffusionCoeffRegions_");
    }


    for (size_t ph = 0; ph < NumVolPhases_; ph++) {
        ThermoPhase* tp = VolPhaseList[ph];
        int iph = getGlobalPhaseIndex(tp);
        if (iph == metalPhaseIndex() || iph == solnPhaseIndex()) {
            //do nothing
        } else {
            phaseIndexSolidPhases_.push_back(iph);
            numSpecInSolidPhases_.push_back(tp->nSpecies());
        }
    }
    developBaseE0();

    atolNLS_.resize(2, 1.0E-10);

    /*
     *  Figure out the normalizing quantity for between relative and absolute moles of FeS2.
     *  What we will do is to calculate the initial number of moles of FeS2 at zero discharge
     *  in the electrode. Then, we will set the moles of A and B to be 1/2 of that no matter
     *  where the electrode is in the DOD.
     */
    int ip_FeS2_A = globalPhaseIndex("FeS2_A(S)");
    int is_FeS2_A = globalSpeciesIndex("FeS2_A(S)");
    int ip_FeS2_B = globalPhaseIndex("FeS2_B(S)");
    int is_FeS2_B = globalSpeciesIndex("FeS2_B(S)");

    /*
     *  Now save the mole numbers as a normalization number
     *     This can be considered as the number of moles of FeS2 in the solid electrode initially
     */
    spMoles_FeS2_Normalization_ = spMoles_final_[is_FeS2_A] + spMoles_final_[is_FeS2_B];

    /*
     *  resize the mole numbers to the geometry.
     */
    resizeMoleNumbersToGeometry();



    electronKmolDischargedToDate_ = -  RelativeExtentRxn_final_ * spMoles_FeS2_Normalization_;

    xRegion_final_ = findRegion(RelativeExtentRxn_final_);
    xRegion_init_ =   xRegion_final_;
    xRegion_init_init_ =   xRegion_final_;
    xRegion_final_final_ = xRegion_final_;
    onRegionBoundary_init_ = -1;
    onRegionBoundary_final_ = -1;
    for (int i = 0; i < (int) RegionBoundaries_ExtentRxn_.size(); i++) {
        if (fabs(RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[i]) < 1.0E-12) {
            onRegionBoundary_init_ = i;
            onRegionBoundary_final_ = i;
            break;
        }
    }

    /*
     *  Transfer solid phase diffusion parameters
     */
    solidDiffusionModel_  =  ei->solidDiffusionModel_;
    locationOfReactingSurface_ = ei->locationOfReactingSurface_;
    diffusionCoeffRegions_ = ei->diffusionCoeffRegions_;


    if (solidDiffusionModel_ == 0) {
        ROPModificationType_ = 0;
    } else {
        if (locationOfReactingSurface_ == 0 ||  locationOfReactingSurface_ == 2) {
            ROPModificationType_ = 2;
        } else {
            ROPModificationType_ = 1;
        }
    }

    /*
     *  We need at least three surfaces for this model:
     *              0 inner surface
     *              1 exterior surface
     *              2 initial surface that just stays constant
     */

    if (surfaceAreaRS_final_.size() <= 2) {
        /*
         * Add a nonreacting outer surface
         */
        indexOfReactingSurface_ = 0;
        numSurfaces_ = 2;
        double sa = surfaceAreaRS_final_[0];
        surfaceAreaRS_final_.resize(3, sa);
        surfaceAreaRS_init_.resize(3, sa);
        surfaceAreaRS_init_init_.resize(3, sa);
        surfaceAreaRS_final_final_.resize(3, sa);
        spNetProdPerArea_List_.resize(m_NumTotSpecies, numSurfaces_, 0.0);
        isExternalSurface_.resize(3, 0);
        isExternalSurface_[0] = false;
        isExternalSurface_[1] = true;
        isExternalSurface_[2] = true;
        RSD_List_.resize(numSurfaces_, 0);
        numRxns_.resize(numSurfaces_, 0);
        ActiveKineticsSurf_.resize(2, false);
    }
    if (locationOfReactingSurface_ == 1) {
        if (RSD_List_[1] == 0) {
            indexOfReactingSurface_ = 1;
            RSD_List_[1] = RSD_List_[0];
            RSD_List_[0] = 0;
            numRxns_[1] = numRxns_[0];
            numRxns_[0] = 0;
            ActiveKineticsSurf_[1] = true;
            ActiveKineticsSurf_[0] = false;
        } else if (locationOfReactingSurface_ == 0 || locationOfReactingSurface_ == 2) {
            if (RSD_List_[0] == 0) {
                indexOfReactingSurface_ = 0;
                RSD_List_[0] = RSD_List_[1];
                RSD_List_[1] = 0;
                ActiveKineticsSurf_[1] = false;
                ActiveKineticsSurf_[0] = true;
                numRxns_[0] = numRxns_[1];
                numRxns_[1] = 0;
            }
        }
    }

    /*
     *
     */
    molarVolume_final_ = molarVolume_relExtentRxn(RelativeExtentRxn_final_);
    /*
     *  resize the mole numbers to the geometry.
     */
    resizeMoleNumbersToGeometry();

    /*
     *  set the moles numbers for A and B to be in the middle.
     */
    spMoles_final_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
    spMoles_init_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
    spMoles_init_init_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
    phaseMoles_final_[ip_FeS2_A] =  spMoles_FeS2_Normalization_ / 2.0;

    spMoles_final_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
    spMoles_init_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
    spMoles_init_init_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
    phaseMoles_final_[ip_FeS2_B] =  spMoles_FeS2_Normalization_ / 2.0;

    std::copy(spMf_final_.begin(), spMf_final_.end(), spMf_init_init_.begin());
    std::copy(spMf_final_.begin(), spMf_final_.end(), spMf_init_.begin());

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        updateState_Phase(iph);
    }
    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_.begin());
    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_init_.begin());

  
    numIntegratedSrc_ = m_NumTotSpecies;
    create_solvers();
    /*
     *  Calculate the solid volume and the external radius. Apparently, this step is needed
     */
    ElectrodeSolidVolume_ = SolidVol();
    double vol =  ElectrodeSolidVolume_ / particleNumberToFollow_;
    Radius_exterior_final_ = pow(vol * 3.0 / (4.0 * Pi), 0.3333333333333333);
    /*
     *  Update the calculation of the internal radius and all of the surface areas
     */
    Electrode_MP_RxnExtent::updateSurfaceAreas();
    for (size_t i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_init_[i]        = surfaceAreaRS_final_[i];
        surfaceAreaRS_final_final_[i] = surfaceAreaRS_final_[i];
        surfaceAreaRS_init_init_[i]   = surfaceAreaRS_final_[i];
    }
    Radius_internal_final_final_ = Radius_internal_final_;
    Radius_internal_init_init_ = Radius_internal_final_;
    Radius_internal_init_ = Radius_internal_final_;

    if (ei->RelativeCapacityDischargedPerMole != -1) {

        setRelativeCapacityDischargedPerMole(ei->RelativeCapacityDischargedPerMole);
    }

    /*
     *  Transfer solid phase diffusion parameters
     */
    solidDiffusionModel_  =  ei->solidDiffusionModel_;
    locationOfReactingSurface_ = ei->locationOfReactingSurface_;

    diffusionCoeffRegions_ = ei->diffusionCoeffRegions_;

   

    return 0;
}
//====================================================================================================================
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
 * @param ei    ELECTRODE_KEY_INPUT pointer object
 *
 *  @return  Returns zero if successful, and -1 if not successful.
 */
int Electrode_MP_RxnExtent::setInitialConditions(ELECTRODE_KEY_INPUT* eibase)
{

    ELECTRODE_MP_RxnExtent_KEY_INPUT* ei = dynamic_cast<ELECTRODE_MP_RxnExtent_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_MP_RxnExtent::setInitialConditions()",
                           " Expecting a child ELECTRODE_KEY_INPUT object,ELECTRODE_MP_RxnExtent_KEY_INPUT,  and didn't get it");
    }



    int flag = Electrode_Integrator::setInitialConditions(ei);
    if (flag != 0) {
        return flag;
    }

    DspMoles_final_.resize(m_NumTotSpecies, 0.0);

    if (ei->RelativeCapacityDischargedPerMole != -1) {
        setRelativeCapacityDischargedPerMole(ei->RelativeCapacityDischargedPerMole);
    }

    for (size_t ph = 0; ph < NumVolPhases_; ph++) {
        ThermoPhase* tp = VolPhaseList[ph];
        int iph = getGlobalPhaseIndex(tp);
        if (iph == metalPhaseIndex() || iph == solnPhaseIndex()) {
            //do nothing
        } else {
            phaseIndexSolidPhases_.push_back(iph);
            numSpecInSolidPhases_.push_back(tp->nSpecies());
        }
    }
    developBaseE0();

    /*
     *  Figure out the normalizing quantity for between relative and absolute moles of FeS2.
     *  What we will do is to calculate the initial number of moles of FeS2 at zero discharge
     *  in the electrode. Then, we will set the moles of A and B to be 1/2 of that no matter
     *  where the electrode is in the DOD.
     */
    int ip_FeS2_A = globalPhaseIndex("FeS2_A(S)");
    int is_FeS2_A = globalSpeciesIndex("FeS2_A(S)");
    int ip_FeS2_B = globalPhaseIndex("FeS2_B(S)");
    int is_FeS2_B = globalSpeciesIndex("FeS2_B(S)");

    /*
     *  Now save the mole numbers as a normalization number
     *     This can be considered as the number of moles of FeS2 in the solid electrode initially.
     *     This is an initial quantity that is overwritten in resizeMoleNumbersToGeometry()
     */
    spMoles_FeS2_Normalization_ = spMoles_final_[is_FeS2_A] + spMoles_final_[is_FeS2_B];

    /*
     *  resize the mole numbers to the geometry.
     */
    resizeMoleNumbersToGeometry();



    electronKmolDischargedToDate_ = -  RelativeExtentRxn_final_ * spMoles_FeS2_Normalization_;

    xRegion_final_ = findRegion(RelativeExtentRxn_final_);
    xRegion_init_ =   xRegion_final_;
    xRegion_init_init_ =   xRegion_final_;
    xRegion_final_final_ = xRegion_final_;
    onRegionBoundary_init_ = -1;
    onRegionBoundary_final_ = -1;
    for (int i = 0; i < (int) RegionBoundaries_ExtentRxn_.size(); i++) {
        if (fabs(RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[i]) < 1.0E-12) {
            onRegionBoundary_init_ = i;
            onRegionBoundary_final_ = i;
            break;
        }
    }

    /*
     *  We need at least two surfaces for this model
     */
    if (surfaceAreaRS_final_.size() <= 1) {
        numSurfaces_ = 2;
        double sa = surfaceAreaRS_final_[0];
        surfaceAreaRS_final_.resize(2, sa);
        surfaceAreaRS_init_.resize(2, sa);
        surfaceAreaRS_init_init_.resize(2, sa);
        surfaceAreaRS_final_final_.resize(2, sa);
        spNetProdPerArea_List_.resize(m_NumTotSpecies, numSurfaces_, 0.0);
        isExternalSurface_.resize(2, 0);
        isExternalSurface_[0] = false;
        isExternalSurface_[1] = true;
        RSD_List_.resize(numSurfaces_, 0);
        numRxns_.resize(numSurfaces_, 0);
        ActiveKineticsSurf_.resize(2, false);
    }
    if (locationOfReactingSurface_ == 1) {
        if (RSD_List_[1] == 0) {
            RSD_List_[1] = RSD_List_[0];
            RSD_List_[0] = 0;
            numRxns_[1] = numRxns_[0];
            numRxns_[0] = 0;
            ActiveKineticsSurf_[1] = true;
            ActiveKineticsSurf_[0] = false;
        } else if (locationOfReactingSurface_ == 0) {
            if (RSD_List_[0] == 0) {
                RSD_List_[0] = RSD_List_[1];
                RSD_List_[1] = 0;
                ActiveKineticsSurf_[1] = false;
                ActiveKineticsSurf_[0] = true;
                numRxns_[0] = numRxns_[1];
                numRxns_[1] = 0;
            }
        }
    }


    /*
     *  set the moles numbers for A and B to be in the middle.
     */
    spMoles_final_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
    spMoles_init_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
    spMoles_init_init_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
    phaseMoles_final_[ip_FeS2_A] =  spMoles_FeS2_Normalization_ / 2.0;

    spMoles_final_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
    spMoles_init_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
    spMoles_init_init_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
    phaseMoles_final_[ip_FeS2_B] =  spMoles_FeS2_Normalization_ / 2.0;


    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        updateState_Phase(iph);
    }
    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_.begin());
    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_init_.begin());

    /*
     *  Calculate the solid volume and the external radius. Apparently, this step is needed
     */
    ElectrodeSolidVolume_ = SolidVol();
    double vol =  ElectrodeSolidVolume_ / particleNumberToFollow_;
    Radius_exterior_final_ = pow(vol * 3.0 / (4.0 * Pi), 0.3333333333333333);
    /*
     *  Update the calculation of the internal radius and all of the surface areas
     */
    Electrode_MP_RxnExtent::updateSurfaceAreas();
    for (size_t i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_init_[i]        = surfaceAreaRS_final_[i];
        surfaceAreaRS_final_final_[i] = surfaceAreaRS_final_[i];
        surfaceAreaRS_init_init_[i]   = surfaceAreaRS_final_[i];
    }
    Radius_internal_final_final_ = Radius_internal_final_;
    Radius_internal_init_init_ = Radius_internal_final_;
    Radius_internal_init_ = Radius_internal_final_;


    /*
     *  Set the initial state and the init_init state from the final state.
     */
    setInitStateFromFinal_Oin(true);



    return 0;
}

//====================================================================================================================
int
Electrode_MP_RxnExtent::electrode_stateSave_create()
{
    eState_final_ = new EState();
    int rr = eState_final_->initialize(this);
    return rr;
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
void Electrode_MP_RxnExtent::setElectrodeSizeParams(double electrodeArea, double electrodeThickness,
        double porosity)
{

    Electrode_Integrator::setElectrodeSizeParams(electrodeArea, electrodeThickness,  porosity);
}
//====================================================================================================================
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
void Electrode_MP_RxnExtent::resizeMoleNumbersToGeometry()
{

    /*
     *  Here we duplicate reizeMoleMnumbersTo Geometry() because SolidVol() is overwritten.
     */

    double partVol = 4.0 * Pi *  Radius_exterior_final_ * Radius_exterior_final_ *  Radius_exterior_final_ / 3.0;
    double targetSolidVol = partVol * particleNumberToFollow_;
    if (spMoles_FeS2_Normalization_ <= 0.0) {

    }
    int is_FeS2_A = globalSpeciesIndex("FeS2_A(S)");
    int is_FeS2_B = globalSpeciesIndex("FeS2_B(S)");
    spMoles_FeS2_Normalization_ = spMoles_final_[is_FeS2_A] + spMoles_final_[is_FeS2_B];
    double currentSolidVol = SolidVol();
    if (currentSolidVol <= 0.0) {
        throw CanteraError("Electrode_MP_RxnExtent::resizeMoleNumbersToGeometry()",
                           "current solid volume is zero or less");
    }
    double ratio = targetSolidVol / currentSolidVol;
    double tMoles = 0.0;
    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        spMoles_final_[k] *= ratio;
        spMoles_init_[k] = spMoles_final_[k];
        spMoles_init_init_[k] = spMoles_final_[k];
        tMoles += spMoles_final_[k];
    }
    molarAtol_ = tMoles * 1.0E-5;
    for (size_t iph = 0; iph < (size_t) m_NumTotPhases; iph++) {
        Electrode::updateState_Phase(iph);
    }
    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_.begin());
    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_init_.begin());


    /*
     *  Now save the mole numbers as a normalization number
     *     This can be considered as the number of moles of FeS2 in the solid electrode initially.
     */
    spMoles_FeS2_Normalization_ = spMoles_final_[is_FeS2_A] + spMoles_final_[is_FeS2_B];

    currentSolidVol = SolidVol();
    if (fabs(currentSolidVol - targetSolidVol) > 1.0E-8*targetSolidVol) {
        throw CanteraError("Electrode_MP_RxnExtent::resizeMoleNumbersToGeometry() Error",
                           "Couldn't set the solid volume correctly: " + fp2str(currentSolidVol) + " vs " + fp2str(targetSolidVol));
    }


    /*
     * Resize the electrolyte so that the total volume of the electrolyte is consistent with the given
     * porosity, porosity_
     */
    // Electrode::resizeSolutionNumbersToPorosity();
    currentSolidVol = SolidVol();
    double grossVol = currentSolidVol / (1.0 - porosity_);
    double targetSolnVol = grossVol - currentSolidVol;
    double totalVol = TotalVol();
    double currentSolnVol = totalVol - currentSolidVol;
    ratio =  targetSolnVol / currentSolnVol;

    int istart = m_PhaseSpeciesStartIndex[solnPhase_];
    int nspSoln = m_PhaseSpeciesStartIndex[solnPhase_ +1] - m_PhaseSpeciesStartIndex[solnPhase_];
    for (int kk = 0; kk < nspSoln; kk++) {
        int k = istart + kk;
        spMoles_final_[k] *= ratio;
        spMoles_init_[k] = spMoles_final_[k];
        spMoles_init_init_[k] = spMoles_final_[k];
        spMoles_final_final_[k] = spMoles_final_[k];
    }
    updateState_Phase(solnPhase_);


    currentSolidVol = SolidVol();
    totalVol = TotalVol();
    double calcPor = (totalVol - currentSolidVol) / totalVol;
    if (fabs(calcPor - porosity_) > 1.0E-6) {
        throw CanteraError("Electrode_MP_RxnExtent::resizeMoleNumbersToGeometry() Error",
                           "Couldn't set the porosity correctly: " + fp2str(calcPor) + " vs " + fp2str(porosity_));
    }





    spMoles_FeS2_Normalization_ = spMoles_final_[is_FeS2_A] + spMoles_final_[is_FeS2_B];
    currentSolidVol = SolidVol();
    if (fabs(currentSolidVol - targetSolidVol) > 1.0E-8*targetSolidVol) {
        throw CanteraError("Electrode_MP_RxnExtent::resizeMoleNumbersToGeometry() Error",
                           "Couldn't set the solid volume correctly: " + fp2str(currentSolidVol) + " vs " + fp2str(targetSolidVol));
    }

    electronKmolDischargedToDate_ = - RelativeExtentRxn_final_ * spMoles_FeS2_Normalization_;
}
//====================================================================================================================
//  Set the initial value of DH_f for B to get a desired value
/*
 *
 */
void Electrode_MP_RxnExtent::developBaseE0()
{

    /*
     * Calculate Hf_Base_
     */
    double um[5];

    int ip_FeS2_A = globalPhaseIndex("FeS2_A(S)");
    int ip_FeS2_B = globalPhaseIndex("FeS2_B(S)");
    ThermoPhase* FeS2_A = &(thermo(ip_FeS2_A));
    ThermoPhase* FeS2_B = &(thermo(ip_FeS2_B));

    double tempO = temperature();
    double presO = pressure_;

    //setState_TP(723.15, OneAtm);
    temperature_ = 723.15;
    pressure_ = OneAtm;
    Electrode::updateState_OnionOut();

    Li_liq_->setState_TP(723.15, OneAtm);

    FeS2_A->getChemPotentials(um);
    double um_FeS2_A = um[0];

    FeS2_B->getChemPotentials(um);
    double um_FeS2_B = um[0];


    Li_liq_->getChemPotentials(um);
    double um_Li_liq = um[0];


    double deltaG = 1.0/4.0 * um_FeS2_B - 1.0/4.0 * um_FeS2_A - um_Li_liq;

    if (enableExtraPrinting_ && (detailedResidPrintFlag_ > 2)) {
        printf("delta G = %g J/kmol \n",  deltaG);
    }
    double voltsBase = - deltaG / Faraday;

    /*
     *  The base volts is 2.05565 for  2/3 FeS2 + Li+ + e- = 1/3 Li3Fe2S4   @723.15 K
     *  (Li -liquid basis)
     */
    double deltaG_plat = - 2.05565 *  Faraday;
    double hf = FeS2_B->Hf298SS(0);
    if (enableExtraPrinting_ && (detailedResidPrintFlag_ > 2)) {
        printf("  Unmodified hf for FeS2_B = %g\n", hf);
    }
    double deltaH = deltaG_plat - deltaG;

    Hf_B_base_ = hf + 4*deltaH;


    FeS2_B->modifyOneHf298SS(0, Hf_B_base_);


    FeS2_B->getChemPotentials(um);
    um_FeS2_B = um[0];


    deltaG = 1.0/4.0 * um_FeS2_B - 1.0/4.0 * um_FeS2_A - um_Li_liq;
    if (enableExtraPrinting_ && (detailedResidPrintFlag_ > 2)) {
        printf("delta G = %g J/kmol \n",  deltaG);
    }
    voltsBase = - deltaG / Faraday;

    hf = FeS2_B->Hf298SS(0);
    if (enableExtraPrinting_ && (detailedResidPrintFlag_ > 2)) {
        printf("    Modified hf for FeS2_B = %g\n", hf);
    }

    deltaG = 1.0/4.0 * um_FeS2_B - 1.0/4.0 * um_FeS2_A - um_Li_liq;


    voltsBase = - deltaG / Faraday;
    if (enableExtraPrinting_ && (detailedResidPrintFlag_ > 2)) {
        printf("delta G = %g J/kmol \n",  deltaG);
        printf("               in volts = %g\n", voltsBase);
        printf("um_li = %g\n", um_Li_liq);
    }
    // setState_TP(tempO, pres0);
    temperature_ = tempO;
    pressure_ = presO;
    Electrode::updateState_OnionOut();
}
//====================================================================================================================
//  Change the Heat of formation of compound B in the reaction in order to generate a given open circuit standard state
//  voltage at the current temperature
/*
 *  @param E0  Value of the standard state OCV that is desired
 */
void Electrode_MP_RxnExtent::changeToE0(double E0, int doPrint)
{
    static int iprint = 3;
    if (!enableExtraPrinting_ || (detailedResidPrintFlag_ <= 1)) {
        iprint = 0;
    }
    if (doPrint) {
        iprint = 1;
    }
    if (iprint > 0) {
        printf("  changeToE0() -------------------------\n");
    }
    int ip_FeS2_A = globalPhaseIndex("FeS2_A(S)");
    int ip_FeS2_B = globalPhaseIndex("FeS2_B(S)");
    ThermoPhase* FeS2_A = &(thermo(ip_FeS2_A));
    ThermoPhase* FeS2_B = &(thermo(ip_FeS2_B));
    Li_liq_->setState_TP(temperature_, pressure_);
    double um[10];
    FeS2_A->getChemPotentials(um);
    double um_FeS2_A = um[0];

    FeS2_B->getChemPotentials(um);
    double um_FeS2_B = um[0];

    Li_liq_->getChemPotentials(um);
    double um_Li_liq = um[0];

    /*
     *  Do more expensive alternate calculation
     */
    ThermoPhase* tElec = & thermo(metalPhase_);
    ThermoPhase* tSalt = & thermo(solnPhase_);

    tElec->getChemPotentials(um);
    double mu_elec = um[0];

    tSalt->getPureGibbs(um);
    int ilt_Lip = tSalt->speciesIndex("Li+");
    double mu_lip = um[ilt_Lip];
    double RT = GasConstant * temperature_;
    // This should be zero
    double delGBase = um_Li_liq -  mu_lip - mu_elec + RT * log(2.0);

    double um_Li_liq_alt = mu_lip + mu_elec - RT * log(2.0);

    /*
     *  This term is needed because the E0 basis for lithium liquid, E-, and Li+ actually
     *  involves this term. It is based on using Li+ mole fraction that is mixed between
     *  the cation and the anion.
     */
    double RTl2 = GasConstant * temperature_ * log(2.0);

    double hf = FeS2_B->Hf298SS(0);
    double deltaG = 1.0/4.0 * um_FeS2_B - 1.0/4.0 * um_FeS2_A - um_Li_liq_alt - RTl2;
    double voltsBase = - deltaG / Faraday;

    if (iprint > 0) {
        printf("\t Current delta G = %22.16E J/kmol \n",  deltaG);
        printf("\t Current volts = %22.16E volts \n", voltsBase);
        printf("\t current Hf_B  = %22.16E J/kmol \n", hf);
    }
    double deltaG_plat = - E0 * Faraday;



    double deltaH = deltaG_plat - deltaG;
    Hf_B_current_ = hf + 4*deltaH;


    FeS2_B->modifyOneHf298SS(0, Hf_B_current_);


    if (iprint > 0) {
        hf = FeS2_B->Hf298SS(0);
        FeS2_B->getChemPotentials(um);
        um_FeS2_B = um[0];
        deltaG = 1.0/4.0 * um_FeS2_B - 1.0/4.0 * um_FeS2_A - um_Li_liq_alt - RTl2;
        double voltsg = - deltaG / Faraday;
        printf("\t New delta G = %22.16E J/kmol \n",  deltaG);
        printf("\t New volts = %22.16E volts \n", voltsg);
        printf("\t New Hf_B  = %22.16E J/kmol \n", hf);

    }
    /*
     *  This is  a pro bono check that  the Li = Li+ + e- reaction is valid and equal to E0 = 0.
     */
    if (iprint > 0) {

        ReactingSurDomain* rsd = reactingSurface(0);
        double deltaG[10];
        rsd->getDeltaSSGibbs(DATA_PTR(deltaG));

        printf("deltaG_kinetics         = %21.16E\n",  deltaG[0]);

        int nStoichElectrons = -1;
        double PhiRxn = deltaG_[0]/Faraday/ nStoichElectrons;
        printf("PhiRxn = %22.16E \n", PhiRxn);
        printf("E0     = %22.16E \n", E0);

        ThermoPhase* tElec = & thermo(metalPhase_);
        ThermoPhase* tSalt = & thermo(solnPhase_);

        tElec->getChemPotentials(um);
        double mu_elec = um[0];

        tSalt->getPureGibbs(um);
        int ilt_Lip = tSalt->speciesIndex("Li+");
        double mu_lip = um[ilt_Lip];
        double RT = GasConstant * temperature_;
        delGBase = um_Li_liq -  mu_lip - mu_elec + RT * log(2.0);

        printf("delGBase = %22.16E -- SHOULD BE ZERO\n", delGBase);

        iprint--;
    }
}
//====================================================================================================================
// Calculate the mole numbers of the solid phases given the relative extent of reaction
/*
 *   Note these operations are not possible, as the mole numbers are reinvented.  I could do this
 *   if I carried the real solid phases around as nonparticipating actors in the calculation. Therefore,
 *   I'll leave this as a hook for later
 */
void  Electrode_MP_RxnExtent::relativeExtentRxnToMoles_final()
{

}
//====================================================================================================================
//     Calculate the mole numbers of the solid phases given the relative extent of reaction
/*
 *  Calculate the relative extent of reaction given the mole numbers of the solid phases.
 *  This is a reversal of the previous operation. The two must be inverses of each other.
 *
 *   Note these operations are not possible, as the mole numbers are reinvented.  I could do this
 *   if I carried the real solid phases around as nonparticipating actors in the calculation. Therefore,
 *   I'll leave this as a hook for later.
 */
void  Electrode_MP_RxnExtent::molesToRelativeExtentRxn_final()
{

}

//====================================================================================================================


//    Return the total volume of solid material
/*
 *  (virtual from Electrode.h)
 *
 *       This is the main routine for calculating the volume of electrode material.
 *       Here we redo the model. The total moles of A + B are used. However,
 *       the molar volume is obtained from the relative extent of reaction within the current plateau
 *       and the mole volumes of the two end points of the current plateau.
 *
 *       We leave out the solnPhase_ volume from the calculation
 *       units = m**3
 */
double Electrode_MP_RxnExtent::SolidVol() const
{
    /*
     *  The total moles are always held constant in the following variable:
     *
     *  spMoles_FeS2_Normalization_  = spMoles_final_[is_FeS2_A] + spMoles_final_[is_FeS2_B];
     */

    //molarVolume_final_ = molarVolume_relExtentRxn(RelativeExtentRxn_final_);


    double solidVol = molarVolume_final_ * spMoles_FeS2_Normalization_;
    return solidVol;
}
//====================================================================================================================
double Electrode_MP_RxnExtent::TotalVol(bool ignoreErrors) const
{

    int iph = solnPhase_;
    double psum = 0.0;
    int kStart = m_PhaseSpeciesStartIndex[iph];
    ThermoPhase& tp = thermo(iph);
    int nspPhase = tp.nSpecies();
    for (int k = 0; k < nspPhase; k++) {
        psum += spMoles_final_[kStart + k] * VolPM_[kStart + k];
    }
    double mv = tp.molarVolume();
    double palt = mv * phaseMoles_final_[iph];
    if (!ignoreErrors) {
        if (palt < -1.0E-15) {
            throw CanteraError(" Electrode_MP_RxnExtent::TotalVol() ",
                               " phase volume is negative " + fp2str(palt));
        }
        if (psum < -1.0E-15) {
            throw CanteraError(" Electrode_MP_RxnExtent::TotalVol() ",
                               " phase volume is negative " + fp2str(psum));
        }
    }
    double denom = palt + psum + 1.0E-9;
    if (tp.eosType() != cLattice) {
        if (!ignoreErrors && 0) {
            if (fabs((palt - psum) / denom) > 1.0E-4) {
                throw CanteraError(" Electrode_MP_RxnExtent::TotalVol() ",
                                   " internal inconsistency " + fp2str(palt) + " " + fp2str(psum));
            }
        }
    }

    double svol = SolidVol();
    double vol = svol + psum;
    return vol;
}
//====================================================================================================================
//  Returns the total moles in the electrode phases of the electrode
/*
 *  @return total moles (kmol)
 */
double  Electrode_MP_RxnExtent::SolidTotalMoles() const
{
    return spMoles_FeS2_Normalization_;
}
//====================================================================================================================
void Electrode_MP_RxnExtent::getPhaseVol(double* const phaseVols) const
{



    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        phaseVols[iph] = 0.0;
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        for (int k = 0; k < nspPhase; k++) {
            phaseVols[iph]  += spMoles_final_[kStart + k] * VolPM_[kStart + k];
        }
    }


    int ip_FeS2_A = globalPhaseIndex("FeS2_A(S)");
    int ip_FeS2_B = globalPhaseIndex("FeS2_B(S)");

    int numRegions = RegionBoundaries_ExtentRxn_.size() - 1;
    int ireg = findRegion(RelativeExtentRxn_final_);
    if (ireg < 0) {
        ireg = 0;
    } else if (ireg > numRegions - 1) {
        ireg = numRegions - 1;
    }
    double deltaReg = RegionBoundaries_ExtentRxn_[ireg+1] - RegionBoundaries_ExtentRxn_[ireg];
    double relRegion = (RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[ireg]) / deltaReg;
    // double molarVol = molarVolumeRegions_[ireg] * relRegion + molarVolumeRegions_[ireg+1] * (1.0 - relRegion);

    phaseVols[ip_FeS2_A] =  molarVolumeRegions_[ireg] * relRegion * spMoles_FeS2_Normalization_;
    phaseVols[ip_FeS2_B] =  molarVolumeRegions_[ireg] * (1.0 - relRegion) * spMoles_FeS2_Normalization_;
}
//====================================================================================================================
// Take the state (i.e., the final state) within the Electrode_Model and push it down
// to the ThermoPhase Objects
/*
 *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
 *  objects in the electrode.  This routine synchronizes all state information to make it consistent.
 */
void Electrode_MP_RxnExtent::updateState()
{
    /*
     *  Update the Electrode base object properties.
     *   We update the solid volume and external radius here as well
     *        Bring the updateState() calc into this routine and examine every aspect of it
     */
    // Electrode_MP_RxnExtent::updateState();
    /*
     * This may be redundant. However, I want to make sure mole fractions are
     * consistent with final moles.
     */
    ReactingSurDomain* rsd = RSD_List_[indexOfReactingSurface_];
    betaF_AB_ = rsd->electrochem_beta(0);
    betaR_AB_ = 1.0 - betaF_AB_;

    /*
     *  Calculate the molar volume from the relative extent of reaction
     */
    molarVolume_final_ = molarVolume_relExtentRxn(RelativeExtentRxn_final_);

    /*
     * Loop over all phases in the object
     */
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        Electrode::updateState_Phase(iph);
    }


    /*
     *  Calculate the real mole numbers given the relative extent of reaction
     *    -> This is an umimplemented hook for the moment.
     */
    relativeExtentRxnToMoles_final();

    /*
     *  Calculate the OCV to be expected given the temperature, the xRegion value, and the RelativeExtentRxn_final_ value
     */
    volts_OCV_SS_final_ = openCircuitVoltageSS_final();


    volts_OCV_final_ = openCircuitVoltage_Region(RelativeExtentRxn_final_, xRegion_final_);


    /*
     *  Change the open circuit standard state voltage by changing the Heat of Formation of B
     */
    changeToE0(volts_OCV_SS_final_);



    /*
     * Calculate the volume of the electrode phase. This is the main routine to do this.
     */
    ElectrodeSolidVolume_ = SolidVol();

    /*
     *  Calculate the volume of a single particle
     */
    double vol =  ElectrodeSolidVolume_ / particleNumberToFollow_;

    /*
     *  Calculate the external radius of the particle (in meters) assuming that all particles are the same
     *  size and are spherical
     */
    Radius_exterior_final_ = pow(vol * 3.0 / (4.0 * Pi), 0.3333333333333333);

    /*
     * Calculate the internal radius of the final state
     */
    Radius_internal_final_ = calculateRadiusInner(RelativeExtentRxn_final_);

    /*
     *  Calculate the surface area of reacting surfaces
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
 *    A value of locationOfReactingSurface_ = 2 indicates that the surface 0 follows the exterior surface of the particle
 *
 *    We also assume that the surface area is equal to the particle surface area multiplied by the numbers of particles.
 *
 *
 *    Dependent StateVariables Used
 *         Radius_exterior_final_;
 *         particleNumberToFollow_
 *
 *    Dependent StateVariables Calculated
 *          surfaceAreaRS_final_[]
 *          Radius_internal_final_
 */
void Electrode_MP_RxnExtent::updateSurfaceAreas()
{
    Radius_internal_final_ = calculateRadiusInner(RelativeExtentRxn_final_);

    double totalSA = 4. * Pi * Radius_internal_final_ * Radius_internal_final_ * particleNumberToFollow_;
    surfaceAreaRS_final_[0] = totalSA;

    totalSA = 4. * Pi * Radius_exterior_final_ * Radius_exterior_final_ * particleNumberToFollow_;
    surfaceAreaRS_final_[1] = totalSA;
}
//====================================================================================================================
//  Extract the ROP of the two reaction fronts from Cantera within this routine
/*
 *  In this routine we calculate the rates of progress from the two surfaces
 *  The vectors are filled in:
 *
 *        ROP_outer_[jRxn]
 *        ROP_inner_[jRxn]
 *        spNetProdPerArea_List_[isk][kIndexKin]
 *        justBornPhase_[jph]
 */
void Electrode_MP_RxnExtent::extractInfo()
{

    int maxNumRxns = 2;
    std::vector<double> netROP(maxNumRxns, 0.0);
    ReactingSurDomain* rsd = RSD_List_[indexOfReactingSurface_];

    int numRxn = rsd->nReactions();
    /*
     *  Set the reaction rate multiplier for the current region
     */
    for (int i = 0; i < numRxn; i++) {

        rsd->setMultiplier(i, rxnPerturbRegions_[xRegion_final_]);
    }
    /*
     *  There is one surface phase and it is active always.
     */


    /*
     *  For each Reacting surface
     *
     *  Get the species production rates for the reacting surface
     */
    const vector<double>& rsSpeciesProductionRates = rsd->calcNetSurfaceProductionRateDensities();

    double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(indexOfReactingSurface_);

    /*
     *  Calculate an effective Damkoeler number for solid phase diffusion and then
     *  calculate a factor that reduces the overall rate of progress of the reaction
     */


    rsd->getNetRatesOfProgress(DATA_PTR(netROP));
    if (goNowhere_) {
        ROP_[0] = 0.0;
    } else {
        ROP_[0] = netROP[0];
    }
    double fac = 1.0;
    if (solidDiffusionModel_) {
        fac = modifyROPForDiffusion();
    }

    /*
     *  loop over the phases in the reacting surface
     *  Get the net production vector
     */
    std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
    int nphRS = rsd->nPhases();
    int jph, kph;
    int kIndexKin = 0;
    for (kph = 0; kph < nphRS; kph++) {
        jph = rsd->kinOrder[kph];
        int istart = m_PhaseSpeciesStartIndex[jph];
        int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
        for (int k = 0; k < nsp; k++) {
            spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin] * fac;
            kIndexKin++;
        }
    }
}
//====================================================================================================================
double Electrode_MP_RxnExtent::modifyROPForDiffusion()
{

    /*
     * Calculate common quantities
     */
    calculatekABForDa();
    double fac = 1.0;
    if (ROPModificationType_ == 1) {
        // Internal to the particle:
        calculateDaInner();
        fac = 1.0 / (1.0 + DaInner_);
    } else if (ROPModificationType_ == 2) {
        // exterior of the particle
        calculateDaOuter();
        fac = 1.0 / (1.0 + DaOuter_);
    }
    ROP_[0] *= fac;
    ROP_AB_ *= fac;
    return fac;
}
//====================================================================================================================
void Electrode_MP_RxnExtent::calculatekABForDa()
{
    ReactingSurDomain* rsd = RSD_List_[indexOfReactingSurface_];
    double kRates_[10];
    double kRRates_[10];
    double ca[10];
    kf_id_ = 0;
    kf_dir_ = -1;
    // -> could start with ropf and ropr. But, you end up at the same place.
    //   double ropf_[10], ropr_[10];
    // rsd->getFwdRatesOfProgress(ropf_);
    // rsd->getRevRatesOfProgress(ropr_);
    if (kf_dir_ > 0) {
        rsd->getFwdRateConstants(kRates_);
        rsd->getRevRateConstants(kRRates_);
        ROP_AB_ = ROP_[kf_id_];
    } else if (kf_dir_ < 0) {
        rsd->getRevRateConstants(kRates_);
        rsd->getFwdRateConstants(kRRates_);
        ROP_AB_ = -ROP_[kf_id_];
    }



    kfAB_ = kRates_[kf_id_];
    krAB_ = kRRates_[kf_id_];

    ThermoPhase& mthermo = thermo(metalPhase_);

    mthermo.getActivityConcentrations(ca);
    //double ca_electron_ = ca[0];

    ThermoPhase& sthermo = thermo(solnPhase_);

    sthermo.getActivityConcentrations(ca);
    int ilt_Lip = sthermo.speciesIndex("Li+");
    ca_Lip_ = ca[ilt_Lip];

    /*
     * Get the open circuit voltage
     */
    Eocv_ = openCircuitVoltage_Region(RelativeExtentRxn_final_, xRegion_final_);
}
//====================================================================================================================
//    Damkoeler number calculated assuming the reaction occurs on the outer surface, and we have interstial diffusion
//    of a neutral diffusing through a region to a new reaction front, whose reactions are fast so that they are
//    in equilibrium
/*
 *  We also can do checks here to make sure that the approximation makes sense. In particular we should
 *  make sure that the mole fraction is between 0 and 1 for the interstitial diffuser.
 */
void Electrode_MP_RxnExtent::calculateDaOuter()
{

    kfExt_ = 1.0E3  * kfAB_;
    krExt_ = krAB_;
    kfInner_ = 1.0E5 * kfExt_;
    krInner_ = 1.0E2 * kfExt_;
    Lin_ = krInner_ / kfInner_;

#ifdef DEBUG_MODE_NOT
    if (counterNumberSubIntegrations_ >= 50497) {
        printf(" calculateDaOuter:\n");
        printf(" thetaR = (ca_Lip_ * krExt_ - kfExt_ * Lin_ ) = %g\n", (ca_Lip_ * krExt_ - kfExt_ * Lin_));
        printf("\t\t krAB_ = %g\n",  kfAB_);
        printf("\t\t ca_Lip_ = %g\n",  ca_Lip_);
        printf("\t\t krExt_ = %g\n",  krExt_);
        printf("\t\t kfExt_ = %g\n",  kfExt_);
        printf("\t\t Lin_ = %g\n",  Lin_);
        printf("\t\t ROP       = %g\n", ROP_[0]);
//	printf("\t\t kRates_[0] = %g\n",kRates_[0]);
        //printf("\t\t kRRates_[0] = %g\n",kRRates_[0]);
    }
#endif

    /*
     *   We increment the region flag by one here to get to the diffusion coefficient and molar concentration because
     *   the inner plateau (with id 0) doesn't actually contribute to the physical problem here.
     */
    double rExt = 0.5 * (Radius_exterior_final_ + Radius_exterior_init_);
    double rInt = 0.5 * (Radius_internal_final_ + Radius_internal_init_);
    DaOuter_Bar_ = (kfExt_ * molarVolumeRegions_[xRegion_final_+1] * (rExt - rInt) * (rExt)) /
                   (diffusionCoeffRegions_[xRegion_final_+1]);

    DaOuter_ = DaOuter_Bar_ / rInt;


    double Kext = kfExt_ / krExt_;

    Lout_ = (DaOuter_ * ca_Lip_ / Kext + Lin_) / (1.0 + DaOuter_);


}
//====================================================================================================================
//    Damkoeler number calculated assuming the reaction occurs on the inner surface, and we have interstial diffusion
//    of a neutral diffusing through a region to a new reaction front, whose reactions are slow so that they are
//    not in equilibrium. The reaction at the edge of the particle, which involves charge transfer is assumed to be fast.
/*
 *  We also can do checks here to make sure that the approximation makes sense. In particular we should
 *  make sure that the mole fraction is between 0 and 1 for the interstitial diffuser.
 */
void Electrode_MP_RxnExtent::calculateDaInner()
{
    double nu = deltaVoltage_ - Eocv_;
    if (innerDaTreatmentType_ == 1) {
        kfInner_ = 1.0E3 * krAB_;
        krInner_ = kfAB_;

        kfExt_ = 1.0E3 * kfInner_;
        krExt_ = 5.0 * 1.0E-3 * kfExt_;

        Lout_ = krExt_ * ca_Lip_ / kfExt_ ;

    } else {

        double tmp =  Faraday * nu / (GasConstant * temperature_);

        kfInner_ = 1.0E3 * krAB_ * exp(betaR_AB_ * tmp);
        krInner_ = kfAB_ * exp(- betaF_AB_ * tmp);

        kfExt_ = 1.0E3 * kfInner_;

        tmp = -betaF_AB_ * Faraday * nu / (GasConstant * temperature_);
        krExt_ = 5.0 * 1.0E3 * krAB_ * exp(tmp);

        Lout_ = krExt_ * ca_Lip_ / kfExt_ ;
    }


    /*
     *   We increment the region flag by one here to get to the diffusion coefficient and molar concentration because
     *   the inner plateau (with id 0) doesn't actually contribute to the physical problem here.
     */
    double rExt = 0.5 * (Radius_exterior_final_ + Radius_exterior_init_);
    double rInt = 0.5 * (Radius_internal_final_ + Radius_internal_init_);
    DaInner_ = (kfInner_ * molarVolumeRegions_[xRegion_final_+1] * (rExt - rInt) * (rInt)) /
               (diffusionCoeffRegions_[xRegion_final_+1]  * rExt);


    double Kinner = kfInner_ / krInner_;

    Lin_ = (DaInner_  / Kinner + Lout_) / (1.0 + DaInner_);

}
//====================================================================================================================
// Update the molar production rates for all annular regions
/*
 *  Requires that SolidInnerKSpecies_IS_, SolidInnerKSpeciesStoichCoeff_IS_, and the like
 *  and ROP_inner_ and ROP_outer_ are already known
 */
void Electrode_MP_RxnExtent::updateMoleRatesFinal()
{
    // Think I don't need this
}
//====================================================================================================================
// Calculate the production rate of species in the electrode at the final time of the time step
/*
 *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
 *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
 *   all species in the electrode.
 */
void Electrode_MP_RxnExtent::updateSpeciesMoleChangeFinal()
{

    std::fill(DspMoles_final_.begin(), DspMoles_final_.end(), 0.);

    if (limitingEquationBehavior_ == 0) {
        if (locationOfReactingSurface_== 0) {
            double mult = (surfaceAreaRS_init_[2] + surfaceAreaRS_final_[2]);
            mult /= 2.0;
            for (size_t i = 0; i < m_totNumVolSpecies; i++) {
                for (int j = 0; j < numRxns_[indexOfReactingSurface_]; j++) {
                    DspMoles_final_[i] +=  mult * productStoichCoeff(indexOfReactingSurface_,i,j) * ROP_[j] ;
                    DspMoles_final_[i] -=  mult * reactantStoichCoeff(indexOfReactingSurface_,i,j) * ROP_[j];
                }
            }
        } else  if (locationOfReactingSurface_== 1) {
            double mult = surfaceAreaRS_init_[0] + 2.0 * surfaceAreaRS_final_[0];
            mult /= 3.0;
            for (size_t i = 0; i < m_totNumVolSpecies; i++) {
                for (int j = 0; j < numRxns_[indexOfReactingSurface_]; j++) {
                    DspMoles_final_[i] +=  mult * productStoichCoeff(indexOfReactingSurface_,i,j) * ROP_[j] ;
                    DspMoles_final_[i] -=  mult * reactantStoichCoeff(indexOfReactingSurface_,i,j) * ROP_[j];
                }
            }
        } else  if (locationOfReactingSurface_== 2) {
            double mult = surfaceAreaRS_init_[1] + surfaceAreaRS_final_[1];
            mult /= 2.0;
            for (size_t i = 0; i < m_totNumVolSpecies; i++) {
                for (int j = 0; j < numRxns_[indexOfReactingSurface_]; j++) {
                    DspMoles_final_[i] +=  mult * productStoichCoeff(indexOfReactingSurface_,i,j)  * ROP_[j];
                    DspMoles_final_[i] -=  mult * reactantStoichCoeff(indexOfReactingSurface_,i,j) * ROP_[j];
                }
            }
        }
        /*
         * We define the srcdot for the extent of reaction as the net electron loss out of the electrode. Thus,
         * we put a negative sign here.
         */
        SrcDot_ExtentRxn_final_  = - DspMoles_final_[kElectron_];
        SrcDot_RelativeExtentRxn_final_  = SrcDot_ExtentRxn_final_ / spMoles_FeS2_Normalization_;

    } else {
        /*
         * special case integration
         */
        if (goNowhere_ == 1) {
            SrcDot_ExtentRxn_final_ = 0.0;
            SrcDot_RelativeExtentRxn_final_ = 0.0;
            for (size_t i = 0; i < m_totNumVolSpecies; i++) {
                DspMoles_final_[i] = 0.0;
            }
        } else {
            double rInt = 0.5 * (Radius_internal_final_ + Radius_internal_init_);
            double theta1 = Radius_internal_init_ * surfaceAreaRS_final_[1] / (rInt + DaOuter_Bar_) / spMoles_FeS2_Normalization_;
            theta1 *= (ca_Lip_ * krExt_ - kfExt_ * Lin_);
#ifdef DEBUG_MODE
            ZZCantera::checkFinite(theta1);
#endif
            double extentLeft1 = 0.0;
            double extentLeft0 = RegionBoundaries_ExtentRxn_[xRegion_final_+1] - RelativeExtentRxn_init_;
            if (extentLeft0 < 1.0E-200) {
                extentLeft1 = 0.0;
            } else {
                double extentLeft0_13 = pow(extentLeft0, 0.3333333333333);
                double extentLeft0_23 = extentLeft0_13 * extentLeft0_13;
                if (theta1 > 0.0) {
                    deltaTdeath_ = extentLeft0 * 3.0 / 2.0 / theta1;
                } else {
                    deltaTdeath_ = 10.0 * deltaTsubcycle_;
                }
                double sgg = 1.0;
                double inside = extentLeft0_23 - 2. * theta1 * deltaTsubcycleCalc_ / (3.00 * extentLeft0_13);
                if (inside < 0.0) {
                    sgg = -1.0;
                    inside = -inside;
                }


                extentLeft1 = sgg * pow(inside, 1.5);
            }

            RelativeExtentRxn_tmp_ = RegionBoundaries_ExtentRxn_[xRegion_final_+1] - extentLeft1;

            SrcDot_ExtentRxn_final_ = (RelativeExtentRxn_tmp_ -  RelativeExtentRxn_init_) /  deltaTsubcycleCalc_ * spMoles_FeS2_Normalization_;

            double fac =  SrcDot_ExtentRxn_final_;

            for (size_t i = 0; i < m_totNumVolSpecies; i++) {
                for (int j = 0; j < numRxns_[indexOfReactingSurface_]; j++) {
                    DspMoles_final_[i] +=  fac * productStoichCoeff(indexOfReactingSurface_,i, j);
                    DspMoles_final_[i] -=  fac * reactantStoichCoeff(indexOfReactingSurface_,i,j);
                }
            }

            if (fabs(SrcDot_ExtentRxn_final_  + DspMoles_final_[kElectron_]) > 1.0E-8 * (fabs(SrcDot_ExtentRxn_final_) + 1.0E-22)) {
                throw CanteraError("", "Logic error");
            }
            SrcDot_RelativeExtentRxn_final_  = SrcDot_ExtentRxn_final_ / spMoles_FeS2_Normalization_;

        }
#ifdef DEBUG_MODE
        ZZCantera::checkFinite(SrcDot_RelativeExtentRxn_final_);
#endif
    }
}

//====================================================================================================================
//
/*
 *  @return   Returns the success of the operation
 *                 1  A predicted solution is achieved
 *                 2  A predicted solution with a multispecies phase pop is acheived
 *                 0  A predicted solution is not achieved, but go ahead anyway
 *                -1  The predictor suggests that the time step be reduced and a retry occur.
 */
int Electrode_MP_RxnExtent::predictSoln_SpeciaType1()
{


    double rInt = 0.5 * (Radius_internal_final_ + Radius_internal_init_);
    double theta1 = Radius_internal_init_ * surfaceAreaRS_final_[1] / (Radius_internal_final_ + DaOuter_Bar_) / spMoles_FeS2_Normalization_;
    double rop = (ca_Lip_ * krExt_ - kfExt_ * Lin_);

    if (goNowhere_) {
        rop = 0.0;
    }
    double rop_orig = rop;
    theta1 *= rop_orig;

    if (rop < 0.0) {
        printf(" we are here at a rop < 0.0 in predictSoln_SpeciaType1 1\n");
    }
    if (goNowhere_) {
        RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
        SrcDot_RelativeExtentRxn_final_ = 0.0;
    } else {

    }

    double chi0 = RegionBoundaries_ExtentRxn_[xRegion_init_ + 1] - RelativeExtentRxn_init_;
    chi0 = fmax(chi0, 0.0);
    double chi0_13 = pow(chi0, 0.333333);
    double chi0_23 = chi0_13 * chi0_13;
    double chi_pred = 0.0;
    if (!goNowhere_) {
        chi_pred = chi0_23 - 2. / 3. * theta1 / chi0_13 *  deltaTsubcycleCalc_;
    }
    if (chi_pred <= 0.0) {
        onRegionBoundary_final_ = xRegion_init_ + 1;
        chi_pred = 0.0;
        //double deltaTsubcycledeath =  chi0_23 * 3. * chi0_13  /  (2. * theta1);
    } else {
        onRegionBoundary_final_ = -1;
    }
    chi_pred = fmax(chi_pred, 0.0);
    chi_pred = pow(chi_pred, 1.5);
    double RelPred = RegionBoundaries_ExtentRxn_[xRegion_init_ + 1] - chi_pred;

    /*
     * calculate an updated inner radius prediction
     */
    double rad_inner_2 = calculateRadiusInner(RelPred);
    Radius_internal_final_ = rad_inner_2;

    extractInfo();

    double rExt = 0.5 * (Radius_exterior_final_ + Radius_exterior_init_);
    rInt = 0.5 * (Radius_internal_init_ + Radius_internal_final_);
    double da_outer_bar_2 = (kfExt_ * molarVolumeRegions_[xRegion_final_+1] * (rExt - rInt) * (rExt)) / (diffusionCoeffRegions_[xRegion_final_+1]);

    double theta2 = Radius_internal_init_ * surfaceAreaRS_final_[1] / (rInt + da_outer_bar_2) / spMoles_FeS2_Normalization_;
    theta2 *= rop_orig;

    if (chi0_13 > 0.0) {
        chi_pred = chi0_23 - 2. / 3. * theta2 / chi0_13 *  deltaTsubcycleCalc_;
    } else {
        chi_pred = 0.0;
    }
    //  Need cutoff above zero here, or else roundoff will cause the an error below
    if (goNowhere_) {
        onRegionBoundary_final_ = xRegion_final_ + 1;
        chi_pred = 0.0;
    } else {
        if (chi_pred <= 1.0E-8) {
            onRegionBoundary_final_ = xRegion_final_ + 1;
            chi_pred = 0.0;
            deltaTsubcycleCalc_ =  chi0_23 * 3. * chi0_13  / (2. * theta2);
        } else {
            onRegionBoundary_final_ = -1;
        }
    }
    chi_pred = fmax(chi_pred, 0.0);
    chi_pred = pow(chi_pred, 1.5);
    double RelPred2 = RegionBoundaries_ExtentRxn_[xRegion_final_ + 1] - chi_pred;


    RelativeExtentRxn_final_ =  RelPred2;

    /*
     *  We just changed the relative extent of reaction. update to possibly get a new open circuit voltage and a radically different
     *  answer.
     */
    updateState();
    /*
     * Calculate ROP_inner and ROP_outer, and justBornPhase_[];
     */
    extractInfo();

    rop = (ca_Lip_ * krExt_ - kfExt_ * Lin_);
    if (goNowhere_) {
        rop = 0.0;
    }
    double theta2_coeff = Radius_internal_init_ * surfaceAreaRS_final_[1] / (rInt + da_outer_bar_2) / spMoles_FeS2_Normalization_;
    theta2 = theta2_coeff * rop;

    double delta_rop = rop - rop_orig;
    if (fabs(delta_rop) > 0.6 * (fabs(rop_orig) + fabs(rop))) {
        double delrel =  deltaTsubcycleCalc_ * theta2_coeff * (delta_rop);
        if (fabs(delrel) > 0.1 *  chi0) {
            /*
             *  we are here when the integration formula has serious problems. Fail the predictor and retry with a smaller time step.
             */
            printf("predictSoln_SpecialType1(): New predictor error condition\n");
            return -1;
        }
    }

    if (onRegionBoundary_final_ == xRegion_final_ + 1) {
        if (rop < 0.0) {
            /*
             *  we are here when the integration formula has serious problems. Fail the predictor and retry with a smaller time step.
             */
            printf("predictSoln_SpecialType1(): New predictor error condition\n");
            return -1;
        }
    }
    if (onRegionBoundary_final_ == xRegion_final_) {
        if (rop > 0.0) {
            /*
             *  we are here when the integration formula has serious problems. Fail the predictor and retry with a smaller time step.
             */
            printf("predictSoln_SpecialType1(): New predictor error condition\n");
            return -1;
        }
    }


    if (RelativeExtentRxn_final_ >=  RegionBoundaries_ExtentRxn_[xRegion_init_ + 1]) {
        if (RelativeExtentRxn_final_ >  RegionBoundaries_ExtentRxn_[xRegion_init_ + 1]) {
            printf(" we are here at an error in predictSoln_SpeciaType1 1\n");
        }
        if (onRegionBoundary_final_ < 0) {
            printf(" we are here at an error in predictSoln_SpeciaType1 2\n");
        }
    }

    /*
     *  If we go backwards by any amount fail the predictor
     */
    if (RelativeExtentRxn_final_ <= RegionBoundaries_ExtentRxn_[xRegion_init_] + 0.3 * chi0) {
        printf(" :predictSoln_SpeciaType1() WARNING: Failing time step because going back too far.\n");
        printf("                         %g  %g %g %g \n", RelativeExtentRxn_final_ , RelativeExtentRxn_init_ ,
               RegionBoundaries_ExtentRxn_[xRegion_init_ + 1],
               RegionBoundaries_ExtentRxn_[xRegion_init_]);
        return  -1;
    }


    updateSurfaceAreas();

    SrcDot_RelativeExtentRxn_final_ = (RelativeExtentRxn_final_ - RelativeExtentRxn_init_) / deltaTsubcycleCalc_;
    SrcDot_ExtentRxn_final_ = SrcDot_RelativeExtentRxn_final_ * spMoles_FeS2_Normalization_;

#ifdef DEBUG_MODE_PREDICTION
    FILE* fp;
    if (counterNumberSubIntegrations_ == 1) {
        fp = fopen("predict.txt", "w");
    } else {
        fp = fopen("predict.txt", "a");
    }
    fprintf(fp, "\n\n");
    fprintf(fp, "predictSoln_SpeciaType1:              counterNumberSubIntegrations_    = %d\n", counterNumberSubIntegrations_);
    fprintf(fp, "                                      RelativeExtentRxn_init_          = %17.9E\n", RelativeExtentRxn_init_);
    fprintf(fp, "                                      RelativeExtentRxn_final_         = %17.9E\n", RelativeExtentRxn_final_);
    fprintf(fp, "                                      theta                            = %17.9E\n", theta2);
    double _rop = (ca_Lip_ * krExt_ - kfExt_ * Lin_);
    double srcDot_ExtentRxn_final = theta1 * spMoles_FeS2_Normalization_;
    fprintf(fp, "                                           SrcDot_ExtentRxn_final          = %17.9E\n",  srcDot_ExtentRxn_final);
    fprintf(fp, "                                           ROP_external                    = %17.9E\n", (ca_Lip_ * krExt_ - kfExt_ * Lin_));
    fprintf(fp, "                                           Rad_internal_init               = %17.9E\n",   Radius_internal_init_);
    fprintf(fp, "                                           surfaceAreaRS_final_[1]         = %17.9E\n",   surfaceAreaRS_final_[1]);
    fprintf(fp, "                                           Denom                           = %17.9E\n", rInt + da_outer_bar_2);
    fprintf(fp, "                                           Radius_internal_final_          = %17.9E\n", rInt);
    fprintf(fp, "                                           Radius_exterior_final_          = %17.9E\n", Radius_exterior_final_);
    fprintf(fp, "                                           deltaR_final_                   = %17.9E\n", rExt - rInt);
    fprintf(fp, "                                           DaOuter_Bar_                    = %17.9E\n", DaOuter_Bar_);
    fprintf(fp, "                                           kfExt_                          = %17.9E\n", kfExt_);
    fprintf(fp, "                                           SrcDot_RelativeExtentExnt_calc  = %17.9E\n",
            _rop * Radius_internal_init_ * surfaceAreaRS_final_[1] / (Radius_internal_final_ + DaOuter_Bar_));

    fprintf(fp, "                                      deltaTsubcycleCalc_              = %17.9E\n", deltaTsubcycleCalc_);
    fprintf(fp, "                                      ExtentLeftInit                   = %17.9E\n", chi0);
    fprintf(fp, "                                      ExtentLeftFinal                  = %17.9E\n", chi_pred);
    fprintf(fp, "                                      onBoundaryFinal                  = %d\n",  onRegionBoundary_final_);
    fclose(fp);

    if (enableExtraPrinting_ && (detailedResidPrintFlag_ > 2)) {
        printf("predictSoln REPORT 1:  deltaTsubcycleCalc_ = %g  RelativeExtentRxn_init_ = %g        RelativeExtentRxn_final_ %g\n",
               deltaTsubcycleCalc_,  RelativeExtentRxn_init_ ,  RelativeExtentRxn_final_);
    }

    if (deltaTsubcycleCalc_ < 1.0E-20) {
        printf("predictSoln(): Small deltaT!\n");
    }


    if (counterNumberSubIntegrations_ >= 0) {
        if (0) {
            printf(" RelativeExtentRxn_final_ = %g\n", 	RelativeExtentRxn_final_);
            printf(" ROP_[0] = %g\n", ROP_[0]);
            printf(" DaOuter_Bar_ = %g\n",  DaOuter_Bar_);
            printf(" DaOuter_Bar_2_ = %g\n",  da_outer_bar_2);
            printf(" DaOuter_ = %g\n",  DaOuter_);
            printf(" DaInner_ = %g\n",  DaInner_);
            printf(" surfaceAreaRS_final_[0] = %g\n", surfaceAreaRS_final_[0]);
            printf(" surfaceAreaRS_final_[1] = %g\n", surfaceAreaRS_final_[1]);
            printf(" Radius_internal_final_ = %g\n",  Radius_internal_final_);
            printf(" Radius_exterior_final_ = %g\n",  Radius_exterior_final_);
            printf(" delta R                = %g\n",  Radius_exterior_final_ - Radius_internal_final_);
            printf(" fac                    = %g\n", 1.0 / (1.0 + DaOuter_));
            printf(" SrcDot_RelativeExtentRxn_final_ = %g\n", SrcDot_RelativeExtentRxn_final_);
        }
        predictSave[0] =   RelativeExtentRxn_final_;
        predictSave[1] =   ROP_[0];
        predictSave[2] =   DaOuter_;
        predictSave[3] =   DaInner_;
        predictSave[4] =   surfaceAreaRS_final_[0];
        predictSave[5] =   surfaceAreaRS_final_[1];
        predictSave[6] =   Radius_internal_final_;
        predictSave[7] =   Radius_exterior_final_;
        predictSave[8] =   1.0 / (1.0 + DaOuter_);
        predictSave[9] =   SrcDot_RelativeExtentRxn_final_ ;
        predictSave[10] = chi_pred;
        predictSave[11] = _rop;
        predictSave[12] = da_outer_bar_2;
        predictSave[13] = deltaTsubcycleCalc_;
        predictSave[14] = onRegionBoundary_final_;
    }


#endif
    return 1;
}
//====================================================================================================================
int  Electrode_MP_RxnExtent::changeRegion(int newRegion)
{
    if (RelativeExtentRxn_init_ != RelativeExtentRxn_final_) {
        throw CanteraError(" Electrode_MP_RxnExtent::changeRegion(int newRegion) ",
                           "function called at the wrong time. extents not the same");
    }
    int tmp =  xRegion_final_;
    xRegion_init_ = newRegion;
    xRegion_final_ = newRegion;
    Radius_internal_final_ = calculateRadiusInner(RelativeExtentRxn_final_);
    Radius_internal_init_ = Radius_internal_final_;
    updateSurfaceAreas();
    surfaceAreaRS_init_[0] = surfaceAreaRS_final_[0];
    surfaceAreaRS_init_[1] = surfaceAreaRS_final_[1];
    return tmp;
}

//==================================================================================================================
// Set the base tolerances for the nonlinear solver within the integrator
/*
 *   The tolerances are based on controlling the integrated electron source term
 *   for the electrode over the integration interval.  The integrated source term
 *   has units of kmol.
 *
 *   Because the electron is only one molar quantity within a bunch of molar quantities,
 *   this requirement will entail that we control the source terms of all species within the
 *   electrode to the tolerance requirements of the electron source term.
 *
 *   @param rtolResid  Relative tolerance allowed for the electron source term over the interval.
 *                     This is a unitless quantity
 *   @param atolResid  Absolute tolerance cutoff for the electron source term over the interval.
 *                     Below this value we do not care about the results.
 *                     atol has units of kmol.
 */
void Electrode_MP_RxnExtent::setNLSGlobalSrcTermTolerances(double rtolResid)
{
    double sum = SolidTotalMoles();
    double val = 1.0E-14 * sum;

    for (int i = 0; i < numIntegratedSrc_; i++) {
        atol_IntegratedSrc_global_[i] = val;
    }
    rtol_IntegratedSrc_global_ = rtolResid;
}
//====================================================================================================================
// Set the Residual absolute error tolerances
/*
 *  (virtual from Electrode_Integrator)
 *
 *   Set the absolute error tolerances fror the nonlinear solvers. This is called at the top
 *   of the integrator() routine.
 *
 *   Calculates atolResidNLS_[]
 *   Calculates atolNLS_[]
 */
void  Electrode_MP_RxnExtent::setResidAtolNLS()
{
    double deltaT = t_final_final_ - t_init_init_;
    deltaT = std::max(deltaT, 1.0E-3);
    /*
     *  We don't care about differences that are 1E-6 less than the global time constant.
     *  residual has the same units as soln.
     */
    atolNLS_[0] = deltaT * 1.0E-6;



    /*
     *  We don't care about differences that are 1E-6 less than the global time constant.
     *  residual has the same units as soln.
     */
    atolResidNLS_[0] = deltaT * 1.0E-6;

    /*
     *  The residual for the relative extent has units of relative extent, which is ~1 quantity. Therefore, set
     *  it to 3 or 4 times the roundoff error, which is the usual treatment.
     */
    atolNLS_[1] = 1.0E-12;
    atolResidNLS_[1] = 1.0E-12;

}
//====================================================================================================================
// Predict the solution
/*
 * Ok at this point we have a time step deltalimiTsubcycle_
 * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
 * We now calculate predicted solution components from these conditions.
 *
 * @return   Returns the success of the operation
 *                 1  A predicted solution is achieved
 *                 2  A predicted solution with a multispecies phase pop is acheived
 *                 0  A predicted solution is not achieved, but go ahead anyway
 *                -1  The predictor suggests that the time step be reduced and a retry occur.
 */
int  Electrode_MP_RxnExtent::predictSoln()
{
    int retn = 1;
    double vtop, vbot;
    double SrcDot_RelativeExtentRxn_final_1;

/*
    bool behaviorChangePossible = false;
    if (solidDiffusionModel_) {
        if (ROPModificationType_ == 2) {
            behaviorChangePossible = true;
        }
    }
*/
    /*
     * Copy initial to final
     */
    copy(surfaceAreaRS_init_.begin(), surfaceAreaRS_init_.end(), surfaceAreaRS_final_.begin());

    // predict that the calculated deltaT is equal to the input deltaT
    deltaTsubcycleCalc_ = deltaTsubcycle_;
#ifdef DEBUG_HKM_NOT
    if (counterNumberSubIntegrations_ >= 55526 && electrodeCellNumber_ == 0) {
        if (counterNumberSubIntegrations_ <= 55540) {
            enableExtraPrinting_ = 10;
            detailedResidPrintFlag_ = 10;
            printLvl_ = 10;
        } else {
            enableExtraPrinting_ = 1;
            detailedResidPrintFlag_ = 0;
            printLvl_ = 3;
        }
    } else {
    }
#endif
#ifdef DEBUG_MODE_NOT
    if (counterNumberSubIntegrations_ >= 1919) {
        enableExtraPrinting_ = 10;
        detailedResidPrintFlag_ = 10;
        printLvl_ = 10;
    }
#endif


    int redoSteps = 0;
    int reDo = 0;
    do {
        /*
         * Copy initial to final
         */
        redoSteps++;
        RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
        copy(spMf_init_.begin(), spMf_init_.end(), spMf_final_.begin());
        copy(spMoles_init_.begin(), spMoles_init_.end(), spMoles_final_.begin());
        copy(phaseMoles_init_.begin(), phaseMoles_init_.end(), phaseMoles_final_.begin());
        molarVolume_final_ = molarVolume_init_;

        /*
         * Handle the case where we start on a Region Boundary, including of course the DoD = 0 boundary
         */
        goNowhere_ = 0;
        if (onRegionBoundary_init_ >= 0 && onRegionBoundary_final_ >= 0) {
            /*
             * Check to see that it is indeed on a boundary
             */
            double relExtentRxn = RegionBoundaries_ExtentRxn_[onRegionBoundary_init_];
            if (fabs(relExtentRxn - RelativeExtentRxn_init_) > 1.0E-9) {
                throw CanteraError("predictSoln():", "ERROR : We are confused\n");
            }

            /*
             * Get the top and bottom voltages
             */
            int nRegions = RegionBoundaries_ExtentRxn_.size() - 1;
            if (onRegionBoundary_init_ == 0) {
                vtop = 10.;
            } else {
                //double vtopSS = openCircuitVoltageSS_Region(relExtentRxn, onRegionBoundary_init_ - 1 );
                vtop = openCircuitVoltage_Region(relExtentRxn, onRegionBoundary_init_ - 1);
            }
            if (onRegionBoundary_init_ == nRegions) {
                vbot =  -10.;
            } else {
                //double vbotSS = openCircuitVoltageSS_Region(relExtentRxn, onRegionBoundary_init_ );
                vbot = openCircuitVoltage_Region(relExtentRxn, onRegionBoundary_init_);
            }
            if (vtop < vbot) {
                throw CanteraError("predictSoln():", "ERROR: Expected but unexpected\n");
            }
            //  If the voltage says we move backwards, we decrement the region before doing the prediction
            //  An extra check is supplied to make sure that we are not at the absolute region 0 boundary, where we can't go backwards
            //  The extra check was necessary because it caused a NaN to appear when deltaVoltage_ is rediculously high
            if ((deltaVoltage_ > (vtop + 1.0E-9)) && (onRegionBoundary_init_ > 0)) {
                changeRegion(onRegionBoundary_init_ - 1);
                onRegionBoundary_final_ = -1;
                goNowhere_ = 0;
            }
            //  If the voltage says we move forwards, we increment the region before doing the prediction
            //  An extra check is supplied to make sure that we are not at the absolute last region boundary, where we can't go forwards
            //  The extra check was necessary because it caused a NaN to appear when deltaVoltage_ is rediculously low
            else if ((deltaVoltage_ < (vbot - 1.0E-9))  && (onRegionBoundary_init_ < numRegions_)) {
                changeRegion(onRegionBoundary_init_);
                onRegionBoundary_final_ = -1;
                goNowhere_ = 0;
            }
            //  If we are at the ends of the calculation or we are in the middle voltage we set a toggle switch to
            //  control the logic
            else {
                goNowhere_ = 1;
                SrcDot_RelativeExtentRxn_final_ = 0.0;
            }
        }

        updateState();
        /*
         * At the start of this step the _final_ values are equal to the _init_ values
         * - we may want to put a check here that this is true. However, for now we will assume it
         *
         * Set the phase existence -> this depends on SpMoles_init_ and SpMoles_final_
         */
        if (redoSteps == 1) {
            setPhaseExistenceForReactingSurfaces(true);
        } else {
            setPhaseExistenceForReactingSurfaces(false);
        }
        /*
         * Calculate ROP_inner and ROP_outer, and justBornPhase_[];
         */
        extractInfo();
        /*
         *  Calculate DspMoles_final_[], the change in moles of each species and  SrcDot_ExtentRxn_final_
         */
        updateSpeciesMoleChangeFinal();
        /*
         *  Check to make sure that if we're starting on a boundary, then we're heading in the right direction
         */
        if (!goNowhere_ && onRegionBoundary_init_ >= 0) {
            /*
             *  Case where we are trying to advance the region
             */
            if (xRegion_final_ >= onRegionBoundary_init_) {
                /*
                 *  If we find ourselves going in the wrong direction, change the calculation to goNowhere
                 */
                if (SrcDot_RelativeExtentRxn_final_ <= 0.0) {
                    printf("New condition with advance in region but negative srcdot caught!\n");
                    goNowhere_ = 1;
                    SrcDot_RelativeExtentRxn_final_ = 0.0;
                    SrcDot_ExtentRxn_final_ = 0.0;
                    onRegionBoundary_final_ = onRegionBoundary_init_;
                }
            }
            /*
             *  Case where we are trying to decrement the region
             */
            if (xRegion_final_ < onRegionBoundary_init_) {
                /*
                 *  If we find ourselves going in the wrong direction, change the calculation to goNowhere
                 */
                if (SrcDot_RelativeExtentRxn_final_ >= 0.0) {
                    printf("New condition with decrement in region and positive srcdot caught!\n");
                    goNowhere_ = 1;
                    SrcDot_RelativeExtentRxn_final_ = 0.0;
                    SrcDot_ExtentRxn_final_ = 0.0;
                    onRegionBoundary_final_ = onRegionBoundary_init_;
                }
            }
            if (goNowhere_) {
                extractInfo();
                updateSpeciesMoleChangeFinal();
            }

        }
	SrcDot_RelativeExtentRxn_final_1 = SrcDot_RelativeExtentRxn_final_;


        double extLeftStart =  RegionBoundaries_ExtentRxn_[xRegion_init_ + 1] - RelativeExtentRxn_init_;

        double relLeftStart = extLeftStart / (RegionBoundaries_ExtentRxn_[xRegion_init_ + 1] - RegionBoundaries_ExtentRxn_[xRegion_init_]);


        limitingEquationBehavior_ = 0;
        if (DaOuter_ > 1.0  && !goNowhere_) {
            if (relLeftStart < 0.5) {
                limitingEquationBehavior_ = 1;
            }
        }




        if (limitingEquationBehavior_ == 1) {

            retn = predictSoln_SpeciaType1();
            //double extLeftFinal =  RegionBoundaries_ExtentRxn_[xRegion_init_ + 1] - RelativeExtentRxn_final_;
            //double relLeftFinal = extLeftFinal / (RegionBoundaries_ExtentRxn_[xRegion_init_ + 1] - RegionBoundaries_ExtentRxn_[xRegion_init_]);
            if (retn != 1) {
                return retn;
            }

            double xBd = RegionBoundaries_ExtentRxn_[xRegion_init_ + 1];
            double deltax;


            if ((SrcDot_ExtentRxn_final_ > 1.0E-200) &&
                    (RelativeExtentRxn_final_ > RegionBoundaries_ExtentRxn_[xRegion_init_+1])) {
                deltax = xBd - RelativeExtentRxn_init_;
                deltaTdeath_ = deltax  * spMoles_FeS2_Normalization_ / SrcDot_ExtentRxn_final_;
                onRegionBoundary_final_ = xRegion_init_ + 1;
                deltaTsubcycleCalc_ =  deltaTdeath_;
                RelativeExtentRxn_final_ = RegionBoundaries_ExtentRxn_[xRegion_init_ + 1];

            } else if ((SrcDot_ExtentRxn_final_ < -1.0E-200) &&
                       (RelativeExtentRxn_final_ < RegionBoundaries_ExtentRxn_[xRegion_init_])) {
                xBd = RegionBoundaries_ExtentRxn_[xRegion_init_];
                deltax = xBd - RelativeExtentRxn_init_;
                deltaTdeath_ = deltax * spMoles_FeS2_Normalization_/ SrcDot_ExtentRxn_final_;
                onRegionBoundary_final_ = xRegion_init_;
                deltaTsubcycleCalc_ =  deltaTdeath_;
                RelativeExtentRxn_final_  = RegionBoundaries_ExtentRxn_[xRegion_init_];
            }

        } else {

            if (goNowhere_) {
                RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
            } else {
                RelativeExtentRxn_final_ = RelativeExtentRxn_init_ + SrcDot_RelativeExtentRxn_final_ * deltaTsubcycleCalc_;
            }
	    SrcDot_RelativeExtentRxn_final_1 = SrcDot_RelativeExtentRxn_final_;

            /*
             *  Calculate the inner radius via one successive approximation to improve the estimate
             *     - HKM - decided to add this code, because it really affects the predicted result at the end of a plateau.
             *            The surface area change is a first order effect.
             */
            if (!goNowhere_) {
                // double tmp = RelativeExtentRxn_final_;

                //  First thing to do is to crop the update to the current plateau
                if ((RelativeExtentRxn_final_ > RegionBoundaries_ExtentRxn_[xRegion_init_+1])) {
                    RelativeExtentRxn_final_ = RegionBoundaries_ExtentRxn_[xRegion_init_ + 1];
                } else if ((RelativeExtentRxn_final_ < RegionBoundaries_ExtentRxn_[xRegion_init_])) {
                    RelativeExtentRxn_final_  = RegionBoundaries_ExtentRxn_[xRegion_init_];
                }

                updateState();

                // Update the surface areas with the new estimate of the final surface area.
                updateSurfaceAreas();

                // needed to update the Da # which depends on the surface area.
                extractInfo();
                // Update the SrcDot_RelativeExtentRxn_final_  value
                updateSpeciesMoleChangeFinal();





#ifdef DEBUG_HKM_NOT
                double tmpp =  RelativeExtentRxn_final_ ;
#endif
                // Calculate a new value of the relative extent.
                RelativeExtentRxn_final_ = RelativeExtentRxn_init_ + SrcDot_RelativeExtentRxn_final_ * deltaTsubcycleCalc_;
#ifdef DEBUG_HKM_NOT
                if (RelativeExtentRxn_final_ < 1.75 && tmpp > 1.75) {
                    printf("we are here\n");
                }
#endif

		if (SrcDot_RelativeExtentRxn_final_1 * SrcDot_RelativeExtentRxn_final_ < 0.0) {
		    if (fabs((SrcDot_RelativeExtentRxn_final_1 - SrcDot_RelativeExtentRxn_final_) * deltaTsubcycleCalc_) > 0.0001) {
			/*
			 *  we are here when the integration formula has serious problems. The source term changed sign
			 *  Fail the step and use a smaller step.
			 */
			// printf("predictSoln(): New predictor error condition\n");
			return -1;
		    }
		} else {
		    // if the source term changed a lot return also.
		    if (fabs((SrcDot_RelativeExtentRxn_final_1 - SrcDot_RelativeExtentRxn_final_) * deltaTsubcycleCalc_) > 0.01) {
			return -1;
		    }
		}

                if (RelativeExtentRxn_final_ >= RegionBoundaries_ExtentRxn_[xRegion_init_ + 1]) {
                    if (RelativeExtentRxn_init_  < RelativeExtentRxn_final_) {
                        if (SrcDot_RelativeExtentRxn_final_ <= 0.0) {
                            /*
                             *  we are here when the integration formula has serious problems. Fail the predictor and retry with a smaller time step.
                             */
                            // printf("predictSoln(): New predictor error condition\n");
                            return -1;
                        }
                    }
                }
                if (RelativeExtentRxn_final_ <= RegionBoundaries_ExtentRxn_[xRegion_init_]) {
                    if (RelativeExtentRxn_init_  > RelativeExtentRxn_final_) {
                        if (SrcDot_RelativeExtentRxn_final_ >= 0.0) {
                            /*
                             *  we are here when the integration formula has serious problems. Fail the predictor and retry with a smaller time step.
                             */
                            // printf("predictSoln(): New predictor error condition\n");
                            return -1;
                        }
                    }
                }




                //  Last thing to do is to crop the update to the current plateau
                if ((RelativeExtentRxn_final_ > RegionBoundaries_ExtentRxn_[xRegion_init_+1])) {
                    RelativeExtentRxn_final_ = RegionBoundaries_ExtentRxn_[xRegion_init_ + 1];
                    onRegionBoundary_final_ = xRegion_init_ + 1;
                    if (SrcDot_RelativeExtentRxn_final_ != 0.0) {
                        deltaTsubcycleCalc_ = (RelativeExtentRxn_final_ - RelativeExtentRxn_init_) / SrcDot_RelativeExtentRxn_final_;
                    }
                } else if ((RelativeExtentRxn_final_ < RegionBoundaries_ExtentRxn_[xRegion_init_])) {
                    RelativeExtentRxn_final_  = RegionBoundaries_ExtentRxn_[xRegion_init_];
                    RelativeExtentRxn_final_ = RegionBoundaries_ExtentRxn_[xRegion_init_];
                    onRegionBoundary_final_ = xRegion_init_;
                    if (SrcDot_RelativeExtentRxn_final_ != 0.0) {
                        deltaTsubcycleCalc_ = (RelativeExtentRxn_final_ - RelativeExtentRxn_init_) / SrcDot_RelativeExtentRxn_final_;

                    }
                }
            }



#if DEBUG_MODE_PREDICTION
            FILE* fp;
            if (counterNumberSubIntegrations_ == 1) {
                fp = fopen("predict.txt", "w");
            } else {
                fp = fopen("predict.txt", "a");
            }
            fprintf(fp, "\n\n");
            fprintf(fp, "predictSoln_Normal:                  counterNumberSubIntegrations_    = %d\n", counterNumberSubIntegrations_);
            fprintf(fp, "                                      RelativeExtentRxn_init_          = %17.9E\n", RelativeExtentRxn_init_);
            fprintf(fp, "                                      RelativeExtentRxn_final_         = %17.9E\n", RelativeExtentRxn_final_);
            fprintf(fp, "                                      SrcDot_RelativeExtentRxn_final_      = %17.9E\n", SrcDot_RelativeExtentRxn_final_);
            fprintf(fp, "                                      deltaTsubcycleCalc_              = %17.9E\n", deltaTsubcycleCalc_);

            fprintf(fp, "                                           Radius_internal_final_ = %17.9E\n",  Radius_internal_final_);
            fprintf(fp, "                                           Radius_exterior_final_ = %17.9E\n",  Radius_exterior_final_);
            fclose(fp);

            if (counterNumberSubIntegrations_ >= 0) {
                if (s_printLvl_PREDICTOR_CORRECTOR) {
                    printf(" RelativeExtentRxn_final_ = %g\n", 	RelativeExtentRxn_final_);
                    printf(" ROP_[0] = %g\n", ROP_[0]);
                    printf(" DaOuter_ = %g\n",  DaOuter_);
                    printf(" DaInner_ = %g\n",  DaInner_);
                    printf(" surfaceAreaRS_final_[0] = %g\n", surfaceAreaRS_final_[0]);
                    printf(" surfaceAreaRS_final_[1] = %g\n", surfaceAreaRS_final_[1]);
                    printf(" Radius_internal_final_ = %g\n",  Radius_internal_final_);
                    printf(" Radius_exterior_final_ = %g\n",  Radius_exterior_final_);
                    printf(" delta R                = %g\n",  Radius_exterior_final_ - Radius_internal_final_);
                    printf(" fac                    = %g\n", 1.0 / (1.0 + DaOuter_));
                    printf(" SrcDot_RelativeExtentRxn_final_ = %g\n", SrcDot_RelativeExtentRxn_final_);
                }
                predictSave[0] =      RelativeExtentRxn_final_;
                predictSave[1] =      ROP_[0];
                predictSave[2] =      DaOuter_;
                predictSave[3] =      DaInner_;
                predictSave[4] =      surfaceAreaRS_final_[0];
                predictSave[5] =      surfaceAreaRS_final_[1];
                predictSave[6] =      Radius_internal_final_;
                predictSave[7] =      Radius_exterior_final_;
                predictSave[8] =      1.0 / (1.0 + DaOuter_);
                predictSave[9] =      SrcDot_RelativeExtentRxn_final_ ;

                predictSave[13] = deltaTsubcycleCalc_;
                predictSave[14] = onRegionBoundary_final_;

            }
#endif

        }

#ifdef DEBUG_MODE
        if (enableExtraPrinting_ && (detailedResidPrintFlag_ > 2)) {
            printf("predictSoln REPORT 2: SubIntegrationCounter = %d t_init = %g t_final = %g deltaTsubcycle = %g\n",
                   counterNumberSubIntegrations_, tinit_,  tfinal_, deltaTsubcycle_);

            printf("                       deltaTsubcycleCalc_ = %g  RelativeExtentRxn_init_ = %17.9E  RelativeExtentRxn_final_ %17.9E,"
                   "onboundary_final = %d\n",
                   deltaTsubcycleCalc_,  RelativeExtentRxn_init_ , RelativeExtentRxn_final_ , onRegionBoundary_final_);
            if (deltaTsubcycleCalc_ < 1.0E-20) {
                printf("predictSoln(): Small deltaT %g\n", deltaTsubcycleCalc_);
            }
        }
#endif

        /*
         *  Make sure we are bounded
         *  Estimate spMoles_final_ so that we can get the final surface area
         */

        if (!goNowhere_) {
            //RelativeExtentRxn_final_ = RelativeExtentRxn_init_ + SrcDot_ExtentRxn_final_ * deltaTsubcycleCalc_ / spMoles_FeS2_Normalization_;

            double deltax =  RelativeExtentRxn_final_ -  RelativeExtentRxn_init_;
            if (xRegion_init_ == 2 && (fabs(deltax) > 0.08)) {
                return -1;
            }

            for (int ph = 0; ph < (int) phaseIndexSolidPhases_.size(); ph++) {
                int iph = phaseIndexSolidPhases_[ph];
                int nsp = numSpecInSolidPhases_[ph];
                for (int sp = 0; sp < nsp; sp++) {
                    int isp = getGlobalSpeciesIndex(iph,sp);
                    spMoles_final_[isp] = spMoles_init_[isp] + DspMoles_final_[isp] * deltaTsubcycleCalc_;
                    if (spMoles_final_[isp] < 0.0) {
                        spMoles_final_[isp] = 0.0;
                    }
                }
            }
        }

    } while (reDo);
#ifdef DEBUG_MODE
    ZZCantera::checkFinite(SrcDot_RelativeExtentRxn_final_);

    if (RelativeExtentRxn_final_ < 1.73 && RelativeExtentRxn_init_ > 1.76) {
        printf("we are here\n");
    }
#endif
    if (deltaTsubcycleCalc_ <= 0.0) {
	throw CanteraError("predictSoln():", "deltaTsubcycleCalc is less than or equal to zero");
    }
    soln_predict_[0] = deltaTsubcycleCalc_ ;
    soln_predict_[1] = RelativeExtentRxn_final_;
    soln_predict_[2] =   onRegionBoundary_final_;
    return 1;
}
//====================================================================================================================
//   Unpack the soln vector
/*
 *  This function unpacks the solution vector into deltaTsubcycleCalc_  and   RelativeExtentRxn_final_
 */
void  Electrode_MP_RxnExtent::unpackNonlinSolnVector(const double* const y)
{
    deltaTsubcycleCalc_ = y[0];
    tfinal_ = tinit_ + deltaTsubcycleCalc_;
    RelativeExtentRxn_final_ = y[1];
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
 *   resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
 *   resid[1] = RxnExtent_final_ - ( RxnExtent_init_ + Icurr * sa * (t_final_ - t_init_) * DiffResistance)
 *
 *  The formulation of the solution vector is as follows. The solution vector will consist of the following form
 *
 *     y =   deltaTsubcycleCalc_
 *          RxnExtent_final_
 *
 *  @param resid    value of the residual
 *
 *  @return  1 Means a good calculation that produces a valid result
 *           0 Bad calculation that means that the current nonlinear iteration should be terminated
 */
int Electrode_MP_RxnExtent::calcResid(double* const resid, const ResidEval_Type_Enum evalType)
{


    resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
    double xBd = RegionBoundaries_ExtentRxn_[xRegion_init_ + 1];
    deltaTdeath_ = 1.0E300;
    double deltax;

    if (limitingEquationBehavior_ == 0) {
        if (SrcDot_ExtentRxn_final_ > 1.0E-200) {
            deltax = xBd - RelativeExtentRxn_init_;
            deltaTdeath_ = deltax  * spMoles_FeS2_Normalization_ / SrcDot_ExtentRxn_final_;
        } else  if (SrcDot_ExtentRxn_final_ < -1.0E-200) {
            xBd = RegionBoundaries_ExtentRxn_[xRegion_init_ ];
            deltax = xBd - RelativeExtentRxn_init_;
            deltaTdeath_ = deltax * spMoles_FeS2_Normalization_/ SrcDot_ExtentRxn_final_;
        } else {
            deltax = xBd - RelativeExtentRxn_init_;
            deltaTdeath_ = 1000. * deltaTsubcycle_;
        }

        if (evalType != JacDelta_ResidEval) {
            if (deltaTdeath_ < deltaTsubcycle_) {
                if (SrcDot_ExtentRxn_final_ > 1.0E-200) {
                    if (onRegionBoundary_final_ != -1) {
                        onRegionBoundary_final_ = xRegion_final_+1;
                    } else {
                        onRegionBoundary_final_ = -1;
                    }
                } else {
                    if (onRegionBoundary_final_ != -1) {
                        onRegionBoundary_final_ = xRegion_final_;
                    } else {
                        onRegionBoundary_final_ = -1;
                    }
                }
            } else {
                if (SrcDot_ExtentRxn_final_ > 1.0E-200 || SrcDot_ExtentRxn_final_ < -1.0E-200) {
                    onRegionBoundary_final_ = -1;
                }
            }
            if (onRegionBoundary_final_ >= 0) {
                if (onRegionBoundary_final_ == onRegionBoundary_init_) {
                    if (!goNowhere_) {
                        //  resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
                        //resid[1] =  RelativeExtentRxn_final_ -  RelativeExtentRxn_init_;
                        //return 0;
                        onRegionBoundary_final_ = -1;
                    }
                }
            }
        }

        if (onRegionBoundary_final_ >= 0 && !goNowhere_) {
            resid[0] = deltaTsubcycleCalc_ - deltax * spMoles_FeS2_Normalization_  / SrcDot_ExtentRxn_final_;
        } else {
            resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
        }

        RelativeExtentRxn_tmp_ = RelativeExtentRxn_init_
                                 + SrcDot_ExtentRxn_final_ * deltaTsubcycleCalc_ / spMoles_FeS2_Normalization_;

#ifdef DEBUG_MODE_PREDICTION
        if (evalType == Base_ShowSolution) {
            double pred_RelativeExtentRxn_final_   =  predictSave[0];
            double pred_ROP                        =  predictSave[1];
            double pred_DaOuter_                   =  predictSave[2];
            double pred_DaInner_                   =  predictSave[3];
            double pred_surfaceAreaRS0_final_      =  predictSave[4];
            double pred_surfaceAreaRS1_final_      =  predictSave[5];
            double pred_Radius_internal_final_     =  predictSave[6];
            double pred_Radius_exterior_final_     =  predictSave[7];
            double pred_fac                        =  predictSave[8];
            double pred_SrcDot_RelativeExtentRxn_final_  = predictSave[9];
            double pred_deltaTsubcycleCalc_        = predictSave[13];
            int pred_onRegionBoundary_final_       = predictSave[14];

            printf("\t\t -------------------------------------------------------------------------------------------------------------------\n");
            printf("\t\t Predicted vs Calculated for limitingEquationBehavior = %d,  counterNumberSubIntegrations_ = %d\n", limitingEquationBehavior_,
                   counterNumberSubIntegrations_);
            //printf("\t\t      RelativeExtentRxn_init_ = %15.6E\n", RelativeExtentRxn_init_);
            printf("\t\t             DeltaVoltage = %19.12E\n", deltaVoltage_);
            printf("\t\t                                          predicted          actual \n");
            printf("\t\t  RelativeExtentRxn_init_        :    %15.6E        \n", RelativeExtentRxn_init_);
            printf("\t\t  RelativeExtentRxn_final_       :    %15.6E  %15.6E\n", pred_RelativeExtentRxn_final_, RelativeExtentRxn_final_);
            printf("\t\t  ROP_[0]                        :    %15.6E  %15.6E\n", pred_ROP, ROP_[0]);
            printf("\t\t  Radius_internal_final_         :    %15.6E  %15.6E\n", pred_Radius_internal_final_,  Radius_internal_final_);
            printf("\t\t  Radius_exterior_final_         :    %15.6E  %15.6E\n", pred_Radius_exterior_final_,  Radius_exterior_final_);
            printf("\t\t  SurfaceAreaRS0_final_          :    %15.6E  %15.6E\n", pred_surfaceAreaRS0_final_,   surfaceAreaRS_final_[0]);
            printf("\t\t  SurfaceAreaRS1_final_          :    %15.6E  %15.6E\n", pred_surfaceAreaRS1_final_,   surfaceAreaRS_final_[1]);
            printf("\t\t  DaOuter                        :    %15.6E  %15.6E\n", pred_DaOuter_, DaOuter_);
            printf("\t\t  DaInner                        :    %15.6E  %15.6E\n", pred_DaInner_, DaInner_);
            printf("\t\t  fac                            :    %15.6E  %15.6E\n", pred_fac, 1.0 / (1.0 + DaOuter_));
            printf("\t\t  SrcDot_RelativeExtentRxn_final_:    %15.6E  %15.6E\n", pred_SrcDot_RelativeExtentRxn_final_, SrcDot_RelativeExtentRxn_final_);
            printf("\t\t  DeltaTsubcycleCalc_            :    %15.6E  %15.6E\n", pred_deltaTsubcycleCalc_ , deltaTsubcycleCalc_);
            printf("\t\t  Onboundary_init_               :    %15d      \n",     onRegionBoundary_init_);
            printf("\t\t  Onboundary_final_              :    %15d  %15d\n",     pred_onRegionBoundary_final_ , onRegionBoundary_final_);
            printf("\t\t -------------------------------------------------------------------------------------------------------------------\n");
        }
#endif

    } else if (limitingEquationBehavior_ == 1) {
        /*
         * special case integration
         */
        //double rExt = 0.5 * (Radius_exterior_final_ + Radius_exterior_init_);
        double rInt = 0.5 * (Radius_internal_final_ + Radius_internal_init_);
        double theta1 = Radius_internal_init_ * surfaceAreaRS_final_[1] / (rInt + DaOuter_Bar_) / spMoles_FeS2_Normalization_;
        theta1 *= (ca_Lip_ * krExt_ - kfExt_ * Lin_);

        double extentLeft0 = RegionBoundaries_ExtentRxn_[xRegion_final_+1] - RelativeExtentRxn_init_;
        double extentLeft1 = extentLeft0;
        if (goNowhere_) {
            extentLeft1 = extentLeft0;
            RelativeExtentRxn_tmp_ = RegionBoundaries_ExtentRxn_[xRegion_final_+1] - extentLeft1;
            SrcDot_RelativeExtentRxn_final_ = 0.0;
        } else {
            if (extentLeft0 == 0.0) {
                printf("we are here\n");
            }
            double extentLeft0_13 = pow(extentLeft0, 0.3333333333333);
            double extentLeft0_23 = extentLeft0_13 * extentLeft0_13;

            double sgg = 1.0;
            double inside = extentLeft0_23 - 2. * theta1 * deltaTsubcycleCalc_ / (3.00 * extentLeft0_13);
            if (inside < 0.0) {
                sgg = -1.0;
                inside = -inside;
            }

            extentLeft1 = sgg * pow(inside, 1.5);
            RelativeExtentRxn_tmp_ = RegionBoundaries_ExtentRxn_[xRegion_final_+1] - extentLeft1;
            SrcDot_RelativeExtentRxn_final_ = (RelativeExtentRxn_tmp_ -  RelativeExtentRxn_init_) / deltaTsubcycleCalc_;
        }
#ifdef DEBUG_MODE
        checkFinite(SrcDot_RelativeExtentRxn_final_);
#endif
        if (theta1 > 0.0) {
            deltaTdeath_ = extentLeft0 * 3.0 / 2.0 / theta1;
        } else if (theta1 < 0.0) {
            if (onRegionBoundary_final_ == -1 || onRegionBoundary_final_ == xRegion_final_) {
                xBd = RegionBoundaries_ExtentRxn_[xRegion_init_];
                deltax = xBd - RelativeExtentRxn_init_;
                deltaTdeath_ = deltax * spMoles_FeS2_Normalization_/ SrcDot_ExtentRxn_final_;
            } else {
                if (RelativeExtentRxn_final_ <=  RegionBoundaries_ExtentRxn_[xRegion_init_+1]) {
                    if (!goNowhere_) {
                        return 1;
                    }
                }
            }
        } else {
            deltaTdeath_ = 10.0 * deltaTsubcycle_;
        }

        if (RelativeExtentRxn_tmp_ < RelativeExtentRxn_init_) {
            /*
             *  It may happen that we go the other way. This section treats the other way.
             */
            deltaTdeath_ = 1000. * deltaTsubcycle_;
            if (RelativeExtentRxn_tmp_ <= RegionBoundaries_ExtentRxn_[xRegion_final_]) {
                // Copy logic from above
                xBd = RegionBoundaries_ExtentRxn_[xRegion_init_];
                deltax = xBd - RelativeExtentRxn_init_;
                deltaTdeath_ = deltax * spMoles_FeS2_Normalization_/ SrcDot_ExtentRxn_final_;
                /*
                 *   This is the only spot where we may allow a change in equation type.
                 *   However, here we do not allow a change in equation type.
                 */
                if (evalType != JacDelta_ResidEval) {
                    if (onRegionBoundary_final_ >= 0) {
                        if (onRegionBoundary_final_ == xRegion_final_ + 1) {
                            //  we're searching for a calculation that ends up at the forward boundary, and we just ended up at the back boundary.
                            //  Terminate and get a better initial guess (this has happened).
                            return 0;
                        }
                    }
                }

                if (onRegionBoundary_final_ >= 0 && !goNowhere_) {
                    resid[0] = deltaTsubcycleCalc_ - deltax * SrcDot_RelativeExtentRxn_final_;
                } else {
                    resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
                }
            }

        }

        if (onRegionBoundary_final_ >= 0 && !goNowhere_) {
            resid[0] = deltaTsubcycleCalc_ - deltaTdeath_;
        } else {
            resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
        }
        if (goNowhere_) {
            extentLeft1 = 0.0;
        }


        RelativeExtentRxn_tmp_ = RegionBoundaries_ExtentRxn_[xRegion_final_+1] - extentLeft1;

#ifdef DEBUG_MODE_PREDICTION
        double srcDot_ExtentRxn_final = theta1 * spMoles_FeS2_Normalization_;
        double _rop = (ca_Lip_ * krExt_ - kfExt_ * Lin_);
        if (evalType == Base_ShowSolution) {
            FILE* fp = fopen("predict.txt", "a");
            fprintf(fp, "calcResid                             counterNumberSubIntegrations_    = %d\n", counterNumberSubIntegrations_);
            fprintf(fp, "                                      RelativeExtentRxn_init_          = %17.9E\n", RelativeExtentRxn_init_);
            fprintf(fp, "                                      RelativeExtentRxn_final_         = %17.9E\n", RelativeExtentRxn_tmp_);
            fprintf(fp, "                                      theta                            = %17.9E\n", theta1);

            fprintf(fp, "                                           SrcDot_ExtentRxn_final_         = %17.9E\n",  srcDot_ExtentRxn_final);
            fprintf(fp, "                                           ROP                             = %17.9E\n",  _rop);
            fprintf(fp, "                                           Rad_internal_init               = %17.9E\n",   Radius_internal_init_);
            fprintf(fp, "                                           surfaceAreaRS_final_[1]         = %17.9E\n",   surfaceAreaRS_final_[1]);
            fprintf(fp, "                                           Denom                           = %17.9E\n", Radius_internal_final_ + DaOuter_Bar_);
            fprintf(fp, "                                           Radius_internal_final_          = %17.9E\n", Radius_internal_final_);
            fprintf(fp, "                                           Radius_exterior_final_          = %17.9E\n", Radius_exterior_final_);
            fprintf(fp, "                                           deltaR_final_                   = %17.9E\n", Radius_exterior_final_ - Radius_internal_final_);
            fprintf(fp, "                                           DaOuter_Bar_                    = %17.9E\n", DaOuter_Bar_);
            fprintf(fp, "                                           kfExt_                          = %17.9E\n", kfExt_);
            fprintf(fp, "                                           SrcDot_RelativeExtentExnt_calc  = %17.9E\n",
                    _rop * Radius_internal_init_ * surfaceAreaRS_final_[1] / (Radius_internal_final_ + DaOuter_Bar_));

            fprintf(fp, "                                      deltaTsubcycleCalc_              = %17.9E\n", deltaTsubcycleCalc_);
            fprintf(fp, "                                      ExtentLeftInit                   = %17.9E\n", extentLeft0);
            fprintf(fp, "                                      ExtentLeftFinal                  = %17.9E\n", extentLeft1);
            fprintf(fp, "                                      onBoundaryFinal                  = %d\n", onRegionBoundary_final_);
            fclose(fp);
        }
#endif

#ifdef DEBUG_MODE_PREDICTION
        if (evalType == Base_ShowSolution) {
            double pred_RelativeExtentRxn_final_   =  predictSave[0];
            double pred_ROP                        =  predictSave[1];
            double pred_DaOuter_                   =  predictSave[2];
            double pred_DaInner_                   =  predictSave[3];
            double pred_surfaceAreaRS0_final_      =  predictSave[4];
            double pred_surfaceAreaRS1_final_      =  predictSave[5];
            double pred_Radius_internal_final_     =  predictSave[6];
            double pred_Radius_exterior_final_     =  predictSave[7];
            //double pred_fac                        =  predictSave[8];
            double pred_SrcDot_RelativeExtentRxn_final_  = predictSave[9];
            double pred_rop_external               =  predictSave[11];
            double pred_DaOuter_Bar_               = predictSave[12];
            double pred_deltaTsubcycleCalc_        = predictSave[13];
            int pred_onRegionBoundary_final_       = predictSave[14];

            pred_DaOuter_ = pred_DaOuter_Bar_ * 2.0 / (pred_Radius_internal_final_ + Radius_internal_init_);

            printf("\t\t -------------------------------------------------------------------------------------------------------------------\n");
            printf("\t\t Predicted vs Calculated for limitingEquationBehavior = %d,  counterNumberSubIntegrations_ = %d\n", limitingEquationBehavior_,
                   counterNumberSubIntegrations_);
            //printf("\t\t      RelativeExtentRxn_init_ = %15.6E\n", RelativeExtentRxn_init_);
            printf("\t\t             DeltaVoltage = %19.12E\n", deltaVoltage_);
            printf("                                               predicted          actual \n");
            printf("\t\t  RelativeExtentRxn_init_        :    %15.6E        \n", RelativeExtentRxn_init_);
            printf("\t\t  RelativeExtentRxn_final_       :    %15.6E  %15.6E\n", pred_RelativeExtentRxn_final_, RelativeExtentRxn_final_);
            printf("\t\t  ROP_external                   :    %15.6E  %15.6E\n", pred_rop_external, _rop);
            printf("\t\t  ROP_[0]                        :    %15.6E  %15.6E\n", pred_ROP, ROP_[0]);
            printf("\t\t  Radius_internal_init_          :    %15.6E        \n", Radius_internal_init_);
            printf("\t\t  Radius_internal_final_         :    %15.6E  %15.6E\n", pred_Radius_internal_final_,  Radius_internal_final_);
            printf("\t\t  Radius_exterior_final_         :    %15.6E  %15.6E\n", pred_Radius_exterior_final_,  Radius_exterior_final_);
            printf("\t\t  SurfaceAreaRS0_final_          :    %15.6E  %15.6E\n", pred_surfaceAreaRS0_final_,   surfaceAreaRS_final_[0]);
            printf("\t\t  SurfaceAreaRS1_final_          :    %15.6E  %15.6E\n", pred_surfaceAreaRS1_final_,   surfaceAreaRS_final_[1]);
            printf("\t\t  DaOuter                        :    %15.6E  %15.6E\n", pred_DaOuter_Bar_, DaOuter_Bar_);
            printf("\t\t  DaOuter                        :    %15.6E  %15.6E\n", pred_DaOuter_, DaOuter_);
            printf("\t\t  DaInner                        :    %15.6E  %15.6E\n", pred_DaInner_, DaInner_);
            printf("\t\t  SrcDot_RelativeExtentRxn_final_:    %15.6E  %15.6E\n", pred_SrcDot_RelativeExtentRxn_final_, SrcDot_RelativeExtentRxn_final_);
            printf("\t\t  Resid_RelativeExtentRxn_final_ :    %15.6E        \n", RelativeExtentRxn_final_ - RelativeExtentRxn_tmp_);
            printf("\t\t  DeltaTsubcycleCalc_            :    %15.6E  %15.6E\n", pred_deltaTsubcycleCalc_ , deltaTsubcycleCalc_);
            printf("\t\t  Onboundary_init_               :    %15d      \n",     onRegionBoundary_init_);
            printf("\t\t  Onboundary_final_              :    %15d  %15d\n",     pred_onRegionBoundary_final_ , onRegionBoundary_final_);
            printf("\t\t -------------------------------------------------------------------------------------------------------------------\n");
            if (pred_onRegionBoundary_final_ != -1 ||  onRegionBoundary_final_ != -1) {
                printf("we are here\n");
            }
        }
#endif


    }
    resid[1] =  RelativeExtentRxn_final_ -  RelativeExtentRxn_tmp_;

    if (goNowhere_) {
        if (SrcDot_ExtentRxn_final_ != 0.0) {
            printf("Electrode_MP_RxnExtent::calcResid() ERROR: CONFUSED\n");
            return 0;
        }
    }
    return 1;
}
//====================================================================================================================
//!  Gather the predicted solution values and the predicted integrated source terms
/*!
 *  (virtual from Electrode_Integrator)
 *
 *  Both the predicted solution values and the predicted integrated source terms are used
 *  in the time step control
 */
void Electrode_MP_RxnExtent::gatherIntegratedSrcPrediction()
{
    extractInfo();
    updateSpeciesMoleChangeFinal();
    IntegratedSrc_Predicted.resize(m_NumTotSpecies);
    for (size_t isp = 0; isp < m_NumTotSpecies; isp++) {
        IntegratedSrc_Predicted[isp] = DspMoles_final_[isp] * deltaTsubcycleCalc_;
    }
}
//====================================================================================================================
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
 *
 * @return 1 Means a good calculation that produces a valid result
 *           0 Bad calculation that means that the current nonlinear iteration should be terminated
 */
int Electrode_MP_RxnExtent::integrateResid(const double t, const double delta_t,
        const double* const y, const double* const ySolnDot,
        double* const resid,
        const ResidEval_Type_Enum evalType, const int id_x,
        const double delta_x)
{
    if ((enableExtraPrinting_ && detailedResidPrintFlag_ > 1) || evalType == Base_ShowSolution) {
        printf("\t\t===============================================================================================================================\n");
        printf("\t\t  EXTRA PRINTING FROM NONLINEAR RESIDUAL: ");
        if (evalType ==  Base_ResidEval) {
            printf(" BASE RESIDUAL");
        } else if (evalType == JacBase_ResidEval) {
            printf(" BASE JAC RESIDUAL");
        } else  if (evalType == JacDelta_ResidEval) {
            printf(" DELTA JAC RESIDUAL");
            printf(" var = %d delta_x = %12.4e Y_del = %12.4e Y_base = %12.4e", id_x, delta_x, y[id_x], y[id_x] - delta_x);
        } else  if (evalType == Base_ShowSolution) {
            printf(" BASE RESIDUAL - SHOW SOLUTION");
        }
        printf(" DomainNumber = %2d , CellNumber = %2d , SubIntegrationCounter = %d\n",
               electrodeDomainNumber_, electrodeCellNumber_, counterNumberSubIntegrations_);
        if (onRegionBoundary_final_ >= 0) {
            printf("\t\t\tSPECIAL SOLVE: Final Point is on Boundary %d\n", onRegionBoundary_final_);
        }
        if (limitingEquationBehavior_ == 1) {
            printf("\t\t\tSPECIAL SOLVE: Da_outer End of plateau equation formulation\n");
        }
        printf("\t\t\tT_init = %12.5g       T_final_ = %12.5g        deltaTsubcycle_ = %12.5e \n", tinit_, tfinal_, deltaTsubcycle_);
    }
    /*
     *  UNPACK THE SOLUTION VECTOR
     */
    unpackNonlinSolnVector(y);

    /*
     *  We either
     */
    //  if (evalType == JacBase_ResidEval) {
    double rTop = RegionBoundaries_ExtentRxn_[xRegion_final_ + 1];
    double rBot = RegionBoundaries_ExtentRxn_[xRegion_final_];
    if (onRegionBoundary_final_ < 0) {
        if (RelativeExtentRxn_final_ > rTop * (1.0000001)) {
            if (evalType == JacBase_ResidEval) {
                onRegionBoundary_final_ = xRegion_init_ + 1;

                bool behaviorChangePossible = false;
                if (solidDiffusionModel_) {
                    if (ROPModificationType_ == 2) {
                        behaviorChangePossible = true;
                    }
                }
                if (behaviorChangePossible) {
                    if (limitingEquationBehavior_ == 0) {
                        return -1;
                    }
                }
            }
        }
        if (onRegionBoundary_final_ == xRegion_init_ + 1) {
            if (deltaTsubcycleCalc_ > 1.5 * deltaTsubcycle_) {
                onRegionBoundary_final_ = -1;
            }
        }
        if (RelativeExtentRxn_final_ < rBot * (0.99999999)) {
            onRegionBoundary_final_ = xRegion_init_;
        }
    }
    //   }


    updateState();

    extractInfo();
    updateMoleRatesFinal();
    updateSpeciesMoleChangeFinal();

    /*
     * Calculate the residual
     */
    int info = calcResid(resid, evalType);

    if ((enableExtraPrinting_ && detailedResidPrintFlag_ > 1) || evalType == Base_ShowSolution) {

        if (info != 1) {
            printf("\t\t    !!!!  FAILED RESIDUAL CALCULATION info = %d !!!!\n", info);
        }
        if (deltaTdeath_ < deltaTsubcycle_) {
            double deltax = 0.0;
            double xBd = RegionBoundaries_ExtentRxn_[xRegion_init_ + 1];
            if (SrcDot_ExtentRxn_final_ > 1.0E-200) {
                deltax = xBd - RelativeExtentRxn_init_;
            } else  if (SrcDot_ExtentRxn_final_ < -1.0E-200) {
                xBd = RegionBoundaries_ExtentRxn_[xRegion_init_ ];
                deltax = xBd - RelativeExtentRxn_init_;
            }
            printf("\t\t regionEnd = %d  deltaTcalc = %16.7E  Pinit/-DelP = %16.7E  Res = %16.7E   \n",
                   xRegion_init_,  deltaTsubcycleCalc_, deltax * spMoles_FeS2_Normalization_/ SrcDot_ExtentRxn_final_,  resid[0]);
        }
        printf("\t\t        PhaseName        Moles_Init    Moles_final     KExists  |   Src_Moles "
               "Calc_Moles_Final  |    Resid     |\n");
        printf("\t\t ---------------------------------------------------------------------"
               "----------------------------------------------------------\n");
        double src =   SrcDot_ExtentRxn_final_ * deltaTsubcycleCalc_ / spMoles_FeS2_Normalization_;
        double res =   RelativeExtentRxn_final_ -(RelativeExtentRxn_init_ + src);
        printf("\t\t %20.20s %13.6e  %13.6e           %2d |  %13.6e %13.6e |  %12.4e  |\n", "RelativeExtentRxn",  RelativeExtentRxn_init_,
               RelativeExtentRxn_final_, 1, src, RelativeExtentRxn_init_ + src, res);
        printf("\t\t====================================================================="
               "==========================================================\n");
    }
    return info;
}
//====================================================================================================================
double Electrode_MP_RxnExtent::l0normM(const std::vector<double>& v1, const std::vector<double>& v2, int num,
                                       const std::vector<double>& atolVec, const double rtol) const
{
    double max0 = 0.0;
    double denom, diff, ee;

    for (int k = 0; k < num; k++) {

        diff = fabs(v1[k] - v2[k]);
        denom = rtol * std::max(fabs(v1[k]), fabs(v2[k]));
        denom = std::max(denom, atolVec[k]);
        ee = diff / denom;
        if (ee > max0) {
            max0 = ee;
        }
    }
    return max0;
}
//====================================================================================================================
void Electrode_MP_RxnExtent::checkRegion(int regionID) const
{
    if (onRegionBoundary_final_ >= 0) {
        int regLow = onRegionBoundary_final_ -1;
        if (regLow < 0) {
            regLow = 0;
        }
        int regHigh = onRegionBoundary_final_;
        if (regHigh >= (int)(RegionBoundaries_ExtentRxn_.size() - 1)) {
            regHigh = onRegionBoundary_final_ - 1;
        }
        if (regionID == regHigh || regionID == regLow) {
            // we are fine
        } else {
            throw CanteraError("checkRegion()", "On a boundary but regionID doesn't agree");
        }
    } else {
        int reg = findRegion(RelativeExtentRxn_final_);
        if (reg != regionID) {
            /*
             *  If we're on the boundary, we allow for a little slop
             */
            if (reg > regionID) {
                double reg_deltalow = RelativeExtentRxn_final_ * (1.0 - 1.0E-5);
                int reg_low = findRegion(reg_deltalow);
                if (reg_low != regionID) {
                    throw CanteraError("checkRegion()", "region ID " + int2str(regionID) +
                                       " doesn't agree with current extent " +
                                       fp2str(RelativeExtentRxn_final_) + ", which is in region " + int2str(reg));
                }
            } else {
                double reg_deltahigh = RelativeExtentRxn_final_ * (1.0 + 1.0E-5);
                int reg_high = findRegion(reg_deltahigh);
                if (reg_high != regionID) {
                    throw CanteraError("checkRegion()", "region ID " + int2str(regionID) +
                                       " doesn't agree with current extent " +
                                       fp2str(RelativeExtentRxn_final_) + ", which is in region " + int2str(reg));
                }
            }
        }
    }
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
 */
void  Electrode_MP_RxnExtent::initialPackSolver_nonlinFunction()
{
    /*
     *  Let's formulate the solution vector. The solution vector will consist of the following form
     *
     *          0    deltaTsubcycleCalc_
     *          1    RxnExtent
     */
    checkRegion(xRegion_final_);
    int nTop = (int) RegionBoundaries_ExtentRxn_.size() - 1;
    double delBounds = 0.01 * fabs(RegionBoundaries_ExtentRxn_[nTop]);

    yvalNLS_[0] = deltaTsubcycleCalc_;
    ylowNLS_[0] = 0.0;
    yhighNLS_[0] = 1.0E300;


    yvalNLS_[1] = RelativeExtentRxn_final_;
    ylowNLS_[1] = -delBounds;
    yhighNLS_[1] = RegionBoundaries_ExtentRxn_[nTop] + delBounds;



}
//==================================================================================================================
//   Calculate the integrated source terms and do other items now that we have a completed time step
/*
 *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
 *  for the current step, and we are calculating post-processed quantities like source terms.
 */
void Electrode_MP_RxnExtent::calcSrcTermsOnCompletedStep()
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
//==================================================================================================================
//   Accumulate src terms and other results from the local step into the global holding bins.
/*
 *  Accumulate source terms on completion of a step. At this point we have solved the nonlinear problem
 *  for the current step and we have satisfied all accuracy requirements.
 *  The step is good. We now accumulate the results before going on to a new local step.
 */
void Electrode_MP_RxnExtent::accumulateSrcTermsOnCompletedStep(bool remove)
{
    if (remove) {
        for (size_t i = 0; i < m_NumTotSpecies; i++) {
            spMoleIntegratedSourceTerm_[i] -= spMoleIntegratedSourceTermLast_[i];
        }
    } else {
        for (size_t i = 0; i < m_NumTotSpecies; i++) {
            spMoleIntegratedSourceTerm_[i] += spMoleIntegratedSourceTermLast_[i];
        }
    }
}
//==================================================================================================================
//! Check to see that the preceding step is a successful one
/*!
 *  (virtual from Electrode_Integrator)
 *
 *   We check to see if the preceding step is a successful one.
 *
 *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
 */
bool  Electrode_MP_RxnExtent::checkSubIntegrationStepAcceptable() const
{
    bool stepAcceptable = true;

    return stepAcceptable;
}
//==================================================================================================================
// Possibly change the solution due to phase births and deaths.
/*
 *   (virtual from Electrode_Integrator)
 *
 *  @return  Returns true if the solution step is good. It returns false if there is a problem.
 */
bool Electrode_MP_RxnExtent::changeSolnForBirthDeaths()
{
    bool stepAcceptable = true;
    if (onRegionBoundary_final_ < 0) {
        int reg = findRegion(RelativeExtentRxn_final_);
        if (xRegion_final_ != reg) {
            if (RelativeExtentRxn_final_ > RegionBoundaries_ExtentRxn_[xRegion_final_+1]) {
                if (fabs(RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[xRegion_final_+1]) < 0.001) {
                    RelativeExtentRxn_final_ = RegionBoundaries_ExtentRxn_[xRegion_final_+1];
                    onRegionBoundary_final_ = xRegion_final_+1;
                    updateState();

                    extractInfo();
                    updateMoleRatesFinal();
                    updateSpeciesMoleChangeFinal();
                    double st =  RelativeExtentRxn_final_ -  RelativeExtentRxn_init_;
                    st *= spMoles_FeS2_Normalization_;
                    double fac = fabs(st / (DspMoles_final_[kElectron_] * deltaTsubcycleCalc_));
                    deltaTsubcycleCalc_ *= fac;
                } else {
                    stepAcceptable = false;
                }
            }
            if (RelativeExtentRxn_final_ < RegionBoundaries_ExtentRxn_[xRegion_final_]) {
                if (fabs(RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[xRegion_final_]) < 0.001) {
                    RelativeExtentRxn_final_ = RegionBoundaries_ExtentRxn_[xRegion_final_];
                    onRegionBoundary_final_ = xRegion_final_;
                    updateState();

                    extractInfo();
                    updateMoleRatesFinal();
                    updateSpeciesMoleChangeFinal();
                    double st =  RelativeExtentRxn_final_ -  RelativeExtentRxn_init_;
                    st *= spMoles_FeS2_Normalization_;
                    double fac = fabs(st / (DspMoles_final_[kElectron_] * deltaTsubcycleCalc_));
                    deltaTsubcycleCalc_ *= fac;


                } else {
                    stepAcceptable = false;
                }
            }
        }
    }
    if (! stepAcceptable) {
        return stepAcceptable;
    }
    if (onRegionBoundary_final_ < 0) {
        /*
         *  If we are close to the boundary, we have to declare that we are on the boundary
         */
        if (fabs(RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[xRegion_final_+1]) < 1.0E-10) {
            if ((RelativeExtentRxn_final_ > RelativeExtentRxn_init_) && (RegionBoundaries_ExtentRxn_[xRegion_final_+1] >  RelativeExtentRxn_init_)) {
                RelativeExtentRxn_final_ = RegionBoundaries_ExtentRxn_[xRegion_final_+1];
                onRegionBoundary_final_ = xRegion_final_ + 1;
                updateState();
                extractInfo();
                updateMoleRatesFinal();
                updateSpeciesMoleChangeFinal();
                double st =  RelativeExtentRxn_final_ -  RelativeExtentRxn_init_;
                st *= spMoles_FeS2_Normalization_;
                double fac = fabs(st / (DspMoles_final_[kElectron_] * deltaTsubcycleCalc_));
                deltaTsubcycleCalc_ *= fac;
            }
        } else  if (fabs(RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[xRegion_final_]) < 1.0E-10) {
            if ((RelativeExtentRxn_final_ < RelativeExtentRxn_init_) && (RegionBoundaries_ExtentRxn_[xRegion_final_] < RelativeExtentRxn_init_)) {
                RelativeExtentRxn_final_ = RegionBoundaries_ExtentRxn_[xRegion_final_];
                onRegionBoundary_final_ = xRegion_final_;
                updateState();
                extractInfo();
                updateMoleRatesFinal();
                updateSpeciesMoleChangeFinal();
                double st =  RelativeExtentRxn_final_ -  RelativeExtentRxn_init_;
                st *= spMoles_FeS2_Normalization_;
                double fac = fabs(st / (DspMoles_final_[kElectron_] * deltaTsubcycleCalc_));
                deltaTsubcycleCalc_ *= fac;
            }
        }
    }
    return stepAcceptable;
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
double Electrode_MP_RxnExtent::predictorCorrectorWeightedSolnNorm(const std::vector<double>& yval)
{
    //double atolVal =  1.0E-8;
    atolNLS_.resize(2);

    atolNLS_[0] = 0.01 * (t_final_final_ - t_init_init_);
    atolNLS_[0] = fmax(atolNLS_[0], 1.0E-6);
    atolNLS_[1] = 1.0E-8;
    double pnorm = l0normM(soln_predict_, yval, 2, atolNLS_, rtolNLS_);
    return pnorm;
}
//====================================================================================================================
// Calculate the vector of predicted errors in the source terms that this integrator is responsible for
/*!
 *  (virtual from Electrode_Integrator)
 *
 *    In the base implementation we assume that the there are just one source term, the electron
 *    source term.
 *    However, this will be wrong in almost all cases.
 *    The number of source terms is unrelated to the number of unknowns in the nonlinear problem.
 *    Source terms will have units associated with them.
 *    For example the integrated source term for electrons will have units of kmol
 */
void Electrode_MP_RxnExtent::predictorCorrectorGlobalSrcTermErrorVector()
{

}
//====================================================================================================================
// Print table representing prediction vs. corrector information
/*
 *  @param yval           Vector of corrector values
 *  @param pnormSrc       Norm of the predictor-corrector comparison for the source vector.
 *  @param pnormSoln      Norm of the predictor-corrector comparison for the solution vector.
 */
void Electrode_MP_RxnExtent::predictorCorrectorPrint(const std::vector<double>& yval, double pnormSrc, double pnormSoln) const
{
    double atolVal =  1.0E-8;
    double denom;
    double tmp;
    int onRegionPredict = soln_predict_[2];
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf(" PREDICTOR_CORRECTOR  SubIntegrationCounter = %7d       t_init = %12.5E,       t_final = %12.5E\n",
           counterNumberSubIntegrations_, tinit_, tfinal_);
    printf("                         IntegrationCounter = %7d  t_init_init = %12.5E, t_final_final = %12.5E\n",
           counterNumberIntegrations_, t_init_init_, t_final_final_);
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("                               Initial       Prediction      Actual          Difference         Tol   Contrib      |\n");
    printf("onRegionPredict          |         %3d            %3d          %3d      |                                          |\n",
           onRegionBoundary_init_, onRegionPredict, onRegionBoundary_final_);
    denom = std::max(fabs(yval[0]), fabs(soln_predict_[0]));
    denom = std::max(denom, atolVal);
    tmp = fabs((yval[0] - soln_predict_[0])/ denom);
    printf("DeltaT                   | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E |\n",
           deltaTsubcycle_, soln_predict_[0],  yval[0], yval[0] - soln_predict_[0], atolVal, tmp);
    denom = std::max(fabs(yval[1]), fabs(soln_predict_[1]));
    denom = std::max(denom, atolVal);
    tmp = fabs((yval[1] - soln_predict_[1])/ denom);
    printf("RelativeExtentRxn_final_ | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E | \n",
           RelativeExtentRxn_init_, soln_predict_[1],  yval[1], yval[1] - soln_predict_[1], atolVal, tmp);
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("                                                                                                        %10.3E\n",
           pnormSoln);

    /*
     *  Do a special call to the residual evaluator. This will print out predicted quantities versus actual quantities.
     *  This is very useful for debugging the predictor.
     */
    //double resid[8];
    //integrateResid(tfinal_, deltaTsubcycle_, &yvalNLS_[0], &ydotNLS_[0], resid, Base_ShowSolution, 0, 0.0);


}

//======================================================================================================================
//!  Calculate the norm of the errors in the global source terms
/*!
 *  (virtual from Electrode_Integrator)
 *
 *   This routine make use of the source term error vector along with rtols and atols for the
 *   individual source terms to calculated a normalized error measure. This is the single number
 *   that the integration routine will try to control as it calculates a time stepping strategy.
 *
 *   @return  Returns a single nondimensional number representing the normalized error
 *            for the calculation of the source term
 */
double Electrode_MP_RxnExtent::predictorCorrectorGlobalSrcTermErrorNorm()
{
    return 0.0;
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
void  Electrode_MP_RxnExtent::resetStartingCondition(double Tinitial, bool doTestsAlways)
{
    bool resetToInitInit = false;
    /*
     * If the initial time input from the parameter list, Tinitial, is the same as the current initial time,
     * Then, we don't advance the time step.
     */
    double tbase = std::max(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-9 * tbase) && !doTestsAlways) {
        resetToInitInit = true;
    }


    tbase = std::max(Tinitial, tbase);
    tbase = std::max(tbase, t_final_final_);
    if (!resetToInitInit) {
	if (fabs(Tinitial - t_final_final_) > (1.0E-9 * tbase)) {
	    throw CanteraError("Electrode_MP_RxnExtent::resetStartingCondition()", "tinit " + fp2str(Tinitial) +" not compat with t_final_final_ "
			       + fp2str(t_final_final_));
	}
    }

    Electrode_Integrator::resetStartingCondition(Tinitial);


    if (!resetToInitInit) {

	Radius_internal_init_init_ =  Radius_internal_final_final_;
	Radius_internal_init_      =  Radius_internal_final_final_;
	Radius_internal_final_     =  Radius_internal_final_final_;
	
	// Copy The final Extent of reaction to the beginning extent
	RelativeExtentRxn_init_init_ =   RelativeExtentRxn_final_final_;
	RelativeExtentRxn_init_      =   RelativeExtentRxn_final_final_;
	
	xRegion_init_init_ = xRegion_final_final_;
	xRegion_init_      = xRegion_final_final_;
	
	molarVolume_init_init_ = molarVolume_final_final_;
	molarVolume_init_      = molarVolume_final_final_;
	molarVolume_final_     = molarVolume_final_final_;
	
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
	int ip_FeS2_A = globalPhaseIndex("FeS2_A(S)");
	int is_FeS2_A = globalSpeciesIndex("FeS2_A(S)");
	int ip_FeS2_B = globalPhaseIndex("FeS2_B(S)");
	int is_FeS2_B = globalSpeciesIndex("FeS2_B(S)");
	spMoles_final_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
	spMoles_init_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
	spMoles_init_init_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
	phaseMoles_final_[ip_FeS2_A] =  spMoles_FeS2_Normalization_ / 2.0;
	
	spMoles_final_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
	spMoles_init_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
	spMoles_init_init_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
	phaseMoles_final_[ip_FeS2_B] =  spMoles_FeS2_Normalization_ / 2.0;
	
	if (eState_final_) {
	    if (!xmlStateData_final_) {
		eState_final_->copyElectrode_intoState(this);
		xmlStateData_final_ = eState_final_->write_electrodeState_ToXML();
	    }
	    delete xmlStateData_init_init_;
	    xmlStateData_init_init_ =   xmlStateData_final_;
	    delete xmlStateData_init_;
	    xmlStateData_init_ = new XML_Node(*xmlStateData_final_);
	    xmlStateData_final_ = 0;
	    delete xmlStateData_final_final_;
	    xmlStateData_final_final_ = 0;
	}
    }
}
//====================================================================================================================
// Set the internal initial intermediate and initial global state from the internal final state
/*
 *  (non-virtual function onionize in-first)
 *
 *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 * @param setInitInit   Boolean indicating whether you should set the init_init state as well
 */
void Electrode_MP_RxnExtent::setInitStateFromFinal_Oin(bool setInitInit)
{
    /*
     * Call the parent object
     */
    Electrode_Integrator::setInitStateFromFinal(setInitInit);


    // Copy The final Extent of reaction to the beginning extent

    RelativeExtentRxn_init_ = RelativeExtentRxn_final_;
    RelativeExtentRxn_final_final_ = RelativeExtentRxn_final_;

    xRegion_init_ = xRegion_final_;
    xRegion_final_final_ = xRegion_final_;

    onRegionBoundary_init_ = onRegionBoundary_final_;

    molarVolume_init_ = molarVolume_final_;
    molarVolume_final_final_ = molarVolume_final_;

    Radius_exterior_init_ = Radius_exterior_final_;

    Radius_internal_final_final_ = Radius_internal_final_;
    Radius_internal_init_ = Radius_internal_final_;

    if (setInitInit) {
        RelativeExtentRxn_init_init_ = RelativeExtentRxn_final_;
        xRegion_init_init_ = xRegion_final_;
        molarVolume_init_init_ = molarVolume_final_;
        Radius_internal_init_init_ = Radius_internal_final_;
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
void Electrode_MP_RxnExtent::setInitStateFromFinal(bool setInitInit)
{
    setInitStateFromFinal_Oin(setInitInit);
}
//====================================================================================================================
// Set the internal final intermediate and from the internal init state
/*
 *  (virtual function from Electrode)
 *
 *  Set the final state from the init state. This is commonly called during a failed time step
 *
 */
void  Electrode_MP_RxnExtent::setFinalStateFromInit()
{
    Electrode_Integrator::setFinalStateFromInit_Oin();
    /*
     * Do stuff not done in base class
     /*/
    RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
    xRegion_final_ = xRegion_init_;
    onRegionBoundary_final_ = onRegionBoundary_init_;
    Radius_internal_final_ = Radius_internal_init_;
    molarVolume_final_ = molarVolume_init_;
    /*
     *  If this object becomes a parent, then we will have to create a setFinalStateFromInit_Oin() to avoid
     *  calling this routine.
     */
    updateState();
}
//====================================================================================================================
// Set the internal initial intermediatefrom the internal initial global state
/*
 *  Set the intial state from the init init state. We also can set the final state from this
 *  routine as well.
 *
 *  The final_final is not touched.
 *
 * @param setFinal   Boolean indicating whether you should set the final as well
 */
void Electrode_MP_RxnExtent::setInitStateFromInitInit(bool setFinal)
{

    Electrode_Integrator::setInitStateFromInitInit(setFinal);

    RelativeExtentRxn_init_ = RelativeExtentRxn_init_init_;
    xRegion_init_ = xRegion_init_init_;
    Radius_internal_init_ = Radius_internal_init_init_;
    molarVolume_init_ = molarVolume_init_init_;

#ifdef DEBUG_MODE_NOT
    if (counterNumberSubIntegrations_ >= 7384 && electrodeCellNumber_ == 0) {
        printf("we are here\n");
    }
#endif
    onRegionBoundary_init_ = -1;
    for (int i = 0; i < (int) RegionBoundaries_ExtentRxn_.size(); i++) {
        if (fabs(RelativeExtentRxn_init_ - RegionBoundaries_ExtentRxn_[i]) < 1.0E-12) {
            onRegionBoundary_init_ = i;
            break;
        }
    }

    if (setFinal) {
        RelativeExtentRxn_final_ = RelativeExtentRxn_init_;
        xRegion_final_ = xRegion_init_;
        onRegionBoundary_final_ = onRegionBoundary_init_;
        Radius_internal_final_ = Radius_internal_init_;
        molarVolume_final_ = molarVolume_init_;
        onRegionBoundary_final_ = -1;
        for (int i = 0; i < (int) RegionBoundaries_ExtentRxn_.size(); i++) {
            if (fabs(RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[i]) < 1.0E-12) {
                onRegionBoundary_final_ = i;
                break;
            }
        }
        updateState();
    }

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
void Electrode_MP_RxnExtent::setInitInitStateFromFinalFinal()
{
    Electrode_Integrator::setInitInitStateFromFinalFinal();

    RelativeExtentRxn_init_init_ = RelativeExtentRxn_final_final_;
    xRegion_init_init_           = xRegion_final_final_;
    Radius_internal_init_init_   = Radius_internal_final_final_;
    molarVolume_init_init_       = molarVolume_final_final_;

    RelativeExtentRxn_init_   = RelativeExtentRxn_final_final_;
    xRegion_init_             = xRegion_final_final_;
    Radius_internal_init_     = Radius_internal_final_final_;
    molarVolume_init_         = molarVolume_final_final_;

    //   onRegionBoundary_init_init_  = onRegionBoundary_final_final_;
    onRegionBoundary_init_ = -1;
    for (int i = 0; i < (int) RegionBoundaries_ExtentRxn_.size(); i++) {
        if (fabs(RelativeExtentRxn_init_ - RegionBoundaries_ExtentRxn_[i]) < 1.0E-12) {
            onRegionBoundary_final_ = i;
            break;
        }
    }

}
//====================================================================================================================
//! Set the internal final global state from the internal final intermediate state
/*!
 *  (virtual function from Electrode)
 *
 *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
 */
void Electrode_MP_RxnExtent::setFinalFinalStateFromFinal()
{
    Electrode_Integrator::setFinalFinalStateFromFinal_Oin();
    RelativeExtentRxn_final_final_ = RelativeExtentRxn_final_;
    xRegion_final_final_ = xRegion_final_;
    Radius_internal_final_final_ = Radius_internal_final_;
    molarVolume_final_final_ = molarVolume_final_;
}


//  -----------------------------------------------------------------------------------------------------------------
int Electrode_MP_RxnExtent::getInitialConditions(const double t0, double* const y, double* const ydot)
{
    for (int k = 0; k < neq_; k++) {
        y[k] = 0.0;
    }
    return 1;
}
//====================================================================================================================
double
Electrode_MP_RxnExtent::filterNewStep(const double timeCurrent,
                                      const double* const ybase,
                                      double* const step0)
{
    return -1;
}
//====================================================================================================================
int Electrode_MP_RxnExtent::calcDeltaSolnVariables(const double t, const double* const ySoln,
        const double* const ySolnDot, double* const deltaYSoln,
        const double* const solnWeights)
{
    if (!solnWeights) {
        for (int i = 0; i < neq_; i++) {
            deltaYSoln[i] = m_atol + fabs(1.0E-6 * ySoln[i]);
        }
    } else {
        for (int i = 0; i < neq_; i++) {
            deltaYSoln[i] = fmax(1.0E-2 * solnWeights[i], 1.0E-6 * fabs(ySoln[i]));
        }
    }
    return 1;
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
int Electrode_MP_RxnExtent::evalResidNJ(const double tdummy, const double delta_t_dummy,
                                        const double* const y,
                                        const double* const ySolnDot,
                                        double* const resid,
                                        const ResidEval_Type_Enum evalType,
                                        const int id_x,
                                        const double delta_x)
{
    int retn =  integrateResid(tdummy, delta_t_dummy, y, ySolnDot, resid, evalType, id_x, delta_x);
    return retn;
}

//====================================================================================================================
void  Electrode_MP_RxnExtent::setState_TP(double temperature, double pressure)
{
    temperature_ = temperature;
    pressure_ = pressure;


    double lb;
    // Have concerns about changing a boundary in the middle of a calculation HKM
    double temperatureTmp = 723.15;
    lb = 1.539 + 0.289E-3 * temperatureTmp;
    RegionBoundaries_ExtentRxn_[2] = lb;
    /*
     *  We also need to change the standard state OCV. This is taken care of within updateState().
     */
    updateState();
}
//====================================================================================================================
int Electrode_MP_RxnExtent::setVoltageVsExtent_FeS2()
{
    RegionBoundaries_ExtentRxn_.resize(5,0.0);
    RegionBoundaries_ExtentRxn_[0] = 0.0;
    RegionBoundaries_ExtentRxn_[1] = 1.5;
    RegionBoundaries_ExtentRxn_[2] = 1.75;
    RegionBoundaries_ExtentRxn_[3] = 2.0;
    RegionBoundaries_ExtentRxn_[4] = 4.0;


    return 0;
}
//====================================================================================================================
//  Get the region
/*
 *     if RelativeExtentRxn is below the low boundary return -1
 *     If RelativeExtentRxn is at the low boundary return region number containing the low boundary
 *     If RelativeExtentRxn is at the high boundary return region number of the next highest region
 */
int Electrode_MP_RxnExtent::findRegion(double RelativeExtentRxn) const
{
    int n = RegionBoundaries_ExtentRxn_.size() - 1;
    for (int i = -1; i < n; i++) {
        if (RelativeExtentRxn < RegionBoundaries_ExtentRxn_[i+1]) {
            return i;
        }
    }
    return n;
}
//====================================================================================================================
double Electrode_MP_RxnExtent::relativeExtentRxn(double time) const
{
    if (fabs(time - tfinal_) < 1.0E-50) {
        return RelativeExtentRxn_final_;
    }
    if (fabs(time - tinit_) < 1.0E-50) {
        return RelativeExtentRxn_init_;
    }
    if (fabs(time - t_init_init_) < 1.0E-50) {
        return RelativeExtentRxn_init_init_;
    }
    if (fabs(time - t_final_final_) < 1.0E-50) {
        return RelativeExtentRxn_final_final_;
    }
    throw CanteraError("Electrode_MP_RxnExtent::relativeExtentRxn(double time)", "unknown time: " + fp2str(time));
    return 0.0;
}
//====================================================================================================================
double Electrode_MP_RxnExtent::molarVolume_relExtentRxn(double relativeExtentRxn) const
{
    int numRegions = RegionBoundaries_ExtentRxn_.size() - 1;
    int ireg = findRegion(relativeExtentRxn);
    if (ireg < 0) {
        ireg = 0;
    } else if (ireg > numRegions - 1) {
        ireg = numRegions - 1;
    }
    double deltaReg = RegionBoundaries_ExtentRxn_[ireg+1] - RegionBoundaries_ExtentRxn_[ireg];
    double relRegion = (relativeExtentRxn - RegionBoundaries_ExtentRxn_[ireg]) / deltaReg;
    double molarVol = molarVolumeRegions_[ireg] * relRegion + molarVolumeRegions_[ireg+1] * (1.0 - relRegion);
    return molarVol;
}
//====================================================================================================================

double Electrode_MP_RxnExtent::openCircuitVoltageSSRxn(int isk, int iReaction) const
{
    if (isk != indexOfReactingSurface_) {
        return 0.0;
        //  printf("Electrode_MP_RxnExtent::openCircuitVoltageSS() ERROR: bad isk\n");
        // exit(-1);
    }
    double voltsSSBase = Electrode::openCircuitVoltageSSRxn(isk, iReaction);
    double voltsSS =  openCircuitVoltageSS_Region(RelativeExtentRxn_final_, xRegion_final_);
    if (fabs(voltsSSBase - voltsSS) > 1.0E-6) {
        throw CanteraError("Electrode_MP_RxnExtent::openCircuitVoltageSSRxn()",
                           "Internal inconsistency: " +  fp2str(voltsSSBase) + "  " + fp2str(voltsSS));
    }
    return voltsSSBase;
}
//====================================================================================================================
double Electrode_MP_RxnExtent::openCircuitVoltageSS_Region(double relativeExtentRxn, int xRegion) const
{

    double lb;
    double voltage, voltage_b, voltage_3;
    double nElectrons = relativeExtentRxn;


    if ((temperature_ < 650) | (temperature_ > 800)) {
        throw CanteraError("Electrode_MP_RxnExtent::openCircuitVoltageSS_Region()",
                           "Temperature outside region of fit");
    }
    int ireg = findRegion(relativeExtentRxn);
    if (ireg != xRegion) {

        double deltaLook = 0.00001 + fabs(relativeExtentRxn) * 0.001;
        if (ireg < xRegion) {
            ireg = findRegion(relativeExtentRxn + deltaLook);
            if (ireg != xRegion) {
                if (enableExtraPrinting_ && (detailedResidPrintFlag_ > 4)) {
                    printf("Electrode_MP_RxnExtent::openCircuitVoltageSS_Region(): we may have some trouble actual = %d, requested= %d\n",
                           ireg, xRegion);
                }
            }
        } else if (ireg > xRegion) {
            ireg = findRegion(relativeExtentRxn - deltaLook);
            if (ireg != xRegion) {
                if (enableExtraPrinting_ && (detailedResidPrintFlag_ > 4)) {
                    printf("Electrode_MP_RxnExtent::openCircuitVoltageSS_Region(): we may have some trouble actual = %d, requested= %d\n",
                           ireg, xRegion);
                }
            }
        }
    }
    if (ireg < 0) {
        ireg = 0;
    } else if (ireg > 3) {
        ireg = 3;
    }


    Electrode_Types_Enum et = electrodeType();

    if (et == MP_RXNEXTENT_ET) {
    switch (xRegion) {
    case 0:
        voltage = 1.889538 + 0.2297e-3*temperature_;
        break;
    case 1:
        voltage = 1.6613762 + 0.4178e-3*temperature_;
        break;
    case 2:
        lb = RegionBoundaries_ExtentRxn_[2];
        voltage_b = 1.6613762 + 0.4178e-3*temperature_;
        voltage_3 = 1.803338 - 0.2355e-3*temperature_;
        voltage = voltage_b + (nElectrons - lb) / (2.0 - lb) * (voltage_3 - voltage_b);
        break;
    case 3:
        voltage = 1.896438 - 0.3958e-3*temperature_;
        break;
    default:
        voltage = 1000.;
        //throw CanteraError("Electrode_MP_RxnExtent::openCircuitVoltageSS_Region()",
        //			 "nElectrons outside valid region, less than 0.0");
        break;
    }

    } else if (et == MP_RXNEXTENT_FES2_ET) {
	voltage = openCircuitVoltageSS_Region_NoCheck(relativeExtentRxn, xRegion);
    } else {
       throw CanteraError("", "Bad eletrodeType");
    }
    return voltage;
}
//=======================================================================================================
double 
Electrode_MP_RxnExtent::openCircuitVoltageSS_Region_NoCheck(double relativeExtentRxn,
							    int xRegion) const
{
    double voltage, voltage_b, voltage_3, lb;
    switch (xRegion) {
    case 0:
        voltage = 1.889538 + 0.2297e-3*temperature_;
        break;
    case 1:
        voltage = 1.6613762 + 0.4178e-3*temperature_;
        break;
    case 2:
        lb = RegionBoundaries_ExtentRxn_[2];
        voltage_b = 1.6613762 + 0.4178e-3*temperature_;
        voltage_3 = 1.803338 - 0.2355e-3*temperature_;
        voltage = voltage_b + (relativeExtentRxn - lb) / (2.0 - lb) * (voltage_3 - voltage_b);
        break;
    case 3:
        voltage = 1.896438 - 0.3958e-3*temperature_;
        break;
    case -1:
	voltage = -1000;
	break;
    default:
        voltage = 1000.;
        break;
    }
    return voltage;
}
//============================================================================================================
double Electrode_MP_RxnExtent::openCircuitVoltageSS_final() const
{
    double volts = openCircuitVoltageSS_Region(RelativeExtentRxn_final_, xRegion_final_);
    return volts;
}
//====================================================================================================================
// Calculate the inner radius at the final state
/*
 *  An inner radius is needed for diffusion approximations. Here we calculate the inner radius
 *  by the extent of reaction carried out in the current region. Everything is relative
 *  to the exterior radius.
 *
 *  When the extent of reaction in the current plateau is zero, we assume that the inner
 *  radius is equal to the external radius (mult by (1 -small)). When the extent of reaction is
 *  nearly equal to the end of the region, we assume the internal radius is nearly zero
 *     (radius_ext * (small))
 *
 *  Here small is defined as 1.0E-4 (will experiment with that number.
 */
double Electrode_MP_RxnExtent::calculateRadiusInner(double relativeExtentRxn_final) const
{
    const double smallFactor = 20.;
    double radius_inner_SmallBound = Radius_exterior_final_ / smallFactor;

    double radius_inner_LargeBound = Radius_exterior_final_ * (1.0 - 1.0/smallFactor);
    if (xRegion_final_ < 0) {
        return  radius_inner_LargeBound;
    } else if (xRegion_final_ == (int) RegionBoundaries_ExtentRxn_.size() - 1) {
        return radius_inner_SmallBound;
    }
    double relExtBot = RegionBoundaries_ExtentRxn_[xRegion_final_];
    double relExtTop = RegionBoundaries_ExtentRxn_[xRegion_final_+1];
    double perCentUnreacted = (relExtTop - relativeExtentRxn_final) / (relExtTop - relExtBot);
    if (perCentUnreacted <0.0) {
        perCentUnreacted = 0.0;
    } else if (perCentUnreacted > 1.0) {
        perCentUnreacted = 1.0;
    }
    double molarVol = molarVolumeRegions_[xRegion_final_];

    double extLeft = relExtTop - relativeExtentRxn_final;
    extLeft = fmax(extLeft, 0.0);
    extLeft = fmin(relExtTop - relExtBot, extLeft);
    double numLeft = extLeft * spMoles_FeS2_Normalization_;

    double r3 = 3.0 * molarVol * numLeft / (4.0 * Pi *  particleNumberToFollow_);
    double r = pow(r3, 0.33333333333333);

    double radius_inner = pow(perCentUnreacted, 0.33333333333333333) * Radius_exterior_final_;
    if (fabs(radius_inner - r) > 0.1) {
        throw CanteraError(" Electrode_MP_RxnExtent::calculateRadiusInner_final()",
                           " of concern, look into ");
    }

    if (radius_inner < radius_inner_SmallBound) {
        return radius_inner_SmallBound;
    }
    if (radius_inner > radius_inner_LargeBound) {
        return radius_inner_LargeBound;
    }
    return  radius_inner;
}
//====================================================================================================================
double Electrode_MP_RxnExtent::openCircuitVoltage(int isk, bool comparedToReferenceElectrode)
{
    return openCircuitVoltageRxn(isk, -1, comparedToReferenceElectrode);
}
//======================================================================================================================
double Electrode_MP_RxnExtent::openCircuitVoltageRxn(int isk, int iReaction, bool comparedToReferenceElectrode) const
{
    if (isk != indexOfReactingSurface_) {
        return 0.0;
    }
    //FIXME: This function is not const and is only used for a debug check here, is that debug check really important?
    //double voltBase = Electrode::openCircuitVoltage(indexOfReactingSurface_);
    double volts =  openCircuitVoltage_Region(RelativeExtentRxn_final_, xRegion_final_);
    /*if (fabs(voltBase - volts) > 1.0E-6) {
        throw CanteraError("Electrode_MP_RxnExtent::openCircuitVoltage()",
                           "Internal inconsistency: " +  fp2str(voltBase) + "  " + fp2str(volts));
    }*/
    return volts;
}
//====================================================================================================================
double Electrode_MP_RxnExtent::openCircuitVoltage_Region(double relativeExtentRxn, int xRegion, bool comparedToReferenceElectrode) const
{
    double voltsSS = openCircuitVoltageSS_Region(relativeExtentRxn, xRegion);
    double deltaGSS = - voltsSS * Faraday;
    double mu[30];
    ThermoPhase* tSalt = &(thermo(solnPhase_));
    //  We assume that the ion in the salt is called Li+
    int ilt_Lip = tSalt->speciesIndex("Li+");
    double deltaG;
    if (comparedToReferenceElectrode) {
        deltaG = deltaGSS; 
    } else {
        tSalt->getChemPotentials(mu);
        double muLip = mu[ilt_Lip];
        tSalt->getPureGibbs(mu);
        double muLipSS = mu[ilt_Lip];
        deltaG = deltaGSS - muLip + muLipSS;
    }
    double volts = - deltaG / Faraday;
    return volts;
}
//====================================================================================================================
// Returns the total capacity of the electrode in Amp seconds
/*
 *  Returns the capacity of the electrode in Amps seconds.
 *  This is the same as the number of coulombs that can be delivered at any voltage.
 *  Note, this number differs from the capacity of electrodes that is usually quoted for
 *  a battery. That number depends on the rate of discharge and also depends on the
 *  specification of a cutoff voltage. Here, we dispense with both of these specifications.
 *  So, it should be considered a theoretical capacity at zero current and minimal cutoff voltage.
 *  It will also include all plateaus that are defined by the electrode object.
 *
 *  @param  platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 *
 *  @return returns the theoretical capacity of the electrode in Amp sec = coulombs.
 */
double Electrode_MP_RxnExtent::capacity(int platNum) const
{
    double fac = 4.0;
    if (platNum >= 0) {
        fac =  RegionBoundaries_ExtentRxn_[platNum+1] -  RegionBoundaries_ExtentRxn_[platNum];
    }
    double capZeroDoD = fac * spMoles_FeS2_Normalization_;
    double tmp = capZeroDoD * Faraday;
    return tmp;
}
//====================================================================================================================
// Initial capacity of the elctrode in Amp seconds
/*
 *  This is the initial capacity of the electrode before any degradation occurs.
 *
 *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 */
double Electrode_MP_RxnExtent::capacityInitial(int platNum) const
{
    return capacity(platNum);
}
//======================================================================================================================
// Amount of charge that the electrode that has available to be discharged
/*
 *   We report the number in terms of Amp seconds = coulombs
 *   TODO -> Put in voltage max and min dependence.
 */
double Electrode_MP_RxnExtent::capacityLeft(int platNum, double voltsMax, double voltsMin) const
{
    if (platNum >= 0) {
        double fac =  RegionBoundaries_ExtentRxn_[platNum+1] -  RegionBoundaries_ExtentRxn_[platNum];
        if (RelativeExtentRxn_final_ < RegionBoundaries_ExtentRxn_[platNum]) {
            return capacity(platNum);
        } else if (RelativeExtentRxn_final_ >= RegionBoundaries_ExtentRxn_[platNum+1]) {
            return 0.0;
        } else {
            fac = RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[platNum];
        }
        return fac * spMoles_FeS2_Normalization_ * Faraday;
    }
    double total = capacity();
    double dis =  - Faraday * electronKmolDischargedToDate_;
    if (pendingIntegratedStep_) {
        dis -= Faraday * spMoleIntegratedSourceTerm_[kElectron_];
    }
    double left = total - dis;
    /*
     *  Do a check
     */
    double left_b = (4.0 - RelativeExtentRxn_final_) * spMoles_FeS2_Normalization_ *  Faraday;
    if (fabs(left_b - left) / total > 0.0001) {
        printf("Electrode_MP_RxnExtent::capacityLeft() ERROR: cell %d, SubIntegrationCounter %d\n",
               electrodeCellNumber_, counterNumberSubIntegrations_);
        printf("                         capacity error at relExtent %17.9E : left_b = %17.9E left = %17.9E\n",
               RelativeExtentRxn_final_, left_b, left);
        double leftR = left / (spMoles_FeS2_Normalization_ *  Faraday);
        printf("                         capacity error at relExtent %17.9E expressed in extent : left_b = %17.9E, left = %17.9E\n",
               RelativeExtentRxn_final_, (4.0 - RelativeExtentRxn_final_), leftR);
        // exit(-1);
        //  throw CanteraError("Electrode_MP_RxnExtent::capacityLeft() ERROR: cell " + int2str(electrodeCellNumber_) +
        //			 ", integrationCounter " + int2str(counterNumberIntegrations_),
        //		 "capacity error at relExtent " + fp2str(RelativeExtentRxn_final_) + " : left_b = "
        //		 + fp2str(left_b) + ", left = " + fp2str(left));

    }
    return left;
}
//====================================================================================================================
//   Set the current capacity discharged in amp seconds
/*
 * This is roughly equal to the total number of electrons that has been discharged
 * from a fully charged state.
 *
 *  @param relDischarged      Current value of the discharged electrons
 *  @param platNum       Plateau number. Default is -1 which treats all plateaus as a single entity.
 */
void  Electrode_MP_RxnExtent::setRelativeCapacityDischargedPerMole(double relDischargedPerMole, int platNum)
{
    double relativeExtentRxn;
    if (platNum < 0) {
        relativeExtentRxn = relDischargedPerMole;
    } else {
        double bb =  RegionBoundaries_ExtentRxn_[platNum+1] -  RegionBoundaries_ExtentRxn_[platNum];
        if (relDischargedPerMole > bb) {
            throw CanteraError(" Electrode_MP_RxnExtent::setCapacityDischarged()",
                               " Too much capacity in plateau " + int2str(platNum) + " assumed");
        }
        relativeExtentRxn =  RegionBoundaries_ExtentRxn_[platNum] + relDischargedPerMole;
    }
    setState_relativeExtentRxn(relativeExtentRxn);
}
//====================================================================================================================
//   Set the current capacity discharged in amp seconds
/*
 * This is roughly equal to the total number of electrons that has been discharged
 * from a fully charged state.
 *
 *  @param  relativeExtentRxn  Relative extent of reaction variable (input)
 */
void  Electrode_MP_RxnExtent::setState_relativeExtentRxn(double relativeExtentRxn)
{
    if (pendingIntegratedStep_) {
        throw CanteraError(" Electrode_MP_RxnExtent::setState_relativeExtentRxn()",
                           " pending integration");
    }
    RelativeExtentRxn_final_ = relativeExtentRxn;

    int ll = RegionBoundaries_ExtentRxn_.size();
    if (RelativeExtentRxn_final_ > RegionBoundaries_ExtentRxn_[ll-1]) {
        throw CanteraError(" Electrode_MP_RxnExtent::setCapacityDischarged()",
                           " Too much capacity in electrode ");
    }
    if (RelativeExtentRxn_final_ < 0.0) {
        throw CanteraError(" Electrode_MP_RxnExtent::setCapacityDischarged()",
                           " Negative capacity in electrode ");
    }
    xRegion_final_ = findRegion(RelativeExtentRxn_final_);
    if (spMoles_FeS2_Normalization_ <= 0.0) {
        int is_FeS2_A = globalSpeciesIndex("FeS2_A(S)");
        int is_FeS2_B = globalSpeciesIndex("FeS2_B(S)");
        spMoles_FeS2_Normalization_ = spMoles_final_[is_FeS2_A] + spMoles_final_[is_FeS2_B];;
    }

    electronKmolDischargedToDate_ = - RelativeExtentRxn_final_ * spMoles_FeS2_Normalization_;

    RelativeExtentRxn_init_  = RelativeExtentRxn_final_;
    RelativeExtentRxn_init_init_ = RelativeExtentRxn_final_;
    RelativeExtentRxn_final_final_ = RelativeExtentRxn_final_;
    xRegion_init_ =  xRegion_final_;
    xRegion_init_init_ = xRegion_final_;
    xRegion_final_final_ = xRegion_final_;

    onRegionBoundary_final_ = -1;
    for (int i = 0; i < (int) RegionBoundaries_ExtentRxn_.size(); i++) {
        if (fabs(RelativeExtentRxn_final_ - RegionBoundaries_ExtentRxn_[i]) < 1.0E-12) {
            onRegionBoundary_init_ = i;
            onRegionBoundary_final_ = i;
            break;
        }
    }
    /*
     *  Now go through an update State
     */
    updateState();

    setInitStateFromFinal(true);

    Radius_internal_init_ = Radius_internal_final_;
    Radius_internal_init_init_ = Radius_internal_final_;
    Radius_internal_final_final_ = Radius_internal_final_;

}
//====================================================================================================================
// This is used to set Phase information that is implicit but not set by a restart or an initialization
/*
 *  (virtual function from Electrode)
 *
 *  Extra informat that may be needed in advance of a successful updateState() call that specifies all of the
 *  information in the state
 *
 *  @param flagErrors If true any changes in the current flags caused by a mismatch between the state
 *                    and the values of the flags will cause an error exit.
 */
bool  Electrode_MP_RxnExtent::stateToPhaseFlagsReconciliation(bool flagErrors)
{
    bool retn = Electrode_Integrator::stateToPhaseFlagsReconciliation(flagErrors);
    xRegion_final_ = findRegion(RelativeExtentRxn_final_);
    setState_relativeExtentRxn(RelativeExtentRxn_final_);
    return retn;
}


//====================================================================================================================
//  Return the number of extra print tables
int Electrode_MP_RxnExtent::getNumPrintTables() const
{
    return 1;
}
//====================================================================================================================
// Get the values that are printed in tables for the 1D code.
/*
 *   @param itable    table id
 *   @param colNames   string names of the header (length is the length of the column)
 *   @param colValues    Value of the columns (length is the length of the column)
 */
void  Electrode_MP_RxnExtent::getPrintTable(int itable, std::vector<std::string>& colNames,
        std::vector<double>& colValues) const
{
    colNames.clear();
    colValues.clear();
    string tmp = "ExtRxn_GI";
    colNames.push_back(tmp);
    tmp = "ExtRxn_GF";
    colNames.push_back(tmp);
    colValues.push_back(RelativeExtentRxn_init_init_);
    colValues.push_back(RelativeExtentRxn_final_final_);
}

//====================================================================================================================
//! Possibly change the solution due to phase births and deaths after phase has been accepted.
/*!
 *   (virtual from Electrode_Integrator)
 *
 *  This routine is carried out after the step is deemed a success. Massaging of the solution
 *  must be carried out within strict tolerances.
 */
void  Electrode_MP_RxnExtent::manageBirthDeathSuccessfulStep()
{
    checkRegion(xRegion_final_);

    if (onRegionBoundary_final_ >= 0) {
        if (fabs(RelativeExtentRxn_final_  - RegionBoundaries_ExtentRxn_[onRegionBoundary_final_]) > 1.0E-5) {
            throw CanteraError("Electrode_MP_RxnExtent::integrate() ERROR: cell " + int2str(electrodeCellNumber_) +
                               ", integrationCounter " + int2str(counterNumberIntegrations_),
                               " RelativeExtentRxn_final_ not on boundary and onRegionBoundary_final_ is");
        } else {
            RelativeExtentRxn_final_  = RegionBoundaries_ExtentRxn_[onRegionBoundary_final_];
        }

    }

    if (deltaTsubcycle_ <= 0.0) {
        throw CanteraError("Electrode_MP_RxnExtent::manageBirthDeathSuccessfulStep() ERROR: cell " + int2str(electrodeCellNumber_) +
                           ", counterNumberIntegrations = "  + int2str(counterNumberIntegrations_) +
                           ", counterNumberSubIntegrations " + int2str(counterNumberSubIntegrations_),
                           "Negative or zero deltaTsubcycle_ = "     + fp2str(deltaTsubcycle_));
    }
    if (deltaTsubcycleCalc_ <= 0.0) {
        throw CanteraError("Electrode_MP_RxnExtent::manageBirthDeathSuccessfulStep() ERROR: cell " + int2str(electrodeCellNumber_) +
                           ", counterNumberIntegrations = "    + int2str(counterNumberIntegrations_) +
                           ", counterNumberSubIntegrations = " + int2str(counterNumberSubIntegrations_),
                           "Negative or zero deltaTsubcycleCalc_ = " + fp2str(deltaTsubcycleCalc_));
    }
    /*
     *   In this routine, we renormalize the fake species moles
     */
    int ip_FeS2_A = globalPhaseIndex("FeS2_A(S)");
    int is_FeS2_A = globalSpeciesIndex("FeS2_A(S)");
    int ip_FeS2_B = globalPhaseIndex("FeS2_B(S)");
    int is_FeS2_B = globalSpeciesIndex("FeS2_B(S)");
    spMoles_final_[is_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
    spMoles_final_[is_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
    phaseMoles_final_[ip_FeS2_A] = spMoles_FeS2_Normalization_ / 2.0;
    phaseMoles_final_[ip_FeS2_B] = spMoles_FeS2_Normalization_ / 2.0;
}

//====================================================================================================================
//   Error check on the routine step
/*
 *    (virtual from Electrode_Integrator)
 *
 *   Error checks go here. All errors are fatal exits.
 */
void Electrode_MP_RxnExtent::check_final_state()
{
}
//====================================================================================================================
//   Check the nonlinear residual equations for completeness and the ability to be solved
/*
 *
 */
int Electrode_MP_RxnExtent::check_nonlinResidConditions()
{
    return 0;
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
void Electrode_MP_RxnExtent::printElectrode(int pSrc, bool subTimeStep)
{

    double egv = TotalVol();
    double tsm = SolidTotalMoles();

    /*
     * We print reaction info here, so need to call this
     */
    updateState();
    extractInfo();

    printf("   ==============================================================================================\n");
    if (subTimeStep) {
        printf("      Electrode at intermediate-step time final = %14.7E\n", tfinal_);
        printf("                   intermediate-step time init  = %14.7E                deltaT = %14.7E\n",
               tinit_, deltaTsubcycleCalc_);
    } else {
        printf("      Electrode at time final = %14.7E\n", t_final_final_);
        printf("                   time init  = %14.7E                         deltaTglobal = %14.7E\n", t_init_init_,
               t_final_final_ - t_init_init_);
    }
    printf("\n");
    printf("                ChemModel Type = %3d , DomainNumber = %2d , CellNumber = %2d , SubIntegrationCounter = %d\n",
           electrodeChemistryModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberSubIntegrations_);
    printf("   ==============================================================================================\n");
    printf("          Voltage = %g volts\n", deltaVoltage_);
    printf("          Number of external surfaces = %d\n", (int) numExternalInterfacialSurfaces_);
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
        printf("          Molar Volume of Solid = %11.4E cm3 gmol-1\n",  ElectrodeSolidVolume_ / tsm * 1.0E3);
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


    double capacd = capacityDischarged();
    printf("          Capacity Discharged (total)  = %12.6g coulombs = %12.6g Ah\n", capacd, capacd / 3600.);
    double dod = depthOfDischarge();
    printf("          Depth of Discharge           = %12.6g coulombs\n", dod);
    double capLeft = capacityLeft();
    double capZero = capacity();
    printf("          Capacity Left                = %12.6g coulombs  \n",   capLeft);
    printf("          Capacity at Zero DOD         = %12.6g coulombs\n",  capZero);
    if (subTimeStep) {
        printf("          Electron Src This int step   = %12.6g coulombs\n",  spMoleIntegratedSourceTermLast_[kElectron_] * Faraday);
    } else {
        printf("          Electron Src This Period     = %12.6g coulombs\n",  spMoleIntegratedSourceTerm_[kElectron_] * Faraday);
    }
    double invDelT = 1.0;
    if (! subTimeStep) {
        if (t_final_final_ > t_init_init_) {
            invDelT = 1.0/ (t_final_final_ - t_init_init_);
        }
    } else {
        if (tfinal_ > tinit_) {
            invDelT = 1.0/ (tfinal_ - tinit_);
        }
    }
    double amps = invDelT * spMoleIntegratedSourceTerm_[kElectron_] * Faraday;
    printf("          Amps         This Period     = %12.6g coulombs sec-1\n", amps);

    printf("\n");
    printf("          followElectrolyteMoles = %d\n", followElectrolyteMoles_);
    printf("          ElectrolytePseudoMoles = %g\n", electrolytePseudoMoles_);

    printf("                                                   Initial       Final\n");
    if (subTimeStep) {
        printf("          RelativeExtent:                           %- 12.5g % -12.5g\n",  RelativeExtentRxn_init_,
               RelativeExtentRxn_final_);
        printf("          Region:                                %5d        %5d\n",  xRegion_init_, xRegion_final_);
        printf("          onRegionBound:                         %5d        %5d\n",  onRegionBoundary_init_, onRegionBoundary_final_);
    } else {
        printf("          RelativeExtent:                     %- 12.5g %- 12.5g\n",  RelativeExtentRxn_init_init_,
               RelativeExtentRxn_final_final_);
        printf("          Region:                                %5d        %5d\n",  xRegion_init_init_, xRegion_final_final_);
    }
    if (onRegionBoundary_final_ == -1) {
        double ooss = openCircuitVoltageSSRxn(indexOfReactingSurface_);
        printf("          OpenCircuitVoltageSS                                   % -12.7g \n", ooss);
        printf("          OpenCircuitVoltage                                     % -12.7g \n", openCircuitVoltage(indexOfReactingSurface_));
    } else {

        if ((onRegionBoundary_final_-1)< 0) {
            printf("          OpenCircuitVoltageSS_top                             NA \n");
            printf("          OpenCircuitVoltage_top                               NA \n");
        } else {
            printf("          OpenCircuitVoltageSS_top                               %- 12.7g \n",
                   openCircuitVoltageSS_Region(RelativeExtentRxn_final_, onRegionBoundary_final_-1));
            printf("          OpenCircuitVoltage_top                                 %- 12.7g \n",
                   openCircuitVoltage_Region(RelativeExtentRxn_final_, onRegionBoundary_final_-1));
        }
        if (onRegionBoundary_final_ >= (int)(RegionBoundaries_ExtentRxn_.size() - 1)) {
            printf("          OpenCircuitVoltageSS_bot                               NA \n");
            printf("          OpenCircuitVoltage_bot                                 NA \n");
        } else {
            printf("          OpenCircuitVoltageSS_bot                               %- 12.7g \n",
                   openCircuitVoltageSS_Region(RelativeExtentRxn_final_, onRegionBoundary_final_));
            printf("          OpenCircuitVoltage_bot                                 %- 12.7g \n",
                   openCircuitVoltage_Region(RelativeExtentRxn_final_, onRegionBoundary_final_));
        }
    }
    printf("          ElectrodeVoltage                                       %- 12.7g\n", deltaVoltage_);


    ReactingSurDomain* rsd = RSD_List_[indexOfReactingSurface_];
    int irxn = 0;
    double dStoich;
    double OCV;
    double io;
    double nu;
    double beta, resist;
#ifdef DONOTREMOVE
    double icurr = rsd->getExchangeCurrentDensityFormulation(irxn, &dStoich, &OCV, &io, &nu, &beta, &resist);
#else
    bool okk = rsd->getExchangeCurrentDensityFormulation(irxn, dStoich, OCV, io, nu, beta, resist);
    double icurr = 0.0;
    if (okk) {
	icurr = rsd->calcCurrentDensity(nu, dStoich, io, beta, temperature_, resist);
    }
#endif


    printf("          OCV from (rsd)                                         %- 12.7g Volts\n", OCV);
    printf("          i_o                                                    %- 12.7E coul/sec/m2\n", io);
    printf("          nStoich                                                %- 12.7g\n", dStoich);
    printf("          nu                                                     %- 12.7g Volts\n", nu);
    printf("          beta                                                   %- 12.7g beta\n", beta);
    printf("          i_curr                                                 %- 12.7E coul/sec/m2\n", icurr);


#ifdef DEBUG_MODE_NOT
    double beta = rsd->electrochem_beta(0);
    double RT = GasConstant * temperature_;

    double ee1 = beta * nu * Faraday *dStoich / RT;
    double ee2 = - (1 - beta) * nu * Faraday * dStoich / RT;

    double iii = io * (exp(ee1) - exp(ee2));
    double iiM = iii / Faraday;

    double netROP[10];
    rsd->getNetRatesOfProgress(DATA_PTR(netROP));

    printf("                                         (compare net ROP from two formulations: %g  %g)\n", -iiM, netROP[0]);

#endif

    if (solidDiffusionModel_) {
        if (ROPModificationType_ == 1) {
            printf("          Electrode reaction modified by an inner region Damkoeler Number:\n");
            printf("                   EqnSystem_limitingEquationBehavior_              %3d\n", limitingEquationBehavior_);
            printf("                   DaInner                            %12.4E\n", DaInner_);
            printf("                   D(xregion)                         %12.4E m2 s-1\n", diffusionCoeffRegions_[xRegion_final_+1]);
            printf("                   MolarVol(xregion)                  %12.4E m3 kmol-1\n", molarVolumeRegions_[xRegion_final_+1]);
            printf("                   kfAB                               %12.4E kmol m-2 s-1\n", kfAB_);
            printf("                   krAB                               %12.4E kmol m-2 s-1\n", krAB_);
            printf("                   ca_Li+                             %12.4E\n", ca_Lip_);
            printf("                   kfExt_                             %12.4E kmol m-2 s-1\n", kfExt_);
            printf("                   krExt_                             %12.4E kmol m-2 s-1\n", krExt_);
            printf("                   kfInner_                           %12.4E kmol m-2 s-1\n", kfInner_);
            printf("                   krInner_                           %12.4E kmol m-2 s-1\n", krInner_);
            printf("                   Lin_                               %12.4E\n", Lin_);
            printf("                   Lout_                              %12.4E\n", Lout_);
            printf("                   Radius_ext_final                   %12.4E m\n", Radius_exterior_final_);
            printf("                   Radius_int_final                   %12.4E m\n", Radius_internal_final_);
            if (kfExt_ < kfInner_) {
                printf("                      WARNING, kfExt is less than kfInner, violating one of the assumptions\n");
            }
            if (Lin_ > 1.0) {
                printf("                      WARNING, Lin is greater than 1.0, violating one of the assumptions\n");
            }
        } else if (ROPModificationType_ == 2) {
            printf("          Electrode reaction modified by an outer region Damkoeler Number:\n");
            printf("                   EqnSystem_limitingEquationBehavior_              %3d\n", limitingEquationBehavior_);
            printf("                   DaOuter                            %12.4E\n", DaOuter_);
            printf("                   D(xregion)                         %12.4E m2 s-1\n", diffusionCoeffRegions_[xRegion_final_+1]);
            printf("                   MolarVol(xregion)                  %12.4E m3 kmol-1\n", molarVolumeRegions_[xRegion_final_+1]);
            printf("                   kfAB                               %12.4E kmol m-2 s-1\n", kfAB_);
            printf("                   krAB                               %12.4E kmol m-2 s-1\n", krAB_);
            printf("                   ca_Li+                             %12.4E\n", ca_Lip_);
            printf("                   kfExt_                             %12.4E kmol m-2 s-1\n", kfExt_);
            printf("                   krExt_                             %12.4E kmol m-2 s-1\n", krExt_);
            printf("                   kfInner_                           %12.4E kmol m-2 s-1\n", kfInner_);
            printf("                   krInner_                           %12.4E kmol m-2 s-1\n", krInner_);
            printf("                   Lin_                               %12.4E\n", Lin_);
            printf("                   Lout_                              %12.4E\n", Lout_);
            printf("                   Radius_ext_final                   %12.4E m\n", Radius_exterior_final_);
            printf("                   Radius_int_final                   %12.4E m\n", Radius_internal_final_);
            if (kfExt_ > kfInner_) {
                printf("                      WARNING, kfExt is greater than kfInner, violating one of the assumptions\n");
            }
            if (Lout_ > 1.0) {
                printf("                      WARNING, Lout is greater than 1.0, violating one of the assumptions\n");
            }

        }
    } else {
        printf("          Solid Phase Diffusion Model:                            NONE\n");
    }

    printElectrodePhaseList(pSrc, subTimeStep);

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        printElectrodePhase(iph, pSrc, subTimeStep);
    }
}

//====================================================================================================================
void Electrode_MP_RxnExtent::printElectrodePhase(int iphI, int pSrc, bool subTimeStep)
{
    size_t iph = iphI;
    int printDetail = 10;
    //int isph;
    int isurfA;
    int isurf=-1;
    double* netROP = new double[m_NumTotSpecies];
    ThermoPhase& tp = thermo(iph);
    string pname = tp.id();
    int istart = m_PhaseSpeciesStartIndex[iph];
    int nsp = tp.nSpecies();
    printf("     ============================================================================================\n");
    printf("          PHASE %d %s \n", iphI, pname.c_str());
    printf("                Total moles = %g\n", phaseMoles_final_[iph]);
    double mv = tp.molarVolume();
    if (iph >= NumVolPhases_) {
        printf("                Molar Volume = %11.5E cm3 gmol-1\n", 0.0);
        printf("                Molar Area   = %11.5E cm2 gmol-1\n", mv * 10.);
    } else {
        printf("                Molar Volume = %11.5E cm3 gmol-1\n", mv * 1.0E3);
    }

    if (iphI == metalPhase_) {
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
    if (iphI == metalPhase_ || iphI == solnPhase_) {
        printf("                Electric Potential = %g volts\n", tp.electricPotential());
    }
    /*
     * Do specific surface phase printouts
     */
    double radius;
    if (iph >= NumVolPhases_) {
        //isph = iph - NumVolPhases_;
        if (indexOfReactingSurface_ == 1) {
            isurf = 1;
        } else {
            isurf = 0;
        }
        if (locationOfReactingSurface_ == 0) {
            isurfA = 2;
            printf("                FOR PURPOSES OF RXN, Using these constant surface areas:\n");
        } else if (locationOfReactingSurface_ == 1) {
            isurfA = 0;
            printf("                FOR PURPOSES OF RXN, Using the internal surface area:\n");
        } else if (locationOfReactingSurface_ == 2) {
            printf("                FOR PURPOSES OF RXN, Using the external surface area:\n");
            isurfA = 1;
        } else {
          throw CanteraError("Electrode_MP_RxnExtent::printElectrodePhase", "locationOfReactingSurface_ != 0,1,2");
        }

        radius = sqrt(surfaceAreaRS_final_[isurfA] / (4.0 * Pi * particleNumberToFollow_));
        radius *= 1.0E6;
        printf("                surface area (final) = %11.5E m**2      Radius(final) = %11.5E um\n",  surfaceAreaRS_final_[isurfA], radius);
        if (subTimeStep) {
            radius = sqrt(surfaceAreaRS_init_[isurfA] / (4.0 * Pi * particleNumberToFollow_));
            radius *= 1.0E6;
            printf("                surface area (init)  = %11.5E m**2      Radius(init)  = %11.5E um\n",  surfaceAreaRS_init_[isurfA], radius);
        } else {
            radius = sqrt(surfaceAreaRS_init_init_[isurfA] / (4.0 * Pi * particleNumberToFollow_));
            radius *= 1.0E6;
            printf("                surface area (init)  = %11.5E m**2      Radius(init)  = %11.5E um\n",  surfaceAreaRS_init_init_[isurfA], radius);
        }

        int ttt = isExternalSurface_[isurfA];
        printf("                IsExternalSurface = %d\n", ttt);
        double oc = openCircuitVoltage(isurf);
        if (oc != 0.0) {
            printf("                 Open Circuit Voltage = %g Volts\n", oc);
        }
    }
    if (printDetail > 3) {
        printf("\n");
        printf("                Name              MoleFrac_final kMoles_final kMoles_init SrcTermIntegrated(kmol)\n");
        for (int k = 0; k < nsp; k++) {
            string sname = tp.speciesName(k);
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
    if (pSrc) {
        if (iph >= NumVolPhases_) {
            const vector<double>& rsSpeciesProductionRates = RSD_List_[isurf]->calcNetSurfaceProductionRateDensities();
            RSD_List_[isurf]->getNetRatesOfProgress(netROP);

            double* spNetProdPerArea = (double*) spNetProdPerArea_List_.ptrColumn(isurf);
            std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
            int nphRS = RSD_List_[isurf]->nPhases();
            int kIndexKin = 0;
            for (int kph = 0; kph < nphRS; kph++) {
                int jph = RSD_List_[isurf]->kinOrder[kph];
                int istart = m_PhaseSpeciesStartIndex[jph];
                int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
                if (goNowhere_) {
                    for (int k = 0; k < nsp; k++) {
                        spNetProdPerArea[istart + k] = 0.0;
                        kIndexKin++;
                    }
                } else {
                    for (int k = 0; k < nsp; k++) {
                        spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                        kIndexKin++;
                    }
                }
            }
            printf("\n");
            printf("                           spName                  SourceRateLastStep (kmol/m2/s) \n");
            for (size_t k = 0; k <  m_NumTotSpecies; k++) {
                string ss = speciesName(k);
                printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
            }
        }
    }
    printf("     ============================================================================================\n");
    delete [] netROP;
}
//==================================================================================================================================
} // End of namespace
//----------------------------------------------------------------------------------------------------------------------------------
