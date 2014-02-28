/*
 * $Id: Electrode_RadialRegion.cpp 298 2012-08-08 20:15:48Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"

#include "Electrode_RadialRegion.h"
#include "cantera/integrators.h"

#include "Electrode_RadialDiffRegions.h"

#include "BlockEntryGlobal.h"

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

namespace Cantera
{

//====================================================================================================================
ELECTRODE_RadialRegion_KEY_INPUT::ELECTRODE_RadialRegion_KEY_INPUT(int printLvl) :
    ELECTRODE_KEY_INPUT(printLvl),
    indexRegion_(0),
    numRadialCells_(3),
    phaseIndeciseKRsolidPhases_(0),
    solidDiffusionModel_(0),
    diffusionCoeffSpecies_(0),
    defaultDiffusionCoeff_(1.0E-12)
{

}
//====================================================================================================================
ELECTRODE_RadialRegion_KEY_INPUT::~ELECTRODE_RadialRegion_KEY_INPUT()
{

}
//====================================================================================================================
ELECTRODE_RadialRegion_KEY_INPUT::ELECTRODE_RadialRegion_KEY_INPUT(const ELECTRODE_RadialRegion_KEY_INPUT& right) :
    ELECTRODE_KEY_INPUT(right),
    indexRegion_(right.indexRegion_),
    numRadialCells_(right.numRadialCells_),
    phaseIndeciseKRsolidPhases_(right.phaseIndeciseKRsolidPhases_),
    solidDiffusionModel_(right.solidDiffusionModel_),
    diffusionCoeffSpecies_(0),
    defaultDiffusionCoeff_(right.defaultDiffusionCoeff_)

{
    diffusionCoeffSpecies_          = right.diffusionCoeffSpecies_;
}
//====================================================================================================================
ELECTRODE_RadialRegion_KEY_INPUT&
ELECTRODE_RadialRegion_KEY_INPUT::operator=(const ELECTRODE_RadialRegion_KEY_INPUT& right)
{
    if (this == &right) {
        return *this;
    }

    ELECTRODE_KEY_INPUT::operator=(right);

    indexRegion_                    = right.indexRegion_;
    numRadialCells_                 = right.numRadialCells_;
    phaseIndeciseKRsolidPhases_     = right.phaseIndeciseKRsolidPhases_;
    solidDiffusionModel_            = right.solidDiffusionModel_;
    diffusionCoeffSpecies_          = right.diffusionCoeffSpecies_;
    defaultDiffusionCoeff_          = right.defaultDiffusionCoeff_;

    return *this;
}
//====================================================================================================================
void ELECTRODE_RadialRegion_KEY_INPUT::setup_input_child1(BEInput::BlockEntry* cf)
{
    /*
     * Obtain the number of regions
     */
    LE_OneInt* s1 = new LE_OneInt("Index of the Region", &(indexRegion_), 0, "indexRegion");
    s1->set_default(0);
    cf->addLineEntry(s1);

    /* --------------------------------------------------------------
     *   numRadialCells
     *  Input the number of cells in the region
     */
    LE_OneInt* nRm = new LE_OneInt("Number of Cells in Region", &(numRadialCells_), 0, "numRadialCells");
    nRm->set_default(3);
    cf->addLineEntry(nRm);

    /*
     * Obtain the number of regions
     */


    LE_OneInt* sdm = new LE_OneInt("Solid Diffusion Model", &(solidDiffusionModel_), 0, "solidDiffusionModel");
    sdm->set_default(0);
    cf->addLineEntry(sdm);
    BaseEntry::set_SkipUnknownEntries(true);

    /*
     * NAME OF THE THERMOPHASE PHASE ASSOCIATED WITH THE REGION
     * NOTE WE WILL EXPAND THIS WHEN WE GO TO MULTIPLE PHASE REGIONS.
     */
    LE_OneStr *pnames = new LE_OneStr("Phase Names within Distributed region", &phaseName_, 1, 1, 0, "phaseNames");

      //  LE_MultiCStr* pnames = new LE_MultiCStr("Phase Names within Distributed region", 1, 10, 1, 1, "phaseNames");

    cf->addLineEntry(pnames);

    BaseEntry::set_SkipUnknownEntries(true);

}
//====================================================================================================================
void ELECTRODE_RadialRegion_KEY_INPUT::setup_input_child2(BEInput::BlockEntry* cf)
{

    /*
     * Obtain the number of regions
     */
/*
    int reqd = 0;
    if (solidDiffusionModel_) {
        reqd = 1;
    }
*/
    /*
     * Make the diffusion coefficients optional unless the solid diffusion model is turned on
     */
/*
    reqd = 0;
    if (solidDiffusionModel_) {
        reqd = 1;
    }
*/
    LE_OneInt* iR = new LE_OneInt("Index of the Region",  &(indexRegion_));
    iR->set_default(0);
    cf->addLineEntry(iR);

  
    BaseEntry::set_SkipUnknownEntries(false);
}

//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_RadialRegion::Electrode_RadialRegion() :
    Electrode_Integrator(),
    numKRSpecies_(1),
    numRCells_(5),
    numSPhases_(1),
    numEqnsCell_(0),

    phaseIndeciseKRsolidPhases_(0),
    KRsolid_speciesList_(0),

    spMoles_KRsolid_Cell_final_(0),
    spMoles_KRsolid_Cell_init_(0),
    spMoles_KRsolid_Cell_final_final_(0),
    spMoles_KRsolid_Cell_init_init_(0),

    spMf_KRsolid_Cell_final_(0),
    spMf_KRsolid_Cell_init_(0),
    spMf_KRsolid_Cell_final_final_(0),
    spMf_KRsolid_Cell_init_init_(0),

    concTot_SPhase_Cell_final_(0),
    concTot_SPhase_Cell_init_(0),
    concTot_SPhase_Cell_final_final_(0),
    concTot_SPhase_Cell_init_init_(0),

    concKRSpecies_Cell_final_(0),
    concKRSpecies_Cell_init_(0),
    concKRSpecies_Cell_final_final_(0),
    concKRSpecies_Cell_init_init_(0),

    molarVolume_refLat_Cell_final_(0),
    molarVolume_refLat_Cell_init_(0),
    molarVolume_refLat_Cell_final_final_(0),
    molarVolume_refLat_Cell_init_init_(0),


    MolarVolume_refLat_Ref_(55.55),

    rnodePos_final_(0),
    rnodePos_init_(0),
    rnodePos_final_final_(0),
    rnodePos_init_init_(0),

    rRefPos_final_(0),
    rRefPos_init_(0),
    rRefPos_final_final_(0),
    rRefPos_init_init_(0),

    cellBoundR_final_(0),
    cellBoundR_init_(0),
    cellBoundR_final_final_(0),
    cellBoundR_init_init_(0),

    partialMolarVolKRSpecies_Cell_final_(0),
    partialMolarVolKRSpecies_Cell_init_(0),
    partialMolarVolKRSpecies_Cell_final_final_(0),
    partialMolarVolKRSpecies_Cell_init_init_(0),

    fractionVolExpansion_Cell_(0),
    DspMoles_RightSurf_final_(0),


    radiusLeft_Ref_(0.0),
    radiusLeft_final_(0.0),
    radiusLeft_init_(0.0),
    radiusLeft_final_final_(0.0),
    radiusLeft_init_init_(0.0),

    Diff_Coeff_KRSolid_(0),
    DphMoles_RightSurf_final_(0),

    surfIndex_RightSurf_(-1),
    actCoeff_Cell_final_(0)
{


}

//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_RadialRegion::Electrode_RadialRegion(const Electrode_RadialRegion& right) :
    Electrode_Integrator(),
    numKRSpecies_(1),
    numRCells_(5),
    numSPhases_(1),
    numEqnsCell_(0),

    phaseIndeciseKRsolidPhases_(0),
    KRsolid_speciesList_(0),

    spMoles_KRsolid_Cell_final_(0),
    spMoles_KRsolid_Cell_init_(0),
    spMoles_KRsolid_Cell_final_final_(0),
    spMoles_KRsolid_Cell_init_init_(0),

    spMf_KRsolid_Cell_final_(0),
    spMf_KRsolid_Cell_init_(0),
    spMf_KRsolid_Cell_final_final_(0),
    spMf_KRsolid_Cell_init_init_(0),


    concTot_SPhase_Cell_final_(0),
    concTot_SPhase_Cell_init_(0),
    concTot_SPhase_Cell_final_final_(0),
    concTot_SPhase_Cell_init_init_(0),

    concKRSpecies_Cell_final_(0),
    concKRSpecies_Cell_init_(0),
    concKRSpecies_Cell_final_final_(0),
    concKRSpecies_Cell_init_init_(0),

    molarVolume_refLat_Cell_final_(0),
    molarVolume_refLat_Cell_init_(0),
    molarVolume_refLat_Cell_final_final_(0),
    molarVolume_refLat_Cell_init_init_(0),

    MolarVolume_refLat_Ref_(55.55),

    rnodePos_final_(0),
    rnodePos_init_(0),
    rnodePos_final_final_(0),
    rnodePos_init_init_(0),

    rRefPos_final_(0),
    rRefPos_init_(0),
    rRefPos_final_final_(0),
    rRefPos_init_init_(0),

    cellBoundR_final_(0),
    cellBoundR_init_(0),
    cellBoundR_final_final_(0),
    cellBoundR_init_init_(0),


    partialMolarVolKRSpecies_Cell_final_(0),
    partialMolarVolKRSpecies_Cell_init_(0),
    partialMolarVolKRSpecies_Cell_final_final_(0),
    partialMolarVolKRSpecies_Cell_init_init_(0),

    fractionVolExpansion_Cell_(0),
    DspMoles_RightSurf_final_(0),



    radiusLeft_Ref_(0.0),
    radiusLeft_final_(0.0),
    radiusLeft_init_(0.0),
    radiusLeft_final_final_(0.0),
    radiusLeft_init_init_(0.0),

    Diff_Coeff_KRSolid_(0),
    DphMoles_RightSurf_final_(0),

    surfIndex_RightSurf_(-1),
    actCoeff_Cell_final_(0)


{
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_RadialRegion&
Electrode_RadialRegion::operator=(const Electrode_RadialRegion& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode_Integrator::operator=(right);

    numKRSpecies_                       = right.numKRSpecies_;
    numRCells_                          = right.numRCells_;
    numSPhases_                         = right.numSPhases_;
    numEqnsCell_                        = right.numEqnsCell_;

    phaseIndeciseKRsolidPhases_         = right.phaseIndeciseKRsolidPhases_;
    KRsolid_speciesList_                = right.KRsolid_speciesList_;

    spMoles_KRsolid_Cell_final_         = right.spMoles_KRsolid_Cell_final_;
    spMoles_KRsolid_Cell_init_          = right.spMoles_KRsolid_Cell_init_;
    spMoles_KRsolid_Cell_final_final_   = right.spMoles_KRsolid_Cell_final_final_;
    spMoles_KRsolid_Cell_init_init_     = right.spMoles_KRsolid_Cell_init_init_;

    spMf_KRsolid_Cell_final_            = right.spMf_KRsolid_Cell_final_;
    spMf_KRsolid_Cell_init_             = right.spMf_KRsolid_Cell_init_;
    spMf_KRsolid_Cell_final_final_      = right.spMf_KRsolid_Cell_final_final_;
    spMf_KRsolid_Cell_init_init_        = right.spMf_KRsolid_Cell_init_init_;



    concTot_SPhase_Cell_final_          = right.concTot_SPhase_Cell_final_;
    concTot_SPhase_Cell_init_           = right.concTot_SPhase_Cell_init_;
    concTot_SPhase_Cell_final_final_    = right.concTot_SPhase_Cell_final_final_;
    concTot_SPhase_Cell_init_init_      = right.concTot_SPhase_Cell_init_init_;

    concKRSpecies_Cell_final_           = right.concKRSpecies_Cell_final_;
    concKRSpecies_Cell_init_            = right.concKRSpecies_Cell_init_;
    concKRSpecies_Cell_final_final_     = right.concKRSpecies_Cell_final_final_;
    concKRSpecies_Cell_init_init_       = right.concKRSpecies_Cell_init_init_;


    molarVolume_refLat_Cell_final_      = right.molarVolume_refLat_Cell_final_;
    molarVolume_refLat_Cell_init_       = right.molarVolume_refLat_Cell_init_;
    molarVolume_refLat_Cell_final_final_= right.molarVolume_refLat_Cell_final_final_;
    molarVolume_refLat_Cell_init_init_  = right.molarVolume_refLat_Cell_init_init_;

    MolarVolume_refLat_Ref_             = right.MolarVolume_refLat_Ref_;

    rnodePos_final_                     = right.rnodePos_final_;
    rnodePos_init_                      = right.rnodePos_init_;
    rnodePos_final_final_               = right.rnodePos_final_final_;
    rnodePos_init_init_                 = right.rnodePos_init_init_;

    rRefPos_final_                      = right.rRefPos_final_;
    rRefPos_init_                       = right.rRefPos_init_;
    rRefPos_final_final_                = right.rRefPos_final_final_;
    rRefPos_init_init_                  = right.rRefPos_init_init_;

    cellBoundR_final_                   = right.cellBoundR_final_;
    cellBoundR_init_                    = right.cellBoundR_init_;
    cellBoundR_final_final_             = right.cellBoundR_final_final_;
    cellBoundR_init_init_               = right.cellBoundR_init_init_;

    partialMolarVolKRSpecies_Cell_final_        = right.partialMolarVolKRSpecies_Cell_final_;
    partialMolarVolKRSpecies_Cell_init_         = right.partialMolarVolKRSpecies_Cell_init_;
    partialMolarVolKRSpecies_Cell_final_final_  = right.partialMolarVolKRSpecies_Cell_final_final_;
    partialMolarVolKRSpecies_Cell_init_init_    = right.partialMolarVolKRSpecies_Cell_init_init_;


    fractionVolExpansion_Cell_             = right.fractionVolExpansion_Cell_;
    DspMoles_RightSurf_final_           = right.DspMoles_RightSurf_final_;



    radiusLeft_Ref_                     = right.radiusLeft_Ref_;
    radiusLeft_final_                   = right.radiusLeft_final_;
    radiusLeft_init_                    = right.radiusLeft_init_;
    radiusLeft_final_final_             = right.radiusLeft_final_final_;
    radiusLeft_init_init_               = right.radiusLeft_init_init_;

    Diff_Coeff_KRSolid_                 = right.Diff_Coeff_KRSolid_;
    DphMoles_RightSurf_final_           = right.DphMoles_RightSurf_final_;

    surfIndex_RightSurf_                = right.surfIndex_RightSurf_;

    KRsolid_speciesList_                = right.KRsolid_speciesList_;



    actCoeff_Cell_final_                = right.actCoeff_Cell_final_;

    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
/*
 *
 * :destructor
 *
 * We need to manually free all of the arrays.
 */
Electrode_RadialRegion::~Electrode_RadialRegion()
{


}
//======================================================================================================================
//    Return the type of electrode
/*
 *  Returns the enum type of the electrode. This is used in the factory routine.
 *
 *  @return Returns an enum type, called   Electrode_Types_Enum
 */
Electrode_Types_Enum Electrode_RadialRegion::electrodeType() const
{
    return SIMPLE_DIFF_ET;
}
//======================================================================================================================
int
Electrode_RadialRegion::electrode_model_create(ELECTRODE_KEY_INPUT* eibase)
{

    //ee_->electrode_model_create(eibase);


    ELECTRODE_RadialRegion_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialRegion_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_RadialRegion::electrode_model_create()",
                           " Expecting a child ELECTRODE_KEY_INPUT object and didn't get it");
    }

    /*
     * Number of cells - hard code for now
     */
    numRCells_ = 5;


    string phaseName = ei->phaseName_;
    /*
     * Find phase name
     */
    int index = ee_->globalPhaseIndex(phaseName);

    phaseIndeciseKRsolidPhases_.resize(1);
    phaseIndeciseKRsolidPhases_[0] = index;

    /*
     *  Determine the number of diestributed phases
     */
    numSPhases_ = phaseIndeciseKRsolidPhases_.size();

    /*
     *  Calculate the number of equations at each node from phaseIndeciseKRsolidPhases_
     *   2 + sum (nsp_each_distrib_phase)
     */
    numKRSpecies_ = 0;
    for (int i = 0; i < (int) phaseIndeciseKRsolidPhases_.size(); i++) {
        int iPh =  phaseIndeciseKRsolidPhases_[i];
        ThermoPhase& th = thermo(iPh);
        int nsp = th.nSpecies();
        numKRSpecies_ += nsp;
    }
    numEqnsCell_ =  numKRSpecies_ + 2;

    /*
     *  Need a value for the molar volume ref
     *      We've set it here to the molar volume of water
     */
    MolarVolume_refLat_Ref_ = 55.55;


    /*
     * Initialize the arrays in this object now that we know the number of equations
     */
    init_sizes();

    /*
     *  Initialize the species
     */
    initializeAsEvenDistribution();




    return 0;
}

//=====================================================================================================================
int
Electrode_RadialRegion::setInitialConditions(ELECTRODE_KEY_INPUT* eibase)
{

    /*
     *  Downcast the Key input to make sure we are being fed the correct child object
     */
    ELECTRODE_RadialRegion_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialRegion_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_RadialRegion::electrode_model_create()",
                           " Expecting a child ELECTRODE_KEY_INPUT object and didn't get it");
    }

    Electrode::setInitialConditions(ei);

    Electrode::setInitStateFromFinal(true);

    setPhaseExistenceForReactingSurfaces();
    // printf("spMoles_final_[2] = %g\n", spMoles_final_[2]);

    neq_ = nEquations();

    /*
     *  Call a routine to
     */
    create_solvers();

    return 0;
}
//====================================================================================================================
void
Electrode_RadialRegion::init_sizes()
{
    int kspCell = numKRSpecies_ *  numRCells_;
    int kphCell = numSPhases_ * numRCells_;

    phaseIndeciseKRsolidPhases_.resize(numSPhases_, -1);
    KRsolid_speciesList_.resize(numKRSpecies_, -1);

    spMoles_KRsolid_Cell_final_.resize(kspCell, 0.0);
    spMoles_KRsolid_Cell_init_.resize(kspCell, 0.0);
    spMoles_KRsolid_Cell_final_final_.resize(kspCell, 0.0);
    spMoles_KRsolid_Cell_init_init_.resize(kspCell, 0.0);

    spMf_KRsolid_Cell_final_.resize(kspCell, 0.0);
    spMf_KRsolid_Cell_init_.resize(kspCell, 0.0);
    spMf_KRsolid_Cell_final_final_.resize(kspCell, 0.0);
    spMf_KRsolid_Cell_init_init_.resize(kspCell, 0.0);

    concTot_SPhase_Cell_final_.resize(kphCell, 0.0);
    concTot_SPhase_Cell_init_.resize(kphCell, 0.0);
    concTot_SPhase_Cell_final_final_.resize(kphCell, 0.0);
    concTot_SPhase_Cell_init_init_.resize(kphCell, 0.0);

    concKRSpecies_Cell_final_.resize(kspCell, 0.0);
    concKRSpecies_Cell_init_.resize(kspCell, 0.0);
    concKRSpecies_Cell_final_final_.resize(kspCell, 0.0);
    concKRSpecies_Cell_init_init_.resize(kspCell, 0.0);

    molarVolume_refLat_Cell_final_.resize(numRCells_, 0.0);
    molarVolume_refLat_Cell_init_.resize(numRCells_, 0.0);
    molarVolume_refLat_Cell_final_final_.resize(numRCells_, 0.0);
    molarVolume_refLat_Cell_init_init_.resize(numRCells_, 0.0);


    rnodePos_final_.resize(numRCells_, 0.0);
    rnodePos_init_.resize(numRCells_, 0.0);
    rnodePos_final_final_.resize(numRCells_, 0.0);
    rnodePos_init_init_.resize(numRCells_, 0.0);

    rRefPos_final_.resize(numRCells_, 0.0);
    rRefPos_init_.resize(numRCells_, 0.0);
    rRefPos_final_final_.resize(numRCells_, 0.0);
    rRefPos_init_init_.resize(numRCells_, 0.0);

    cellBoundR_final_.resize(numRCells_, 0.0);
    cellBoundR_init_.resize(numRCells_, 0.0);
    cellBoundR_final_final_.resize(numRCells_, 0.0);
    cellBoundR_init_init_.resize(numRCells_, 0.0);

    partialMolarVolKRSpecies_Cell_final_.resize(kspCell, 0.0);
    partialMolarVolKRSpecies_Cell_init_.resize(kspCell, 0.0);
    partialMolarVolKRSpecies_Cell_final_final_.resize(kspCell, 0.0);
    partialMolarVolKRSpecies_Cell_init_init_.resize(kspCell, 0.0);

    fractionVolExpansion_Cell_.resize(numRCells_, 0.0);
    DspMoles_RightSurf_final_.resize(m_NumTotSpecies, 0.0);

    Diff_Coeff_KRSolid_.resize(numKRSpecies_, 0.0);
    DphMoles_RightSurf_final_.resize(m_NumTotPhases, 0.0);

    actCoeff_Cell_final_.resize(kspCell, 1.0);

}
//====================================================================================================================
void
Electrode_RadialRegion::initializeAsEvenDistribution()
{
    int jRPh;
    int iPh;
    int kspStart;
    int nSpecies;
    int kSp;
    int iCell;

    /*
     *  We assume that the non-distributed quantities are pertinent and up to date.
     *  Molar volumes should be good here too, since we are assuming that everything is well mixed.
     */

    /*
     *  We calculate the volumes of the phases that we are to distribute.
     *    -> Note could use SolidVol() here but am thinking that we will be more specific about which
     *       phases are included here in the future. Therefore, it's not appropriate
     */
    double currentPhaseVols = 0.0;
    for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
        iPh = phaseIndeciseKRsolidPhases_[jRPh];
        double phaseVol = 0.0;
        int kStart = m_PhaseSpeciesStartIndex[iPh];
        ThermoPhase& tp = thermo(iPh);
        nSpecies = tp.nSpecies();
        if (iPh != solnPhase_) {
            for (int k = 0; k < nSpecies; k++) {
                phaseVol += spMoles_final_[kStart + k] * VolPM_[kStart + k];
            }
        }
        currentPhaseVols += phaseVol;
    }

    /*
     * Calculate the volume of a single particle occupied by the phases
     */
    double particleVol =  currentPhaseVols / particleNumberToFollow_;

    /*
     *  Here we calculate the inner and outer radius' of the domain to be initialized.
     *    We will assume here that the current domain is the inner core until otherwise needed
     */
    double radiusInner = 0.0;
    double radiusOuter = pow((particleVol * 3.0 / (4.0 * Pi)), 0.3333333333333333333333);

    /*
     *  Generate an initial grid between radiusInner and radiusOuter containing numRCells_.
     *      This grid is storred in rnodePos_final_[]
     */
    radialGridGenerate(numRCells_ + 1, radiusInner, radiusOuter, 1.0, rnodePos_final_);

    for (iCell = 0; iCell < numRCells_; iCell++) {
        fractionVolExpansion_Cell_[iCell] = (double) iCell / (numRCells_ - 1.0);
    }

    /*
     *  Calculate the cell boundaries in a pre - loop
     */
    double cbR3_final = 0.0;
    double cbL3_final = 0.0;
    for (iCell = 0; iCell < numRCells_; iCell++) {
        cellBoundR_final_[iCell]  = 0.5 * (rnodePos_final_[iCell] + rnodePos_final_[iCell+1]);
    }

    for (iCell = 0; iCell < numRCells_; iCell++) {
        int indexMidKRSpecies =  iCell * numKRSpecies_;
        int indexMidKRPhases =  iCell * numSPhases_;
        int kstart = 0;
        cbL3_final =  cbR3_final;
        cbR3_final = cellBoundR_final_[iCell] *  cellBoundR_final_[iCell] *  cellBoundR_final_[iCell];

        /*
         *  Calculate the volume of the cell
         */
        double volTotalCell = 4. * Pi / 3. * (cbR3_final - cbL3_final) * particleNumberToFollow_;
        /*
         *  Calculate the fraction of the total volume that is in the cell
         */
        double fracTotalVol = volTotalCell / particleVol;

        /*
         *  Loop over distributed phases
         */
        for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jRPh];
            kspStart = m_PhaseSpeciesStartIndex[iPh];
            ThermoPhase* th = & thermo(iPh);
            nSpecies = th->nSpecies();
            /*
             *  Loop over species
             */
            for (kSp = 0; kSp < nSpecies; kSp++) {
                int iKRSpecies = kstart + kSp;

                /*
                 * Calculate the total moles of the species in the current cell
                 */
                spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] = fracTotalVol * spMoles_final_[kspStart + kSp];
                /*
                 *  Calculate the concentration of the species in the cell
                 */
                concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] = spMoles_final_[kspStart + kSp] / particleVol;
                /*
                 * Calculate the mole fraction of the species in the cell
                 */
                spMf_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] = spMf_final_[kspStart + kSp];
                /*
                 * Calculate the total concentration of phases in this cell
                 */
                concTot_SPhase_Cell_final_[indexMidKRPhases + jRPh] += phaseMoles_final_[iPh] /  particleVol;



            }

            /*
             *  Set the state of the thermophase object
             */
            th->setState_TPX(temperature_, pressure_, &spMf_KRsolid_Cell_final_[indexMidKRSpecies + kstart]);
            /*
             *  Find the partial molar volumes of the species and the activity coefficients
             */
            th->getPartialMolarVolumes(& partialMolarVolKRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
            th->getActivityCoefficients(& actCoeff_Cell_final_[indexMidKRSpecies + kstart]);
            /*
             * increment the species index
             */
            kstart += nSpecies;

        }


    }



    /*
     *  Generate an initial value of the lattice position refence field.
     *  This will be identically equal to the radius initially.
     */
    MolarVolume_refLat_Ref_ =  1.0 / concTot_SPhase_Cell_final_[0];
    for (iCell = 0; iCell < numRCells_; iCell++) {
        molarVolume_refLat_Cell_final_[iCell] =  MolarVolume_refLat_Ref_;
        rRefPos_final_[iCell] = rnodePos_final_[iCell];
    }


}
//====================================================================================================================
//    The internal state of the electrode must be kept for the initial and final times of an integration step.
/*
 *  This function advances the initial state to the final state that was calculated
 *  in the last integration step.
 *
 * @param Tinitial   This is the New initial time. This time is compared against the "old"
 *                   final time, to see if there is any problem.
 */
void  Electrode_RadialRegion::resetStartingCondition(double Tinitial, bool doResetAlways)
{
    /*
    * If the initial time is input, then the code doesn't advance
    */
    double tbase = MAX(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-9 * tbase) && !doResetAlways) {
        return;
    }
    Electrode_Integrator::resetStartingCondition(Tinitial, doResetAlways);
}
//====================================================================================================================
//! Take the state (i.e., the final state) within the Electrode_Model and push it down
//! to the ThermoPhase objects and propogate it to all other aspects of the final state
/*!
 *  This routine does not calculate the global summed state. It does add contributions
 *  into spMoles_final_ and that's it.
 *
 *  We take the unknowns from the nonlinear problem and calculate all the values
 *
 *  All of these properties are defined for the _final_ state.
 *
 *  Thus, this is the main routine that reconciles all of the state information within the object.
 *  At the end of this routine, all aspects of the final state are consistent with each other.
 *
 *  Fundamental Variables:
 *       concKRSpecies_Cell_final_[]
 *
 *  Quantitities filled by this routine
 *     Distributed Quantities:
 *        cellBoundR_final_[iCell]
 *        spMoles_KRsolid_Cell_final_
 *        spMf_KRsolid_Cell_final_
 *        concTot_SPhase_Cell_final_
 *        partialMolarVolKRSpecies_Cell_final_
 *     Global Quantities:
 *         spMoles_final_   (additions)
 *
 *
 *  @param zeroGlobals   Zero the globals before adding to spMoles.
 *
 */
void Electrode_RadialRegion::updateStateDistrib(bool zeroGlobals)
{
    // Indexes
    int iCell, iPh, jRPh, kSp, kspStart;
    // Cubes of the cell boundaries, Right and Left, for the initial and final times
    double cbR3_final = 0.0;
    double cbL3_final = 0.0;

    /*
     *   We now have a vector of cells.
     *   Each of the cells must have their own conditions.
     *   We need to loop over the conditions and have their activities calculated
     */
    /*
     *  Calculate the cell boundaries in a pre - loop
     */
    for (iCell = 0; iCell < numRCells_; iCell++) {
        cellBoundR_final_[iCell]  = 0.5 * (rnodePos_final_[iCell] + rnodePos_final_[iCell+1]);
    }

    /*
     *  Zero out the distributed total species moles
     */
    if (zeroGlobals) {
        for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jRPh];
            kspStart = m_PhaseSpeciesStartIndex[iPh];
            ThermoPhase* th = & thermo(iPh);
            int nSpecies = th->nSpecies();
            for (kSp = 0; kSp < nSpecies; kSp++) {
                spMoles_final_[kspStart + kSp] = 0.0;
            }
        }
    }

    for (int iCell = 0; iCell < numRCells_; iCell++) {


        /*
         *  Calculate the cell size
         */
        cbL3_final  = cbR3_final;
        cbR3_final = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
        double volCell = 4./3. * (cbR3_final  - cbL3_final);

        double volTotalCell = volCell * particleNumberToFollow_;

        int indexMidKRSpecies =  iCell * numKRSpecies_;
        int indexMidKRPhases =  iCell * numSPhases_;
        int kstart = 0;
        for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jRPh];
            kspStart = m_PhaseSpeciesStartIndex[iPh];
            ThermoPhase* th = & thermo(iPh);
            int nSpecies = th->nSpecies();

            double moleSum = 0.0;
            concTot_SPhase_Cell_final_[indexMidKRPhases + jRPh] = 0.0;
            for (int kSp = 0; kSp < nSpecies; kSp++) {
                int iKRSpecies = kstart + kSp;
                concTot_SPhase_Cell_final_[indexMidKRPhases + jRPh] += concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];

                /*
                 * Find the mole numbers of species in the cell, spMolesKRSpecies_Cell_final_[indexTopKRSpecies + iKRSpecies]
                 *     from concKRsolid_Cell_final_;
                 */
                spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] = volTotalCell * concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];
                moleSum += spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies];
                spMoles_final_[kspStart + kSp] += spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies];
            }
            /*
             * Find the mole fractions
             *     from spMoles_KRsolid_Cell_final_;
             *     If the phase is not present, don't update the mole fractions. We will go with the mole fractions saved in spMf_KRsolid_Cell_final_
             */
            if (concTot_SPhase_Cell_final_[indexMidKRPhases + jRPh] > 1.0E-200) {
                for (int kSp = 0; kSp < nSpecies; kSp++) {
                    int iKRSpecies = kstart + kSp;
                    spMf_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] = spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] / moleSum;
                }
                th->setState_TPX(temperature_, pressure_,  &spMf_KRsolid_Cell_final_[indexMidKRSpecies + kstart]);
            }

            th->getPartialMolarVolumes(& partialMolarVolKRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
            th->getActivityCoefficients(& actCoeff_Cell_final_[indexMidKRSpecies + kstart]);
            kstart += nSpecies;
        }

        /*
         * Calculate the reference lattice molar volume for the current cell. Here we assume that the phase 0 concentration is
         * inverse of the molar volume of the lattice.
         *  The calculation will have to be generalized more if this is not the case
         */
        molarVolume_refLat_Cell_final_[iCell] =  1.0 / concTot_SPhase_Cell_final_[indexMidKRPhases + 0];


    }


}
//================================================================================================
/*
 * There is a small dependence on mf_external and mf_internal exhibited by this function
 */
void  Electrode_RadialRegion::extractInfoJustBorn(std::vector<int>& justBornMultiSpecies)
{

    updateState();
}

//====================================================================================================================
// Main routine to calculate the residual
/*
 *
 *  Format of the residual equations
 *                                                             Unknown
 * --------------------------------------------------------------------------------------------------------------
 *         Residual (Time)                                     deltaT
 *
 *         Loop over cells
 *
 *           Residual (Reference/lattice Position)            rRefPos_final_[iCell];
 *           Residual (Mesh Position)
 *           Residual (Concentration _ k=0)
 *             . . .
 *           Residual (Concentration _ k=Ns-1)
 *  --------------------------------------------------------------------------------------------------------------
 *
 */
int Electrode_RadialRegion::calcResid(double* const resid, const ResidEval_Type_Enum evalType)
{
    // Indexes
    int iCell, iPh, jPh;
    // Cubes of the cell boundary radii, Right and Left, for the initial and final times
    double cbRadius3_RHS_init = 0.0;
    double cbRadius3_LHS_init = 0.0;
    double cbRadius3_RHS_final = 0.0;
    double cbRadius3_LHS_final = 0.0;

    // Reference lattice squared, Right and Left, at final time
    double refLat2_RHS_final = 0.0;
    // Reference lattice cubed, Right and Left, at final time
    double refLat3_RHS_final = 0.0;
    double refLat3_LHS_final = 0.0;


    double deltaX, fluxTC, fluxMoles;
    double dcadx[10];

    // Diffusive flux
    double flux[10];
    // Velocity of the refrence radius on the right side of the cell
    double vrefR = 0.0;
    // Temporary pointers for accessing the phase total concentrations (kmol m-3), init and final
    double* concTotalVec_SPhase_init = 0;
    double* concTotalVec_SPhase_final = 0;
    // index of the rref equation
    int rindex;
    // index of the mesh motion equation
    int xindex;
    // index of the start of the concentration equations for the first phase at the current cell
    int cindex;
    // solution index at the start of the current phase at the current cell
    int cIndexPhStart;



    // Velocity of the right cell boundary during the time step;
    std::vector<doublereal> cellBoundRVeloc(numRCells_);


    // Velocity of the left cell boundary during the time step;
    std::vector<doublereal> cellBoundLVeloc(numRCells_);

    // Velocity of the lattice at the right cell boundary during the time step
    std::vector<doublereal> refLatVeloc_RHS(numRCells_);

    // Velocity of the lattice at the right cell boundary during the time step
    std::vector<doublereal> refLatVeloc_LHS(numRCells_);


    // Node velocity during the time step
    std::vector<doublereal> rnodeVeloc(numRCells_);

    // refLatticeInterpolatedField cubed
    std::vector<doublereal> refLattice3InterpolatedField(numRCells_);



    /*
     *    Residual equation for the time step -> Right now we don't have a model
     */
    resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;

    /*
     *  Calculate the cell boundaries in a pre - loop
     */
    for (iCell = 0; iCell < numRCells_; iCell++) {
        rnodeVeloc[iCell] = (rnodePos_final_[iCell] - rnodePos_init_[iCell]) / deltaTsubcycleCalc_;

        cellBoundR_init_[iCell]   = 0.5 * (rnodePos_init_[iCell] + rnodePos_init_[iCell+1]);
        cellBoundRVeloc[iCell] = (cellBoundR_final_[iCell] - cellBoundR_init_[iCell]) / deltaTsubcycleCalc_;
        // cellBoundLVeloc[iCell] = (cellBoundL_final_[iCell] - cellBoundL_init_[iCell] ) / deltaTsubcycleCalc_;
    }





    // -------------------------- Identify Total Lattice Production on Right Side ----------------------------------------------------

    iCell = numRCells_ - 1;

    /*
     *  Select the first phase, so that we can calculate the lattice production rates on the right side
     */
    jPh = 0;
    iPh = phaseIndeciseKRsolidPhases_[jPh];
    int iStart = getGlobalSpeciesIndex(iPh, 0);
    ThermoPhase* th = & thermo(iPh);
    int nSpecies = th->nSpecies();

    /*
     *   Find the index start for the current cell
     */
    int kCellStart = iCell *  numKRSpecies_;

    /*
     * Calculate the Solid volume creation rate for the first phase
     */
    double SolidVolCreationRatePhase0 = 0.0;
    for (int kSp = 0; kSp < nSpecies; kSp++) {
        SolidVolCreationRatePhase0 += partialMolarVolKRSpecies_Cell_final_[kCellStart + kSp] * DspMoles_RightSurf_final_[iStart + kSp];
    }

    /*
     *  Calcualte the molar production rate of the first Phase.
     */
    double latticeProductionRate = 0.0;
    for (int kSp = 0; kSp < nSpecies; kSp++) {
        latticeProductionRate += DspMoles_RightSurf_final_[iStart + kSp];
    }

    double molarLatticeVolume_final_ = 0.0;
    for (int kSp = 0; kSp < nSpecies; kSp++) {
        molarLatticeVolume_final_ += spMf_KRsolid_Cell_final_[kCellStart + kSp] * partialMolarVolKRSpecies_Cell_final_[kSp];
    }


    double cellBoundR_star2 = (cellBoundR_final_[iCell] * cellBoundR_final_[iCell]
                               + cellBoundR_final_[iCell] * cellBoundR_init_[iCell] + cellBoundR_init_[iCell] * cellBoundR_init_[iCell])/3.0;
    double areaR_star = 4.0 * Pi * cellBoundR_star2 * particleNumberToFollow_;
    /*
     * Calculation of the movement at the top of the
     *    C dXdt * Area - SolidCreationRate = 0
     *
     *  this is an expression of the conservation of total solid molar volume at the
     *  r_exterior.
     */
    xindex = 1 + numEqnsCell_ * iCell;
    resid[xindex] = rnodeVeloc[iCell]  * areaR_star - SolidVolCreationRatePhase0;




    double deltaRStar = 0.0;



    // --------------------------- First Main Loop Over Cells ----------------------------------------------------
    /*
     *  Calculate lattice and mesh movement
     *     -> Need to calculate convection
     */
    int currCell = 0;
    double r3L_curr = 0.0;
    double refLat3_LHS_curr = 0.0;

    for (int iCell = 0; iCell < numRCells_; iCell++) {
        /*
         *  Copy right side to left side
         */
        cbRadius3_LHS_final  = cbRadius3_RHS_final;
        cbRadius3_LHS_init   = cbRadius3_RHS_init;
        refLat3_LHS_final  = refLat3_RHS_final;

        /*
         *  Calculate indexes for accessing the residual
         */
        rindex = 1 + numEqnsCell_ * iCell;
        xindex = 2 + numEqnsCell_ * iCell;
        cindex = 3 + numEqnsCell_ * iCell;
        cIndexPhStart = cindex;

        /*
         *  Temporary variables for total concentrations of each of the phases
         */
        concTotalVec_SPhase_init  = &(concTot_SPhase_Cell_init_[iCell*numSPhases_]);
        concTotalVec_SPhase_final = &(concTot_SPhase_Cell_final_[iCell*numSPhases_]);

        /*
         * Calculate the cubes of the cell boundary radii on the rhs
         */
        cbRadius3_RHS_init  = cellBoundR_init_[iCell]  * cellBoundR_init_[iCell]  * cellBoundR_init_[iCell];
        cbRadius3_RHS_final = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];

        /*
         * Calculate powers of the reference velocity at the right boundary
         */
        double refLat_RHS_final = rRefPos_final_[iCell];


        refLat2_RHS_final = refLat_RHS_final * refLat_RHS_final;
        refLat3_RHS_final = refLat2_RHS_final * refLat_RHS_final;
        /*
         * Calculate the molar volume of the first phase, final and init
         */
        double vbar_final = 1.0 / concTotalVec_SPhase_final[0];
        //double vbar_init  = 1.0 / concTotalVec_SPhase_init[0];

        /*
         *  Calculate the area of the outer cell which conserves constant functions under mesh movement wrt the Reynolds transport theorum
         */
        cellBoundR_star2 = (cellBoundR_final_[iCell] * cellBoundR_final_[iCell]
                            + cellBoundR_final_[iCell] * cellBoundR_init_[iCell] + cellBoundR_init_[iCell] * cellBoundR_init_[iCell])/3.0;
        areaR_star = 4.0 * Pi * cellBoundR_star2;
        /*
         * Residual calculation - Value of the reference/lattice radius at the right cell boundary
         */
        double rhs = refLat3_LHS_final + MolarVolume_refLat_Ref_ / vbar_final * (cbRadius3_RHS_final - cbRadius3_LHS_final);
        resid[rindex] = refLat3_RHS_final - rhs;

        /*
         *  Calculate the time derivative of the molar volume
         */
        //double vbarDot = (vbar_final - vbar_init) / deltaTsubcycleCalc_;

        for (jPh = 0; jPh < numSPhases_; jPh++) {
            for (int kSp = 0; kSp < nSpecies; kSp++) {

            }
        }


        /*
         * Value to find:
         *       This is the value of the init RHS refLat
         */
        double refLatFind = rRefPos_init_[iCell];
        double refLatFind3 = refLatFind * refLatFind  * refLatFind;
        /*
         * Search the current cell
         */
        if (refLatFind <= rRefPos_final_[currCell]) {
            double r3 =  r3L_curr + vbar_final / MolarVolume_refLat_Ref_ * (refLatFind3 -  refLat3_LHS_curr);
            refLattice3InterpolatedField[iCell] = r3;
        } else {
            do {
                if (currCell < numRCells_) {
                    currCell++;
                }
            } while (refLatFind > rRefPos_final_[currCell] && currCell <= numRCells_);

            double r3 =  r3L_curr + vbar_final / MolarVolume_refLat_Ref_ * (refLatFind3 - refLat3_LHS_curr);
            refLattice3InterpolatedField[iCell] = r3;
        }



    }

    //
    // --------------------------- Second Main Loop Over Cells ----------------------------------------------------
    //
    /*
     *  Calculate speciec conservation
     */

    for (int iCell = 0; iCell < numRCells_; iCell++) {

        /*
         *  Copy left side to right side
         */
        cbRadius3_LHS_final  = cbRadius3_RHS_final;
        cbRadius3_LHS_init   = cbRadius3_RHS_init;
        refLat3_LHS_final  = refLat3_RHS_final;

        /*
         *  Calculate indexes for accessing the residual
         */
        rindex = 1 + numEqnsCell_ * iCell;
        xindex = 2 + numEqnsCell_ * iCell;
        cindex = 3 + numEqnsCell_ * iCell;
        cIndexPhStart = cindex;

        /*
         *  Temporary variables for total concentrations of each of the phases
         */
        concTotalVec_SPhase_init  = &(concTot_SPhase_Cell_init_[iCell*numSPhases_]);
        concTotalVec_SPhase_final = &(concTot_SPhase_Cell_final_[iCell*numSPhases_]);
        /*
         * Calculate the cubes of the cell boundary radii
         */
        cbRadius3_RHS_init  = cellBoundR_init_[iCell]  * cellBoundR_init_[iCell]  * cellBoundR_init_[iCell];
        cbRadius3_RHS_final = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
        /*
         * Calculate powers of the reference velocity at the right boundary
         */
        double r0R_final = rRefPos_final_[iCell];
        refLat2_RHS_final = r0R_final * r0R_final;
        /*
         * Calculate the molar volume of the first phase, final and init
         */
        // double vbar_final = 1.0 / concTotalVec_SPhase_final[0];
        //double vbar_init  = 1.0 / concTotalVec_SPhase_init[0];

        /*
         *  Calculate the area of the outer cell which conserves constant functions under mesh movement wrt the Reynolds transport theorum
         */
        cellBoundR_star2 = (cellBoundR_final_[iCell] * cellBoundR_final_[iCell]
                            + cellBoundR_final_[iCell] * cellBoundR_init_[iCell] + cellBoundR_init_[iCell] * cellBoundR_init_[iCell])/3.0;
        areaR_star = 4.0 * Pi * cellBoundR_star2;

        /*
         *  Calculate the time derivative of the molar volume
         */
        //double vbarDot = (vbar_final - vbar_init) / deltaTsubcycleCalc_;

        for (jPh = 0; jPh < numSPhases_; jPh++) {
            for (int kSp = 0; kSp < nSpecies; kSp++) {

            }
        }




        /*
         * Calculate the molar volume of the lattice at the beginning of the step.
         */
        kCellStart = iCell *  numKRSpecies_;
        double molarLatticeVolume_init_ = 0.0;
        for (int kSp = 0; kSp < nSpecies; kSp++) {
            molarLatticeVolume_init_ += spMf_KRsolid_Cell_init_[kCellStart + kSp] * partialMolarVolKRSpecies_Cell_init_[kCellStart + kSp];
        }

        /*
         *  Calculate the total lattice moles in the initial cell.
         */
        double TotalLatticeMoles_init = 4. / 3. * Pi * (cbRadius3_RHS_init -  cbRadius3_LHS_init) / molarLatticeVolume_init_;

        /*
         *  Calculate the molar lattice volume
         *       -> Assume that the first solid phase is the lattice phase
         */
        jPh = 0;
        molarLatticeVolume_final_ = 0.0;
        for (int kSp = 0; kSp < nSpecies; kSp++) {
            molarLatticeVolume_final_ += spMf_KRsolid_Cell_final_[kSp] * partialMolarVolKRSpecies_Cell_final_[kSp];
        }

        double cellVolume_star = TotalLatticeMoles_init * molarLatticeVolume_final_;

        double cbRadius3_RHS_star = 3. / 4. / Pi * (cellVolume_star + 4./3.* Pi * cbRadius3_LHS_final);

        double cbRadius_RHS_star = pow(cbRadius3_RHS_star, 0.33333333333333);

        rRefPos_final_[iCell+1] =  rnodePos_final_[iCell] + 2.0 * (cbRadius_RHS_star - rnodePos_final_[iCell]);

        rnodePos_final_[iCell+1] =  rRefPos_final_[iCell+1] - rnodePos_final_[iCell] + 0.5 * (fractionVolExpansion_Cell_[iCell+1] + fractionVolExpansion_Cell_[iCell]) * deltaRStar;


        /*
         *  Figure out the volume of lattice that has moved on the rhs of the cell volume
         */

        //bool gainedLatticeSites = false;
        double volGained_RHS;
        if (cbRadius3_RHS_final >= refLattice3InterpolatedField[iCell]) {
            //gainedLatticeSites = true;
            volGained_RHS = 4. * Pi / 3.0 * (cbRadius3_RHS_final - refLattice3InterpolatedField[iCell]);
        } else {
            volGained_RHS = 4. * Pi / 3.0 * (cbRadius3_RHS_final - refLattice3InterpolatedField[iCell]);
        }

        /*
         * Find the velocity of the reference radius at the right cell boundary
         */
        // vrefR = -1.0 / (3. * r0R2_final) * (3.0 *r0L2_final * vrefL - MolarVolume_Ref_ / (vbar_final * vbar_final) * vbarDot * (cbRadius3_RHS_final - cbRadius3_LHS_final));





        /*
         * Node position residual - spline equations, it all depends on the top node's equation formulation.
         * Everything else get's dragged along with it
         *
         *      resid(Rnode_IC) =  Rnode_iC - RLeft_final - frac_iC/frac_iCP1 * (Rnode_iCP1 - RLeft_final) = 0
         */
        resid[xindex] = rnodePos_final_[iCell] - radiusLeft_final_
                        - fractionVolExpansion_Cell_[iCell] * (rnodePos_final_[iCell+1] - radiusLeft_final_);




        int indexMidKRSpecies =  iCell    * numKRSpecies_;
        int indexTopKRSpecies = (iCell+1) * numKRSpecies_;
        int kstart = 0;




        /*
         *  Loop Over phases
         *
         *      new_moles =   4/3pi * (R_iCRHS **3 - R_icLHS **3) C_iPhase_final
         *      old_moles =   4/3pi * (R_iCRHS **3 - R_icLHS **3) C_iPhase_init
         *
         *  mesh movement with no material movement -> idea here is to ensure that constant functions can advect unchanged when the mesh moves.
         *
         *      Resid(phaseMoles_iC) = (new_moles - old_moles) /deltaT
         *                              - ((V_iCRHS -vrefR) *  C_iPhase_final * A_iCRHS)  +  ((V_iCLHS -vrefR) *  C_iPhase_final * A_iCLHS)
         */
        for (jPh = 0; jPh < numSPhases_; jPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jPh];
            ThermoPhase* th = & thermo(iPh);
            int nSpecies = th->nSpecies();
            double old_stuff = concTotalVec_SPhase_init[jPh]  * 4.0 / 3.0 * Pi * (cbRadius3_RHS_init -  cbRadius3_LHS_init) * particleNumberToFollow_;
            double new_stuff = concTotalVec_SPhase_final[jPh] * 4.0 / 3.0 * Pi * (cbRadius3_RHS_final - cbRadius3_LHS_final) * particleNumberToFollow_;

            /*
             *  Change in the total amount of moles of phase jPh
             */
            resid[cIndexPhStart] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;
            /*
             * Convective flux - mesh movement with NO material movement
             */
            double vtotalR = cellBoundRVeloc[iCell] -  vrefR;
            if (iCell < (numRCells_- 1)) {
                if (vtotalR >= 0.0) {
                    fluxTC = vtotalR * concTot_SPhase_Cell_final_[numSPhases_*(iCell+1) + jPh];
                } else {
                    fluxTC = vtotalR * concTot_SPhase_Cell_final_[numSPhases_*(iCell)   + jPh];
                }
                resid[cIndexPhStart]                -= fluxTC * areaR_star;
                resid[cIndexPhStart + numEqnsCell_] += fluxTC * areaR_star;
            }


            /*
             *  Residual for the species in the phase
             *   -  Skip the first species
             *
             *
             *      new_moles =   4/3pi * (R_iCRHS **3 - R_icLHS **3) C_iPhase_final
             *      old_moles =   4/3pi * (R_iCRHS **3 - R_icLHS **3) C_iPhase_init
             *
             *  mesh movement with no material movement -> idea here is to ensure that constant functions can advect unchanged when the mesh moves.
             *
             *      Resid(phaseMoles_iC) = (new_moles - old_moles) /deltaT
             *                              - ((V_iCRHS -vrefR) *  C_iPhase_final * A_iCRHS)  +  (( V_iCLHS -vrefR) *  C_iPhase_final * A_iCLHS)
             *                              - fluxDiffusion_iCRHS * A_iCRHS                 + fluxDiffusion_iCLHS * A_iCLHS
             */
            for (int kSp = 1; kSp < nSpecies; kSp++) {
                int iKRSpecies = kstart + kSp;
                /*
                 * Diffusive flux at the outer boundary for each species
                 */
                if (iCell < (numRCells_ - 1)) {
                    deltaX = rnodePos_final_final_[iCell+1] - rnodePos_final_final_[iCell];
                    double caR = concKRSpecies_Cell_final_[indexTopKRSpecies + iKRSpecies] * actCoeff_Cell_final_[indexTopKRSpecies + iKRSpecies];
                    double caL = concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] * actCoeff_Cell_final_[indexMidKRSpecies + iKRSpecies];
                    dcadx[kSp] = (caR - caL) / deltaX;
                    flux[kSp] = - Diff_Coeff_KRSolid_[iKRSpecies] * dcadx[kSp];

                    resid[cIndexPhStart + kSp]                += flux[kSp] * areaR_star;
                    resid[cIndexPhStart + kSp + numEqnsCell_] -= flux[kSp] * areaR_star;
                }

                /*
                 * Convective flux - mesh movement relative to the material movement
                 *          We don't do it for the last cell
                 *
                 *
                 */
                if (iCell < (numRCells_- 1)) {
                    if (volGained_RHS >= 0.0) {
                        fluxMoles = volGained_RHS / deltaTsubcycleCalc_ * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell+1) + iKRSpecies];
                        resid[cIndexPhStart + kSp]                -= fluxMoles;
                        resid[cIndexPhStart + kSp + numEqnsCell_] += fluxMoles;
                    } else {
                        fluxMoles = volGained_RHS / deltaTsubcycleCalc_ * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell)  + iKRSpecies];
                        resid[cIndexPhStart + kSp]                -= fluxMoles;
                        resid[cIndexPhStart + kSp + numEqnsCell_] += fluxMoles;
                    }
                }

                /*
                 *  Change in the total amount of moles of species jPh
                 */
                old_stuff = concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies] * 4.0 / 3.0 * Pi *
                            (cbRadius3_RHS_init  - cbRadius3_LHS_init) * particleNumberToFollow_;

                new_stuff = concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies] * 4.0 / 3.0 * Pi * (cbRadius3_RHS_final - cbRadius3_LHS_final)* particleNumberToFollow_;

                resid[cIndexPhStart + kSp] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;
            }

            // end of phase loop
            cIndexPhStart += nSpecies;
            kstart += nSpecies;
        }

        if (iCell == numRCells_ - 1) {
            /*
             *  Do the exterior
             */
            iCell =  numRCells_ - 1;
            double SolidVolCreationRate = 0.0;

            for (jPh = 0; jPh < numSPhases_; jPh++) {
                iPh = phaseIndeciseKRsolidPhases_[jPh];
                int iStart = getGlobalSpeciesIndex(iPh, 0);
                ThermoPhase* th = & thermo(iPh);
                int nSpecies = th->nSpecies();

                /*
                 *  Add in the molar production term for the phase due to the reaction on the right surface
                 */
                resid[cIndexPhStart] -= DphMoles_RightSurf_final_[iPh];

                /*
                 *  Add in the molar production rate for each species due to the reaction on the right surface
                 *  Also calculate the volume creation rate.
                 */
                SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart] * DspMoles_RightSurf_final_[iStart];
                for (int kSp = 1; kSp < nSpecies; kSp++) {
                    resid[cIndexPhStart + kSp] -= DspMoles_RightSurf_final_[iStart + kSp];
                    SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart + kSp] * DspMoles_RightSurf_final_[iStart + kSp];
                }
            }
            /*
             * Calculation of the movement at the top of the
             *    C dXdt * Area - SolidCreationRate = 0
             *
             *  this is an expression of the conservation of total solid molar volume at the
             *  r_exterior.
             */
            xindex = 1 + numEqnsCell_ * iCell;
            resid[xindex] = rnodeVeloc[iCell]  * areaR_star - SolidVolCreationRate;
        }
    }
    return 1;
}

//====================================================================================================================
//! Set the internal initial intermediate and initial global state from the internal final state
/*!
 *  (virtual function)
 *
 *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 * @param setInitInit   Boolean indicating whether you should set the init_init state as well
 *                      (default is false)
 */
void  Electrode_RadialRegion::setInitStateFromFinal(bool setInitInit)
{
    int i;
    int kspCell =  numKRSpecies_ *  numRCells_;
    int kphCell = numSPhases_ * numRCells_;

    radiusLeft_init_ = radiusLeft_final_;

    for (i = 0; i <  numRCells_; i++) {
        cellBoundR_init_[i] =  cellBoundR_final_[i];
        rnodePos_init_[i] = rnodePos_final_[i];
        rRefPos_init_[i] = rRefPos_final_[i];
        molarVolume_refLat_Cell_init_[i] = molarVolume_refLat_Cell_final_[i];
    }

    for (i = 0; i < kspCell; i++) {
        spMoles_KRsolid_Cell_init_[i] = spMoles_KRsolid_Cell_final_[i];
        spMf_KRsolid_Cell_init_[i] = spMf_KRsolid_Cell_final_[i];
        concKRSpecies_Cell_init_[i] = concKRSpecies_Cell_final_[i];
        partialMolarVolKRSpecies_Cell_init_[i] = partialMolarVolKRSpecies_Cell_final_[i];
    }
    for (i = 0; i < kphCell; i++) {
        concTot_SPhase_Cell_init_[i] = concTot_SPhase_Cell_final_[i];
    }

    if (setInitInit) {

        radiusLeft_init_init_ = radiusLeft_final_;

        for (i = 0; i <  numRCells_; i++) {
            cellBoundR_init_init_[i] =  cellBoundR_final_[i];
            rnodePos_init_init_[i] = rnodePos_final_[i];
            rRefPos_init_init_[i] = rRefPos_final_[i];
            molarVolume_refLat_Cell_init_init_[i] = molarVolume_refLat_Cell_final_[i];
        }
        for (i = 0; i < kspCell; i++) {
            spMoles_KRsolid_Cell_init_init_[i] = spMoles_KRsolid_Cell_final_[i];
            spMf_KRsolid_Cell_init_init_[i] = spMf_KRsolid_Cell_final_[i];
            concKRSpecies_Cell_init_init_[i] = concKRSpecies_Cell_final_[i];
            partialMolarVolKRSpecies_Cell_init_init_[i] = partialMolarVolKRSpecies_Cell_final_[i];
        }
        for (i = 0; i < kphCell; i++) {
            concTot_SPhase_Cell_init_init_[i] = concTot_SPhase_Cell_final_[i];
        }
    }
}
//====================================================================================================================
// Set the internal final intermediate state from the internal init state
/*
 *  (virtual function from Electrode)
 *
 *  Set the final state from the init state. This is commonly called during a failed time step
 */
void  Electrode_RadialRegion::setFinalStateFromInit()
{
    int i;
    int kspCell =  numKRSpecies_ *  numRCells_;
    int kphCell = numSPhases_ * numRCells_;

    radiusLeft_final_ = radiusLeft_init_;

    for (i = 0; i <  numRCells_; i++) {
        cellBoundR_final_[i] =  cellBoundR_init_[i];
        rnodePos_final_[i] = rnodePos_init_[i];
        rRefPos_final_[i] = rRefPos_init_[i];
        molarVolume_refLat_Cell_final_[i] = molarVolume_refLat_Cell_init_[i];
    }
    for (i = 0; i < kspCell; i++) {
        spMoles_KRsolid_Cell_final_[i] = spMoles_KRsolid_Cell_init_[i];
        spMf_KRsolid_Cell_final_[i] = spMf_KRsolid_Cell_init_[i];
        concKRSpecies_Cell_final_[i] = concKRSpecies_Cell_init_[i];
        partialMolarVolKRSpecies_Cell_final_[i] = partialMolarVolKRSpecies_Cell_init_[i];
    }
    for (i = 0; i < kphCell; i++) {
        concTot_SPhase_Cell_final_[i] = concTot_SPhase_Cell_init_[i];
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
void  Electrode_RadialRegion::setInitStateFromInitInit(bool setFinal)
{
    int i;
    int kspCell =  numKRSpecies_ *  numRCells_;
    int kphCell = numSPhases_ * numRCells_;

    radiusLeft_init_ = radiusLeft_init_init_;

    for (i = 0; i <  numRCells_; i++) {
        cellBoundR_init_[i] =  cellBoundR_init_init_[i];
        rnodePos_init_[i] = rnodePos_init_init_[i];
        rRefPos_init_[i] = rRefPos_init_init_[i];
        molarVolume_refLat_Cell_init_[i] = molarVolume_refLat_Cell_init_init_[i];
    }

    for (i = 0; i < kspCell; i++) {
        spMoles_KRsolid_Cell_init_[i] = spMoles_KRsolid_Cell_init_init_[i];
        spMf_KRsolid_Cell_init_[i] = spMf_KRsolid_Cell_init_init_[i];
        concKRSpecies_Cell_init_[i] = concKRSpecies_Cell_init_init_[i];
        partialMolarVolKRSpecies_Cell_init_[i] = partialMolarVolKRSpecies_Cell_init_init_[i];
    }
    for (i = 0; i < kphCell; i++) {
        concTot_SPhase_Cell_init_[i] = concTot_SPhase_Cell_init_init_[i];
    }


    if (setFinal) {

        radiusLeft_final_ = radiusLeft_init_init_;

        for (i = 0; i <  numRCells_; i++) {
            cellBoundR_final_[i] =  cellBoundR_init_init_[i];
            rnodePos_final_[i] = rnodePos_init_init_[i];
            rRefPos_final_[i] = rRefPos_init_init_[i];
        }
        for (i = 0; i < kspCell; i++) {
            spMoles_KRsolid_Cell_final_[i] = spMoles_KRsolid_Cell_init_init_[i];
            spMf_KRsolid_Cell_final_[i] = spMf_KRsolid_Cell_init_init_[i];
            concKRSpecies_Cell_final_[i] = concKRSpecies_Cell_init_init_[i];
            partialMolarVolKRSpecies_Cell_final_[i] = partialMolarVolKRSpecies_Cell_init_init_[i];
        }
        for (i = 0; i < kphCell; i++) {
            concTot_SPhase_Cell_final_[i] = concTot_SPhase_Cell_init_init_[i];
        }
    }

}
//====================================================================================================================
// Set the internal initial intermediate and initial global state from the internal final state
/*
 *  (virtual function from Electrode)
 *
 *  Set the intial state from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 * @param setInitInit   Boolean indicating whether you should set the init_init state as well. When
 *                      we do this we set the final state as well.
 */
void  Electrode_RadialRegion::setInitStateFromFinalFinal(bool setInitInit)
{
    int i;
    int kspCell =  numKRSpecies_ *  numRCells_;
    int kphCell = numSPhases_ * numRCells_;

    radiusLeft_init_ = radiusLeft_final_final_;

    for (i = 0; i <  numRCells_; i++) {
        cellBoundR_init_[i] =  cellBoundR_final_final_[i];
        rnodePos_init_[i] = rnodePos_final_final_[i];
        rRefPos_init_[i] = rRefPos_final_final_[i];
        molarVolume_refLat_Cell_init_[i] = molarVolume_refLat_Cell_final_final_[i];
    }

    for (i = 0; i < kspCell; i++) {
        spMoles_KRsolid_Cell_init_[i] = spMoles_KRsolid_Cell_final_final_[i];
        spMf_KRsolid_Cell_init_[i] = spMf_KRsolid_Cell_final_final_[i];
        concKRSpecies_Cell_init_[i] = concKRSpecies_Cell_final_final_[i];
        partialMolarVolKRSpecies_Cell_init_[i] = partialMolarVolKRSpecies_Cell_final_final_[i];
    }
    for (i = 0; i < kphCell; i++) {
        concTot_SPhase_Cell_init_[i] = concTot_SPhase_Cell_final_final_[i];
    }


    if (setInitInit) {

        radiusLeft_init_init_ = radiusLeft_final_final_;

        for (i = 0; i <  numRCells_; i++) {
            cellBoundR_init_init_[i] =  cellBoundR_final_final_[i];
            rnodePos_init_init_[i] = rnodePos_final_final_[i];
            rRefPos_init_init_[i] = rRefPos_final_final_[i];
            molarVolume_refLat_Cell_init_init_[i] = molarVolume_refLat_Cell_final_final_[i];
        }
        for (i = 0; i < kspCell; i++) {
            spMoles_KRsolid_Cell_init_init_[i] = spMoles_KRsolid_Cell_final_final_[i];
            spMf_KRsolid_Cell_init_init_[i] = spMf_KRsolid_Cell_final_final_[i];
            concKRSpecies_Cell_init_init_[i] = concKRSpecies_Cell_final_final_[i];
            partialMolarVolKRSpecies_Cell_init_init_[i] = partialMolarVolKRSpecies_Cell_final_final_[i];
        }
        for (i = 0; i < kphCell; i++) {
            concTot_SPhase_Cell_init_init_[i] = concTot_SPhase_Cell_final_final_[i];
        }
    }

    Electrode_RadialRegion::setFinalStateFromInit();
}
//====================================================================================================================
// Set the internal final global state from the internal final intermediate state
/*
 *  (virtual function from Electrode)
 *
 *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
 */
void  Electrode_RadialRegion::setFinalFinalStateFromFinal()
{
    int i;
    int kspCell =  numKRSpecies_ *  numRCells_;
    int kphCell = numSPhases_ * numRCells_;

    radiusLeft_final_final_ = radiusLeft_final_;

    for (i = 0; i <  numRCells_; i++) {
        cellBoundR_final_final_[i] =  cellBoundR_final_[i];
        rnodePos_final_final_[i] = rnodePos_final_[i];
        rRefPos_final_final_[i] = rRefPos_final_[i];
        molarVolume_refLat_Cell_final_final_[i] = molarVolume_refLat_Cell_final_[i];
    }

    for (i = 0; i < kspCell; i++) {
        spMoles_KRsolid_Cell_final_final_[i] = spMoles_KRsolid_Cell_final_[i];
        spMf_KRsolid_Cell_final_final_[i] = spMf_KRsolid_Cell_final_[i];
        concKRSpecies_Cell_final_final_[i] = concKRSpecies_Cell_final_[i];
        partialMolarVolKRSpecies_Cell_final_final_[i] = partialMolarVolKRSpecies_Cell_final_[i];
    }
    for (i = 0; i < kphCell; i++) {
        concTot_SPhase_Cell_final_final_[i] = concTot_SPhase_Cell_final_[i];
    }
}
//====================================================================================================================
void Electrode_RadialRegion::printElectrode(int pSrc, bool subTimeStep)
{
    int iph;
    double* netROP = new double[m_NumTotSpecies];
    double egv = TotalVol();
    printf("   ===============================================================\n");
    if (subTimeStep) {
        printf("      Electrode at intermediate-step time final = %g\n", tfinal_);
        printf("                   intermediate-step time final_final  = %g\n", tinit_);
    } else {
        printf("      Electrode at time final = %g\n", t_final_final_);
        printf("                   time init  = %g\n", t_init_init_);
    }
    printf("   ===============================================================\n");
    printf("          Number of external surfaces = %d\n", numExternalInterfacialSurfaces_);
    printf("          Solid Volume = %10.3E\n", ElectrodeSolidVolume_);
    printf("          Total Volume = %10.3E\n", egv);
    printf("          Temperature = %g\n", temperature_);
    printf("          Pressure = %g\n", pressure_);
    double capacd = capacityDischarged();
    printf("          Capacity Discharged = %g coulombs = %g Ah\n", capacd, capacd / 3600.);
    printf("\n");
    printf("          followElectrolyteMoles = %d\n", followElectrolyteMoles_);
    printf("          ElectrolytePseudoMoles = %g\n",  electrolytePseudoMoles_);

    for (iph = 0; iph < m_NumTotPhases; iph++) {
        printElectrodePhase(iph, pSrc);
        printf("     ===============================================================\n");
    }
    delete [] netROP;
}
//===================================================================================================================

void Electrode_RadialRegion::printElectrodePhase(int iph, int pSrc, bool subTimeStep)
{
    int isph = -1;
    double* netROP = new double[m_NumTotSpecies];
    ThermoPhase& tp = thermo(iph);
    string pname = tp.id();
    int istart = m_PhaseSpeciesStartIndex[iph];
    int nsp = tp.nSpecies();
    printf("     ===============================================================\n");
    printf("          Phase %d %s \n", iph,pname.c_str());
    printf("                Total moles = %g\n", phaseMoles_final_[iph]);
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
        printf("                  Voltage = %g\n", tp.electricPotential());
    }
    if (iph >= NumVolPhases_) {
        isph = iph - NumVolPhases_;
        printf("                surface area (final) = %g\n",  surfaceAreaRS_final_[isph]);
        printf("                surface area (init)  = %g\n",  surfaceAreaRS_init_[isph]);
        int ddd =  isExternalSurface_[isph];
        printf("                IsExternalSurface = %d\n", ddd);
        double oc = openCircuitVoltage(isph);
        if (oc != 0.0) {
            printf("                 Open Circuit Voltage = %g\n", oc);
        }
    }
    printf("\n");
    printf("                Name               MoleFrac_final  kMoles_final kMoles_init SrcTermLastStep(kMoles)\n");
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
                       spMoles_final_[istart + k],   spMoles_init_[istart + k]);
            } else {
                printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                       spMoles_final_[istart + k],   spMoles_init_init_[istart + k]);
            }
        }
    }
    if (iph >= NumVolPhases_) {
        const vector<double>& rsSpeciesProductionRates = RSD_List_[isph]->calcNetProductionRates();
        RSD_List_[isph]->getNetRatesOfProgress(netROP);

        doublereal* spNetProdPerArea = (doublereal*) spNetProdPerArea_List_.ptrColumn(isph);
        std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
        int nphRS = RSD_List_[isph]->nPhases();
        int kIndexKin = 0;
        for (int kph = 0; kph < nphRS; kph++) {
            int jph = RSD_List_[isph]->kinOrder[kph];
            int istart = m_PhaseSpeciesStartIndex[jph];
            int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
            for (int k = 0; k < nsp; k++) {
                spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                kIndexKin++;
            }
        }
        printf("\n");
        printf("                           spName                  Source (kmol/m2/s) \n");
        for (int k = 0; k <  m_NumTotSpecies; k++) {
            string ss = speciesName(k);
            printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
        }
    }
    delete [] netROP;

}
//======================================================================================================================
// Generate a grid of radial points
/*
 *  @param numPoints   Number of grid points
 *  @param radiusInner  Inner radius
 *  @param radiusOuter  Outer Radius
 *  @param geomFactor   Geometric variation of the radii
 *                          Values greater than one increase density at the outer
 *                          Values less than one increase density at the inner
 *  @param radialMesh   Returned mesh
 */
void Electrode_RadialRegion::radialGridGenerate(int numPoints, doublereal radiusInner, doublereal radiusOuter,
        doublereal geomFactor,
        std::vector<double>& radialMesh) const
{
    if (numPoints < 2) {
        throw CanteraError("radialgridGenerate() ERROR",
                           "number of points must be 2 or greater");
    }
    if (radiusInner < 0.0) {
        throw CanteraError("radialgridGenerate() ERROR",
                           "raidusInner must be greater than zero");
    }
    if (radiusOuter < 0.0) {
        throw CanteraError("radialgridGenerate() ERROR",
                           "raidusOuter must be greater than zero");
    }
    if (radiusOuter <= radiusInner) {
        throw CanteraError("radialgridGenerate() ERROR",
                           "raidusOuter must be greater than radiusInner");
    }
    if (geomFactor <= 0.0) {
        throw CanteraError("radialgridGenerate() ERROR",
                           "geomFactor  must be greater than zero");
    }
    /*
     *  General algorithm is to calcuate the radius' based on equalizing the volumes of the cells.
     *  Then we use geomFactor to make a skew of the mesh based on a geometric factor
     */

    int iCell;
    int numCells = numPoints - 1;
    double vol = (4. * Pi / 3. * radiusOuter * radiusOuter * radiusOuter -
                  4. * Pi / 3. * radiusInner * radiusInner * radiusInner);
    double avgCellVol = vol / numCells;


    std::vector<doublereal> cellVol(numCells, avgCellVol);

    int jCell = numCells / 2;


    for (iCell = jCell+1; iCell < numCells; iCell++) {
        cellVol[iCell] = cellVol[iCell-1] * geomFactor;
    }
    for (iCell = jCell-1; iCell >= 0; iCell++) {
        cellVol[iCell] = cellVol[iCell+1] / geomFactor;
    }
    double totalVol = 0.0;
    for (iCell = 0; iCell < numCells; iCell ++) {
        totalVol += cellVol[iCell];
    }

    double fac = vol / totalVol;
    for (iCell = 0; iCell < numCells; iCell ++) {
        cellVol[iCell] *= fac;
    }

    radialMesh.resize(numPoints);
    radialMesh[0] = radiusInner;
    double vInner = 4. * Pi / 3. *  radialMesh[0] *  radialMesh[0] *  radialMesh[0];
    for (iCell = 0; iCell < numCells; iCell ++) {
        double vOuter = vInner + cellVol[iCell];
        double cubed = (vOuter) * 3. / (4.0 * Pi);
        double rad = pow(cubed, 0.3333333333333333333);
        radialMesh[iCell+1] = rad;
        vInner = vOuter;
    }
    radialMesh[numPoints - 1] = radiusOuter;

}

} // End of namespace Cantera
//======================================================================================================================
