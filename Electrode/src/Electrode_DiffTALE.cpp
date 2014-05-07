/*
 * $Id: Electrode_DiffTALE.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include "tok_input_util.h"

#include "Electrode_DiffTALE.h"
#include "cantera/integrators.h"
#include "Electrode_RadialDiffRegions.h"
#include "EState_RadialDistrib.h"

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


static const double ONE_THIRD = 1.0 / 3.0;
namespace Cantera
{
//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_DiffTALE::Electrode_DiffTALE() :
    Electrode_Integrator(),

    electrodeType_(ELECTRODETYPE_ANODE),
    numKRSpecies_(1),
    numRCells_(5),
    numSPhases_(1),
    numEqnsCell_(0),
    thermoSPhase_List_(0),
    spMoles_KRsolid_Cell_final_(0),
    spMoles_KRsolid_Cell_init_(0),
    spMoles_KRsolid_Cell_final_final_(0),
    spMoles_KRsolid_Cell_init_init_(0),
    phaseMoles_KRsolid_Cell_final_(0),

    KRsolid_speciesList_(0),
    KRsolid_speciesNames_(0),
    phaseIndeciseKRsolidPhases_(0),
    distribPhIndexKRsolidPhases_(0),
    numSpeciesInKRSolidPhases_(0),
    KRsolid_phaseNames_(0),
    phaseIndeciseNonKRsolidPhases_(0),
    numNonSPhases_(0),
    concTot_SPhase_Cell_final_final_(0),
    concTot_SPhase_Cell_final_(0),
    concTot_SPhase_Cell_init_(0),
    concTot_SPhase_Cell_init_init_(0),
    molarDensity_SPhase_Cell_final_(0),
    concKRSpecies_Cell_init_(0),
    concKRSpecies_Cell_final_(0),
    concKRSpecies_Cell_init_init_(0),
    concKRSpecies_Cell_final_final_(0),

    spMf_KRSpecies_Cell_final_(0),
    spMf_KRSpecies_Cell_init_(0),
    MolarVolume_Ref_(0),
    rnodePos_final_final_(0),
    rnodePos_final_(0),
    rnodePos_init_(0),
    rnodePos_init_init_(0),
    rLatticeCBR_final_final_(0),
    rLatticeCBR_final_(0),
    rLatticeCBR_init_(0),
    rLatticeCBR_init_init_(0),
    rLatticeCBR_ref_(0),
    vLatticeCBR_cell_(0),
    cellBoundR_final_(0),
    cellBoundR_init_(0),
    cellBoundL_final_(0),
    cellBoundL_init_(0),

    volPP_Cell_final_(0),
    fracVolNodePos_(0),
    partialMolarVolKRSpecies_Cell_final_(0),
    partialMolarCpKRSpecies_Cell_final_(0),
    partialMolarEnthKRSpecies_Cell_final_(0),
    partialMolarEnthKRSpecies_Cell_init_(0),
    ROP_(0),
    DspMoles_final_(0),
    m_rbot0_(0.0),
    Diff_Coeff_KRSolid_(0),
    DphMolesSrc_final_(0),


    surfIndexExteriorSurface_(0),
    NTflux_final_(0.0),
    DiffCoeff_(1.0E-12),
    DiffCoeff_default_(1.0E-12),
    actCoeff_Cell_final_(0),
    actCoeff_Cell_init_(0),
    phaseID_TimeDeathMin_(-1),
    cellID_TimeDeathMin_(-1),
    onRegionBoundary_init_(-1),
    onRegionBoundary_final_(-1),
    onRegionBoundary_init_init_(-1),
    onRegionBoundary_final_final_(-1),
    atolBaseResid_(1.0E-12),
    goNowhere_(0),
    molarDensity_(0.0),
    formulationType_(0),
    formulationTypeTotalConc_(1)
{


}

//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_DiffTALE::Electrode_DiffTALE(const Electrode_DiffTALE& right) :
    Electrode_Integrator(),

    electrodeType_(ELECTRODETYPE_ANODE),
    numKRSpecies_(1),
    numRCells_(5),
    numSPhases_(1),
    numEqnsCell_(0),
    thermoSPhase_List_(0),
    spMoles_KRsolid_Cell_final_(0),
    spMoles_KRsolid_Cell_init_(0),
    spMoles_KRsolid_Cell_final_final_(0),
    spMoles_KRsolid_Cell_init_init_(0),
    phaseMoles_KRsolid_Cell_final_(0),

    KRsolid_speciesList_(0),
    KRsolid_speciesNames_(0),
    phaseIndeciseKRsolidPhases_(0),
    distribPhIndexKRsolidPhases_(0),
    numSpeciesInKRSolidPhases_(0),
    KRsolid_phaseNames_(0),
    phaseIndeciseNonKRsolidPhases_(0),
    numNonSPhases_(0),
    concTot_SPhase_Cell_final_final_(0),
    concTot_SPhase_Cell_final_(0),
    concTot_SPhase_Cell_init_(0),
    concTot_SPhase_Cell_init_init_(0),
    molarDensity_SPhase_Cell_final_(0),

    concKRSpecies_Cell_init_(0),
    concKRSpecies_Cell_final_(0),
    concKRSpecies_Cell_init_init_(0),
    concKRSpecies_Cell_final_final_(0),
  
    spMf_KRSpecies_Cell_final_(0),
    spMf_KRSpecies_Cell_init_(0),
    MolarVolume_Ref_(0),
    rnodePos_final_final_(0),
    rnodePos_final_(0),
    rnodePos_init_(0),
    rnodePos_init_init_(0),
    rLatticeCBR_final_final_(0),
    rLatticeCBR_final_(0),
    rLatticeCBR_init_(0),
    rLatticeCBR_init_init_(0),
    rLatticeCBR_ref_(0),
    vLatticeCBR_cell_(0),
    cellBoundR_final_(0),
    cellBoundR_init_(0),
    cellBoundL_final_(0),
    cellBoundL_init_(0),

    volPP_Cell_final_(0),
    fracVolNodePos_(0),
    partialMolarVolKRSpecies_Cell_final_(0),
    partialMolarCpKRSpecies_Cell_final_(0),
    partialMolarEnthKRSpecies_Cell_final_(0),
    partialMolarEnthKRSpecies_Cell_init_(0),
    ROP_(0),
    DspMoles_final_(0),
    m_rbot0_(0.0),
    Diff_Coeff_KRSolid_(0),
    DphMolesSrc_final_(0),

    surfIndexExteriorSurface_(0),
    NTflux_final_(0.0),
    DiffCoeff_(1.0E-12),
    DiffCoeff_default_(1.0E-12),
    actCoeff_Cell_final_(0),
    actCoeff_Cell_init_(0),
    phaseID_TimeDeathMin_(-1),
    cellID_TimeDeathMin_(-1),
    onRegionBoundary_init_(-1),
    onRegionBoundary_final_(-1),
    onRegionBoundary_init_init_(-1),
    onRegionBoundary_final_final_(-1),
    atolBaseResid_(1.0E-12),
    goNowhere_(0),
    molarDensity_(0.0),
    formulationType_(0),
    formulationTypeTotalConc_(1) 
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
Electrode_DiffTALE&
Electrode_DiffTALE::operator=(const Electrode_DiffTALE& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode_Integrator::operator=(right);

    electrodeType_                      = right.electrodeType_;
    numKRSpecies_                       = right.numKRSpecies_;
    numRCells_                          = right.numRCells_;
    numSPhases_                         = right.numSPhases_;
    numEqnsCell_                        = right.numEqnsCell_;
    thermoSPhase_List_                  = right.thermoSPhase_List_;
    spMoles_KRsolid_Cell_final_         = right.spMoles_KRsolid_Cell_final_;
    spMoles_KRsolid_Cell_init_          = right.spMoles_KRsolid_Cell_init_;
    spMoles_KRsolid_Cell_final_final_   = right.spMoles_KRsolid_Cell_final_final_;
    spMoles_KRsolid_Cell_init_init_     = right.spMoles_KRsolid_Cell_init_init_;
    phaseMoles_KRsolid_Cell_final_      = right.phaseMoles_KRsolid_Cell_final_;
    KRsolid_speciesList_                = right.KRsolid_speciesList_;
    KRsolid_speciesNames_               = right.KRsolid_speciesNames_;
    phaseIndeciseKRsolidPhases_         = right.phaseIndeciseKRsolidPhases_;
    distribPhIndexKRsolidPhases_        = right.distribPhIndexKRsolidPhases_;
    numSpeciesInKRSolidPhases_          = right.numSpeciesInKRSolidPhases_;
    KRsolid_phaseNames_                 = right.KRsolid_phaseNames_;
    concTot_SPhase_Cell_final_final_    = right.concTot_SPhase_Cell_final_final_;
    concTot_SPhase_Cell_final_          = right.concTot_SPhase_Cell_final_;
    concTot_SPhase_Cell_init_           = right.concTot_SPhase_Cell_init_;  
    concTot_SPhase_Cell_init_init_      = right.concTot_SPhase_Cell_init_init_;
    molarDensity_SPhase_Cell_final_     = right.molarDensity_SPhase_Cell_final_;
    concKRSpecies_Cell_init_            = right.concKRSpecies_Cell_init_;
    concKRSpecies_Cell_final_           = right.concKRSpecies_Cell_final_;
    concKRSpecies_Cell_init_init_       = right.concKRSpecies_Cell_init_init_;
    concKRSpecies_Cell_final_final_     = right.concKRSpecies_Cell_final_final_;
    spMf_KRSpecies_Cell_final_          = right.spMf_KRSpecies_Cell_final_;
    spMf_KRSpecies_Cell_init_           = right.spMf_KRSpecies_Cell_init_;
    MolarVolume_Ref_                    = right.MolarVolume_Ref_;

    rnodePos_final_final_               = right.rnodePos_final_final_;
    rnodePos_final_                     = right.rnodePos_final_;
    rnodePos_init_                      = right.rnodePos_init_;
    rnodePos_init_init_                 = right.rnodePos_init_init_;
    rLatticeCBR_final_final_            = right.rLatticeCBR_final_final_;
    rLatticeCBR_final_                  = right.rLatticeCBR_final_;
    rLatticeCBR_init_                   = right.rLatticeCBR_init_;
    rLatticeCBR_init_init_              = right.rLatticeCBR_init_init_;
    rLatticeCBR_ref_                    = right.rLatticeCBR_ref_;
    vLatticeCBR_cell_                   = right.vLatticeCBR_cell_;
    cellBoundR_final_                   = right.cellBoundR_final_;
    cellBoundR_init_                    = right.cellBoundR_init_;
    cellBoundL_final_                   = right.cellBoundL_final_;
    cellBoundL_init_                    = right.cellBoundL_init_;
    volPP_Cell_final_                   = right.volPP_Cell_final_;

 
    fracVolNodePos_                     = right.fracVolNodePos_;
    partialMolarVolKRSpecies_Cell_final_= right.partialMolarVolKRSpecies_Cell_final_;
    partialMolarCpKRSpecies_Cell_final_ = right.partialMolarCpKRSpecies_Cell_final_;
    partialMolarEnthKRSpecies_Cell_final_ = right.partialMolarEnthKRSpecies_Cell_final_;
    partialMolarEnthKRSpecies_Cell_init_ = right.partialMolarEnthKRSpecies_Cell_init_;

    ROP_                                = right.ROP_;
    DspMoles_final_                     = right.DspMoles_final_;
    m_rbot0_                            = right.m_rbot0_;
    Diff_Coeff_KRSolid_                 = right.Diff_Coeff_KRSolid_;
    DphMolesSrc_final_                  = right.DphMolesSrc_final_;
    surfIndexExteriorSurface_           = right.surfIndexExteriorSurface_;

    NTflux_final_                       = right.NTflux_final_;
    DiffCoeff_                          = right.DiffCoeff_;
    DiffCoeff_default_                  = right.DiffCoeff_default_;
    actCoeff_Cell_final_                = right.actCoeff_Cell_final_;
    actCoeff_Cell_init_                 = right.actCoeff_Cell_init_;
    phaseID_TimeDeathMin_               = right.phaseID_TimeDeathMin_;
    cellID_TimeDeathMin_                = right.cellID_TimeDeathMin_;
    onRegionBoundary_init_              = right.onRegionBoundary_init_;
    onRegionBoundary_final_             = right.onRegionBoundary_final_;
    onRegionBoundary_init_init_         = right.onRegionBoundary_init_init_;
    onRegionBoundary_final_final_       = right.onRegionBoundary_final_final_;
    atolBaseResid_                      = right.atolBaseResid_;
    goNowhere_                          = right.goNowhere_;
    molarDensity_                       = right.molarDensity_;
    formulationType_                    = right.formulationType_;
    formulationTypeTotalConc_           = right.formulationTypeTotalConc_;

    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
//   destructor
Electrode_DiffTALE::~Electrode_DiffTALE()
{
}
//======================================================================================================================
//    Return the type of electrode
/*
 *  Returns the enum type of the electrode. This is used in the factory routine.
 *
 *  @return Returns an enum type, called   Electrode_Types_Enum
 */
Electrode_Types_Enum Electrode_DiffTALE::electrodeType() const
{
    return DIFF_TALE_ET;
}
//====================================================================================================================
int Electrode_DiffTALE::electrode_input_child(ELECTRODE_KEY_INPUT** ei_ptr)
{
    /*
     *  Malloc an expanded child input
     */
    ELECTRODE_RadialRegion_KEY_INPUT * ei_radial = new ELECTRODE_RadialRegion_KEY_INPUT();
    /*
     *  Find the command file
     */
    ELECTRODE_KEY_INPUT* ei = *(ei_ptr);
    string commandFile = ei->commandFile_;
    BlockEntry* cf = ei->lastBlockEntryPtr_;
    ei_radial->printLvl_ = ei->printLvl_;
    /*
     *  Parse the complete child input file
     */
    int ok = ei_radial->electrode_input_child(commandFile, cf);
    if (ok == -1) {
        return -1;
    }
    /*
     * Switch the pointers around so that the child input file is returned.
     * Delete the original pointer.
     */
    delete ei;
    *ei_ptr = ei_radial;
    return 0;
}
//======================================================================================================================
int
Electrode_DiffTALE::electrode_model_create(ELECTRODE_KEY_INPUT* eibase)
{
    int iPh, jPh;
    /*
     *  Downcast the Key input to make sure we are being fed the correct child object
     */
    ELECTRODE_RadialDiffRegions_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialDiffRegions_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_DiffTALE::electrode_model_create()",
                           " Expecting a child ELECTRODE_RadialDiffRegions_KEY_INPUT object and didn't get it");
    }

    /*
     *  Create the base model
     */
    Electrode_Integrator::electrode_model_create(eibase);

    /*
     * Get the Number of cells from the input file
     */
    numRCells_ = ei->numRadialCellsRegions_[0];
    /*
     *  Get the volume phases which are distributed across the radial region
     */
    ELECTRODE_RadialRegion_KEY_INPUT& r0 = ei->rregions_[0];
    phaseIndeciseKRsolidPhases_ = r0.phaseIndeciseKRsolidPhases_;
    /*
     *  Get a count of the distributed phases
     */
    numSPhases_ = phaseIndeciseKRsolidPhases_.size();
   
    KRsolid_phaseNames_.resize(numSPhases_);
    for (iPh = 0; iPh < numSPhases_; iPh++) {
	KRsolid_phaseNames_[iPh] = phaseName(phaseIndeciseKRsolidPhases_[iPh]);
    }

    /*
     *  Construct the inverse mapping between regular phase indeces and distributed phases
     */
    distribPhIndexKRsolidPhases_.resize(m_NumTotPhases, -1);
    for (iPh = 0; iPh <  m_NumTotPhases; iPh++) {
	distribPhIndexKRsolidPhases_[iPh] = -1;
    }
    for  (jPh = 0; jPh < numSPhases_; jPh++) {
	iPh = phaseIndeciseKRsolidPhases_[jPh];
	distribPhIndexKRsolidPhases_[iPh] = jPh;
    }

    /*
     *  Construct the null mapping of phaseIndeciseKRsolidPhases_;
     */
    numNonSPhases_ = 0;
    phaseIndeciseNonKRsolidPhases_.clear();
    for (iPh = 0; iPh <  m_NumTotPhases; iPh++) {
	if (std::find(phaseIndeciseKRsolidPhases_.begin(), phaseIndeciseKRsolidPhases_.end(), iPh)
	    == phaseIndeciseKRsolidPhases_.end()) {
	    phaseIndeciseNonKRsolidPhases_.push_back(iPh);
	    numNonSPhases_++;
	}
    }


    /*
     *  Calculate the number of equations at each node from phaseIndeciseKRsolidPhases_
     *   2 + sum (nsp_each_distrib_phase)
     */
    numKRSpecies_ = 0;
    numSpeciesInKRSolidPhases_.clear();
    kstartKRSolidPhases_.clear();
    for (int i = 0; i < numSPhases_; i++) {
        iPh =  phaseIndeciseKRsolidPhases_[i];
        ThermoPhase& th = thermo(iPh);
        int nsp = th.nSpecies();
	numSpeciesInKRSolidPhases_.push_back(nsp);
	kstartKRSolidPhases_.push_back(numKRSpecies_);
        numKRSpecies_ += nsp;
    }
    numEqnsCell_ =  numKRSpecies_ + 2;

    /*
     * Set the size of the diffusion coefficient array and set it to a default value.
     *   -> lots of room for embelishment later.
     */
    Diff_Coeff_KRSolid_.resize(numKRSpecies_, ei->diffusionCoeffRegions_[0]);
    
    /*
     *  Create the mapping array KRsolid_speciesList_[]
     */
    KRsolid_speciesList_.resize(numKRSpecies_, -1);
    thermoSPhase_List_.resize(numSPhases_, 0);
    KRsolid_speciesNames_.clear();
    int KRsolid = 0;
    for  (int i = 0; i < (int) phaseIndeciseKRsolidPhases_.size(); i++) {
	iPh =  phaseIndeciseKRsolidPhases_[i];
	ThermoPhase& th = thermo(iPh);
	thermoSPhase_List_[i] = &th;
	int nsp = th.nSpecies();
	int kstart = getGlobalSpeciesIndex(iPh);
	for (int kk = 0; kk < nsp; kk++) {
	    KRsolid_speciesList_[KRsolid] = kstart + kk;
	    std::string kname = th.speciesName(kk);
	    KRsolid_speciesNames_.push_back(kname);
	    KRsolid++;
	}
    }

    /*
     *  Calculate the number of equations to be solved in the nonlinear system
     */
    neq_ = 1 +  numEqnsCell_ * numRCells_;
    
    /*
     *  We will do a cursory check of surface phases here. The assumption for this object
     *  until further work is that there is one surface, and that it is the exterior
     *  surface of the particle
     */
    if (m_NumSurPhases != 1) {
	printf("Unhandled situation about surface phases\n");
	exit(-1);
    }

    /*
     *  Identity of the surface phase for the exterior of the particle in the arrays,
     *     surfaceAreaRS_final_;
     */
    surfIndexExteriorSurface_ = 0;
 
    /*
     * Initialize the arrays in this object now that we know the number of equations
     */
    init_sizes();

    /*
     *  Do a check on total volume
     */
    double currentSolidVol = Electrode::SolidVol();
    double volPerParticle = currentSolidVol / particleNumberToFollow_;
    Radius_exterior_final_ = pow(volPerParticle * 3 / (4. * Pi), ONE_THIRD);

    /*
     *  Initialize the grid
     *      Take the gross dimensions of the domain and create a grid.
     */
    init_grid();

    /*
     *  Work on the Electrode_Integration quantities
     *     numIntegratedSrc_ = number of species in the PhaseList_
     *                         (note, this means internal diffusion is not tracked here
     *      (already set up)
     */ 
  
    create_solvers();

    /*
     *  Initialize all of the variables on the domain
     *    We take the grid size and the mole numbers from the base Electrode representation
     *    and create a distributed representation.
     *    Any disparity in mole numbers creates an error.
     */
    initializeAsEvenDistribution();

    resetCapacityDischargedToDate(); 

    return 0;
}
//====================================================================================================================
//  Set the electrode initial conditions from the input file.
/*!
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
int Electrode_DiffTALE::setInitialConditions(ELECTRODE_KEY_INPUT* eibase)
{
    ELECTRODE_RadialDiffRegions_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialDiffRegions_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_DiffTALE::electrode_model_create()",
                           " Expecting a child ELECTRODE_RadialDiffRegions_KEY_INPUT object and didn't get it");
    }
   
    int flag = Electrode_Integrator::setInitialConditions(ei);
    if (flag != 0) {
        return flag;
    }

    /*
     *  Initialize all of the variables on the domain
     *    We take the grid size and the mole numbers from the base Electrode representation
     *    and create a distributed representation.
     *    Any disparity in mole numbers creates an error.
     */
    initializeAsEvenDistribution();

    /*
     *  Set the initial state and the init_init state from the final state.
     */
    setInitStateFromFinal(true);

    return 0;
}
//====================================================================================================================
int
Electrode_DiffTALE::electrode_stateSave_create()
{
    eState_final_ = new EState_RadialDistrib();
    int rr = eState_final_->initialize(this);
    if (rr >= 0) {
        rr = 0;
    }
    return rr;
}
//====================================================================================================================
void
Electrode_DiffTALE::init_sizes()
{
    int kspCell =  numKRSpecies_ *  numRCells_;
    int nPhCell = numSPhases_ * numRCells_;

    spMoles_KRsolid_Cell_final_.resize(kspCell, 0.0);
    spMoles_KRsolid_Cell_init_.resize(kspCell, 0.0);
    spMoles_KRsolid_Cell_final_final_.resize(kspCell, 0.0);
    spMoles_KRsolid_Cell_init_init_.resize(kspCell, 0.0);

    phaseMoles_KRsolid_Cell_final_.resize(nPhCell, 0.0);

    concTot_SPhase_Cell_final_final_.resize(nPhCell, 0.0);
    concTot_SPhase_Cell_final_.resize(nPhCell, 0.0);
    concTot_SPhase_Cell_init_.resize(nPhCell, 0.0);
    concTot_SPhase_Cell_init_init_.resize(nPhCell, 0.0);
    molarDensity_SPhase_Cell_final_.resize(nPhCell, 0.0);

    concKRSpecies_Cell_init_.resize(kspCell, 0.0);
    concKRSpecies_Cell_final_.resize(kspCell, 0.0);
    concKRSpecies_Cell_init_init_.resize(kspCell, 0.0);
    concKRSpecies_Cell_final_final_.resize(kspCell, 0.0);
    
    concTot_SPhase_Cell_final_final_.resize(nPhCell, 0.0);
    concTot_SPhase_Cell_final_.resize(nPhCell, 0.0);
    concTot_SPhase_Cell_init_.resize(nPhCell, 0.0);
    concTot_SPhase_Cell_init_init_.resize(nPhCell, 0.0);


    spMf_KRSpecies_Cell_final_.resize(kspCell, 0.0);
    spMf_KRSpecies_Cell_init_.resize(kspCell, 0.0);

    MolarVolume_Ref_ = 1.0;

    rnodePos_final_final_.resize(numRCells_, 0.0);
    rnodePos_final_.resize(numRCells_, 0.0);
    rnodePos_init_.resize(numRCells_, 0.0);
    rnodePos_init_init_.resize(numRCells_, 0.0);

    rLatticeCBR_final_final_.resize(numRCells_, 0.0);
    rLatticeCBR_final_.resize(numRCells_, 0.0);
    rLatticeCBR_init_.resize(numRCells_, 0.0);
    rLatticeCBR_init_init_.resize(numRCells_, 0.0);
    rLatticeCBR_ref_.resize(numRCells_, 0.0);

    vLatticeCBR_cell_.resize(numRCells_, 0.0);

    cellBoundR_final_.resize(numRCells_, 0.0);
    cellBoundR_init_.resize(numRCells_, 0.0);
    cellBoundL_final_.resize(numRCells_, 0.0);
    cellBoundL_init_.resize(numRCells_, 0.0);

    volPP_Cell_final_.resize(numRCells_, 0.0);
    //  fracNodePos_.resize(numRCells_, 0.0);
    fracVolNodePos_.resize(numRCells_, 0.0);

    partialMolarVolKRSpecies_Cell_final_.resize(kspCell, 0.0);
    partialMolarCpKRSpecies_Cell_final_.resize(kspCell, 0.0);
    partialMolarEnthKRSpecies_Cell_final_.resize(kspCell, 0.0);
    partialMolarEnthKRSpecies_Cell_init_.resize(kspCell, 0.0);

    DspMoles_final_.resize(m_NumTotSpecies, 0.0);

    Diff_Coeff_KRSolid_.resize(numKRSpecies_, 0.0);

    DphMolesSrc_final_.resize(m_NumTotPhases, 0.0);


    actCoeff_Cell_final_.resize(kspCell, 1.0);
    actCoeff_Cell_init_.resize(kspCell, 1.0);

    numLatticeCBR_init_.resize(numRCells_, 0.0);
    numLatticeCBR_final_.resize(numRCells_, 0.0);

    int maxNumRxns = RSD_List_[0]->nReactions();
    ROP_.resize(maxNumRxns, 0.0);

}
//====================================================================================================================
void
Electrode_DiffTALE::init_grid()
{
    // inner and outer radius
    double  volContainedCell;
     double r_in = m_rbot0_;
     double rnode3;
     double r_out =  Radius_exterior_final_;
     //double partVol = Pi * inputParticleDiameter_ * inputParticleDiameter_ * inputParticleDiameter_ / 6.0;
     double nn = numRCells_ - 1.0;

     for (int iCell = 0; iCell < numRCells_; iCell++) {
	 fracVolNodePos_[iCell] = 1.0 / nn * iCell;
     }

  
     double vol_total_part = Pi * 4. / 3.0 * ((r_out * r_out * r_out) - (r_in * r_in * r_in));

     double rbot3 = r_in * r_in * r_in;
     //double vbot = rbot3 * 4.0 * Pi / 3.0;
   
     cellBoundL_final_[0] = r_in; 
     rnodePos_final_[0] = r_in;
     for (int iCell = 0; iCell < numRCells_-1; iCell++) { 
	 volContainedCell =  fracVolNodePos_[iCell+1] * vol_total_part;
	 rnode3 = rbot3 + volContainedCell * 3.0 / (4.0 * Pi);
	 rnodePos_final_[iCell+1] = pow(rnode3, ONE_THIRD);
	 if (iCell+1 == numRCells_) {
	     rnodePos_final_[iCell+1] = r_out;
	 }
	 cellBoundR_final_[iCell] = 0.5 * (rnodePos_final_[iCell+1] + rnodePos_final_[iCell]);
	 cellBoundL_final_[iCell+1] = cellBoundR_final_[iCell];
     }
     cellBoundR_final_[numRCells_-1] = rnodePos_final_[numRCells_-1];

     for (int iCell = 0; iCell < numRCells_; iCell++) {
	 rnodePos_final_final_[iCell] = rnodePos_final_[iCell];
	 rnodePos_init_[iCell] = rnodePos_final_[iCell];
	 rnodePos_init_init_[iCell] = rnodePos_final_[iCell];
	 rLatticeCBR_final_[iCell] = cellBoundR_final_[iCell];
	 rLatticeCBR_final_final_[iCell] = cellBoundR_final_[iCell];
	 rLatticeCBR_init_[iCell] = cellBoundR_final_[iCell];
	 rLatticeCBR_init_init_[iCell] = cellBoundR_final_[iCell];

	 double cbl3 = cellBoundL_final_[iCell] * cellBoundL_final_[iCell] * cellBoundL_final_[iCell];
	 double cbR3 = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
	 volPP_Cell_final_[iCell] = 4. * Pi / 3. * (cbR3 - cbl3);

	 if (iCell ==  0) {
	     rLatticeCBR_ref_[iCell] =  cellBoundR_final_[iCell];
	 } else {
	     rLatticeCBR_ref_[iCell] =  cellBoundR_final_[iCell];
	 }
     }
}
//====================================================================================================================
void
Electrode_DiffTALE::initializeAsEvenDistribution()
{
    /*
     *  Overall algorithm is to fill out the final state variables. Then populate the other times
     */
    int iCell, i, k, KRSolid,  kspCell, iphCell;
    /*
     *  First, get the mole fractions at all cell points
     */
    for (iCell = 0; iCell < numRCells_; ++iCell) {
	for (KRSolid = 0; KRSolid <  numKRSpecies_; KRSolid++) {
	    k = KRsolid_speciesList_[KRSolid];
	    i = KRSolid + numKRSpecies_ * iCell;
	    spMf_KRSpecies_Cell_final_[i] = spMf_final_[k];
	    spMf_KRSpecies_Cell_init_[i]  = spMf_final_[k];
	}
    }
    /*
     *  Calculate the total concentrations of phases at all cell points -  concTot_SPhase_Cell_final_
     *  Calculate the species concentrations at the cell points, -  concKRSpecies_Cell_final_
     */
    for (iCell = 0; iCell < numRCells_; ++iCell) {
	kspCell = iCell * numKRSpecies_;
	iphCell = iCell * numSPhases_;
	for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    // iPh = phaseIndeciseKRsolidPhases_[jRPh];
	    ThermoPhase* tp =  thermoSPhase_List_[jRPh];
	    tp->setState_TPX(temperature_, pressure_, &(spMf_KRSpecies_Cell_final_[kspCell]));
	    concTot_SPhase_Cell_final_[iphCell + jRPh] = tp->molarDensity();
	    tp->getConcentrations(&(concKRSpecies_Cell_final_[kspCell]));

	    int nsp = tp->nSpecies();
      
	    kspCell += nsp;
	}
    }
    /*
     *  Calculate species moles in each cell, spMoles_KRsolid_Cell_final_
     *  Calculate the partial molar volumes of each species, partialMolarVolKRSpecies_Cell_final_
     *  Calculate the activity coefficients.
     */
    for (iCell = 0; iCell < numRCells_; ++iCell) {
	kspCell = iCell * numKRSpecies_;
	iphCell = iCell * numSPhases_;
	double tmp = volPP_Cell_final_[iCell] * particleNumberToFollow_;

	for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    //iPh = phaseIndeciseKRsolidPhases_[jRPh];
	    ThermoPhase* tp =  thermoSPhase_List_[jRPh];
	    tp->setState_TPX(temperature_, pressure_, &(spMf_KRSpecies_Cell_final_[kspCell]));
	    int nsp = tp->nSpecies();
     	   
	    for (int k = 0; k < nsp; ++k) {
		spMoles_KRsolid_Cell_final_[kspCell + k] = tmp * concKRSpecies_Cell_final_[kspCell + k];
	    }
	    tp->getPartialMolarVolumes(&(partialMolarVolKRSpecies_Cell_final_[kspCell]));
	    tp->getActivityCoefficients(&(actCoeff_Cell_final_[kspCell]));

	    if (doThermalPropertyCalculations_) {
                tp->getPartialMolarCp(&(partialMolarCpKRSpecies_Cell_final_[kspCell]));
                tp->getPartialMolarEnthalpies(&(partialMolarEnthKRSpecies_Cell_final_[kspCell]));
            }

	    kspCell += nsp;
	}
    }

    /*
     *  Now transfer that to other states
     */
    setFinalFinalStateFromFinal();
    setInitStateFromFinal(true);
    
}
//====================================================================================================================
//    Return the total volume of solid material
/*
 *
 */
double Electrode_DiffTALE::SolidVol() const
{
    double v0_3 = 4. * Pi * m_rbot0_ * m_rbot0_ * m_rbot0_ / 3.;
    double r_ext = rLatticeCBR_final_[numRCells_ - 1];
    double Vext_3 = 4. * Pi * r_ext *  r_ext *  r_ext / 3.;
    double svol = (Vext_3 -  v0_3) * particleNumberToFollow_;
    return svol;
}
//====================================================================================================================
//  Returns the total Heat Capacity of the Material in the Solid Electrode at constant volume
/*
 *  This is an extensive quantity.
 *  (virtual from Electrode)
 *
 *  @return Joule K-1
 */
double Electrode_DiffTALE::SolidHeatCapacityCV() const
{
    int jRPh, iPh;
    if (!doThermalPropertyCalculations_) {
        for (int iph = 0; iph < m_NumTotPhases; iph++) {
            for (jRPh = 0; jRPh < numNonSPhases_; jRPh++) {
                iPh = phaseIndeciseNonKRsolidPhases_[jRPh];
                int istart = m_PhaseSpeciesStartIndex[iph];
                ThermoPhase& tp = thermo(iph);
                tp.getPartialMolarCp(&(CvPM_[istart]));
            }
        }

        for (int iCell = 0; iCell < numRCells_; iCell++) {
            int indexMidKRSpecies =  iCell * numKRSpecies_;
            int kstart = 0;
            for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
                ThermoPhase* th = thermoSPhase_List_[jRPh];
                int nSpecies = th->nSpecies();
                th->setState_TPX(temperature_, pressure_, &spMf_KRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
                th->getPartialMolarCp(&partialMolarCpKRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
                kstart += nSpecies;
            }
        }
    }
    double heatCapacity = 0.0;

    //
    // For no-distributed phases, we do the normal heatCapacity calc.
    //
    for (jRPh = 0; jRPh < numNonSPhases_; jRPh++) {
        iPh = phaseIndeciseNonKRsolidPhases_[jRPh];
        int kStart = m_PhaseSpeciesStartIndex[iPh];
        ThermoPhase& tp = thermo(iPh);
        int nspPhase = tp.nSpecies();
        if (iPh != solnPhase_) {
            for (int k = 0; k < nspPhase; k++) {
                heatCapacity += spMoles_final_[kStart + k] * CvPM_[kStart + k];
            }
        }
    }

    for (int iCell = 0; iCell < numRCells_; iCell++) {
        int kstart = iCell * numKRSpecies_;
        for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
            ThermoPhase* th = thermoSPhase_List_[jRPh];
            int nSpecies = th->nSpecies();
            for (int kSp = 0; kSp < nSpecies; kSp++) {
                heatCapacity += spMoles_KRsolid_Cell_final_[kstart + kSp] * partialMolarCpKRSpecies_Cell_final_[kstart + kSp];
            }
            kstart += nSpecies;
        }
    }
    return heatCapacity;
}
//====================================================================================================================
//  Returns the total enthalpy of the solid electrode
/*
 *  This is an extensive quantity.
 *  (virtual from Electrode)
 *
 *  @return Joule
 */
double Electrode_DiffTALE::SolidEnthalpy() const
{
    size_t jRPh, iPh;
    if (!doThermalPropertyCalculations_) {
	for (int iph = 0; iph < m_NumTotPhases; iph++) {
	    for (jRPh = 0; jRPh < (size_t) numNonSPhases_; jRPh++) {
		iPh = phaseIndeciseNonKRsolidPhases_[jRPh];
		size_t istart = m_PhaseSpeciesStartIndex[iph];
		ThermoPhase& tp = thermo(iph);
		tp.getPartialMolarEnthalpies(&(enthalpyMolar_final_[istart]));
	    }
	}

	for (int iCell = 0; iCell < numRCells_; iCell++) {
	    int indexMidKRSpecies =  iCell * numKRSpecies_;
	    int kstart = 0;
	    for (jRPh = 0; jRPh < (size_t) numSPhases_; jRPh++) {
		ThermoPhase* th = thermoSPhase_List_[jRPh];
		size_t nSpecies = th->nSpecies();
		th->setState_TPX(temperature_, pressure_, &spMf_KRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
		th->getPartialMolarCp(&partialMolarEnthKRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
		kstart += nSpecies;
	    }	
	}
    }
    double enthalpy = 0.0;

    //
    // For no-distributed phases, we do the normal enthalpy calc.
    //
    for (jRPh = 0; jRPh < (size_t) numNonSPhases_; jRPh++) {
	iPh = phaseIndeciseNonKRsolidPhases_[jRPh];
        int kStart = m_PhaseSpeciesStartIndex[iPh];
        ThermoPhase& tp = thermo(iPh);
        int nspPhase = tp.nSpecies();
        if (iPh != (size_t) solnPhase_) {
            for (size_t k = 0; k < (size_t) nspPhase; k++) {
		enthalpy += spMoles_final_[kStart + k] * enthalpyMolar_final_[kStart + k];
            }
        }
    }

    for (int iCell = 0; iCell < numRCells_; iCell++) {
        size_t kstart = iCell * numKRSpecies_;
        for (jRPh = 0; jRPh < (size_t) numSPhases_; jRPh++) {
            ThermoPhase* th = thermoSPhase_List_[jRPh];
            size_t nSpecies = th->nSpecies();
            for (size_t kSp = 0; kSp < nSpecies; kSp++) {
                enthalpy += spMoles_KRsolid_Cell_final_[kstart + kSp] * partialMolarEnthKRSpecies_Cell_final_[kstart + kSp];
            }
	    kstart += nSpecies;
        }
    }
    return enthalpy;
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
void Electrode_DiffTALE::resetStartingCondition(double Tinitial, bool doResetAlways)
{
    /*
     * If the initial time is input, then the code doesn't advance
     */
    double tbase = std::max(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-9 * tbase) && !doResetAlways) {
        return;
    }
    Electrode_Integrator::resetStartingCondition(Tinitial, doResetAlways);

    int iCell, i;
    int ntotal = numRCells_ * numKRSpecies_;
    for (i = 0; i < ntotal; ++i) {
        spMoles_KRsolid_Cell_init_init_[i] = spMoles_KRsolid_Cell_final_final_[i];
	spMoles_KRsolid_Cell_init_[i]      = spMoles_KRsolid_Cell_final_final_[i];
        concKRSpecies_Cell_init_init_[i]   = concKRSpecies_Cell_final_final_[i];
	concKRSpecies_Cell_init_[i]        = concKRSpecies_Cell_final_final_[i];
	spMf_KRSpecies_Cell_init_[i]       = spMf_KRSpecies_Cell_final_[i];
	actCoeff_Cell_init_[i]             = actCoeff_Cell_final_[i];
    }

    int iTotal =  numSPhases_ * numRCells_;
    for (i = 0; i < iTotal; ++i) {
        concTot_SPhase_Cell_init_init_[i] = concTot_SPhase_Cell_final_final_[i];
        concTot_SPhase_Cell_init_[i]      = concTot_SPhase_Cell_final_[i];
    }

    for (iCell = 0; iCell < numRCells_; ++iCell) {
        rLatticeCBR_init_init_[iCell] = rLatticeCBR_final_final_[iCell];
	rLatticeCBR_init_[iCell]      = rLatticeCBR_final_final_[iCell];
        rnodePos_init_init_[iCell]    = rnodePos_final_final_[iCell];
        rnodePos_init_[iCell]         = rnodePos_final_final_[iCell];
	cellBoundR_init_[iCell]       = cellBoundR_final_[iCell];
	cellBoundL_init_[iCell]       = cellBoundL_final_[iCell];
    }

    onRegionBoundary_init_init_ =  onRegionBoundary_final_final_;
    onRegionBoundary_init_      =  onRegionBoundary_final_final_;
}
//================================================================================================================
//  update the global phase numbers 
/*
 *     This is for distributed phases
 *    We don't calculate a mole fraction vector here. It doesn't make sense to do so.
 *    Instead we take the value of the exterior cell's mole fraction vector.
 *
 *  update phaseMoles_final_[]
 */
void Electrode_DiffTALE::updatePhaseNumbers(int iph)
{ 
    int istart = m_PhaseSpeciesStartIndex[iph];
    ThermoPhase& tp = thermo(iph);
    int nsp = m_PhaseSpeciesStartIndex[iph + 1] - istart;
    
    double tmp = 0.0;
    for (int k = 0; k < nsp; k++) {
        tmp += spMoles_final_[istart + k];
    }
    /*
     *   We will experiment with ways of treating negative values here
     */
    if (tmp > 1.0E-200) {
	phaseMoles_final_[iph] = tmp;
    } else  if (tmp < -1.0E-200) {
	phaseMoles_final_[iph] = tmp;
    } else {
	for (int k = 0; k < nsp; k++) {
	    spMoles_final_[istart + k] = 0.0;
	    phaseMoles_final_[iph] = 0.0;
	}
    }

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
//====================================================================================================================
// Take the state (i.e., the final state) within the Electrode_DiffTALE and push it up
// to the zero-dimensional parent object
/*
 *  update the following variables:
 *
 *          spMoles_final_ [] -> sum solid phase species
 *          spMf_final_[]  -> Use exterior cell values
 *          deltaVoltage_
 *          ElectrodeSolidVolume_ 
 *          Radius_exterior_final_
 */
void Electrode_DiffTALE::updateState_OneToZeroDimensions()
{
    int jRPh, iPh, kspStart;
    /*
     *  Zero out the distributed fields -> before summing up the values
     */
    for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
	iPh = phaseIndeciseKRsolidPhases_[jRPh];
	kspStart = m_PhaseSpeciesStartIndex[iPh];
	ThermoPhase* th = thermoSPhase_List_[jRPh];
	int nSpecies = th->nSpecies();
	for (int kSp = 0; kSp < nSpecies; kSp++) {
	    spMoles_final_[kspStart + kSp] = 0.0;
	}
    }
    /*
     *  Loop over cells summing into spMoles_final_[]
     */
    for (int iCell = 0; iCell < numRCells_; iCell++) {
        int kstart = iCell * numKRSpecies_;
        for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jRPh];
            kspStart = m_PhaseSpeciesStartIndex[iPh];
            ThermoPhase* th = thermoSPhase_List_[jRPh];
            int nSpecies = th->nSpecies();
            for (int kSp = 0; kSp < nSpecies; kSp++) {
                spMoles_final_[kspStart + kSp] += spMoles_KRsolid_Cell_final_[kstart + kSp];
            }
	    kstart += nSpecies;
        }
    }

    /*
     *  Calculate the total mole fractions, note this won't correspond to any
     *  particular cell
     */
    int cellSpecial_ = numRCells_ - 1;
    int kstart = cellSpecial_ * numKRSpecies_;
    for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
	iPh = phaseIndeciseKRsolidPhases_[jRPh];
	kspStart = m_PhaseSpeciesStartIndex[iPh];
	ThermoPhase* th = thermoSPhase_List_[jRPh];
	int nSpecies = th->nSpecies();
	double tmp = 0.0;
	for (int kSp = 0; kSp < nSpecies; kSp++) {
	    tmp += spMoles_final_[kspStart + kSp];
	}
	if (tmp > 1.0E-200) {
	    for (int kSp = 0; kSp < nSpecies; kSp++) {
		spMf_final_[kspStart + kSp] = spMoles_final_[kspStart + kSp] / tmp;
	    }
	}
	kstart += nSpecies;
    }
    /*
     *  For non-distributed phases, do the parent class updatePhaseNumbers() routine
     */
    for (jRPh = 0; jRPh < numNonSPhases_; jRPh++) {
	iPh = phaseIndeciseNonKRsolidPhases_[jRPh];
	Electrode::updatePhaseNumbers(iPh);
    }
    for (jRPh = 0; jRPh < numNonSPhases_; jRPh++) {
	iPh = phaseIndeciseNonKRsolidPhases_[jRPh];
	ThermoPhase& tp = thermo(iPh);
	int istart = m_PhaseSpeciesStartIndex[iPh];
	int nsp = m_PhaseSpeciesStartIndex[iPh + 1] - istart;
	phaseMoles_final_[iPh] = 0.0;
	for (int k = 0; k < nsp; k++) {
	    phaseMoles_final_[iPh] += spMoles_final_[istart + k];
	}
	// Here we set the state within the phase object for nondistributed phases
	tp.setState_TPX(temperature_, pressure_, &spMf_final_[istart]);
	tp.setElectricPotential(phaseVoltages_[iPh]);
	tp.getPartialMolarVolumes(&(VolPM_[istart]));
	tp.getElectrochemPotentials(&(spElectroChemPot_[istart]));
    }
    /*
     *  Calculate the voltage field
     */
    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];
    /*
     * Calculate the volume of the electrode phase. This is the main routine to do this.
     */
    ElectrodeSolidVolume_ = SolidVol();
    /*
     *  Calculate the exterior radius
     */
    Radius_exterior_final_ = rnodePos_final_[numRCells_ - 1];
    /*
     *  Update the surface areas
     */
    updateSurfaceAreas();
}
//====================================================================================================================
// Take the state (i.e., the final state) within the Electrode_Model and push it down
// to the ThermoPhase objects and propogate it to all other aspects of the final state
/*
 *  (virtual function from Electrode)
 *  This virtual function should be written so that child routines are not called by parent routines.
 *
 *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
 *  objects in the electrode.
 *
 *  We also take the state of the electrode as described by the mole numbers and mole fractions
 *  and calculate the geometrical details involved with the model. This includes the radii and
 *  thicknesses of the regions and the existence of the annular regions with their associated boolean flags.
 *
 *  All of these properties are defined for the _final_ state.
 *
 *  Thus, this is the main routine that reconciles all of the state information within the object.
 *  At the end of this routine, all aspects of the final state are consistent with each other.
 *
 *  prerequisites: The object must have been already created.
 *
 *  Fundamental Variables:
 *        concTot_SPhase_Cell_final_[]
 *        concKRSpecies_Cell_final_[1::N-1]
 *        rLatticeCBR_final_[iCell]  
 *        rnodePos_final_[iCell]
 *
 *  Variables to be calculated:  (all of the rest)
 *
 *          spMoles_KRsolid_Cell_final_[]
 *          concTot_SPhase_Cell_final_[]
 *          spMf_KRSpecies_Cell_final_[]
 *          partialMolarVolKRSpecies_Cell_final_[]
 *          actCoeff_Cell_final_[]
 *
 *          cellBoundR_final_[]
 *          cellBoundL_final_;
 *          volPP_Cell_final_;
 *
 * This routine imposes the condition that the cell boundaries are 1/2 between nodes.           
 */
void Electrode_DiffTALE::updateState()
{
    // Indexes
    int iCell, jRPh;
    double tmp;
    // Cubes of the cell boundaries, Right and Left, for the initial and final times
    double cbR3_final = 0.0;
    double cbL3_final = 0.0;

    /*
     *   We now have a vector of cells.
     *   Each of the cells must have their own conditions.
     *   We need to loop over the conditions and have their activivities calculated
     *
     *  Calculate the cell boundaries in a pre - loop
     */
    cellBoundL_final_[0] = rnodePos_final_[0];
    for (iCell = 0; iCell < numRCells_-1; iCell++) {
        cellBoundR_final_[iCell]  = 0.5 * (rnodePos_final_[iCell] + rnodePos_final_[iCell+1]);
	cellBoundL_final_[iCell+1] = cellBoundR_final_[iCell];
    }
    cellBoundR_final_[numRCells_-1] = rnodePos_final_[numRCells_-1];
    for (iCell = 0; iCell < numRCells_-1; iCell++) {
	double l3 =  cellBoundL_final_[iCell] * cellBoundL_final_[iCell] *  cellBoundL_final_[iCell];
	double r3 =  cellBoundR_final_[iCell] * cellBoundR_final_[iCell] *  cellBoundR_final_[iCell];
	volPP_Cell_final_[iCell]  = 4.0 * Pi / 3.0 * (r3 - l3);
    }


    for (int iCell = 0; iCell < numRCells_; iCell++) {
        /*
         *  Calculate the cell size
         */
        cbL3_final  = cbR3_final;
        cbR3_final = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
        double volCell = 4./3. * Pi * (cbR3_final  - cbL3_final);

	// Calculate the total volume of the cell
        double volTotalCell = volCell * particleNumberToFollow_;
	//printf("volTotalCell = %g\n", volTotalCell);

        int indexMidKRSpecies =  iCell * numKRSpecies_;
        int kstart = 0;
	int indexCellPhase = iCell * numSPhases_;
	/*
	 *  Loop over distributed phases
	 */
        for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
         
            ThermoPhase* th = thermoSPhase_List_[jRPh];
            int nSpecies = th->nSpecies();
	   
	    double& concTotPhase = concTot_SPhase_Cell_final_[indexCellPhase  + jRPh];
	    double concKRSpecies0Phase = concTotPhase;
	    double phaseMoles = 0.0;
            for (int kSp = 1; kSp < nSpecies; kSp++) {
                int iKRSpecies = kstart + kSp;
                /*
                 * Find the mole numbers of species in the cell, spMolesKRSpecies_Cell_final_[indexTopKRSpecies + iKRSpecies]
                 *     from concKRsolid_Cell_final_;
                 */
		concKRSpecies0Phase -= concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];
                spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] =
		    volTotalCell * concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];
		phaseMoles += spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies];
            }
	    /*
	     *  Fill in the quantities for species 0, which aren't formally part of the solution vector
	     */
	    concKRSpecies_Cell_final_[indexMidKRSpecies + 0] = concKRSpecies0Phase;
	    spMoles_KRsolid_Cell_final_[indexMidKRSpecies + 0] = volTotalCell * concKRSpecies_Cell_final_[indexMidKRSpecies + 0];
	    phaseMoles += spMoles_KRsolid_Cell_final_[indexMidKRSpecies + 0];

	    phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] = phaseMoles;
            /*
             * Find the mole fractions
             *     from spMoles_KRsolid_Cell_final_;
             */
	    if (concTotPhase > 1.0E-200) {
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMf_KRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] = concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] / concTotPhase;
		}
	    } else if (concTotPhase > 1.0E-200) {
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    tmp = spMf_KRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];
		    if (tmp < 0.0 || tmp > 1.0) {
			throw CanteraError("Electrode::updatePhaseNumbers()",
					   "Mole fractions out of bounds:" + int2str(kSp) + " " + fp2str(tmp));
		    }
		}
	    } else {
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    tmp = spMf_KRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];
		    if (tmp < 0.0 || tmp > 1.0) {
			throw CanteraError("Electrode::updatePhaseNumbers()",
					   "Mole fractions out of bounds:" + int2str(kSp) + " " + fp2str(tmp));
		    }
		    spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] = 0.0;
		}
		concTotPhase = 0.0;
	    }
	    /*
	     * Calculate the activities of the species
	     */
	    th->setState_TPX(temperature_, pressure_, &spMf_KRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
	    molarDensity_SPhase_Cell_final_[indexCellPhase  + jRPh] = th->molarDensity();
	    th->getActivityCoefficients(&actCoeff_Cell_final_[indexMidKRSpecies + kstart]);
	    th->getPartialMolarVolumes(&partialMolarVolKRSpecies_Cell_final_[indexMidKRSpecies + kstart]);

            if (doThermalPropertyCalculations_) {
                th->getPartialMolarCp(&partialMolarCpKRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
            }

	    kstart += nSpecies;
	}
    }

    /*
     *  Update the state to zero-dimensional parent fields
     */
    updateState_OneToZeroDimensions();
}

//=========================================================================
// under construction
/*
 *   This routine does several checks
 *     - Checks to see if the total phase Moles is consistent with the species moles in each cell
 *     - Checks to see if the mesh geometry for each cell is consistent with the 
 *       molar volume multiplied by phase moles.
 */
void Electrode_DiffTALE::checkGeometry() const
{
    double  cbL3_final , cbR3_final = 0.0, CBR = 0.0;
    int iCell, iPh;
    double rdel = 0.0;
    double phaseTot = 0.0;
    for (iCell = 0; iCell < numRCells_; iCell++) {

	/*
	 *  Check value of cellBoundR_final_ against rnodePos_final_
	 *   They must be at the midpoints wrt the radial coordinate.
	 */
	if (iCell <  numRCells_ - 1) {
	    CBR = 0.5 * (rnodePos_final_[iCell] + rnodePos_final_[iCell+1]);
	    rdel = fabs(cellBoundR_final_[iCell] - CBR) / CBR;
	    if (rdel > 1.0E-10) {
		printf("CBR is bad\n");
		exit(-1);
	    }
	    cbL3_final = cbR3_final;
	    cbR3_final = CBR * CBR * CBR;
	} else {
	    CBR = rnodePos_final_[iCell];
	    cbL3_final = cbR3_final;
	    cbR3_final = CBR * CBR * CBR;
	}
	/*
	 *  Calculate the total volume in a cell via the geometry. We will then compare that volume
	 *  to the volume calculated from the volume calculated using the mole numbers and partial molar volumes
	 */
	double totalVolGeom = 4. * Pi / 3.0 * (cbR3_final - cbL3_final) * particleNumberToFollow_;
	
	int indexMidKRSpecies =  iCell    * numKRSpecies_;
        int kstart = 0;
	/*
         *  Copy right side previous to left side current quantities
         */
        cbL3_final  = cbR3_final;
        //r0L3_final  = r0R3_final;
	// This is the chemical volume
	double totalCellVol = 0.0;
	for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    iPh = phaseIndeciseKRsolidPhases_[jRPh];
	    ThermoPhase* th = & thermo(iPh);
	    int nSpecies = th->nSpecies();
	    phaseTot = 0.0;
	    for (int kSp = 0; kSp < nSpecies; kSp++) {
		int iKRSpecies = kstart + kSp;
		phaseTot += spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies];
	    }
	    if (phaseTot > 1.0E-200) {
		rdel = fabs(phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] - phaseTot);
		if (rdel > 1.0E-10) {
		    printf("Electrode_DiffTALE::checkGeometry(): phaseMoles don't agree with spMoles\n");
		    exit(-1);
		}
	    }
	    // th->setState_TPX(temperature_, pressure_, &(spMf_KRSpecies_Cell_final_[indexMidKRSpecies +  kstart]));
	    //concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh] = th->molarDensity();
	    //th->getConcentrations(&(concKRSpecies_Cell_final_[indexMidKRSpecies +  kstart]));
	    totalCellVol += phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] / concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh];
	}
	rdel = fabs(totalCellVol -  totalVolGeom);
	if (rdel > 1.0E-10) {
	    printf("Electrode_DiffTALE::checkGeometry(): Moles and geometry don't match\n");
	    printf("         icell = %d,    volGeom = %g   volMoles = %g   rdelta = %g\n", 
		   iCell, totalVolGeom, totalCellVol, rdel);
	    exit(-1);
	}
    }

}
//========================================================================================================================
//  Here we check to see if we can account for the 
/*
 *
 */
void Electrode_DiffTALE::checkMoles_final_init() const
{
    double sum, sum_f, sum_i;
    bool doErr = false;
    for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	int iPh = phaseIndeciseKRsolidPhases_[jRPh];
	int iStart = getGlobalSpeciesIndex(iPh, 0);
	ThermoPhase* th = & thermo(iPh);
	int nSpecies = th->nSpecies();
	    
	for (int kSp = 0; kSp < nSpecies; kSp++) {
	    int isp = iStart + kSp;
	    sum_f = spMoles_final_[isp] - spMoleIntegratedSourceTermLast_[isp];
	    sum_i = spMoles_init_[isp];
	    sum = sum_f - sum_i;
	    if (abs(sum) > sum_i * 1.0E-6) {
		printf("Electrode_DiffTALE::checkMoles_final_init ERROR: sum = % 19.12E\n", sum);
		printf("                                isp = %2d    sum_f  = % 19.12E sum_i = % 19.12E\n",
		       isp, sum_f, sum_i);
		printf("                             spMoles_final = % 19.12E EletrodeSrc = % 19.12E spMoles_init = % 19.12E\n",
		       spMoles_final_[isp],  spMoleIntegratedSourceTermLast_[isp], spMoles_init_[isp]);
		   
		doErr = true;
	    }		
	}
    
	if (doErr) {
	    for (int kSp = 0; kSp < nSpecies; kSp++) {
		int isp = iStart + kSp;
		sum_f = spMoles_final_[isp] - spMoleIntegratedSourceTermLast_[isp];
		sum_i = spMoles_init_[isp];
		sum = sum_f - sum_i;
	
		printf("Electrode_DiffTALE::checkMoles_final_init ERROR: sum = % 19.12E\n", sum);
		printf("                                isp = %2d    sum_f  = % 19.12E sum_i = % 19.12E\n",
		       isp, sum_f, sum_i);
		printf("                             spMoles_final = % 19.12E EletrodeSrc = % 19.12E spMoles_init = % 19.12E\n",
		       spMoles_final_[isp],  spMoleIntegratedSourceTermLast_[isp], spMoles_init_[isp]);
	    
	    }		
	    exit(-1);
	}
    }
}
//========================================================================================================================

void Electrode_DiffTALE::check_final_state()
{
#ifdef DEBUG_NEWMODELS
    // While in debug mode, check the inventory of capacity to account for all electrons
    checkCapacityBalances_final();
#endif
#ifdef DEBUG_NEWMODELS
    checkMoles_final_init();
#endif
}
//========================================================================================================================
//  Calculate the diffusive flux of all distributed species at the right cell boundary of cell iCell.
/*
 *
 *  Algorithm assumes that species 0 is special. It's usually called the vacancy species. Think of it as the vacency
 *  species. We sum up the diffusive fluxes of all the other species. Then, the diffusive flux of the vacency is calculated
 *  as the negative of that sum. What we are doing is ensuring that the sum of the diffusive flux of all species is equal
 *  to zero.
 *
 *   The diffusive flux is the based on the gradient of the activity concentration rather than the concentration. 
 *   This difference is significant in many battery systems.
 *
 *  1/29/14  Confirmed that this routine makes no assumption that cell boundaries are halfway between nodes 
 *
 */
void Electrode_DiffTALE::diffusiveFluxRCB(double * const fluxRCB, int iCell, bool finalState) const  
{ 
    double caC, caR, dcadxR;
    int jPh, iPh, kSp;
    int indexMidKRSpecies =  iCell    * numKRSpecies_;
    int indexRightKRSpecies = (iCell+1) * numKRSpecies_;
    int kstart = 0;
    if (finalState) {
	double deltaX = rnodePos_final_[iCell+1] - rnodePos_final_[iCell];
	for (jPh = 0; jPh < numSPhases_; jPh++) {
	    iPh = phaseIndeciseKRsolidPhases_[jPh];
	    ThermoPhase* th = & thermo(iPh);
	    int nSpecies = th->nSpecies();
	    double fluxT = 0.0;

	    for (kSp = 1; kSp < nSpecies; kSp++) {
		int iKRSpecies = kstart + kSp;
                double caAvg = 2.0 / 
			       (actCoeff_Cell_final_[indexRightKRSpecies + iKRSpecies] + actCoeff_Cell_final_[indexMidKRSpecies + iKRSpecies]);
                caR = concKRSpecies_Cell_final_[indexRightKRSpecies + iKRSpecies] * actCoeff_Cell_final_[indexRightKRSpecies + iKRSpecies];
                caC = concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] * actCoeff_Cell_final_[indexMidKRSpecies + iKRSpecies];
		dcadxR = caAvg * (caR - caC) / deltaX; 

		// Calculate the flux of moles at the ride side cell boundary out of the cell
		fluxRCB[iKRSpecies] = - Diff_Coeff_KRSolid_[iKRSpecies] * dcadxR;
		fluxT += fluxRCB[iKRSpecies];
	    }
       
	    fluxRCB[kstart] = -fluxT;
	    kstart += nSpecies;
	}
    } else {
	double deltaX = rnodePos_init_[iCell+1] - rnodePos_init_[iCell];
	for (jPh = 0; jPh < numSPhases_; jPh++) {
	    iPh = phaseIndeciseKRsolidPhases_[jPh];
	    ThermoPhase* th = & thermo(iPh);
	    int nSpecies = th->nSpecies();

	    double fluxT = 0.0;
	    for (kSp = 1; kSp < nSpecies; kSp++) {
		int iKRSpecies = kstart + kSp;
	
                double caAvg = 2.0 / (actCoeff_Cell_init_[indexRightKRSpecies + iKRSpecies] + actCoeff_Cell_init_[indexMidKRSpecies + iKRSpecies]);
                caR = concKRSpecies_Cell_init_[indexRightKRSpecies + iKRSpecies] * actCoeff_Cell_init_[indexRightKRSpecies + iKRSpecies];
                caC = concKRSpecies_Cell_init_[indexMidKRSpecies + iKRSpecies] * actCoeff_Cell_init_[indexMidKRSpecies + iKRSpecies];
		dcadxR = caAvg * (caR - caC) / deltaX;

		// Calculate the flux of moles at the ride side cell boundary out of the cell
		fluxRCB[iKRSpecies] = - Diff_Coeff_KRSolid_[iKRSpecies] * dcadxR;
		fluxT += fluxRCB[iKRSpecies];
	    }
       
	    fluxRCB[kstart] = -fluxT;
	    kstart += nSpecies;

	}
    }
}
//=========================================================================================================================
// Predict the solution
/*
 * Ok at this point we have a time step deltalimiTsubcycle_
 * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
 * We now calculate predicted solution components from these conditions.
 *
 *  Essentially this routine has to come up with data on  deltaSubcycleCalc_  ,  rLatticeCBR_final_[iCell];,
 *    concTot_SPhase_Cell_final_[iCell * numSPhases_ + jPh] and
 *     concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies]
 *
 *         Residual (Time)                                     deltaSubcycleCalc_                   0
 *                                                                                            1
 *         Loop over cells                                                            0 <=  iCell < numRCells_
 *                                                                                     j = numEqnsCell_ * iCell
 *                    
 *            Residual (Reference/lattice Position)           rLatticeCBR_final_[iCell];        (1+j)
 *            Residual (Mesh Position)                        rnodePos_final_[iCell]            (j+1) + 1
 *            Loop over distributed Phases
 *            Residual (Concentration _ k=0)                  concTot_SPhase_Cell_final_[iCell * numSPhases_ + jPh]
 *              . . .
 *            Residual (Concentration _ k=Ns-1)               concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies]
 *
 * @return   Returns the success of the operation
 *                 1  A predicted solution is achieved
 *                 2  A predicted solution with a multispecies phase pop is achieved
 *                 0  A predicted solution is not achieved, but go ahead anyway
 *                -1  The predictor suggests that the time step be reduced and a retry occur.
 */
int Electrode_DiffTALE::predictSolnResid()
{
    // Indexes
    int iCell, iPh, jPh;
 
    // Cubes of the cell boundaries, Right and Left, for the initial and final times
    double cbR3_final = 0.0;
    double cbL3_final = 0.0;
    double newRtop, volToBeMoved, newCBR;
 

    // Diffusive flux
    double fluxR[10];
 
    // solution index at the start of the current phase at the current cell
    //int cIndexPhStart;

    // Location of the right cell boundary at the beginning of the step
    std::vector<doublereal> cellBoundR_init(numRCells_);
    // Velocity of the cell boundary during the time step;
    std::vector<doublereal> cellBoundRVeloc(numRCells_);
    // Node velocity during the time step
    std::vector<doublereal> rnodeVeloc(numRCells_);
    numLattices_pred_ = 0.0;
    numLattices_init_ = 0.0;
    checkGeometry();

    // predict that the calculated deltaT is equal to the input deltaT
    deltaTsubcycleCalc_ = deltaTsubcycle_;

#ifdef DEBUG_MODE
    if ( counterNumberSubIntegrations_ > 13000) {
//	printf("we are here\n");
    }
#endif

    /*
     *  Calculate the cell boundaries in a pre - loop
     */
    for (iCell = 0; iCell < numRCells_; iCell++) {
	rnodePos_final_[iCell] = rnodePos_init_[iCell];
        rnodeVeloc[iCell] = (rnodePos_final_[iCell] - rnodePos_init_[iCell]) / deltaTsubcycleCalc_;
	if (iCell == numRCells_ - 1) {
	    cellBoundR_init[iCell]   = rnodePos_init_[iCell];
	} else {
	    cellBoundR_init[iCell]   = 0.5 * (rnodePos_init_[iCell] + rnodePos_init_[iCell+1]);
	}
	cellBoundR_final_[iCell] = cellBoundR_init[iCell];
        cellBoundRVeloc[iCell] = (cellBoundR_final_[iCell] - cellBoundR_init[iCell]) / deltaTsubcycleCalc_;

	rLatticeCBR_ref_[iCell] = cellBoundR_init[iCell];	
    }

    // ---------------------------  Main Loop Over Cells ----------------------------------------------------

    /*
     *  First     1) Move species due to diffusion
     *            2) Add species due to reaction at the edge
     */
    for (int iCell = 0; iCell < numRCells_; iCell++) {
        /*
         *  Copy right side previous to left side current quantities
         */
        cbL3_final  = cbR3_final;
	// r0L3_final  = r0R3_final;
        /*
         *  Calculate the area of the outer cell which conserves constant functions under mesh movement wrt the Reynolds transport theorum
         */
        double cellBoundR_star2 = (cellBoundR_final_[iCell] * cellBoundR_final_[iCell]
                                   - cellBoundR_final_[iCell] * cellBoundR_init[iCell] + cellBoundR_final_[iCell] * cellBoundR_final_[iCell])/3.0;
        double areaR_star = 4.0 * Pi * cellBoundR_star2 * particleNumberToFollow_;

        int indexMidKRSpecies =  iCell    * numKRSpecies_;
        int indexRightKRSpecies = (iCell+1) * numKRSpecies_;
        int kstart = 0;

	/*
	 *  We calculate the diffusion step only on the rhs of each cell. So we skip the last cell.
	 */
	if (iCell != numRCells_ - 1) {
	    /*
	     *  Find the diffusive flux based on the initial conditions at the right cell boundary
	     */
	    diffusiveFluxRCB(fluxR, iCell, false);
	    /*
	     *  Then distribute the diffusive fluxes 
	     */
	    for (jPh = 0; jPh < numSPhases_; jPh++) {
		iPh = phaseIndeciseKRsolidPhases_[jPh];
		ThermoPhase* th = & thermo(iPh);
		int nSpecies = th->nSpecies();
		/*
		 *  Calculation of an explicit diffusion step
		 *  Here we calculate the 0th species as well
		 */
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    // Distribute the flux in an explicit step
		    spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies]   -= fluxR[iKRSpecies] * areaR_star * deltaTsubcycleCalc_;
		    spMoles_KRsolid_Cell_final_[indexRightKRSpecies + iKRSpecies] += fluxR[iKRSpecies] * areaR_star * deltaTsubcycleCalc_;
		    if (spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] < 0.0) {
			return -1;
		    }
		    if (spMoles_KRsolid_Cell_final_[indexRightKRSpecies + iKRSpecies] < 0.0) {
			return -1;
		    }
		}
		kstart += nSpecies;
	    }
	}

	/*
	 *  The end cell has a special treatment.
	 *  There is a reaction there that injects species here.
	 */
        if (iCell == numRCells_ - 1) {
   
            double SolidVolCreation = 0.0;
            for (jPh = 0; jPh < numSPhases_; jPh++) {
                iPh = phaseIndeciseKRsolidPhases_[jPh];
                int iStart = getGlobalSpeciesIndex(iPh, 0);
                ThermoPhase* th = & thermo(iPh);
                int nSpecies = th->nSpecies();
                for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] += DspMoles_final_[iStart + kSp] * deltaTsubcycleCalc_;
		    if (spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] < 0.0) {
			return -1;
		    }
                    SolidVolCreation += partialMolarVolKRSpecies_Cell_final_[iStart + kSp] * DspMoles_final_[iStart + kSp] * deltaTsubcycleCalc_;
                }
		kstart += nSpecies;
            }	    
	    Radius_exterior_final_ = SolidVolCreation / areaR_star + Radius_exterior_init_;
	}
    }

    /*
     *  Previously we had fixed spMoles. Now, we figure out all of the other information
     *  Calculate all other quantities
     *              Calculate control boundaries based on real values
     *              Calculate external radius based on total stuff.
     */
    for (int iCell = 0; iCell < numRCells_; iCell++) {
	int indexMidKRSpecies =  iCell * numKRSpecies_;
        int kstart = 0;
	/*
         *  Copy right side previous to left side current quantities
         */
        cbL3_final  = cbR3_final;
        //r0L3_final  = r0R3_final;
	double totalCellVol = 0.0;
	for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    iPh = phaseIndeciseKRsolidPhases_[jRPh];
	    ThermoPhase* th = & thermo(iPh);
	    int nSpecies = th->nSpecies();
	    phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] = 0.0;
	    for (int kSp = 0; kSp < nSpecies; kSp++) {
		int iKRSpecies = kstart + kSp;
		phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] += spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies];
	    }
	    double tot = phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh];
	    if (tot > 1.0E-200) {
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMf_KRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] = spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] / tot;
		}
	    }
	    th->setState_TPX(temperature_, pressure_, &(spMf_KRSpecies_Cell_final_[indexMidKRSpecies +  kstart]));
	    /*
	     * Now that we have different mole fractions, we have changed the area of the cell
	     */
	    concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh] = th->molarDensity();
	    //molarDensity_ = th->molarDensity();
	    th->getConcentrations(&(concKRSpecies_Cell_final_[indexMidKRSpecies +  kstart]));
	
	    th->getPartialMolarVolumes(&(partialMolarVolKRSpecies_Cell_final_[indexMidKRSpecies]));

	    totalCellVol += phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] / concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh];
	}

	double totalLeftVol = 4.0 * Pi / 3.0 * cbL3_final * particleNumberToFollow_;
	double totalEnclosedVol = totalLeftVol + totalCellVol;
	cbR3_final = totalEnclosedVol * 3.0 / ( 4.0 * Pi * particleNumberToFollow_);
	cellBoundR_final_[iCell] = pow(cbR3_final, ONE_THIRD);
	rLatticeCBR_final_[iCell] = cellBoundR_init_[iCell];
	if (iCell <  numRCells_ - 1) {
	    rnodePos_final_[iCell+1] = rnodePos_final_[iCell] + 2.0 * (cellBoundR_final_[iCell] - rnodePos_final_[iCell]);
	} else {
	    double oldRtop = rnodePos_final_[iCell]; 
	    double newRtop = cellBoundR_final_[iCell];
	    if (fabs(oldRtop - newRtop) > 1.0E-14) {
		// This ignores the possibility of creation/destruction of lattice sites
		rnodePos_final_[iCell] = cellBoundR_final_[iCell];
		Radius_exterior_final_ = rnodePos_final_[iCell];
		//cellBoundR_final_[iCell-1] =  0.5 * (rnodePos_final_[iCell-1] + rnodePos_final_[iCell]);
		//todo fix the moles in the two cells so that it is conservative.
	    }
	}
    }

    /*
     *  Now, we have a new exterior dimension. rnodePos_final_[numRCells_ - 1] = Radius_exterior_final_.
     *  We now have to move the mesh to conform to the spline treatment.
     *
     */
    double rtop = rnodePos_final_[numRCells_ - 1];
    for (int iCell = 0; iCell < numRCells_-1; iCell++) {
	
	int indexMidKRSpecies   =  iCell      * numKRSpecies_;
	int indexRightKRSpecies = (iCell + 1) * numKRSpecies_;
	int kstart = 0;

	if (iCell <  numRCells_-2) {
	    double I_j   = rnodePos_final_[iCell] * rnodePos_final_[iCell] * rnodePos_final_[iCell] / (rtop * rtop * rtop);
	    double tot = I_j + fracVolNodePos_[iCell+1];
	    double rbar = pow(tot, ONE_THIRD);
	    newRtop = rtop * rbar;
	} else {
	    newRtop = rnodePos_final_[numRCells_-1];
	}
	// new cell boundary betwen iCell and iCell+1
	newCBR = 0.5 * (rnodePos_final_[iCell] + newRtop);
	/*
	 *  Move stuff from one cell to another
	 */
	if (newCBR > cellBoundR_final_[iCell]) {
	    volToBeMoved = 4. * Pi / 3.0 * particleNumberToFollow_ *
			   (newCBR * newCBR * newCBR - cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell]);

	    for (jPh = 0; jPh < numSPhases_; jPh++) {
		iPh = phaseIndeciseKRsolidPhases_[jPh];
		ThermoPhase* th = & thermo(iPh);
		int nSpecies = th->nSpecies();
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    double spMolesToBeMoved = volToBeMoved * concKRSpecies_Cell_final_[indexRightKRSpecies + iKRSpecies];
		    spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies  ] += spMolesToBeMoved;
		    spMoles_KRsolid_Cell_final_[indexRightKRSpecies + iKRSpecies] -= spMolesToBeMoved;
		}
		kstart += nSpecies;
	    }

	} else if (newCBR < cellBoundR_final_[iCell]) {
	    volToBeMoved = 4. * Pi / 3.0 * particleNumberToFollow_ *
			   (cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell] - newCBR * newCBR * newCBR);

	    for (jPh = 0; jPh < numSPhases_; jPh++) {
		iPh = phaseIndeciseKRsolidPhases_[jPh];
		ThermoPhase* th = & thermo(iPh);
		int nSpecies = th->nSpecies();
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    double spMolesToBeMoved = volToBeMoved * concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];
		    spMoles_KRsolid_Cell_final_[indexMidKRSpecies   + iKRSpecies] -= spMolesToBeMoved;
		    spMoles_KRsolid_Cell_final_[indexRightKRSpecies + iKRSpecies] += spMolesToBeMoved;
		}
		kstart += nSpecies;
	    }

	}

	cellBoundR_final_[iCell] = newCBR;
	rnodePos_final_[iCell+1] = newRtop;
    }

    /*
     *  After we have changed the moles, we calculate all of the other quantities.
     *  Mole fractions, concentrations, and partial molar volumes.
     */
    cbR3_final = 0.0;
    double cbR3_init = 0.0;
    double cbL3_init = 0.0;
    double CBR = 0.0, rdel, ratio, CBR_init;
    for (int iCell = 0; iCell < numRCells_; iCell++) {

	int indexMidKRSpecies   =  iCell * numKRSpecies_;
	int kstart = 0;
	for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    iPh = phaseIndeciseKRsolidPhases_[jRPh];
	    ThermoPhase* th = & thermo(iPh);
	    int nSpecies = th->nSpecies();
	    phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] = 0.0;
	    for (int kSp = 0; kSp < nSpecies; kSp++) {
		int iKRSpecies = kstart + kSp;
		phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] += spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies];
	    }
	    double tot = phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh];
	    if (tot > 1.0E-200) {
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMf_KRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] = spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] / tot;
		}
	    }
	    th->setState_TPX(temperature_, pressure_, &(spMf_KRSpecies_Cell_final_[indexMidKRSpecies + kstart]));
	    concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh] = th->molarDensity();
	    th->getConcentrations(&(concKRSpecies_Cell_final_[indexMidKRSpecies +  kstart]));
	    th->getPartialMolarVolumes(&(partialMolarVolKRSpecies_Cell_final_[indexMidKRSpecies]));
	}


	if (iCell <  numRCells_ - 1) {
	    CBR = 0.5 * (rnodePos_final_[iCell] + rnodePos_final_[iCell+1]);
	    CBR_init = 0.5 * (rnodePos_init_[iCell] + rnodePos_init_[iCell+1]);
	    rdel = fabs(cellBoundR_final_[iCell] - CBR) / CBR;
	    if (rdel > 1.0E-10) {
		printf("CBR is bad\n");
		exit(-1);
	    }
	    cbL3_final = cbR3_final;
	    cbR3_final = CBR * CBR * CBR;
	    cbL3_init = cbR3_init;
	    cbR3_init = CBR_init * CBR_init * CBR_init;
	} else {
	    CBR = rnodePos_final_[iCell];
	    CBR_init = rnodePos_init_[iCell];
	    rdel = fabs(cellBoundR_final_[iCell] - CBR) / CBR;
	    if (rdel > 1.0E-10) {
		printf("CBR is bad\n");
		exit(-1);
	    }
	    cbL3_final = cbR3_final;
	    cbR3_final = CBR * CBR * CBR;
	    cbL3_init = cbR3_init;
	    cbR3_init = CBR_init * CBR_init * CBR_init;
	}


      	double totalVolGeom = 4. * Pi / 3.0 * (cbR3_final - cbL3_final) * particleNumberToFollow_;
	double totalCellVol = 0.0;
	for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    iPh = phaseIndeciseKRsolidPhases_[jRPh];	   
	    totalCellVol += phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] / concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh];
	}
	ratio = totalVolGeom / totalCellVol;
	if (fabs(ratio - 1.0) > 1.0E-6) {
	    printf("predictor is causing an error\n");
	}

	if (ratio != 1.0) {
	    for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
		iPh = phaseIndeciseKRsolidPhases_[jRPh];
		ThermoPhase* th = & thermo(iPh);
		int nSpecies = th->nSpecies();
		phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] *= ratio;
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] *= ratio;
		}
	    }
	}

	double vol_final =  4.0 / 3.0 * Pi * (cbR3_final - cbL3_final) * particleNumberToFollow_;
	double vol_init =  4.0 / 3.0 * Pi * (cbR3_init - cbL3_init) * particleNumberToFollow_;
	numLattices_pred_ += vol_final * concTot_SPhase_Cell_final_[iCell];
	numLattices_init_ += vol_init * concTot_SPhase_Cell_init_[iCell];
	
    }

    
    
    // Check the internal consistency
    checkGeometry();

    return 1;

}
//=========================================================================================================================
// Predict the solution
/*
 * Ok at this point we have a time step deltalimiTsubcycle_
 * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
 * We now calculate predicted solution components from these conditions.
 *
 *  Essentially this routine has to come up with data on  deltaSubcycleCalc_  ,  rLatticeCBR_final_[iCell];,
 *    concTot_SPhase_Cell_final_[iCell * numSPhases_ + jPh] and
 *     concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies]
 *
 *         Residual (Time)                                     deltaSubcycleCalc_                   0
 *                                                                                            1
 *         Loop over cells                                                            0 <=  iCell < numRCells_
 *                                                                                     j = numEqnsCell_ * iCell
 *                    
 *            Residual (Reference/lattice Position)           rLatticeCBR_final_[iCell];        (1+j)
 *            Residual (Mesh Position)                        rnodePos_final_[iCell]            (j+1) + 1
 *            Loop over distributed Phases
 *            Residual (Concentration _ k=0)                  concTot_SPhase_Cell_final_[iCell * numSPhases_ + jPh]
 *              . . .
 *            Residual (Concentration _ k=Ns-1)               concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies]
 *
 * @return   Returns the success of the operation
 *                 1  A predicted solution is achieved
 *                 2  A predicted solution with a multispecies phase pop is achieved
 *                 0  A predicted solution is not achieved, but go ahead anyway
 *                -1  The predictor suggests that the time step be reduced and a retry occur.
 *
 *          EXPERIMENTAL CODE
 */
int Electrode_DiffTALE::predictSolnResid2()
{
    // Indexes
    int iCell, iPh, jPh;
 
    // Cubes of the cell boundaries, Right and Left, for the initial and final times
    double cbR3_final = 0.0;
    double cbL3_final = 0.0;
 

    // Diffusive flux
    double fluxR[10];
 
    // solution index at the start of the current phase at the current cell
    //int cIndexPhStart;

    // Location of the right cell boundary at the beginning of the step
    std::vector<doublereal> cellBoundR_init(numRCells_);
    // Velocity of the cell boundary during the time step;
    std::vector<doublereal> cellBoundRVeloc(numRCells_);
    // Node velocity during the time step
    std::vector<doublereal> rnodeVeloc(numRCells_);
    numLattices_pred_ = 0.0;
    numLattices_init_ = 0.0;
    checkGeometry();

    // predict that the calculated deltaT is equal to the input deltaT
    deltaTsubcycleCalc_ = deltaTsubcycle_;


    /*
     *  Calculate the cell boundaries in a pre - loop
     */
    for (iCell = 0; iCell < numRCells_; iCell++) {
	rnodePos_final_[iCell] = rnodePos_init_[iCell];
        rnodeVeloc[iCell] = 0.0;
	if (iCell == numRCells_ - 1) {
	    cellBoundR_init[iCell]   = rnodePos_init_[iCell];
	} else {
	    cellBoundR_init[iCell]   = 0.5 * (rnodePos_init_[iCell] + rnodePos_init_[iCell+1]);
	}
	cellBoundR_final_[iCell] = cellBoundR_init[iCell];
        cellBoundRVeloc[iCell] = 0.0;

	rLatticeCBR_ref_[iCell] = cellBoundR_init[iCell];
    }

    // ---------------------------  Main Loop Over Cells ----------------------------------------------------

    /*
     *  First     1) Move species due to diffusion
     *            2) Add species due to reaction at the edge
     *
     *   JUST CHANGE SPMOLES
     */
    for (int iCell = 0; iCell < numRCells_; iCell++) {
        /*
         *  Copy right side previous to left side current quantities
         */
        cbL3_final  = cbR3_final;
	// r0L3_final  = r0R3_final;
        /*
         *  Calculate the area of the outer cell which conserves constant functions under mesh movement wrt the Reynolds transport theorum
         */
        double cellBoundR_star2 = (cellBoundR_final_[iCell] * cellBoundR_final_[iCell]
                                   - cellBoundR_final_[iCell] * cellBoundR_init[iCell] + cellBoundR_final_[iCell] * cellBoundR_final_[iCell])/3.0;
        double areaR_star = 4.0 * Pi * cellBoundR_star2 * particleNumberToFollow_;

        int indexMidKRSpecies =  iCell    * numKRSpecies_;
        int indexRightKRSpecies = (iCell+1) * numKRSpecies_;
        int kstart = 0;

	/*
	 *  We calculate the diffusion step only on the rhs of each cell. So we skip the last cell.
	 */
	if (iCell != numRCells_ - 1) {
	    /*
	     *  Find the diffusive flux based on the initial conditions at the right cell boundary
	     */
	    diffusiveFluxRCB(fluxR, iCell, false);
	    /*
	     *  Then distribute the diffusive fluxes 
	     */
	    for (jPh = 0; jPh < numSPhases_; jPh++) {
		iPh = phaseIndeciseKRsolidPhases_[jPh];
		ThermoPhase* th = & thermo(iPh);
		int nSpecies = th->nSpecies();
		/*
		 *  Calculation of an explicit diffusion step
		 *  Here we calculate the 0th species as well
		 */
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    // Distribute the flux in an explicit step
		    spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies]   -= fluxR[iKRSpecies] * areaR_star * deltaTsubcycleCalc_;
		    spMoles_KRsolid_Cell_final_[indexRightKRSpecies + iKRSpecies] += fluxR[iKRSpecies] * areaR_star * deltaTsubcycleCalc_;
		    if (spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] < 0.0) {
			return -1;
		    }
		    if (spMoles_KRsolid_Cell_final_[indexRightKRSpecies + iKRSpecies] < 0.0) {
			return -1;
		    }
		}
		kstart += nSpecies;
	    }
	}

	/*
	 *  The end cell has a special treatment.
	 *  There is a reaction there that injects species here.
	 */
        if (iCell == numRCells_ - 1) {
            for (jPh = 0; jPh < numSPhases_; jPh++) {
                iPh = phaseIndeciseKRsolidPhases_[jPh];
                int iStart = getGlobalSpeciesIndex(iPh, 0);
                ThermoPhase* th = & thermo(iPh);
                int nSpecies = th->nSpecies();
                for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] += DspMoles_final_[iStart + kSp] * deltaTsubcycleCalc_;
		    if (spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] < 0.0) {
			return -1;
		    }
               
                }
		kstart += nSpecies;
            }	    

	}
    }

    






    /*
     *  Previously we had fixed spMoles. Now, we figure out all of the other information
     *  Calculate all other quantities
     *              Calculate control boundaries based on real values
     *              Calculate external radius based on total stuff.
     */
    for (int iCell = 0; iCell < numRCells_; iCell++) {
	int indexMidKRSpecies =  iCell * numKRSpecies_;
        int kstart = 0;
	/*
         *  Copy right side previous to left side current quantities
         */
        cbL3_final  = cbR3_final;
        //r0L3_final  = r0R3_final;
	double totalCellVol = 0.0;
	for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    iPh = phaseIndeciseKRsolidPhases_[jRPh];
	    ThermoPhase* th = & thermo(iPh);
	    int nSpecies = th->nSpecies();
	    phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] = 0.0;
	    for (int kSp = 0; kSp < nSpecies; kSp++) {
		int iKRSpecies = kstart + kSp;
		phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] += spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies];
	    }
	    double tot = phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh];
	    if (tot > 1.0E-200) {
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMf_KRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] = spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] / tot;
		}
	    }
	    th->setState_TPX(temperature_, pressure_, &(spMf_KRSpecies_Cell_final_[indexMidKRSpecies +  kstart]));
	    /*
	     * Now that we have different mole fractions, we have changed the area of the cell
	     */
	    concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh] = th->molarDensity();
	    //molarDensity_ = th->molarDensity();
	    th->getConcentrations(&(concKRSpecies_Cell_final_[indexMidKRSpecies +  kstart]));
	
	    th->getPartialMolarVolumes(&(partialMolarVolKRSpecies_Cell_final_[indexMidKRSpecies]));

	    totalCellVol += phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] / concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh];
	}

	double totalLeftVol = 4.0 * Pi / 3.0 * cbL3_final * particleNumberToFollow_;
	double totalEnclosedVol = totalLeftVol + totalCellVol;
	cbR3_final = totalEnclosedVol * 3.0 / ( 4.0 * Pi * particleNumberToFollow_);
	cellBoundR_final_[iCell] = pow(cbR3_final, ONE_THIRD);
	rLatticeCBR_final_[iCell] = cellBoundR_final_[iCell];
	if (iCell <  numRCells_ - 1) {
	    rnodePos_final_[iCell+1] = rnodePos_final_[iCell] + 2.0 * (cellBoundR_final_[iCell] - rnodePos_final_[iCell]);
	} else {
	    double oldRtop = rnodePos_final_[iCell]; 
	    double newRtop = cellBoundR_final_[iCell];
	    if (fabs(oldRtop - newRtop) > 1.0E-14) {
		// This ignores the possibility of creation/destruction of lattice sites
		rLatticeCBR_final_[iCell] = newRtop;
		rnodePos_final_[iCell] = cellBoundR_final_[iCell];
		Radius_exterior_final_ = rnodePos_final_[iCell];
		//cellBoundR_final_[iCell-1] =  0.5 * (rnodePos_final_[iCell-1] + rnodePos_final_[iCell]);
		//todo fix the moles in the two cells so that it is conservative.
	    }
	}
    }


  
    /*
     *  After we have changed the moles, we calculate all of the other quantities.
     *  Mole fractions, concentrations, and partial molar volumes.
     */
    cbR3_final = 0.0;
    double cbR3_init = 0.0;
    double cbL3_init = 0.0;
    double CBR = 0.0, rdel, CBR_init;
    for (int iCell = 0; iCell < numRCells_; iCell++) {

	int indexMidKRSpecies   =  iCell * numKRSpecies_;
	int kstart = 0;
	for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    iPh = phaseIndeciseKRsolidPhases_[jRPh];
	    ThermoPhase* th = & thermo(iPh);
	    int nSpecies = th->nSpecies();
	    phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] = 0.0;
	    for (int kSp = 0; kSp < nSpecies; kSp++) {
		int iKRSpecies = kstart + kSp;
		phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh] += spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies];
	    }
	    double tot = phaseMoles_KRsolid_Cell_final_[iCell * numSPhases_ + jRPh];
	    if (tot > 1.0E-200) {
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMf_KRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] = spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] / tot;
		}
	    }
	    th->setState_TPX(temperature_, pressure_, &(spMf_KRSpecies_Cell_final_[indexMidKRSpecies + kstart]));
	    concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh] = th->molarDensity();
	    th->getConcentrations(&(concKRSpecies_Cell_final_[indexMidKRSpecies +  kstart]));
	    th->getPartialMolarVolumes(&(partialMolarVolKRSpecies_Cell_final_[indexMidKRSpecies]));
	}


	if (iCell <  numRCells_ - 1) {
	    CBR = 0.5 * (rnodePos_final_[iCell] +  rnodePos_final_[iCell+1]);
	    CBR_init = cellBoundR_init_[iCell];
	    rdel = fabs(cellBoundR_final_[iCell] - CBR) / CBR;
	    if (rdel > 1.0E-10) {
		printf("CBR %d is bad, %g, %g\n",iCell, CBR, cellBoundR_final_[iCell] );

	    }
	    CBR = cellBoundR_final_[iCell];
	    cbL3_final = cbR3_final;
	    cbR3_final = CBR * CBR * CBR;
	    cbL3_init = cbR3_init;
	    cbR3_init = CBR_init * CBR_init * CBR_init;
	} else {
	    CBR = rnodePos_final_[iCell];
	    CBR_init = rnodePos_init_[iCell];
	    rdel = fabs(cellBoundR_final_[iCell] - CBR) / CBR;
	    if (rdel > 1.0E-10) {
		printf("Last CBR is bad\n");
		exit(-1);
	    }
	    cbL3_final = cbR3_final;
	    cbR3_final = CBR * CBR * CBR;
	    cbL3_init = cbR3_init;
	    cbR3_init = CBR_init * CBR_init * CBR_init;
	}


	double vol_final =  4.0 / 3.0 * Pi * (cbR3_final - cbL3_final) * particleNumberToFollow_;
	double vol_init =  4.0 / 3.0 * Pi * (cbR3_init - cbL3_init) * particleNumberToFollow_;
	numLattices_pred_ += vol_final * concTot_SPhase_Cell_final_[iCell];
	numLattices_init_ += vol_init * concTot_SPhase_Cell_init_[iCell];
	
    }

    
    
    // Check the internal consistency
    checkGeometry();

    return 1;

}
//==================================================================================================================
int Electrode_DiffTALE::predictSoln()
{
    int retn = -1;

    // predict that the calculated deltaT is equal to the input deltaT
    deltaTsubcycleCalc_ = deltaTsubcycle_;

    int redoSteps = 0;
    int reDo;

    do {
	reDo = 0;
	redoSteps++;
	/*
         * Copy initial to final
         */
	setFinalStateFromInit();
	/*
	 *  Update the Internal State of ThermoPhase Objects
	 */
	updateState();
	/*
	 *    Calculate ROP_, and justBornPhase_[];
	 */
	extractInfo();
	/*
	 *    Calculate DspMoles_final_[]
	 */
	updateSpeciesMoleChangeFinal();
	/*
	 *    Now do a prediction    
	 */    
	retn = predictSolnResid();
	//retn = predictSolnResid2();
	if (retn < 0) {
	    reDo = 1;
	}

	if (reDo) {
	    deltaTsubcycle_ *= 0.25;
	}

    } while (reDo);

    /*
     *  Keep a copy of the estimated soln vector to do a predictor corrector step
     */
    packNonlinSolnVector(DATA_PTR(soln_predict_));

    yvalNLS_init_[0] =  deltaTsubcycleCalc_;

    return retn;
}
//==================================================================================================================
void Electrode_DiffTALE::check_yvalNLS_init(bool doOthers)
{
    packNonlinSolnVector(DATA_PTR(yvalNLS_init_));

    if (doOthers) {
	int neqNLS = nEquations();
	for (int i = 0; i < neqNLS; i++) {
	    yvalNLS_[i] = yvalNLS_init_[i];
	    yvalNLS_final_final_[i] = yvalNLS_init_[i];
	    yvalNLS_init_init_[i]           = yvalNLS_init_[i];
	}
    }
}

//====================================================================================================================
// Unpack the soln vector
/*
 *  (virtual from Electrode_Integrator)
 *
 * --------------------------------------------------------------------------------------------------------------
 *         Residual (Time)                                     deltaSubcycleCalc_                   0
 *                                                                                            1
 *         Loop over cells                                                            0 <=  iCell < numRCells_
 *                                                                                     j = numEqnsCell_ * iCell
 *                    
 *            Residual (Reference/lattice Position)          rLatticeCBR_final_[iCell];      (1+j)
 *            Residual (Mesh Position)                                                       (j+1) + 1
 *          Loop over distributed Phases
 *            Residual (Concentration _ k=0)                  concTot_SPhase_Cell_final_[iCell * numSPhase_ + jPh]
 *              . . .
 *            Residual (Concentration _ k=Ns-1)               concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies]
 *  --------------------------------------------------------------------------------------------------------------
 */
void Electrode_DiffTALE::unpackNonlinSolnVector(const double* const y)
{ 
    int index = 0;
    int kstart, jRPh, iKRSpecies;
    deltaTsubcycleCalc_ = y[0];
    tfinal_ = tinit_ + deltaTsubcycleCalc_;
    index++;

    for (int iCell = 0; iCell < numRCells_; iCell++) {
	kstart = 0;

	rLatticeCBR_final_[iCell] = y[index];
	index++;

	rnodePos_final_[iCell] = y[index];
	index++;

	for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    // int iPh = phaseIndeciseKRsolidPhases_[jRPh];
	    int nsp = numSpeciesInKRSolidPhases_[jRPh];
	    concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh] = y[index];
	    concKRSpecies_Cell_final_[iCell * numKRSpecies_ + 0] = concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh];
	    for (int kSp = 1; kSp < nsp; kSp++) {
		iKRSpecies = kstart + kSp;
		concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies] = y[index + kSp];
		concKRSpecies_Cell_final_[iCell * numKRSpecies_ + 0] -=
		    concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies];
	    }
	    kstart += nsp;
	    index += nsp;
	}
    }
}
//==================================================================================================================
/*
 * --------------------------------------------------------------------------------------------------------------
 *         Residual (Time)                                     deltaSubcycleCalc_                   0
 *                                                                                            1
 *         Loop over cells                                                            0 <=  iCell < numRCells_
 *                                                                                     j = numEqnsCell_ * iCell
 *                    
 *            Residual (Reference/lattice Position)          rLatticeCBR_final_[iCell];      (1+j)
 *            Residual (Mesh Position)                                                       (j+1) + 1
 *          Loop over distributed Phases
 *            Residual (Concentration _ k=0)                  concTot_SPhase_Cell_final_[iCell * numSPhase_ + jPh]
 *              . . .
 *            Residual (Concentration _ k=Ns-1)               concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies]
 *  --------------------------------------------------------------------------------------------------------------
*/
void Electrode_DiffTALE::packNonlinSolnVector(double* const y) const
{
    int index = 0;
    /*
     *  Set up the detaTsubcycle variable
     */
    y[0] = deltaTsubcycleCalc_;
    index++;
    int kstart, jRPh, iKRSpecies;
    for (int iCell = 0; iCell < numRCells_; iCell++) {
	kstart = 0;

	y[index] = rLatticeCBR_final_[iCell];
	index++;

	y[index] = rnodePos_final_[iCell];
	index++;

	for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    int nsp = numSpeciesInKRSolidPhases_[jRPh];
	    y[index] = concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh];	  
	    for (int kSp = 1; kSp < nsp; kSp++) {
		iKRSpecies = kstart + kSp;
		y[index + kSp] = concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies];
	    }
	    kstart += nsp;
	    index += nsp;
	}
    }
}
//==================================================================================================================
//   Calculate the integrated source terms and do other items now that we have a completed time step
/*
 *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
 *  for the current step, and we are calculating post-processed quantities like source terms.
 */
void Electrode_DiffTALE::calcSrcTermsOnCompletedStep()
{
    bool doOneWay = false;
    if (doOneWay) {
        /*
         *  Calculate the integrated source term
         *       An alternative would be to redo the residual calculation. However, here we assume that
         *       the residual calculation has been done and the results are in _final_
         */
        for (int i = 0; i < m_NumTotSpecies; i++) {
            spMoleIntegratedSourceTermLast_[i] = spMoles_final_[i] - spMoles_init_[i];
        }
    } else {
        extractInfo();
        updateSpeciesMoleChangeFinal();
        for (int isp = 0; isp < m_NumTotSpecies; isp++) {
            spMoleIntegratedSourceTermLast_[isp] = DspMoles_final_[isp] * deltaTsubcycleCalc_;
        }
    }
    if (doThermalPropertyCalculations_) {
        integratedThermalEnergySourceTermLast_ = thermalEnergySourceTerm_EnthalpyFormulation_SingleStep();
        integratedThermalEnergySourceTerm_overpotential_Last_ = thermalEnergySourceTerm_Overpotential_SingleStep();
        integratedThermalEnergySourceTerm_reversibleEntropy_Last_ = thermalEnergySourceTerm_ReversibleEntropy_SingleStep();
    }
}
//==================================================================================================================
//  Gather the predicted solution values and the predicted integrated source terms
/*
 *  (virtual from Electrode_Integrator)
 *
 *  Both the predicted solution values and the predicted integrated source terms are used
 *  in the time step control
 */
void  Electrode_DiffTALE::gatherIntegratedSrcPrediction()
{
    /*
     *  Here we don't recalculate anything. It was previously calculated in predictSoln().
     */
    for (int isp = 0; isp < m_NumTotSpecies; isp++) {
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
double Electrode_DiffTALE::predictorCorrectorWeightedSolnNorm(const std::vector<double>& yval)
{
    double pnorm = l0normM(soln_predict_, yval, neq_, atolNLS_, rtolNLS_);

    double pnorm_dot = l0normM(soln_predict_fromDot_, yval, neq_, atolNLS_, rtolNLS_);
    if (pnorm_dot < pnorm) {
#ifdef DEBUG_MODE
	if (printLvl_ > 2) {
	    printf("Electrode_DiffTALE::predictorCorrectorWeightedSolnNorm(): pnorm_dot %g beat out pnorm %g\n",
		   pnorm_dot, pnorm);
	}
#endif
	pnorm = pnorm_dot;
    }

    return pnorm;
}
//====================================================================================================================
// Calculate the vector of predicted errors in the source terms that this integrator is responsible for
/*
 *  (virtual from Electrode_Integrator)
 *
 *    In the base implementation we assume that the there are just one source term, the electron
 *    source term.
 *    However, this will be wrong in almost all cases.
 *    The number of source terms is unrelated to the number of unknowns in the nonlinear problem.
 *    Source terms will have units associated with them.
 *    For example the integrated source term for electrons will have units of kmol
 */
void Electrode_DiffTALE::predictorCorrectorGlobalSrcTermErrorVector()
{

}
//====================================================================================================================
static double relv(double a, double b, double atol)
{
    if (a == 0.0 && b == 0.0) {
	return 0.0;
    }
    double denom = max(fabs(a), fabs(b));
    if (denom < atol) {
	denom = atol;
    }
    return fabs((a - b)/ denom);
}

//====================================================================================================================
void Electrode_DiffTALE::predictorCorrectorPrint(const std::vector<double>& yval,
						   double pnormSrc, double pnormSoln) const
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



    denom = MAX(fabs(yval[0]), fabs(soln_predict_[0]));
    denom = MAX(denom, atolVal);
    tmp = fabs((yval[0] - soln_predict_[0])/ denom);
    printf("DeltaT                   | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E |\n",
           deltaTsubcycle_, soln_predict_[0],  yval[0], yval[0] - soln_predict_[0], atolVal, tmp);
   
    int index = 1;


    for (int iCell = 0; iCell < numRCells_; iCell++) {
	int kstart = 0;

	double rLatticeCBR_final_val = yval[index];
	double rLatticeCBR_final_predict = soln_predict_[index];
	tmp =  relv(rLatticeCBR_final_predict, rLatticeCBR_final_val,  atolNLS_[index] / rtolNLS_);
	printf("rLatticeCBR     %3d      | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E |\n", iCell,
	       rLatticeCBR_init_[iCell],
	       rLatticeCBR_final_predict,
	       rLatticeCBR_final_val,
	       rLatticeCBR_final_val -  rLatticeCBR_final_predict,
	       atolNLS_[index], tmp / rtolNLS_);
	index++;

	double rnodePos_final_val = yval[index];
	double rnodePos_final_predict = soln_predict_[index];
	tmp =  relv(rnodePos_final_predict, rnodePos_final_val,  atolNLS_[index]/rtolNLS_);
	printf("rnodePos        %3d      | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E |\n", iCell,
	       rnodePos_init_[iCell],
	       rnodePos_final_predict,
	       rnodePos_final_val,
	       rnodePos_final_val -  rnodePos_final_predict,
	       atolNLS_[index], tmp / rtolNLS_);
	index++;


	for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    // int iPh = phaseIndeciseKRsolidPhases_[jRPh];
	    int nsp = numSpeciesInKRSolidPhases_[jRPh];
	    double concTot_SPhase_Cell_final_val = yval[index];
	    double concTot_SPhase_Cell_final_predict = soln_predict_[index];
	    tmp =  relv(concTot_SPhase_Cell_final_predict,  concTot_SPhase_Cell_final_val,  atolNLS_[index] / rtolNLS_);
	    printf("concTot_SPhase  %3d      | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E |\n", iCell,
		   concTot_SPhase_Cell_init_[iCell * numSPhases_ + jRPh],
		   concTot_SPhase_Cell_final_predict,
		   concTot_SPhase_Cell_final_val ,
		   concTot_SPhase_Cell_final_val -  concTot_SPhase_Cell_final_predict,
		   atolNLS_[index], tmp / rtolNLS_);
	   
	    for (int kSp = 1; kSp < nsp; kSp++) {
		//int iKRSpecies = kstart + kSp;
		double concKRSpecies_Cell_final_val = yval[index + kSp];
		double concKRSpecies_Cell_final_predict = soln_predict_[index + kSp];

		tmp = relv(concKRSpecies_Cell_final_val, concKRSpecies_Cell_final_predict, atolNLS_[index + kSp]/rtolNLS_);
		printf("ConcKRSpeci     %3d %3d  | %14.7E %14.7E %14.7E | %14.7E | %10.3E | %10.3E |\n", iCell, kSp,
		       concKRSpecies_Cell_init_[iCell * numKRSpecies_ + kSp],
		       concKRSpecies_Cell_final_predict,
		       concKRSpecies_Cell_final_val,
		       concKRSpecies_Cell_final_val - concKRSpecies_Cell_final_predict,
		       atolNLS_[index + kSp], tmp / rtolNLS_);
		       
	    }
	    kstart += nsp;
	    index += nsp;
	}

    }
    printf("total Lattices           | %14.7E %14.7E %14.7E | %14.7E |            |            |\n", numLattices_init_,
	   numLattices_pred_, numLattices_final_,  numLattices_final_ -  numLattices_init_);




   
    printf(" -------------------------------------------------------------------------------------------------------------------\n");
    printf("                                                                                                        %10.3E\n",
           pnormSoln);

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
void Electrode_DiffTALE::setNLSGlobalSrcTermTolerances(double rtolResid)
{
    double sum = SolidTotalMoles();
    double val = 1.0E-14 * sum;

    for (int i = 0; i < numIntegratedSrc_; i++) {
        atol_IntegratedSrc_global_[i] = val;
    }
    rtol_IntegratedSrc_global_ = rtolResid;
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
void  Electrode_DiffTALE::setResidAtolNLS()
{ 
    double atolMF = 1.0E-12;
    double deltaT = t_final_final_ - t_init_init_;
    deltaT = MAX(deltaT, 1.0E-3);
    /*
     *  calculate the total moles of solid. This is used to dimensionalize the atol mole numbers
     *  -> Might rethink this, because we are using concentrations here
     */
    //double solidMoles = SolidTotalMoles();
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
    //double atolVal = solidMoles * atolBaseResid_;

    /* 
     *  set the atol on the time step
     */
    atolNLS_[0] = (1.0E-50);
    int index = 0;
    index++;

    for (int iCell = 0; iCell < numRCells_; iCell++) {
	int rindex = 1 + numEqnsCell_ * iCell;
        int xindex = 2 + numEqnsCell_ * iCell;
        int cindex = 3 + numEqnsCell_ * iCell;
        int cIndexPhStart = cindex;
	/*
	 *  Set atol on the radial equation
	 */
	atolNLS_[rindex] = atolMF * Radius_exterior_final_;
	atolResidNLS_[rindex] = atolNLS_[rindex];

	/*
	 * set atol on the lattice radial equation
	 */
	atolNLS_[xindex] = atolMF * Radius_exterior_final_;
	atolResidNLS_[xindex] = atolNLS_[xindex];


	for (int jPh = 0; jPh < numSPhases_; jPh++) {
            int iPh = phaseIndeciseKRsolidPhases_[jPh];
	    /*
	     *  set the tolerance on the phase concentration
	     */ 
	    ThermoPhase* th = & thermo(iPh);
	    int nSpecies = th->nSpecies();
	    double molarVol = th->molarVolume();
	    double phaseConc = 1.0 / molarVol;
	    atolNLS_[cIndexPhStart] = phaseConc * atolMF;
	    atolResidNLS_[xindex] = atolNLS_[cIndexPhStart];

	    for (int kSp = 1; kSp < nSpecies; kSp++) {
		atolNLS_[cIndexPhStart + kSp] = phaseConc * atolMF;
		atolResidNLS_[cIndexPhStart + kSp] = atolNLS_[cIndexPhStart + kSp];
	    }
	    cIndexPhStart += nSpecies;
	}

    }

  
    if ((int) atolNLS_.size() !=  neq_) {
        printf("ERROR\n");
        exit(-1);
    }

    determineBigMoleFractions();
}
//=================================================================================================================================
//   Evaluate the residual function
/*
 * Driver for evaluating the residual
 * (virtual from NonlinearSolver)
 *
 * @param t             Time                    (input)
 * @param delta_t       The current value of the time step (input)
 * @param y             Solution vector (input, do not modify)
 * @param ySolnDot      Rate of change of solution vector. (input, do not modify)
 * @param resid         Value of the residual that is computed (output)
 * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
 * @param id_x          Index of the variable that is being numerically differenced to find
 *                      the jacobian (defaults to -1, which indicates that no variable is being
 *                      differenced or that the residual doesn't take this issue into account)
 * @param delta_x       Value of the delta used in the numerical differencing
 *
 * @return
 */
int  Electrode_DiffTALE::evalResidNJ(const doublereal t, const doublereal delta_t,
                                       const doublereal* const y,
                                       const doublereal* const ySolnDot,
                                       doublereal* const resid,
                                       const ResidEval_Type_Enum evalType,
                                       const int id_x,
                                       const doublereal delta_x)
{
    int retn = integrateResid(t, delta_t, y, ySolnDot, resid, evalType, id_x, delta_x);
    return retn;
}
//=================================================================================================================================
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
int Electrode_DiffTALE::integrateResid(const doublereal t, const doublereal delta_t,
					 const doublereal* const y, const doublereal* const ySolnDot,
					 doublereal* const resid,
					 const ResidEval_Type_Enum evalType, const int id_x,
					 const doublereal delta_x)
{
    
    //int neq = nResidEquations();

    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
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
    }

    /*
     *  UNPACK THE SOLUTION VECTOR
     */
    unpackNonlinSolnVector(y);


    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        if (phaseID_TimeDeathMin_ >= 0) {
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

    /*
     *  Update all of the final state variables.
     */
    updateState();

    /*
     *  Extract information from reaction mechanisms
     */
    extractInfo();

    /*
     *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *   all species in the electrode. (special treatment of the surface area
     */
    updateSpeciesMoleChangeFinal();

    for (int i = 0; i < neq_; i++) {
	resid[i] = 0.0;
    }

    /*
     * Calculate the residual
     */
    //calcResid_2(resid, evalType);
    calcResid(resid, evalType);

  
    int index = 1;
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {

        if (phaseID_TimeDeathMin_ >= 0) {

            printf("\t\t minPH = %3d minCell = %3d delTcalc = %16.7E    Res = %16.7E   \n",
                   phaseID_TimeDeathMin_, cellID_TimeDeathMin_, deltaTsubcycleCalc_,  resid[0]);
        }
        printf("\t\t        PhaseName        Moles_Init    Moles_final              |   Src_Moles  Pred_Moles_Final  |    Resid     |\n");
        for (int iph = 0; iph < m_NumTotPhases; iph++) {
            double src =   DphMolesSrc_final_[iph] * deltaTsubcycleCalc_;
            printf("\t\t %20.20s  %12.4e  %12.4e               | %12.4e %12.4e", PhaseNames_[iph].c_str(), phaseMoles_init_[iph],
                   phaseMoles_final_[iph], src,  phaseMoles_init_[iph] + src);
            bool found = false;
            for (int ph = 0; ph <  numSPhases_; ph++) {
                int jph = phaseIndeciseKRsolidPhases_[ph];
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
            }
            printf("\n");
        }
	printf("\t\t      SpeciesName        Moles_Init    Moles_final   SpMF       |   Src_Moles  Pred_Moles_Final\n");
        for (int k = 0; k < m_NumTotSpecies; k++) {
            string ss = speciesName(k);
            double src =  DspMoles_final_[k] * deltaTsubcycleCalc_;
            printf("\t\t %20s  %12.4e  %12.4e  %12.4e | %12.4e %12.4e", ss.c_str(), spMoles_init_[k], spMoles_final_[k],
                   spMf_final_[k], src, spMoles_init_[k] + src);
            bool found = false;
            for (int ph = 0; ph < numSPhases_; ph++) {
                int iph = phaseIndeciseKRsolidPhases_[ph];
                if (numSpeciesInKRSolidPhases_[ph] > 1) {
                    int iStart = getGlobalSpeciesIndex(iph, 0);
                    for (int sp = 0; sp < numSpeciesInKRSolidPhases_[ph]; sp++) {
			int i = iStart + sp;
			if (i == k) {
			    found = true;
			    double res =  spMoles_final_[k] - (spMoles_init_[k] + src);
			    printf("      |%12.4e  | ", res);
			}
                    }
                }
            }
            if (!found) {
                printf("      |              | ");
            }
            printf("\n");

        }

	showResidual(8, resid);

    }



    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printf("\t\t===============================================================================================================================\n");
    }
    return 1;
}
//====================================================================================================================
// Main routine to calculate the residual
/*
 *
 *  Format of the residual equations
 *                                                             Unknown                           Index
 * --------------------------------------------------------------------------------------------------------------
 *         Residual (Time)                                     deltaSubcycleCalc_                   0
 *                                                                                            1
 *         Loop over cells                                                            0 <=  iCell < numRCells_
 *                                                                                     j = numEqnsCell_ * iCell
 *                    
 *            Residual (Reference/lattice Position)           rLatticeCBR_final_[iCell];       (1+j)
 *            Residual (Mesh Position)                                                       (j+1) + 1
 *            Residual (Concentration _ k=0)
 *              . . .
 *            Residual (Concentration _ k=Ns-1)
 *  --------------------------------------------------------------------------------------------------------------
 *
 *   Notes:
 *        1/31/14 confirmed that this routine makes no inherent assumption that the cell boundaries are halfway between nodes
 *
 */
int Electrode_DiffTALE::calcResid(double* const resid, const ResidEval_Type_Enum evalType)
{
    int iCell, iPh, jPh, jCell;
    double I_j;
    double I_jp1;
    double Lleft, L3left, rjCBL3;
    double vbarLattice_final_jcell, vbarLattice_init_jcell;
    double rnow3, rnow;
    //double Lright;
    // double rLtarget;

    // Cubes of the cell boundaries, Right and Left, for the initial and final times
    double cbR3_init = 0.0;
    double cbL3_init = 0.0;
    double cbR3_final = 0.0;
    double cbL3_final = 0.0;
    // Reference radius squared, Right and Left, at final time
    //double r0R2_final = 0.0;
    //double r0L2_final = 0.0;
    // Reference radius cubed, Right and Left, at final time

    double rRefCBR3 = 0.0;

    double fluxTC, fluxM;

    // Diffusive flux
    double fluxRCB[10];
   
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

    numLattices_final_ = 0.0;
    numLattices_init_ = 0.0;

    // Location of the right cell boundary at the beginning of the step
    //  std::vector<doublereal> cellBoundR_init(numRCells_);
    // Velocity of the cell boundary during the time step;
    std::vector<doublereal> cellBoundRVeloc(numRCells_);
    // Node velocity during the time step
    std::vector<doublereal> rnodeVeloc(numRCells_);

    /*
     *    Residual equation for the time step -> Right now we don't have a model
     */
    resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;

    /*
     *  Determine the lattice velocity on the left side domain boundary
     */
    //    vLatticeR is currently set to zero -> this is the appropriate default until something more is needed
    double volFac = 4. * Pi / 3.0 * particleNumberToFollow_;
    /*
     *  Calculate the cell quantitites in a pre-loop
     */
    double cbl, cbr = 0.0;
    for (iCell = 0; iCell < numRCells_; iCell++) {
	cbl = cbr;

        cbr = cellBoundR_init_[iCell];
	cbL3_init = cbl * cbl * cbl;
	cbR3_init = cbr * cbr * cbr;
	double vol_init_ = volFac * (cbR3_init - cbL3_init);
	/*
	 *  Calculate the nodal velocity
	 */
        rnodeVeloc[iCell] = (rnodePos_final_[iCell] - rnodePos_init_[iCell]) / deltaTsubcycleCalc_;
	/*
	 *  Calculate the the velocity of the right cell boundary
	 */
        cellBoundRVeloc[iCell] = (cellBoundR_final_[iCell] - cellBoundR_init_[iCell]) / deltaTsubcycleCalc_;
	/*
	 *  The reference lattice coordinate is rebaselined so that it is linked with the radius coordinate
	 *  at every time step.
	 */
	rLatticeCBR_ref_[iCell] = cellBoundR_init_[iCell];

	if (iCell == 0) {
	    numLatticeCBR_init_[iCell] = concTot_SPhase_Cell_init_[iCell*numSPhases_] * vol_init_;
	} else {
	    numLatticeCBR_init_[iCell] =  numLatticeCBR_init_[iCell-1] + concTot_SPhase_Cell_init_[iCell*numSPhases_] * vol_init_;
	}
    }

    // ---------------------------  Main Loop Over Cells ----------------------------------------------------

    cbR3_final = 0.0;
    cbR3_init = 0.0;
    for (int iCell = 0; iCell < numRCells_; iCell++) {

        /*
         *  Copy right side previous to left side current quantities
         */
        //vLatticeCBL = vLatticeCBR;
        cbL3_final  = cbR3_final;
        cbL3_init   = cbR3_init;
        //r0L2_final  = r0R2_final;

        /*
         *  Calculate indexes for accessing the residual
         */
        rindex = 1 + numEqnsCell_ * iCell;
        xindex = 2 + numEqnsCell_ * iCell;
        cindex = 3 + numEqnsCell_ * iCell;
        cIndexPhStart = cindex;

        /*
         *  Temporary pointer variables for total concentrations of each of the phases within the cell
         */
        concTotalVec_SPhase_init  = &(concTot_SPhase_Cell_init_[iCell*numSPhases_]);
        concTotalVec_SPhase_final = &(concTot_SPhase_Cell_final_[iCell*numSPhases_]);

        /*
         * Calculate the cubes of the cell boundary radii
         */
        cbR3_init  = cellBoundR_init_[iCell]  * cellBoundR_init_[iCell]  * cellBoundR_init_[iCell];
        cbR3_final = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
        double vol_final_ = volFac * (cbR3_final - cbL3_final);
        /*
         * Calculate powers of the lattice coordinate at the right cell boundary
         */
        /*
         * Calculate the molar volume of the first phase, final and init values
         * We will assume for now that this is the lattice molar volume
         */
        if (iCell == 0) {
            numLatticeCBR_final_[iCell] = concTotalVec_SPhase_final[0] * vol_final_;
        } else {
            numLatticeCBR_final_[iCell] = numLatticeCBR_final_[iCell-1] + concTotalVec_SPhase_final[0] * vol_final_;
        }
        /*
         * Residual calculation - Value of the Lattice radius at the right cell boundary
	 *  This is Eqn. 21 in the Memo writeup
	 *      it's the residual equation for determining  rLatticeCBR_final_[iCell]
	 *
	 *      When there is no expansion of the lattice,  rLatticeCBR_final_[iCell] = cellBoundR_final_[iCell]
         */
	double deltaNum = 0.0;
	if (numLatticeCBR_final_[iCell] == numLatticeCBR_init_[iCell]) {
	    resid[rindex] = rLatticeCBR_final_[iCell] - cellBoundR_init_[iCell];
	} else if (numLatticeCBR_final_[iCell] > numLatticeCBR_init_[iCell]) {
	    // jCell is the cell to look for within init to find the solution
	    jCell = iCell + 1;
	    if (jCell > numRCells_ -1) {
		jCell = numRCells_ - 1;
		deltaNum = numLatticeCBR_final_[iCell] - numLatticeCBR_init_[iCell-1];
	    } else {
		deltaNum = numLatticeCBR_final_[iCell] - numLatticeCBR_init_[iCell];
	    }
	    cbl = cellBoundR_init_[jCell-1];
	    double cbL3_init_j = cbl * cbl * cbl;
	    double r3 = deltaNum / volFac / concTot_SPhase_Cell_init_[jCell*numSPhases_] + cbL3_init_j;
	    double r = pow(r3, ONE_THIRD);
	    resid[rindex] = rLatticeCBR_final_[iCell] - r;
	} else {
	    if (iCell == 0) {
		deltaNum = numLatticeCBR_final_[iCell];
		cbl = 0.0;
	    } else {
		deltaNum = numLatticeCBR_final_[iCell] - numLatticeCBR_init_[iCell-1];
		cbl = cellBoundR_init_[iCell-1];
	    }
	    double cbL3_init_j = cbl * cbl * cbl;
	    double r3 = deltaNum / volFac / concTotalVec_SPhase_init[0] + cbL3_init_j;
	    double r = pow(r3, ONE_THIRD);
	    resid[rindex] = rLatticeCBR_final_[iCell] - r;
	}


	//double rhs = r0L3_final + vbarLattice_init / vbarLattice_final * (cbR3_final - cbL3_final);
        //resid[rindex] = r0R_final - pow(rhs, ONE_THIRD);

	if (formulationType_ > 0) {
	    resid[rindex] = rLatticeCBR_final_[iCell] - rLatticeCBR_init_[iCell];
	}

        /*
         *  Calculate the time derivative of the molar volume
         */
        //double vbarDot = (vbarLattice_final - vbarLattice_init) / deltaTsubcycleCalc_;

        /*
         * Find the lattice velocity of the reference radius at the right cell boundary
         */
        //double vLatticeCBR = -1.0 / (3. * r0R2_final) *
	//		     (3.0 *r0L2_final * vLatticeCBL - MolarVolume_Ref_ / (vbarLattice_final * vbarLattice_final) *
	//		      vbarDot * (cbR3_final - cbL3_final));
	double vLatticeCBR = 0.0;

	//rLtarget =  rLatticeCBR_ref_[iCell];
	if ((rLatticeCBR_final_[iCell] > rLatticeCBR_ref_[iCell]) || (iCell == (numRCells_-1))) {
	    // we are here to when the lattice velocity has put the previous lattice CBR into the current cell
	    // Assign the cell id the current cell.
	    //  (In other words the current cell is expanding
	    jCell = iCell;
	    // The CBL value is the CBR of the next lowest cell.
	    // except for the lowest cell
	    if (jCell == 0) {
		Lleft =  m_rbot0_;
	    } else {
		Lleft =  rLatticeCBR_final_[jCell-1];
	    }
	    //Lright =  rLatticeCBR_final_[jCell];
	    // If the lattice has moved more than one cell, then throw an error condition for now.
	    // We'll come back to fix this later
	    if (Lleft > rLatticeCBR_ref_[iCell]) {
		printf("ERROR: lattice movement not handled\n");
		exit(-1);
	    }
	    // Take the cube of that
	    L3left = Lleft *  Lleft *  Lleft;
	    // Find the cubes of the right CBR and left CBR values
	    //rjCBR3 = cellBoundR_final_[jCell] * cellBoundR_final_[jCell] *  cellBoundR_final_[jCell];
	    rRefCBR3 = rLatticeCBR_ref_[iCell] *rLatticeCBR_ref_[iCell] * rLatticeCBR_ref_[iCell];
	    rjCBL3 = cellBoundL_final_[jCell] * cellBoundL_final_[jCell] *  cellBoundL_final_[jCell];
	    vbarLattice_final_jcell = 1.0 / concTot_SPhase_Cell_final_[iCell*numSPhases_];
	    vbarLattice_init_jcell = 1.0 /  concTot_SPhase_Cell_init_[iCell*numSPhases_];
	    rnow3 =  L3left +  vbarLattice_init_jcell / vbarLattice_final_jcell * (rRefCBR3 - rjCBL3);
	    rnow = pow(rnow3, ONE_THIRD);
	    if (rnow >  cellBoundR_final_[jCell]) {
		if (rnow > cellBoundR_final_[jCell]*(1.0 + 1.0E-14)) {
		    // printf("analyse condition 1 \n");
		}
	    }
	    vLatticeCBR = (rnow - rLatticeCBR_ref_[iCell]) / deltaTsubcycleCalc_;
	} else {
	    // we are here to when the lattice velocity has put the previous lattice CBR into the next cell
	    // Assign the cell id the current cell.
	    jCell = iCell+1;
	    // The CBL value is the CBR of the next lowest cell.
	    Lleft =  rLatticeCBR_final_[jCell-1];
	    //Lright =  rLatticeCBR_final_[jCell];
	    // If the lattice has moved more than one cell, then throw an error condition for now.
	    // We'll come back to fix this later
	    if (Lleft > rLatticeCBR_ref_[iCell]) {
		printf("ERROR: lattice movement not handled\n");
		exit(-1);
	    }
	    // Take the cube of that
	    L3left = Lleft *  Lleft *  Lleft;
	    // Find the cubes of the right CBR and left CBR values
	    // rjCBR3 = cellBoundR_final_[jCell] * cellBoundR_final_[jCell] *  cellBoundR_final_[jCell];
	    rRefCBR3 = rLatticeCBR_ref_[iCell] * rLatticeCBR_ref_[iCell] * rLatticeCBR_ref_[iCell];
	    rjCBL3 = cellBoundL_final_[jCell] * cellBoundL_final_[jCell] * cellBoundL_final_[jCell];
	    vbarLattice_final_jcell = 1.0 / concTot_SPhase_Cell_final_[iCell*numSPhases_];
	    vbarLattice_init_jcell = 1.0 /  concTot_SPhase_Cell_init_[iCell*numSPhases_];
	    rnow3 =  L3left +  vbarLattice_init_jcell / vbarLattice_final_jcell * (rRefCBR3 - rjCBL3);
	    rnow = pow(rnow3, ONE_THIRD);
	    vLatticeCBR = (rnow - rLatticeCBR_ref_[iCell]) / deltaTsubcycleCalc_;
	}
	/*
	 *  Store the calculated lattice velocity
	 */
	vLatticeCBR_cell_[iCell] = vLatticeCBR;

        /*
         * Node position residual - spline equations, it all depends on the top node's equation formulation.
         * Everything else get's dragged along with it
	 * We will calculate the bc when we've calculated the surface reaction rates below.
	 *         -> special case iCell= 0 set to m_rbot0_
	 *         ->              iCell= numRCells_-1 -> distinguishing condition set to boundary condition
         */

	if (iCell == 0) {
	    resid[xindex] = rnodePos_final_[iCell] - m_rbot0_;
	} else if (iCell == numRCells_-1) {
	    // This condition is handled below
	} else {
	    double rtop = rnodePos_final_[numRCells_ - 1];
	    I_j   =  rnodePos_final_[iCell] * rnodePos_final_[iCell] * rnodePos_final_[iCell] / (rtop * rtop * rtop);
	    I_jp1 =  rnodePos_final_[iCell+1] * rnodePos_final_[iCell+1] * rnodePos_final_[iCell+1] / (rtop * rtop * rtop);
	    resid[xindex] = I_jp1 - I_j - fracVolNodePos_[iCell];
	}
	if (formulationType_ > 0) {
	    resid[xindex] = rnodePos_final_[iCell] - rnodePos_init_[iCell];
	}
    

        /*
         *  Calculate the area of the outer cell which conserves constant functions under mesh movement wrt the Reynolds transport theorum
         */
        double cellBoundR_star2 = (cellBoundR_init_[iCell] * cellBoundR_init_[iCell]
                                   + cellBoundR_final_[iCell] * cellBoundR_init_[iCell] + cellBoundR_final_[iCell] * cellBoundR_final_[iCell])/3.0;
        double areaR_star = 4.0 * Pi * cellBoundR_star2 * particleNumberToFollow_;

   
        int kstart = 0;
	/*
	 * Calculate the molar diffusive flux at the right cell boundary. Note the sum of the molar fluxes over all species
	 * is zero. This is crucial. 
	 */
	if (iCell < (numRCells_-1)) {
	    diffusiveFluxRCB(fluxRCB, iCell, true);
	}

        for (jPh = 0; jPh < numSPhases_; jPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jPh];
            ThermoPhase* th = & thermo(iPh);
            int nSpecies = th->nSpecies();

	    /*
	     *  Residual for the total concentration of the phase in the cell
	     *   -  Skip the first species
	     */
            double old_stuff = concTotalVec_SPhase_init[jPh]  * 4.0 / 3.0 * Pi * (cbR3_init -  cbL3_init)  * particleNumberToFollow_;
            double new_stuff = concTotalVec_SPhase_final[jPh] * 4.0 / 3.0 * Pi * (cbR3_final - cbL3_final) * particleNumberToFollow_;
            /*
             *  Change in the total amount of moles of phase jPh
             */
            resid[cIndexPhStart + kstart] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;

            /*
             * Convective flux - mesh movement with NO material movement
             */
            double vtotalR = cellBoundRVeloc[iCell] - vLatticeCBR;
	    if (formulationType_ > 0) {
		vtotalR = 0.0;
	    }
            if (iCell < (numRCells_- 1)) {
                if (vtotalR >= 0.0) {
                    fluxTC = vtotalR * concTot_SPhase_Cell_final_[numSPhases_*(iCell+1) + jPh];
                } else {
                    fluxTC = vtotalR * concTot_SPhase_Cell_final_[numSPhases_*(iCell)   + jPh];
                }
                resid[cIndexPhStart + kstart]                -= fluxTC * areaR_star;
                resid[cIndexPhStart + kstart + numEqnsCell_] += fluxTC * areaR_star;
            }

	    if (formulationTypeTotalConc_ == 1) {
		resid[cIndexPhStart + kstart] =  concTotalVec_SPhase_final[jPh] - molarDensity_SPhase_Cell_final_[iCell*numSPhases_ + jPh];
	    }

            /*
             *  Residual for the species in the phase
             *   -  Skip the first species
             */
            for (int kSp = 1; kSp < nSpecies; kSp++) {
                int iKRSpecies = kstart + kSp;
                /*
                 * Diffusive flux at the outer boundary for each species
                 */
                if (iCell < (numRCells_ - 1)) {
                    resid[cIndexPhStart + iKRSpecies]                += fluxRCB[iKRSpecies] * areaR_star;
                    resid[cIndexPhStart + iKRSpecies + numEqnsCell_] -= fluxRCB[iKRSpecies] * areaR_star;
                }

                /*
                 * Convective flux - mesh movement with NO material movement
                 */
                if (iCell < (numRCells_- 1)) {
                    if (vtotalR >= 0.0) {
                        fluxM = vtotalR * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell+1) + iKRSpecies];
                        resid[cIndexPhStart + iKRSpecies]                -= fluxM * areaR_star;
                        resid[cIndexPhStart + iKRSpecies + numEqnsCell_] += fluxM * areaR_star;
                    } else {
                        fluxM = vtotalR * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell)   + iKRSpecies];
                        resid[cIndexPhStart + iKRSpecies]                -= fluxM * areaR_star;
                        resid[cIndexPhStart + iKRSpecies + numEqnsCell_] += fluxM * areaR_star;
                    }
                }

                /*
                 *  Change in the total amount of moles of phase jPh
                 */
                old_stuff = concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]  * 4.0 / 3.0 * Pi * (cbR3_init  - cbL3_init) * particleNumberToFollow_;
                new_stuff = concKRSpecies_Cell_final_[numKRSpecies_ * iCell + iKRSpecies] * 4.0 / 3.0 * Pi * (cbR3_final - cbL3_final)* particleNumberToFollow_;
                resid[cIndexPhStart + iKRSpecies] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;
            }

            // end of phase loop
            kstart += nSpecies;
        }

        if (iCell == numRCells_ - 1) {
            /*
             *  Do the exterior
             */
            double SolidVolCreationRate = 0.0;
	    kstart = 0;
            for (jPh = 0; jPh < numSPhases_; jPh++) {
                iPh = phaseIndeciseKRsolidPhases_[jPh];
                int iStart = getGlobalSpeciesIndex(iPh, 0);
                ThermoPhase* th = & thermo(iPh);
                int nSpecies = th->nSpecies();
                resid[cIndexPhStart + kstart] -= DphMolesSrc_final_[iPh];
                SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart] * DspMoles_final_[iStart];
                for (int kSp = 1; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
                    resid[cIndexPhStart + iKRSpecies] -= DspMoles_final_[iStart + kSp];
                    SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart + kSp] *  DspMoles_final_[iStart + kSp];
                }


		if (formulationTypeTotalConc_ == 1) {
		    resid[cIndexPhStart + kstart] = concTotalVec_SPhase_final[jPh] - molarDensity_SPhase_Cell_final_[iCell*numSPhases_ + jPh];
		}

		kstart += nSpecies;
            }
            /*
             * Calculation of the movement at the top of the
             *    C dXdt * Area - SolidCreationRate = 0
             *
             *  this is an expression of the conservation of total solid molar volume at the
             *  r_exterior.
             */

	    if (formulationType_ > 0) {
		resid[xindex] += 0.0;
	    } else {
		resid[xindex] = rnodeVeloc[iCell] * areaR_star - SolidVolCreationRate;
	    }

        }

	double vol_final =  4.0 / 3.0 * Pi * (cbR3_final - cbL3_final) * particleNumberToFollow_;
	double vol_init =  4.0 / 3.0 * Pi * (cbR3_init - cbL3_init) * particleNumberToFollow_;
	numLattices_final_ += vol_final * concTot_SPhase_Cell_final_[iCell];
	numLattices_init_ += vol_init * concTot_SPhase_Cell_init_[iCell];

	
    }
    return 1;
}

//====================================================================================================================
// Main routine to calculate the residual
/*
 *
 *  Format of the residual equations
 *                                                             Unknown                           Index
 * --------------------------------------------------------------------------------------------------------------
 *         Residual (Time)                                     deltaSubcycleCalc_                   0
 *                                                                                            1
 *         Loop over cells                                                            0 <=  iCell < numRCells_
 *                                                                                     j = numEqnsCell_ * iCell
 *                    
 *            Residual (Reference/lattice Position)           rLatticeCBR_final_[iCell];       (1+j)
 *            Residual (Mesh Position)                                                       (j+1) + 1
 *            Residual (Concentration _ k=0)
 *              . . .
 *            Residual (Concentration _ k=Ns-1)
 *  --------------------------------------------------------------------------------------------------------------
 *  **********************************DEAD  BRANCH *****************************************************************
 */
int Electrode_DiffTALE::calcResid_2(double* const resid, const ResidEval_Type_Enum evalType)
{
    int iCell, iPh, jPh, jCell;
    double I_j;
    double I_jp1;
    double Lleft, L3left, rjCBL3;
    double vbarLattice_final_jcell, vbarLattice_init_jcell;
    double rnow3, rnow;
    //double Lright;
    // double rLtarget;

    // Cubes of the cell boundaries, Right and Left, for the initial and final times
    double cbR3_init = 0.0;
    double cbL3_init = 0.0;
    double cbR3_final = 0.0;
    double cbL3_final = 0.0;
    // Reference radius squared, Right and Left, at final time
    //double r0R2_final = 0.0;
    //double r0L2_final = 0.0;
    // Reference radius cubed, Right and Left, at final time

    double r0R3_final = 0.0;
    double r0L3_final = 0.0;
    double rRefCBR3 = 0.0;

    double fluxTC, fluxM;

    // Diffusive flux
    double fluxRCB[10];
   
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


    // Location of the right cell boundary at the beginning of the step
    //  std::vector<doublereal> cellBoundR_init(numRCells_);
    // Velocity of the cell boundary during the time step;
    std::vector<doublereal> cellBoundRVeloc(numRCells_);
    // Node velocity during the time step
    std::vector<doublereal> rnodeVeloc(numRCells_);

    /*
     *    Residual equation for the time step -> Right now we don't have a model
     */
    resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;

    /*
     *  Determine the lattice velocity on the left side domain boundary
     */
    //    vLatticeR is currently set to zero -> this is the appropriate default until something more is needed
    double volFac = 4. * Pi / 3.0 * particleNumberToFollow_;
    /*
     *  Calculate the cell quantitites in a pre-loop
     */
    double cbl, cbr = 0.0;
    for (iCell = 0; iCell < numRCells_; iCell++) {
	cbl = cbr;

        cbr = cellBoundR_init_[iCell];
	cbL3_init = cbl * cbl * cbl;
	cbR3_init = cbr * cbr * cbr;
	double vol_init_ = volFac * (cbR3_init - cbL3_init);
	/*
	 *  Calculate the nodal velocity
	 */
        rnodeVeloc[iCell] = (rnodePos_final_[iCell] - rnodePos_init_[iCell]) / deltaTsubcycleCalc_;
	/*
	 *  Calculate the the velocity of the right cell boundary
	 */
        cellBoundRVeloc[iCell] = (cellBoundR_final_[iCell] - cellBoundR_init_[iCell]) / deltaTsubcycleCalc_;
	/*
	 *  The reference lattice coordinate is rebaselined so that it is linked with the radius coordinate
	 *  at every time step.
	 */
	rLatticeCBR_ref_[iCell] =  cellBoundR_init_[iCell];
	if (iCell == 0) {
	    numLatticeCBR_init_[iCell] = concTot_SPhase_Cell_init_[iCell*numSPhases_] * vol_init_;
	} else {
	    numLatticeCBR_init_[iCell] =  numLatticeCBR_init_[iCell-1] + concTot_SPhase_Cell_init_[iCell*numSPhases_] * vol_init_;
	}
    }

    // ---------------------------  Main Loop Over Cells ----------------------------------------------------

    cbR3_final = 0.0;
    cbR3_init = 0.0;
    r0R3_final = 0.0;
    for (int iCell = 0; iCell < numRCells_; iCell++) {

        /*
         *  Copy right side previous to left side current quantities
         */
        //vLatticeCBL = vLatticeCBR;
        cbL3_final  = cbR3_final;
        cbL3_init   = cbR3_init;
        //r0L2_final  = r0R2_final;
        r0L3_final  = r0R3_final;

        /*
         *  Calculate indexes for accessing the residual
         */
        rindex = 1 + numEqnsCell_ * iCell;
        xindex = 2 + numEqnsCell_ * iCell;
        cindex = 3 + numEqnsCell_ * iCell;
        cIndexPhStart = cindex;

        /*
         *  Temporary pointer variables for total concentrations of each of the phases within the cell
         */
        concTotalVec_SPhase_init  = &(concTot_SPhase_Cell_init_[iCell*numSPhases_]);
        concTotalVec_SPhase_final = &(concTot_SPhase_Cell_final_[iCell*numSPhases_]);

        /*
         * Calculate the cubes of the cell boundary radii
         */
        cbR3_init  = cellBoundR_init_[iCell]  * cellBoundR_init_[iCell]  * cellBoundR_init_[iCell];
        cbR3_final = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
        /*
         * Calculate powers of the lattice coordinate at the right cell boundary
         */
        double r0R_final = rLatticeCBR_final_[iCell];
        r0R3_final = r0R_final * r0R_final * r0R_final;
        /*
         * Calculate the molar volume of the first phase, final and init values
	 *  We will assume for now that this is the lattice molar volume
         */
        double vbarLattice_final = 1.0 / concTotalVec_SPhase_final[0];
        double vbarLattice_init  = 1.0 / concTotalVec_SPhase_init[0];

        /*
         * Residual calculation - Value of the Lattice radius at the right cell boundary
	 *  This is Eqn. 21 in the Memo writeup
	 *      it's the residual equation for determining  rLatticeCBR_final_[iCell]
	 *
	 *      When there is no expansion of the lattice,  rLatticeCBR_final_[iCell] = cellBoundR_final_[iCell]
         */
        double rhs = r0L3_final + vbarLattice_init / vbarLattice_final * (cbR3_final - cbL3_final);
        resid[rindex] = r0R_final - pow(rhs, ONE_THIRD);

	if (formulationType_ > 0) {
	    resid[rindex] = rLatticeCBR_final_[iCell] - rLatticeCBR_init_[iCell];
	}

        /*
         *  Calculate the time derivative of the molar volume
         */
        //double vbarDot = (vbarLattice_final - vbarLattice_init) / deltaTsubcycleCalc_;

        /*
         * Find the lattice velocity of the reference radius at the right cell boundary
         */
        //double vLatticeCBR = -1.0 / (3. * r0R2_final) *
	//		     (3.0 *r0L2_final * vLatticeCBL - MolarVolume_Ref_ / (vbarLattice_final * vbarLattice_final) *
	//		      vbarDot * (cbR3_final - cbL3_final));
	double vLatticeCBR = 0.0;

	//rLtarget =  rLatticeCBR_ref_[iCell];
	if ((rLatticeCBR_final_[iCell] > rLatticeCBR_ref_[iCell]) || (iCell == (numRCells_-1))) {
	    // we are here to when the lattice velocity has put the previous lattice CBR into the current cell
	    // Assign the cell id the current cell.
	    //  (In other words the current cell is expanding
	    jCell = iCell;
	    // The CBL value is the CBR of the next lowest cell.
	    // except for the lowest cell
	    if (jCell == 0) {
		Lleft =  m_rbot0_;
	    } else {
		Lleft =  rLatticeCBR_final_[jCell-1];
	    }
	    //Lright =  rLatticeCBR_final_[jCell];
	    // If the lattice has moved more than one cell, then throw an error condition for now.
	    // We'll come back to fix this later
	    if (Lleft > rLatticeCBR_ref_[iCell]) {
		printf("ERROR: lattice movement not handled\n");
		exit(-1);
	    }
	    // Take the cube of that
	    L3left = Lleft *  Lleft *  Lleft;
	    // Find the cubes of the right CBR and left CBR values
	    //rjCBR3 = cellBoundR_final_[jCell] * cellBoundR_final_[jCell] *  cellBoundR_final_[jCell];
	    rRefCBR3 = rLatticeCBR_ref_[iCell] *rLatticeCBR_ref_[iCell] * rLatticeCBR_ref_[iCell];
	    rjCBL3 = cellBoundL_final_[jCell] * cellBoundL_final_[jCell] *  cellBoundL_final_[jCell];
	    vbarLattice_final_jcell = 1.0 / concTot_SPhase_Cell_final_[iCell*numSPhases_];
	    vbarLattice_init_jcell = 1.0 /  concTot_SPhase_Cell_init_[iCell*numSPhases_];
	    rnow3 =  L3left +  vbarLattice_init_jcell / vbarLattice_final_jcell * (rRefCBR3 - rjCBL3);
	    rnow = pow(rnow3, ONE_THIRD);
	    if (rnow >  cellBoundR_final_[jCell]) {
		if (rnow > cellBoundR_final_[jCell]*(1.0 + 1.0E-14)) {
		    // printf("analyse condition 1 \n");
		}
	    }
	    vLatticeCBR = (rnow - rLatticeCBR_ref_[iCell]) / deltaTsubcycleCalc_;
	} else {
	    // we are here to when the lattice velocity has put the previous lattice CBR into the next cell
	    // Assign the cell id the current cell.
	    jCell = iCell+1;
	    // The CBL value is the CBR of the next lowest cell.
	    Lleft =  rLatticeCBR_final_[jCell-1];
	    //Lright =  rLatticeCBR_final_[jCell];
	    // If the lattice has moved more than one cell, then throw an error condition for now.
	    // We'll come back to fix this later
	    if (Lleft > rLatticeCBR_ref_[iCell]) {
		printf("ERROR: lattice movement not handled\n");
		exit(-1);
	    }
	    // Take the cube of that
	    L3left = Lleft *  Lleft *  Lleft;
	    // Find the cubes of the right CBR and left CBR values
	    // rjCBR3 = cellBoundR_final_[jCell] * cellBoundR_final_[jCell] *  cellBoundR_final_[jCell];
	    rRefCBR3 = rLatticeCBR_ref_[iCell] * rLatticeCBR_ref_[iCell] * rLatticeCBR_ref_[iCell];
	    rjCBL3 = cellBoundL_final_[jCell] * cellBoundL_final_[jCell] * cellBoundL_final_[jCell];
	    vbarLattice_final_jcell = 1.0 / concTot_SPhase_Cell_final_[iCell*numSPhases_];
	    vbarLattice_init_jcell = 1.0 /  concTot_SPhase_Cell_init_[iCell*numSPhases_];
	    rnow3 =  L3left +  vbarLattice_init_jcell / vbarLattice_final_jcell * (rRefCBR3 - rjCBL3);
	    rnow = pow(rnow3, ONE_THIRD);
	    vLatticeCBR = (rnow - rLatticeCBR_ref_[iCell]) / deltaTsubcycleCalc_;
	}
	/*
	 *  Store the calculated lattice velocity
	 */
	vLatticeCBR_cell_[iCell] = vLatticeCBR;

        /*
         * Node position residual - spline equations, it all depends on the top node's equation formulation.
         * Everything else get's dragged along with it
	 * We will calculate the bc when we've calculated the surface reaction rates below.
	 *         -> special case iCell= 0 set to m_rbot0_
	 *         ->              iCell= numRCells_-1 -> distinguishing condition set to boundary condition
         */

	if (iCell == 0) {
	    resid[xindex] = rnodePos_final_[iCell] - m_rbot0_;
	} else if (iCell == numRCells_-1) {
	    // This condition is handled below
	} else {
	    double rtop = rnodePos_final_[numRCells_ - 1];
	    I_j   =  rnodePos_final_[iCell] * rnodePos_final_[iCell] * rnodePos_final_[iCell] / (rtop * rtop * rtop);
	    I_jp1 =  rnodePos_final_[iCell+1] * rnodePos_final_[iCell+1] * rnodePos_final_[iCell+1] / (rtop * rtop * rtop);
	    resid[xindex] = I_jp1 - I_j - fracVolNodePos_[iCell];
	}
	if (formulationType_ > 0) {
	    resid[xindex] = rnodePos_final_[iCell] -  rnodePos_init_[iCell];
	}
    

        /*
         *  Calculate the area of the outer cell which conserves constant functions under mesh movement wrt the Reynolds transport theorum
         */
        double cellBoundR_star2 = (cellBoundR_init_[iCell] * cellBoundR_init_[iCell]
                                   + cellBoundR_final_[iCell] * cellBoundR_init_[iCell] + cellBoundR_final_[iCell] * cellBoundR_final_[iCell])/3.0;
        double areaR_star = 4.0 * Pi * cellBoundR_star2 * particleNumberToFollow_;

   
        int kstart = 0;
	/*
	 * Calculate the molar diffusive flux at the right cell boundary. Note the sum of the molar fluxes over all species
	 * is zero. This is crucial. 
	 */
	if (iCell < (numRCells_-1)) {
	    diffusiveFluxRCB(fluxRCB, iCell, true);
	}

        for (jPh = 0; jPh < numSPhases_; jPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jPh];
            ThermoPhase* th = & thermo(iPh);
            int nSpecies = th->nSpecies();

	    /*
	     *  Residual for the total concentration of the phase in the cell
	     *   -  Skip the first species
	     */
            double old_stuff = concTotalVec_SPhase_init[jPh]  * 4.0 / 3.0 * Pi * (cbR3_init -  cbL3_init)  * particleNumberToFollow_;
            double new_stuff = concTotalVec_SPhase_final[jPh] * 4.0 / 3.0 * Pi * (cbR3_final - cbL3_final) * particleNumberToFollow_;
            /*
             *  Change in the total amount of moles of phase jPh
             */
            resid[cIndexPhStart + kstart] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;

            /*
             * Convective flux - mesh movement with NO material movement
             */
            double vtotalR = cellBoundRVeloc[iCell] - vLatticeCBR;
	    if (formulationType_ > 0) {
		vtotalR = 0.0;
	    }
            if (iCell < (numRCells_- 1)) {
                if (vtotalR >= 0.0) {
                    fluxTC = vtotalR * concTot_SPhase_Cell_final_[numSPhases_*(iCell+1) + jPh];
                } else {
                    fluxTC = vtotalR * concTot_SPhase_Cell_final_[numSPhases_*(iCell)   + jPh];
                }
                resid[cIndexPhStart + kstart]                -= fluxTC * areaR_star;
                resid[cIndexPhStart + kstart + numEqnsCell_] += fluxTC * areaR_star;
            }

	    if (formulationTypeTotalConc_ == 1) {
		resid[cIndexPhStart + kstart] =  concTotalVec_SPhase_final[jPh] - molarDensity_SPhase_Cell_final_[iCell*numSPhases_ + jPh];
	    }

            /*
             *  Residual for the species in the phase
             *   -  Skip the first species
             */
            for (int kSp = 1; kSp < nSpecies; kSp++) {
                int iKRSpecies = kstart + kSp;
                /*
                 * Diffusive flux at the outer boundary for each species
                 */
                if (iCell < (numRCells_ - 1)) {
                    resid[cIndexPhStart + iKRSpecies]                += fluxRCB[iKRSpecies] * areaR_star;
                    resid[cIndexPhStart + iKRSpecies + numEqnsCell_] -= fluxRCB[iKRSpecies] * areaR_star;
                }

                /*
                 * Convective flux - mesh movement with NO material movement
                 */
                if (iCell < (numRCells_- 1)) {
                    if (vtotalR >= 0.0) {
                        fluxM = vtotalR * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell+1) + iKRSpecies];
                        resid[cIndexPhStart + iKRSpecies]                -= fluxM * areaR_star;
                        resid[cIndexPhStart + iKRSpecies + numEqnsCell_] += fluxM * areaR_star;
                    } else {
                        fluxM = vtotalR * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell)   + iKRSpecies];
                        resid[cIndexPhStart + iKRSpecies]                -= fluxM * areaR_star;
                        resid[cIndexPhStart + iKRSpecies + numEqnsCell_] += fluxM * areaR_star;
                    }
                }

                /*
                 *  Change in the total amount of moles of phase jPh
                 */
                old_stuff = concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies]  * 4.0 / 3.0 * Pi * (cbR3_init  - cbL3_init) * particleNumberToFollow_;
                new_stuff = concKRSpecies_Cell_final_[numKRSpecies_ * iCell + iKRSpecies] * 4.0 / 3.0 * Pi * (cbR3_final - cbL3_final)* particleNumberToFollow_;
                resid[cIndexPhStart + iKRSpecies] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;
            }

            // end of phase loop
            kstart += nSpecies;
        }

        if (iCell == numRCells_ - 1) {
            /*
             *  Do the exterior
             */
            double SolidVolCreationRate = 0.0;
	    kstart = 0;
            for (jPh = 0; jPh < numSPhases_; jPh++) {
                iPh = phaseIndeciseKRsolidPhases_[jPh];
                int iStart = getGlobalSpeciesIndex(iPh, 0);
                ThermoPhase* th = & thermo(iPh);
                int nSpecies = th->nSpecies();
                resid[cIndexPhStart + kstart] -= DphMolesSrc_final_[iPh];
                SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart] * DspMoles_final_[iStart];
                for (int kSp = 1; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
                    resid[cIndexPhStart + iKRSpecies] -= DspMoles_final_[iStart + kSp];
                    SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart + kSp] *  DspMoles_final_[iStart + kSp];
                }


		if (formulationTypeTotalConc_ == 1) {
		    resid[cIndexPhStart + kstart] = concTotalVec_SPhase_final[jPh] - molarDensity_SPhase_Cell_final_[iCell*numSPhases_ + jPh];
		}

		kstart += nSpecies;
            }
            /*
             * Calculation of the movement at the top of the
             *    C dXdt * Area - SolidCreationRate = 0
             *
             *  this is an expression of the conservation of total solid molar volume at the
             *  r_exterior.
             */

	    if (formulationType_ > 0) {
		resid[xindex] += 0.0;
	    } else {
		resid[xindex] = rnodeVeloc[iCell] * areaR_star - SolidVolCreationRate;
	    }

        }
    }
    return 1;
}
//=======================================================================================================================
void  Electrode_DiffTALE::showSolution(int indentSpaces)
{
    std::string title = "Lattice Radius CBR (m)";
    vector<std::string> colNames;
    colNames.push_back("LatticeRadius");
    showOneField(title, indentSpaces, &rnodePos_final_[0], numRCells_, &rLatticeCBR_final_[0], colNames, 1);
     

    title = "Phase Concentrations (kmol m-3)";
    showOneField(title, indentSpaces, &rnodePos_final_[0], numRCells_, &concTot_SPhase_Cell_final_[0], KRsolid_phaseNames_,
		 numSPhases_);

    title = "Species Concentrations (kmol /m3)";
    showOneField(title, indentSpaces, &rnodePos_final_[0], numRCells_, &concKRSpecies_Cell_final_[0], KRsolid_speciesNames_,
		 numKRSpecies_);
}
//====================================================================================================================
static void drawline(int sp, int ll)
{
  for (int i = 0; i < sp; i++) {
      printf(" ");
  }
  for (int i = 0; i < ll; i++) {
      printf("-");
  }
  printf("\n");
}
//====================================================================================================================
void Electrode_DiffTALE::showResidual(int indentSpaces, const double * const residual) 
{

    drawline(124, indentSpaces);

    std::string title = "Residual for Lattice Radius CBR (m)";

    const double * const res = residual + 1;
    int iTerm = 0;

    showOneResid(title, indentSpaces, &cellBoundR_final_[0], numRCells_, 1, 0, &rLatticeCBR_init_[0],
		 &rLatticeCBR_final_[0], numEqnsCell_, iTerm, &errorLocalNLS_[1], &atolResidNLS_[1],
		 res);

    title = "Residual for position (m)";
    iTerm = 1;
    showOneResid(title, indentSpaces, &rnodePos_final_[0], numRCells_, 1, 0, &rnodePos_init_[0],
		 &rnodePos_final_[0], numEqnsCell_, iTerm, &errorLocalNLS_[1], &atolResidNLS_[1],
		 res);

    int kstart = 0;
    for (int jPh = 0; jPh < numSPhases_; jPh++) {
	int iPh = phaseIndeciseKRsolidPhases_[jPh];
	//int iStart = getGlobalSpeciesIndex(iPh, 0);
	ThermoPhase* th = & thermo(iPh);
	int nSpecies = th->nSpecies();

	title = "Residual for total phase concentration " + th->name();
	iTerm = 2 + kstart;
	showOneResid(title, indentSpaces, &rnodePos_final_[0], numRCells_, 1, 0, &concTot_SPhase_Cell_init_[0],
		     &concTot_SPhase_Cell_final_[0], numEqnsCell_, iTerm, &errorLocalNLS_[1], &atolResidNLS_[1],
		     res);


	for (int kSp = 1; kSp < nSpecies; kSp++) {
	    int iKRSpecies = kstart + kSp;

	    title = "Residual for total species concentration of " + th->speciesName(kSp);
	    iTerm = 2 + iKRSpecies;
	    showOneResid(title, indentSpaces, &rnodePos_final_[0], numRCells_, numKRSpecies_, iKRSpecies, &concKRSpecies_Cell_init_[0],
			 &concKRSpecies_Cell_final_[0], numEqnsCell_, iTerm, &errorLocalNLS_[1], &atolResidNLS_[1],
			 res);

	}
	kstart += nSpecies;
    }


    drawline(124, indentSpaces);

}

//====================================================================================================================
void  Electrode_DiffTALE::showOneField(const std::string &title, int indentSpaces,
					 const double * const radialValues, int numRadialVals, 
					 const double * const vals, const std::vector<std::string> &varNames, int numFields)
{
    int n, iCell;
    double v;
    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
	indent += " ";
    }
    int numBlockRows = numFields / 5;
    if (title.size() > 0) {
	drawline(indentSpaces, 80);
	printf("%s  %s\n", indent.c_str(), title.c_str());
    }
    for (int iBlock = 0; iBlock < numBlockRows; iBlock++) {
	drawline(indentSpaces, 80);
	printf("%s        z   ", indent.c_str());
	for (n = 0; n < 5; n++) {
	    int ivar = iBlock * 5 + n;
	    string name = varNames[ivar];
	    printf(" %15s", name.c_str());
	}
	printf("\n");
	drawline(indentSpaces, 80);

	for (iCell = 0; iCell < numRadialVals; iCell++) {
	    doublereal r = radialValues[iCell];
	    printf("\n%s    %- 10.4E ", indent.c_str(), r);
	    int istart = iCell * numFields;
	    for (n = 0; n < 5; n++) {
		v = vals[istart + iBlock * 5 + n];
		printf(" %- 10.4E ", v);
	    }
	}
	printf("\n");
    }
    int nrem = numFields - 5 * numBlockRows;
    if (nrem > 0) {
	drawline(indentSpaces, 80);
	printf("%s        z   ", indent.c_str());
	for (n = 0; n < nrem; n++) {
	    int ivar = numBlockRows * 5 + n;
	    string name = varNames[ivar];
	    printf(" %15s", name.c_str());
	}
	printf("\n");
	drawline(indentSpaces, 80);

	for (iCell = 0; iCell < numRadialVals; iCell++) {
	    doublereal r = radialValues[iCell];
	    printf("\n%s    %- 10.4E ", indent.c_str(), r);
	    int istart = iCell * numFields;
	    for (n = 0; n < nrem; n++) {
		v = vals[istart + numBlockRows * 5 + n];
		printf(" %- 10.4E ", v);
	    }
	}
	printf("\n");
    }
}
//====================================================================================================================
void  Electrode_DiffTALE::showOneFieldInitFinal(const std::string &title, int indentSpaces, 
						  const double * const radialValues, int numRadialVals, 
						  const double * const vals_init,  const double * const vals_final,
						  const std::vector<std::string> &varNames, int numFields)
{
    int n, iCell;
    double v_init, v_final;
    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
	indent += " ";
    }
    int numBlockRows = numFields / 4;
    if (title.size() > 0) {
	drawline(indentSpaces, 80);
	printf("%s  %s\n", indent.c_str(), title.c_str());
    }
    for (int iBlock = 0; iBlock < numBlockRows; iBlock++) {
	drawline(indentSpaces, 80);
	printf("%s        z      | ", indent.c_str());
	for (n = 0; n < 4; n++) {
	    int ivar = iBlock * 4 + n;
	    string name = varNames[ivar];
	    printf("(f) %-13.13s (i) | ", name.c_str());
	}
	printf("\n");
	drawline(indentSpaces, 80);

	for (iCell = 0; iCell < numRadialVals; iCell++) {
	    doublereal r = radialValues[iCell];
	    printf("%s    %- 10.4E |", indent.c_str(), r);
	    int istart = iCell * numFields;
	    for (n = 0; n < 4; n++) {
		v_init = vals_init[istart + iBlock * 4 + n];
		v_final = vals_final[istart + iBlock * 4 + n];
		printf(" %- 10.4E %-10.4E |", v_final, v_init);
	    }
	    printf("\n");
	}
    }
    int nrem = numFields - 4 * numBlockRows;
    if (nrem > 0) {
	drawline(indentSpaces, 80);
	printf("%s        z      | ", indent.c_str());
	for (n = 0; n < nrem; n++) {
	    int ivar = numBlockRows * 4 + n;
	    string name = varNames[ivar];
	    printf("(f) %-13.13s (i) | ", name.c_str());
	}
	printf("\n");
	drawline(indentSpaces, 80);

	for (iCell = 0; iCell < numRadialVals; iCell++) {
	    doublereal r = radialValues[iCell];
	    printf("%s    %- 10.4E |", indent.c_str(), r);
	    int istart = iCell * numFields;
	    for (n = 0; n < nrem; n++) {
		v_init = vals_init[istart + numBlockRows * 4 + n];
		v_final = vals_final[istart + numBlockRows * 4 + n];
		printf(" %- 10.4E %- 10.4E |", v_final, v_init);
	    }
	    printf("\n");
	}
    }
    printf("\n");
}

//====================================================================================================================
void  Electrode_DiffTALE::showOneResid(const std::string &title, int indentSpaces,
					 const double * const radialValues, int numRadialVals, 
					 int numFields, int iField, const double * const val_init,  
					 const double * const val_final,
					 int numEqnsCell, int iEqn, const double * const resid_error, 
					 const double * const solnError_tol, 
					 const double * const residual)
{
    int iCell;
    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
	indent += " ";
    }
  
    drawline(indentSpaces, 80);
    printf("%s  %s\n", indent.c_str(), title.c_str());
 
    drawline(indentSpaces, 80);
    printf("%s        z      ", indent.c_str());
    
    printf(" %-11s", "Init_Value");
    printf(" %-11s", "Final_Value");
    if (resid_error) {
	printf(" %-11s", "Error");
    }
    if (solnError_tol) {
	printf(" %-11s |  ", "Toler");
    }
    printf(" %-11s\n", "Resid");
    
    drawline(indentSpaces, 80);

    for (iCell = 0; iCell < numRadialVals; iCell++) {
	doublereal r = radialValues[iCell];
	printf("%s    %- 10.4E ", indent.c_str(), r);
	int istart = iCell * numFields + iField;
	printf(" %- 10.4E ", val_init[istart]);
	printf(" %- 10.4E ", val_final[istart]);
	istart = iCell *  numEqnsCell + iEqn;
	if (resid_error) {
	    printf(" %- 10.4E ", resid_error[istart]);
	}
	if (solnError_tol) {
	    printf(" %- 10.4E ", solnError_tol[istart]);
	}
	printf(" | ");
	printf(" %- 10.4E ", residual[istart]);
	printf("\n");
    }
    drawline(indentSpaces, 80);
}
//====================================================================================================================
//  Extract the ROP of the reaction fronts from Cantera within this routine
/*
 *  In this routine we calculate the rates of progress of reactions and species on all active reacting surfaces.
 *
 *        ROP_[jRxn]
 *        spNetProdPerArea_List_[isk][kIndexKin]
 */
void  Electrode_DiffTALE::extractInfo()
{
    //double mfSig = 0.0;
    /*
     *  this is an initial copy from Electrode_CSTR, since it should be nearly the same
     */
    int maxNumRxns = RSD_List_[0]->nReactions();
    std::vector<double> netROP(maxNumRxns, 0.0);

    /*
     *  Load the conditions of the last cell into the ThermoPhase object
     */
    int iCell = numRCells_ - 1;
    int kspCell = iCell * numKRSpecies_;
    for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	ThermoPhase* tp =  thermoSPhase_List_[jRPh];
	double* spMf_ptr =  &(spMf_KRSpecies_Cell_final_[kspCell]);
	tp->setState_TPX(temperature_, pressure_, spMf_ptr);
	int nsp = tp->nSpecies();
	kspCell += nsp;
    }

    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    for (int isk = 0; isk < numSurfaces_; isk++) {
        double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
        std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.0);
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
                 *  Get the species production rates for the reacting surface, spNetProdPerArea
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

		//printf(" extractInfo : MoleFraction = %20.10E  ROP = %20.10E \n", mfSig, ROP_[0]);
		       
            }
        }
    }
}
//====================================================================================================================
//! Collect mole change information
/*!
 *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
 *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
 *   all species in the electrode due to surface reactions.
 *
 *   We use sa_star calculation that combines r_final and r_init to make a conserved system.
 *
 *   (inherited from Electrode_Integrator)
 */
void Electrode_DiffTALE::updateSpeciesMoleChangeFinal()
{
    double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(0);
    std::fill(DspMoles_final_.begin(), DspMoles_final_.end(), 0.0);

    double r_init  = Radius_exterior_init_;
    double r_final = Radius_exterior_final_;

    double surfaceArea_star =  4. * Pi / 3. * (r_init * r_init + r_init * r_final + r_final * r_final);
    double mult = surfaceArea_star * particleNumberToFollow_;
    for (int i = 0; i < m_totNumVolSpecies; i++) {
        DspMoles_final_[i] += mult * spNetProdPerArea[i];
    }
}

//====================================================================================================================
// Pack the nonlinear solver proplem
/*
 *  formulate the nonlinear solver problem to be solved.
 *     Fields to be filled in
 *             yvalNLS_
 *             ylowNLS_
 *             yhighNLS_
 *             deltaBoundsMagnitudesNLS_
 *  (virtual from Electrode_Integrator)
 *
 * --------------------------------------------------------------------------------------------------------------
 *         Residual (Time)                                     deltaSubcycleCalc_                   0
 *                                                                                            1
 *         Loop over cells                                                            0 <=  iCell < numRCells_
 *                                                                                     j = numEqnsCell_ * iCell
 *                    
 *            Residual (Reference/lattice Position)          rLatticeCBR_final_[iCell];      (1+j)
 *            Residual (Mesh Position)                                                       (j+1) + 1
 *          Loop over distributed Phases
 *            Residual (Concentration _ k=0)                  concTot_SPhase_Cell_final_[iCell * numSPhase_ + jPh]
 *              . . .
 *            Residual (Concentration _ k=Ns-1)               concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies]
 *  --------------------------------------------------------------------------------------------------------------
 */
void Electrode_DiffTALE::initialPackSolver_nonlinFunction()
{
    int index = 0;
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
     *  electrons. So, we will set atolBaseResid_ to 1.0E-25 with possible additions later.
     */
    double solidMoles = SolidTotalMoles();
    double atolVal = solidMoles * atolBaseResid_;



    int kstart, jRPh, iKRSpecies;
  

    for (int iCell = 0; iCell < numRCells_; iCell++) {
	kstart = 0;

	yvalNLS_[index] = rLatticeCBR_final_[iCell];
	ylowNLS_[index] = -1.0E300;
	yhighNLS_[index] = 1.0E300;
	deltaBoundsMagnitudesNLS_[index] = 0.1 * Radius_exterior_init_;
	index++;

	yvalNLS_[index] = rnodePos_final_[iCell];
	ylowNLS_[index] = -1.0E300;
	yhighNLS_[index] = 1.0E300;
	deltaBoundsMagnitudesNLS_[index] = 0.1 * Radius_exterior_init_;
	index++;

	for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    int nsp = numSpeciesInKRSolidPhases_[jRPh];
	    yvalNLS_[index] = concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh];
             /*
	      * If we are birthing a phase, we decide to do this in the predictor. We search for
	      * a solution where this is the case. If this doesn't happen, we produce a hard error.
	      *  if (justBornPhase_[iph]) {
	      *     ylowNLS_[index] = atolVal;
	      *  } else {
	      *     ylowNLS_[index] = (-0.0001 * solidMoles);
	      *  }
	      *
	      */
	    // we allow the total phase moles to go negative, initially.  For now, we won't worry about
	    // phase fronts within the radial direction.
	    ylowNLS_[index] = (-0.0001 * yvalNLS_[index]);
	    yhighNLS_[index] = 1.0E3 * yvalNLS_[index];
 
	  
	    for (int kSp = 1; kSp < nsp; kSp++) {
		iKRSpecies = kstart + kSp;
		yvalNLS_[index + kSp] = concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies];
		ylowNLS_[index + kSp]  = 0.0;
                yhighNLS_[index + kSp] = 1.0E3 * yvalNLS_[index];
                deltaBoundsMagnitudesNLS_[index + kSp] = 1000. * atolVal;
	    }
	    kstart += nsp;
	    index += nsp;
	}
    }
}
//====================================================================================================================
// Set the internal initial intermediate and initial global state from the internal final state
/*
 *  (virtual function from Electrode.h)
 *
 *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 * @param setInitInit   Boolean indicating whether you should set the init_init state as well
 */
void  Electrode_DiffTALE::setInitStateFromFinal(bool setInitInit)
{
    /*
     * Call the parent object
     */
    Electrode_Integrator::setInitStateFromFinal(setInitInit);

    int iCell, i;
 
    int ntotal = numRCells_ * numKRSpecies_;
    for (i = 0; i < ntotal; ++i) {
	spMoles_KRsolid_Cell_init_[i] = spMoles_KRsolid_Cell_final_[i];
	concKRSpecies_Cell_init_[i]   = concKRSpecies_Cell_final_[i];
	spMf_KRSpecies_Cell_init_[i]  = spMf_KRSpecies_Cell_final_[i];
	actCoeff_Cell_init_[i]        = actCoeff_Cell_final_[i];
    }

    int iTotal =  numSPhases_ * numRCells_;
    for (i = 0; i < iTotal; ++i) {
	concTot_SPhase_Cell_init_[i] = concTot_SPhase_Cell_final_[i];
    }

    for (iCell = 0; iCell < numRCells_; ++iCell) {
	rLatticeCBR_init_[iCell] = rLatticeCBR_final_[iCell];
	rnodePos_init_[iCell]    = rnodePos_final_[iCell];
	cellBoundR_init_[iCell]  = cellBoundR_final_[iCell];
	cellBoundL_init_[iCell]  = cellBoundL_final_[iCell];
    }

    /*
     *  Now transfer that to other states
     */
    //setFinalFinalStateFromFinal();

    if (setInitInit) {
	for (i = 0; i < ntotal; ++i) {
	    spMoles_KRsolid_Cell_init_init_[i] = spMoles_KRsolid_Cell_final_[i];
	    concKRSpecies_Cell_init_init_[i]   = concKRSpecies_Cell_final_[i];
	}
	
	for (i = 0; i < iTotal; ++i) {
	    concTot_SPhase_Cell_init_init_[i] = concTot_SPhase_Cell_final_[i];
	}

	for (iCell = 0; iCell < numRCells_; ++iCell) {
	    rLatticeCBR_init_init_[iCell] = rLatticeCBR_final_[iCell];
	    rnodePos_init_init_[iCell]    = rnodePos_final_[iCell];
	}
    }
}
//====================================================================================================================
// Set the internal initial intermediate and initial global state from the internal final_final state
/*
 *  (virtual function from Electrode.h)
 *
 *  Set the initial  and init_init state from the final state.
 */
void Electrode_DiffTALE::setInitInitStateFromFinalFinal()
{
    Electrode_Integrator::setInitInitStateFromFinalFinal();

    int iCell, i;
    int ntotal = numRCells_ * numKRSpecies_;
    for (i = 0; i < ntotal; ++i) {
        spMoles_KRsolid_Cell_init_init_[i] = spMoles_KRsolid_Cell_final_final_[i];
	spMoles_KRsolid_Cell_init_[i]      = spMoles_KRsolid_Cell_final_final_[i];
        concKRSpecies_Cell_init_init_[i]   = concKRSpecies_Cell_final_final_[i];
	concKRSpecies_Cell_init_[i]        = concKRSpecies_Cell_final_final_[i];
	spMf_KRSpecies_Cell_init_[i]       = spMf_KRSpecies_Cell_final_[i];
	actCoeff_Cell_init_[i]             = actCoeff_Cell_final_[i];
    }

    int iTotal =  numSPhases_ * numRCells_;
    for (i = 0; i < iTotal; ++i) {
        concTot_SPhase_Cell_init_init_[i] = concTot_SPhase_Cell_final_final_[i];
        concTot_SPhase_Cell_init_[i]      = concTot_SPhase_Cell_final_[i];
    }

    for (iCell = 0; iCell < numRCells_; ++iCell) {
        rLatticeCBR_init_init_[iCell] = rLatticeCBR_final_final_[iCell];
	rLatticeCBR_init_[iCell]      = rLatticeCBR_final_final_[iCell];
        rnodePos_init_init_[iCell]    = rnodePos_final_final_[iCell];
        rnodePos_init_[iCell]         = rnodePos_final_final_[iCell];
	cellBoundR_init_[iCell]       = cellBoundR_final_[iCell];
    }

    onRegionBoundary_init_init_ =  onRegionBoundary_final_final_;
    onRegionBoundary_init_      =  onRegionBoundary_final_final_;
}
//====================================================================================================================
void Electrode_DiffTALE::setFinalStateFromInit()
{   
    /*
     * Call the parent object
     */
    Electrode_Integrator::setFinalStateFromInit();
    
    int iCell, i;
    int iTotal =  numSPhases_ * numRCells_;
    int ntotal = numRCells_ * numKRSpecies_;
    for (i = 0; i < ntotal; ++i) {
	spMoles_KRsolid_Cell_final_[i] = spMoles_KRsolid_Cell_init_[i];
	concKRSpecies_Cell_final_[i]   = concKRSpecies_Cell_init_[i];
	actCoeff_Cell_final_[i]        = actCoeff_Cell_init_[i];
	spMf_KRSpecies_Cell_final_[i]  = spMf_KRSpecies_Cell_init_[i];
    }
    
    for (i = 0; i < iTotal; ++i) {
	concTot_SPhase_Cell_final_[i] = concTot_SPhase_Cell_init_[i];
    }
    
    for (iCell = 0; iCell < numRCells_; ++iCell) {
	rLatticeCBR_final_[iCell] = rLatticeCBR_init_[iCell];
	rnodePos_final_[iCell]    = rnodePos_init_[iCell];
	cellBoundR_final_[iCell]  = cellBoundR_init_[iCell];
	cellBoundL_final_[iCell]  = cellBoundL_init_[iCell];
    }
    onRegionBoundary_final_ =  onRegionBoundary_init_;
}
//====================================================================================================================
//    Set the internal initial intermediate from the internal initial global state
/*
 *  Set the intial state from the init init state. We also can set the final state from this
 *  routine as well.
 *
 *  The final_final is not touched.
 *
 * @param setFinal   Boolean indicating whether you should set the final as well
 */
void Electrode_DiffTALE::setInitStateFromInitInit(bool setFinal)
{
    /*
     * Call the parent object
     */
    Electrode_Integrator::setInitStateFromInitInit(setFinal);

    int iCell, i; 
    int ntotal = numRCells_ * numKRSpecies_;
    for (i = 0; i < ntotal; ++i) {
	spMoles_KRsolid_Cell_init_[i] = spMoles_KRsolid_Cell_init_init_[i];
	concKRSpecies_Cell_init_[i] =  concKRSpecies_Cell_init_init_[i];
    }
    int iTotal =  numSPhases_ * numRCells_;
    for (i = 0; i < iTotal; ++i) {
	concTot_SPhase_Cell_init_[i] = concTot_SPhase_Cell_init_init_[i];
    }
    for (iCell = 0; iCell < numRCells_; ++iCell) {
	rLatticeCBR_init_[iCell] = rLatticeCBR_init_init_[iCell];
	rnodePos_init_[iCell] = rnodePos_init_init_[iCell];
    }

    if (setFinal) {
	for (i = 0; i < ntotal; ++i) {
	    spMoles_KRsolid_Cell_final_[i] = spMoles_KRsolid_Cell_init_init_[i];
	    concKRSpecies_Cell_final_[i] =  concKRSpecies_Cell_init_init_[i];
	}
	for (i = 0; i < iTotal; ++i) {
	    concTot_SPhase_Cell_final_[i] = concTot_SPhase_Cell_init_init_[i];
	}
	for (iCell = 0; iCell < numRCells_; ++iCell) {
	    rLatticeCBR_final_[iCell] = rLatticeCBR_init_init_[iCell];
	    rnodePos_final_[iCell] = rnodePos_init_init_[iCell];
	}
    }
}
//====================================================================================================================
// Set the internal final global state from the internal final intermediate state
/*
 *  (virtual function from Electrode.h)
 *
 *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
 */
void Electrode_DiffTALE::setFinalFinalStateFromFinal() 
{
    /*
     * Call the parent object
     */
    Electrode_Integrator::setFinalFinalStateFromFinal();

    int iCell, i;
 
    int ntotal = numRCells_ * numKRSpecies_;
    for (i = 0; i < ntotal; ++i) {
	spMoles_KRsolid_Cell_final_final_[i] = spMoles_KRsolid_Cell_final_[i];
	concKRSpecies_Cell_final_final_[i] =  concKRSpecies_Cell_final_[i];
    }

    int iTotal =  numSPhases_ * numRCells_;
    for (i = 0; i < iTotal; ++i) {
	concTot_SPhase_Cell_final_final_[i] = concTot_SPhase_Cell_final_[i];
    }

    for (iCell = 0; iCell < numRCells_; ++iCell) {
	rLatticeCBR_final_final_[iCell] = rLatticeCBR_final_[iCell];
	rnodePos_final_final_[iCell] = rnodePos_final_[iCell];
    }
}
//====================================================================================================================
// Returns the equilibrium OCV for the selected ReactingSurfaceDomain and current conditions (virtual)
/*
 *  This routine uses a root finder to find the voltage at which there
 *  is zero net electron production.  It leaves the object unchanged. However, it
 *  does change the voltage of the phases during the calculation, so this is a non const function.
 *
 * @param isk  Reacting surface domain id
 */
 double  Electrode_DiffTALE::openCircuitVoltage(int isk)
 {
     /*
      *  Load the conditions of the last cell into the ThermoPhase object
      */
     int iCell = numRCells_ - 1;
     int kspCell = iCell * numKRSpecies_;
     for (int jRPh = 0; jRPh < numSPhases_; jRPh++) {
	ThermoPhase* tp =  thermoSPhase_List_[jRPh];
	double* spMf_ptr =  &(spMf_KRSpecies_Cell_final_[kspCell]);
	tp->setState_TPX(temperature_, pressure_, spMf_ptr);
	int nsp = tp->nSpecies();
	kspCell += nsp;
     }
     
     double val = Electrode::openCircuitVoltage(isk);
     return val;
 }
//====================================================================================================================
void Electrode_DiffTALE::printElectrode(int pSrc, bool subTimeStep)
{
    int iph;
    vector<std::string> colNames;
    double* netROP = new double[m_NumTotSpecies];
    double egv = TotalVol();
    double tsm = SolidTotalMoles();
    printf("   ==============================================================================================\n");
    if (subTimeStep) {
        printf("      Electrode_DiffTALE at intermediate-step time final = %12.5E\n", tfinal_);
        printf("                              intermediate-step time init  = %12.5E\n", tinit_);
	printf("                ChemModel Type = %3d , DomainNumber = %2d , CellNumber = %2d , SubIntegrationCounter = %d\n",
               electrodeChemistryModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberSubIntegrations_);
    } else {
        printf("      Electrode_DiffTALE at time final = %12.5E\n", t_final_final_);
        printf("                              time init  = %12.5E\n", t_init_init_);
	printf("                Model Type = %3d , DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n",
               electrodeChemistryModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberIntegrations_);
    }
    printf("   ==============================================================================================\n");
    printf("          Voltage = %g volts\n", deltaVoltage_);
    if (subTimeStep) {
        printf("          Current = %g amps\n", integratedLocalCurrent());
    } else {
        printf("          Current = %g amps\n", integratedCurrent());
    }

    printf("          Number of external surfaces = %d\n", numExternalInterfacialSurfaces_);
    printf("          Solid Volume = %10.3E\n", ElectrodeSolidVolume_);
    printf("          Total Volume = %10.3E\n", egv);
    if (egv > 0.0) {
        printf("          Porosity     = %10.3E\n", (egv - ElectrodeSolidVolume_) / egv);
    } else {
        printf("          Porosity     =       NA\n");
    }
    printf("          Temperature = %g\n", temperature_);
    printf("          Pressure = %g\n", pressure_);
    printf("          Total Solid Moles = %11.4E kmol\n", tsm);
    if (tsm > 0.0) {
        printf("          Molar Volume of Solid = %11.4E cm3 gmol-1\n", ElectrodeSolidVolume_ / tsm * 1.0E3);
    } else {
    }
    printf("          Particle Number to Follow = %11.4E\n", particleNumberToFollow_);
    printf("          followElectrolyteMoles = %d\n", followElectrolyteMoles_);
    printf("          ElectrolytePseudoMoles = %g\n",  electrolytePseudoMoles_);

    printElectrodeCapacityInfo(pSrc, subTimeStep);

    colNames.push_back("LatticeRadius");
    std::string title = "";
    int indentSpaces = 10;
    showOneField(title, indentSpaces, &rnodePos_final_[0], numRCells_, &rLatticeCBR_final_[0], colNames, 1);

    for (iph = 0; iph < m_NumTotPhases; iph++) {
        printElectrodePhase(iph, pSrc, subTimeStep);
    }
    printf("   ==============================================================================================\n");
    delete [] netROP;
}
//===================================================================================================================

void Electrode_DiffTALE::printElectrodePhase(int iph, int pSrc, bool subTimeStep)
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
        printf("                Electric Potential = %g\n", tp.electricPotential());
    }
    if (iph >= NumVolPhases_) {
        isph = iph - NumVolPhases_;
        printf("                surface area (final) = %g\n",  surfaceAreaRS_final_[isph]);
        printf("                surface area (init)  = %g\n",  surfaceAreaRS_init_[isph]);
        int ddd =  isExternalSurface_[isph];
        printf("                IsExternalSurface = %d\n", ddd);
        double oc = openCircuitVoltage(isph);
        if (oc != 0.0) {
            printf("                Open Circuit Voltage = %g\n", oc);
        }
    }
    printf("\n");
    printf("                Name               MoleFrac_final MoleFrac_init kMoles_final kMoles_init SrcTermLastStep(kMoles)\n");
    for (int k = 0; k < nsp; k++) {
        string sname = tp.speciesName(k);
        if (pSrc) {
            if (subTimeStep) {
                printf("                %-22s %10.3E %10.3E %10.3E %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k], spMf_init_[istart + k],
                       spMoles_final_[istart + k], spMoles_init_[istart + k],
                       spMoleIntegratedSourceTermLast_[istart + k]);
            } else {
                printf("                %-22s %10.3E %10.3E %10.3E %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k], spMf_init_[istart + k],
                       spMoles_final_[istart + k], spMoles_init_init_[istart + k],
                       spMoleIntegratedSourceTerm_[istart + k]);
            }
        } else {
            if (subTimeStep) {
                printf("                %-22s %10.3E %10.3E %10.3E %10.3E\n", sname.c_str(), spMf_final_[istart + k], spMf_init_[istart + k],
                       spMoles_final_[istart + k],   spMoles_init_[istart + k]);
            } else {
                printf("                %-22s %10.3E %10.3E %10.3E %10.3E\n", sname.c_str(), spMf_final_[istart + k], spMf_init_[istart + k],
                       spMoles_final_[istart + k],   spMoles_init_init_[istart + k]);
            }
        }
    }
    if (iph >= NumVolPhases_) {
        const vector<double>& rsSpeciesProductionRates = RSD_List_[isph]->calcNetSurfaceProductionRateDensities();
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
    /*
     * Add distributed printouts.
     */
    if (distribPhIndexKRsolidPhases_[iph] >= 0) {
	int jPh = distribPhIndexKRsolidPhases_[iph];
	std::vector<double> concKRSpecies_iph_init(numRCells_ * nsp);
	std::vector<double> concKRSpecies_iph_final(numRCells_ * nsp);
	std::vector<double> mf_iph_init(numRCells_ * nsp);
	std::vector<double> mf_iph_final(numRCells_ * nsp);
	std::vector<double> spMoles_iph_init(numRCells_ * nsp);
	std::vector<double> spMoles_iph_final(numRCells_ * nsp);
	std::vector<double> concTot_iph_final(numRCells_);
	std::vector<double> concTot_iph_init(numRCells_);
	for (int iCell = 0; iCell < numRCells_; iCell++) {
	    int istart = iCell * nsp;
	    int jstart = iCell * numKRSpecies_;
	    concTot_iph_final[iCell] = concTot_SPhase_Cell_final_[iCell * numSPhases_ + jPh];
	    concTot_iph_init[iCell] = concTot_SPhase_Cell_init_[iCell * numSPhases_ + jPh];

	    for (int kSp = 0; kSp < nsp; kSp++) {
		int iKRSpecies = kstartKRSolidPhases_[jPh] + kSp;
		concKRSpecies_iph_init[istart + kSp] = concKRSpecies_Cell_init_[jstart + iKRSpecies];
		concKRSpecies_iph_final[istart + kSp] = concKRSpecies_Cell_final_[jstart + iKRSpecies];
		mf_iph_init[istart + kSp]  = spMf_KRSpecies_Cell_init_[jstart + iKRSpecies];
		mf_iph_final[istart + kSp] = spMf_KRSpecies_Cell_final_[jstart + iKRSpecies];

		spMoles_iph_init[istart + kSp]  = spMoles_KRsolid_Cell_init_[jstart + iKRSpecies];
	        spMoles_iph_final[istart + kSp] = spMoles_KRsolid_Cell_final_[jstart + iKRSpecies];
	    }
	}
	std::vector<string> speciesNames;
	for (int kSp = 0; kSp < nsp; kSp++) {
	    speciesNames.push_back( tp.speciesName(kSp));
	}
        string pName = tp.name();
	std::vector<string> pNames;
	pNames.push_back(pName);

	string title = "    Species Cell Moles (final and init)";
	showOneFieldInitFinal(title, 14, &rnodePos_final_[0], numRCells_, &spMoles_iph_init[0], &spMoles_iph_final[0],
			      speciesNames, nsp);

	title = "    Phase Concentration (kmol /m3) (final and init)";
	showOneFieldInitFinal(title, 14, &rnodePos_final_[0], numRCells_, &concTot_iph_init[0], &concTot_iph_final[0],
			      pNames, 1);


	title = "   Species Concentrations (kmol /m3) (final and init)";

	showOneFieldInitFinal(title, 14, &rnodePos_final_[0], numRCells_, &concKRSpecies_iph_init[0], &concKRSpecies_iph_final[0],
			      speciesNames, nsp);

        title = "   Mole Fractions (final and init)";

	showOneFieldInitFinal(title, 14, &rnodePos_final_[0], numRCells_, &mf_iph_init[0], &mf_iph_final[0],
			      speciesNames, nsp);


    }


    delete [] netROP;

}


} // End of namespace Cantera
//======================================================================================================================
