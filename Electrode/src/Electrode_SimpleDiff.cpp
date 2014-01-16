/*
 * $Id: Electrode_SimpleDiff.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"

#include "Electrode_SimpleDiff.h"
#include "cantera/integrators.h"
#include "Electrode_RadialDiffRegions.h"

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
//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_SimpleDiff::Electrode_SimpleDiff() :
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
    cellBoundL_final_(0),

    volPP_Cell_final_(0),
    fracVolNodePos_(0),
    partialMolarVolKRSpecies_Cell_final_(0),
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
    goNowhere_(0)
{


}

//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_SimpleDiff::Electrode_SimpleDiff(const Electrode_SimpleDiff& right) :
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
    cellBoundL_final_(0),

    volPP_Cell_final_(0),
    fracVolNodePos_(0),
    partialMolarVolKRSpecies_Cell_final_(0),
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
    goNowhere_(0)  
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
Electrode_SimpleDiff&
Electrode_SimpleDiff::operator=(const Electrode_SimpleDiff& right)
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
    cellBoundL_final_                   = right.cellBoundL_final_;
    volPP_Cell_final_                   = right.volPP_Cell_final_;

 
    fracVolNodePos_                     = right.fracVolNodePos_;
    partialMolarVolKRSpecies_Cell_final_= right.partialMolarVolKRSpecies_Cell_final_;
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

    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
//   destructor
Electrode_SimpleDiff::~Electrode_SimpleDiff()
{
}
//======================================================================================================================
//    Return the type of electrode
/*
 *  Returns the enum type of the electrode. This is used in the factory routine.
 *
 *  @return Returns an enum type, called   Electrode_Types_Enum
 */
Electrode_Types_Enum Electrode_SimpleDiff::electrodeType() const
{
    return SIMPLE_DIFF_ET;
}
//====================================================================================================================
int Electrode_SimpleDiff::electrode_input_child(ELECTRODE_KEY_INPUT** ei_ptr)
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
Electrode_SimpleDiff::electrode_model_create(ELECTRODE_KEY_INPUT* eibase)
{
    int iPh, jPh;
    /*
     *  Downcast the Key input to make sure we are being fed the correct child object
     */
    ELECTRODE_RadialDiffRegions_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialDiffRegions_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_SimpleDiff::electrode_model_create()",
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
    RadialDiffRegionSpec& r0 = ei->rregions_[0];
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
    Radius_exterior_final_ = pow(volPerParticle * 3 / (4. * Pi), 0.333333333333333333333);

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
  

    /*
     *  Initialize all of the variables on the domain
     *    We take the grid size and the mole numbers from the base Electrode representation
     *    and create a distributed representation.
     *    Any disparity in mole numbers creates an error.
     */
    initializeAsEvenDistribution();

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
int Electrode_SimpleDiff::setInitialConditions(ELECTRODE_KEY_INPUT* eibase)
{
    ELECTRODE_RadialDiffRegions_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialDiffRegions_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_SimpleDiff::electrode_model_create()",
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
void
Electrode_SimpleDiff::init_sizes()
{
    int kspCell =  numKRSpecies_ *  numRCells_;
    int nPhCell = numSPhases_ * numRCells_;

    spMoles_KRsolid_Cell_final_.resize(kspCell, 0.0);
    spMoles_KRsolid_Cell_init_.resize(kspCell, 0.0);
    spMoles_KRsolid_Cell_final_final_.resize(kspCell, 0.0);
    spMoles_KRsolid_Cell_init_init_.resize(kspCell, 0.0);

    concTot_SPhase_Cell_final_final_.resize(nPhCell, 0.0);
    concTot_SPhase_Cell_final_.resize(nPhCell, 0.0);
    concTot_SPhase_Cell_init_.resize(nPhCell, 0.0);
    concTot_SPhase_Cell_init_init_.resize(nPhCell, 0.0);

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
    cellBoundL_final_.resize(numRCells_, 0.0);

    volPP_Cell_final_.resize(numRCells_, 0.0);
    //  fracNodePos_.resize(numRCells_, 0.0);
    fracVolNodePos_.resize(numRCells_, 0.0);

    partialMolarVolKRSpecies_Cell_final_.resize(kspCell, 0.0);

    DspMoles_final_.resize(m_NumTotSpecies, 0.0);

    Diff_Coeff_KRSolid_.resize(numKRSpecies_, 0.0);

    DphMolesSrc_final_.resize(m_NumTotPhases, 0.0);


    actCoeff_Cell_final_.resize(kspCell, 1.0);
    actCoeff_Cell_init_.resize(kspCell, 1.0);

    int maxNumRxns = RSD_List_[0]->nReactions();
    ROP_.resize(maxNumRxns, 0.0);

}
//====================================================================================================================
void
Electrode_SimpleDiff::init_grid()
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
	 rnodePos_final_[iCell+1] = pow(rnode3, 0.33333333333333333);
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
	 rLatticeCBR_final_[iCell] = rnodePos_final_[iCell];
	 rLatticeCBR_final_final_[iCell] = rnodePos_final_[iCell];
	 rLatticeCBR_init_[iCell] = rnodePos_final_[iCell];
	 rLatticeCBR_init_init_[iCell] = rnodePos_final_[iCell];

	 double cbl3 = cellBoundL_final_[iCell] * cellBoundL_final_[iCell] * cellBoundL_final_[iCell];
	 double cbR3 = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
	 volPP_Cell_final_[iCell] = 4. * Pi / 3. * (cbR3 - cbl3);

	 if (iCell ==  0) {
	     rLatticeCBR_ref_[iCell] =  m_rbot0_;
	 } else {
	     rLatticeCBR_ref_[iCell] =  cellBoundR_final_[iCell-1];
	 }
     }
}
//====================================================================================================================
void
Electrode_SimpleDiff::initializeAsEvenDistribution()
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
double Electrode_SimpleDiff::SolidVol() const
{
    double v0_3 = 4. * Pi * m_rbot0_ * m_rbot0_ * m_rbot0_ / 3.;
    double r_ext = rLatticeCBR_final_[numRCells_ - 1];
    double Vext_3 = 4. * Pi * r_ext *  r_ext *  r_ext / 3.;
    double svol = (Vext_3 -  v0_3) * particleNumberToFollow_;
    return svol;
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
void Electrode_SimpleDiff::resetStartingCondition(double Tinitial, bool doResetAlways)
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
//================================================================================================================
//  update the global phase numbers 
/*
 *     This is for distributed phases
 *    We don't calculate a mole fraction vector here. It doesn't make sense to do so.
 *    Instead we take the value of the exterior cell's mole fraction vector.
 *
 *  update phaseMoles_final_[]
 */
void Electrode_SimpleDiff::updatePhaseNumbers(int iph)
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
// Take the state (i.e., the final state) within the Electrode_SimpleDiff and push it up
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
void Electrode_SimpleDiff::updateState_OneToZeroDimensions()
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
     *  Pick a cell to calculate the mole fractions from
     */
    int cellSpecial_ = numRCells_ - 1;
    int kstart = cellSpecial_ * numKRSpecies_;
    for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
	iPh = phaseIndeciseKRsolidPhases_[jRPh];
	kspStart = m_PhaseSpeciesStartIndex[iPh];
	ThermoPhase* th = thermoSPhase_List_[jRPh];
	int nSpecies = th->nSpecies();
	for (int kSp = 0; kSp < nSpecies; kSp++) {
	    spMf_final_[kspStart + kSp] = spMf_KRSpecies_Cell_final_[kstart + kSp];
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
    
	// Here we set the state within the phase object
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
 */
void Electrode_SimpleDiff::updateState()
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

        double volTotalCell = volCell * particleNumberToFollow_;

        int indexMidKRSpecies =  iCell * numKRSpecies_;
        int kstart = 0;
	int indexCellPhase = iCell * numSPhases_;
	/*
	 *  Loop over distributed phases
	 */
        for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
            //iPh = phaseIndeciseKRsolidPhases_[jRPh];
            //kspStart = m_PhaseSpeciesStartIndex[iPh];
            ThermoPhase* th = thermoSPhase_List_[jRPh];
            int nSpecies = th->nSpecies();
	   
            //double total = 0.0;
	    double& concPhase = concTot_SPhase_Cell_final_[indexCellPhase  + jRPh];
	    concPhase = 0.0;
            for (int kSp = 0; kSp < nSpecies; kSp++) {
                int iKRSpecies = kstart + kSp;
                /*
                 * Find the mole numbers of species in the cell, spMolesKRSpecies_Cell_final_[indexTopKRSpecies + iKRSpecies]
                 *     from concKRsolid_Cell_final_;
                 */
		concPhase += concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];
                spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] =
		    volTotalCell * concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];
            }
            /*
             * Find the mole fractions
             *     from spMoles_KRsolid_Cell_final_;
             */
	    if (concPhase > 1.0E-200) {
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMf_KRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] = 
			concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] / concPhase;
		}
	    } else if (concPhase > 1.0E-200) {
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
		concPhase = 0.0;
	    }
	    /*
	     * Calculate the activities of the species
	     */
	    th->setState_TPX(temperature_, pressure_, &spMf_KRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
	    th->getActivityCoefficients(&actCoeff_Cell_final_[indexMidKRSpecies + kstart]);
	    th->getPartialMolarVolumes(&partialMolarVolKRSpecies_Cell_final_[indexMidKRSpecies + kstart]);
	}
    }

    /*
     *  Update the state to zero-dimensional parent fields
     */
    updateState_OneToZeroDimensions();
}
//=========================================================================================================================
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
int Electrode_SimpleDiff::predictSolnResid()
{
   // Indexes
    int iCell, iPh, jPh, jCell;
    double I_j;
    double I_jp1;
    double Lleft, L3left, rjCBR3, rjCBL3;
    double vbarLattice_final_jcell, vbarLattice_init_jcell;
    double caR, caC, caL;
    double rnow3, rnow;
    //double Lright;
    // double rLtarget;

    // Cubes of the cell boundaries, Right and Left, for the initial and final times
    double cbR3_init = 0.0;
    double cbL3_init = 0.0;
    double cbR3_final = 0.0;
    double cbL3_final = 0.0;
    // Reference radius squared, Right and Left, at final time
    double r0R2_final = 0.0;
    double r0L2_final = 0.0;
    // Reference radius cubed, Right and Left, at final time

    double r0R3_final = 0.0;
    double r0L3_final = 0.0;


    double deltaX, fluxTC, fluxM;
    double dcadxR[10];
    double dcadxL[10];

    // Diffusive flux
    double fluxR[10];
    double fluxL[10];
    // Lattice Velocity on the left side of the cell
    double vLatticeCBL = 0.0;
    // Lattice Velocity on the right side of the cell
    double vLatticeCBR = 0.0;
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
    std::vector<doublereal> cellBoundR_init(numRCells_);
    // Velocity of the cell boundary during the time step;
    std::vector<doublereal> cellBoundRVeloc(numRCells_);
    // Node velocity during the time step
    std::vector<doublereal> rnodeVeloc(numRCells_);

    double resid[70];
  

    // predict that the calculated deltaT is equal to the input deltaT
    deltaTsubcycleCalc_ = deltaTsubcycle_;

    /*
     *  Calculate the cell boundaries in a pre - loop. First assume no movement of the mesh
     *  Deal with it later.
     */
    for (iCell = 0; iCell < numRCells_; iCell++) {
	rnodePos_final_[iCell] = rnodePos_init_[iCell];
	rnodeVeloc[iCell] = 0.0;
	cellBoundR_final_[iCell] = 0.0;
	cellBoundRVeloc[iCell] = 0.0;
	rLatticeCBR_ref_[iCell] = cellBoundR_init[iCell];
    }

    // ---------------------------  Main Loop Over Cells ----------------------------------------------------

    for (int iCell = 0; iCell < numRCells_; iCell++) {
        /*
         *  Copy right side previous to left side current quantities
         */
        vLatticeCBL = vLatticeCBR;
        cbL3_final  = cbR3_final;
        cbL3_init   = cbR3_init;
        r0L2_final  = r0R2_final;
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
        cbR3_init  = cellBoundR_init[iCell]   * cellBoundR_init[iCell]   * cellBoundR_init[iCell];
        cbR3_final = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
        /*
         * Calculate powers of the lattice velocity at the right boundary
         */
        double r0R_final = rLatticeCBR_final_[iCell];
        r0R3_final = r0R_final * r0R_final * r0R_final;
        /*
         *  Calculate the molar volume of the first phase, final and init values.
	 *  We will assume for now that this is the lattice molar volume.
         */
        double vbarLattice_final = 1.0 / concTotalVec_SPhase_final[0];
        double vbarLattice_init  = 1.0 / concTotalVec_SPhase_init[0];
        /*
         *  Value of the Lattice radius at the right cell boundary
         */
        double rhs = r0L3_final + vbarLattice_init / vbarLattice_final * (cbR3_final - cbL3_final);
	r0R_final = pow(rhs, 0.333333333333333333);
	//  Find an estimate of the position of the node.
	rLatticeCBR_final_[iCell] = r0R_final;

        /*
         *  Calculate the time derivative of the molar volume
         */
        double vbarDot = (vbarLattice_final - vbarLattice_init) / deltaTsubcycleCalc_;
        /*
         * Find the lattice velocity of the reference radius at the right cell boundary
         */
        double vLatticeCBR = -1.0 / (3. * r0R2_final) *
			     (3.0 *r0L2_final * vLatticeCBL - MolarVolume_Ref_ / (vbarLattice_final * vbarLattice_final) *
			      vbarDot * (cbR3_final - cbL3_final));

	//rLtarget =  rLatticeCBR_ref_[iCell];
	if (rLatticeCBR_final_[iCell] > rLatticeCBR_ref_[iCell]) {
	    // we are here to when the lattice velocity has put the previous lattice CBR into the current cell
	    // Assign the cell id the current cell.
	    jCell = iCell;
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
	    rjCBR3 = cellBoundR_final_[jCell] * cellBoundR_final_[jCell] *  cellBoundR_final_[jCell];
	    rjCBL3 = cellBoundL_final_[jCell] * cellBoundL_final_[jCell] *  cellBoundL_final_[jCell];
	    vbarLattice_final_jcell = 1.0 / concTot_SPhase_Cell_final_[iCell*numSPhases_];
	    vbarLattice_init_jcell  = 1.0 / concTot_SPhase_Cell_init_[iCell*numSPhases_];
	    rnow3 =  L3left + vbarLattice_init_jcell / vbarLattice_final_jcell * (rjCBR3 - rjCBL3);
	    rnow = pow(rnow3, 0.333333333333333);
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
	    rjCBR3 = cellBoundR_final_[jCell] * cellBoundR_final_[jCell] *  cellBoundR_final_[jCell];
	    rjCBL3 = cellBoundL_final_[jCell] * cellBoundL_final_[jCell] *  cellBoundL_final_[jCell];
	    vbarLattice_final_jcell = 1.0 / concTot_SPhase_Cell_final_[iCell*numSPhases_];
	    vbarLattice_init_jcell  = 1.0 /  concTot_SPhase_Cell_init_[iCell*numSPhases_];
	    rnow3 =  L3left +  vbarLattice_init_jcell / vbarLattice_final_jcell * ( rjCBR3 - rjCBL3);
	    rnow = pow(rnow3, 0.333333333333333);
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
	I_j = rnodePos_final_[iCell] * rnodePos_final_[iCell] * rnodePos_final_[iCell];
	if (iCell == numRCells_-1) {
	} else {
	    if (iCell == 0) {
		rnodePos_final_[iCell] = 0;
	    }
	    I_jp1 =  rnodePos_final_[iCell+1] * rnodePos_final_[iCell+1] * rnodePos_final_[iCell+1];
	}
	I_jp1 = I_j - fracVolNodePos_[iCell+1];
	rnodePos_final_[iCell+1] = pow(I_jp1, 0.3333333);

        /*
         *  Calculate the area of the outer cell which conserves constant functions under mesh movement wrt the Reynolds transport theorum
         */
        double cellBoundR_star2 = (cellBoundR_final_[iCell] * cellBoundR_final_[iCell]
                                   - cellBoundR_final_[iCell] * cellBoundR_init[iCell] + cellBoundR_final_[iCell] * cellBoundR_final_[iCell])/3.0;
        double areaR_star = 4.0 * Pi * cellBoundR_star2 * particleNumberToFollow_;

        int indexMidKRSpecies =  iCell    * numKRSpecies_;
        int indexRightKRSpecies = (iCell+1) * numKRSpecies_;
	int indexLeftKRSpecies = (iCell-1) * numKRSpecies_;
        int kstart = 0;

	/*
	 *  We assume that we have the mesh movement calculated before these steps.
	 *           
	 */
        for (jPh = 0; jPh < numSPhases_; jPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jPh];
            ThermoPhase* th = & thermo(iPh);
            int nSpecies = th->nSpecies();
            double old_stuff = concTotalVec_SPhase_init[jPh]  * 4.0 / 3.0 * Pi * (cbR3_init -  cbL3_init)  * particleNumberToFollow_;
            double new_stuff = concTotalVec_SPhase_final[jPh] * 4.0 / 3.0 * Pi * (cbR3_final - cbL3_final) * particleNumberToFollow_;

	    /*
	     * Relative movement
	     */
	    double vtotalR = cellBoundRVeloc[iCell] - vLatticeCBR_cell_[iCell];

            /*
             *  Change in the total amount of moles of phase jPh
             */
	    double volBefore = 4.0 / 3.0 * Pi * (cbR3_init -  cbL3_init)* particleNumberToFollow_;
	    double volAfter =  4.0 / 3.0 * Pi * (cbR3_final - cbL3_final) * particleNumberToFollow_;
	    double molesBefore = concTotalVec_SPhase_init[jPh] * volAfter;


	    if (iCell < (numRCells_- 1)) {
                if (vtotalR >= 0.0) {
                    fluxTC = vtotalR * concTot_SPhase_Cell_final_[numSPhases_*(iCell+1) + jPh];
                } else {
                    fluxTC = vtotalR * concTot_SPhase_Cell_final_[numSPhases_*(iCell)   + jPh];
                }
                resid[cIndexPhStart]                -= fluxTC * areaR_star;
                resid[cIndexPhStart + numEqnsCell_] += fluxTC * areaR_star;
            }

	    double molesAfter = molesBefore + 	vLatticeCBR_cell_[iCell];


            resid[cIndexPhStart] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;
            /*
             * Convective flux - mesh movement with NO material movement
             */
            vtotalR = cellBoundRVeloc[iCell] - vLatticeCBR;
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
             */
            for (int kSp = 1; kSp < nSpecies; kSp++) {
                int iKRSpecies = kstart + kSp;
                /*
                 * Diffusive flux at the outer boundary for each species
                 */
                if (iCell < (numRCells_ - 1)) {
                    deltaX = rnodePos_init_[iCell+1] - rnodePos_init_[iCell];
                    caR = concKRSpecies_Cell_init_[indexRightKRSpecies + iKRSpecies] * actCoeff_Cell_init_[indexRightKRSpecies + iKRSpecies];
                    caC = concKRSpecies_Cell_init_[indexMidKRSpecies + iKRSpecies] * actCoeff_Cell_init_[indexMidKRSpecies + iKRSpecies];
                    dcadxR[kSp] = (caR - caL) / deltaX;
     
                    fluxR[kSp] = - Diff_Coeff_KRSolid_[iKRSpecies] * dcadxR[kSp];
		    

                    resid[cIndexPhStart + kSp]                += fluxR[kSp] * areaR_star;
                    resid[cIndexPhStart + kSp + numEqnsCell_] -= fluxR[kSp] * areaR_star;
                }

                /*
                 * Convective flux - mesh movement with NO material movement
                 */
                if (iCell < (numRCells_- 1)) {
                    if (vtotalR >= 0.0) {
                        fluxM = vtotalR * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell+1) + iKRSpecies];
                        resid[cIndexPhStart + kSp]                -= fluxM * areaR_star;
                        resid[cIndexPhStart + kSp + numEqnsCell_] += fluxM * areaR_star;
                    } else {
                        fluxM = vtotalR * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell)   + iKRSpecies];
                        resid[cIndexPhStart + kSp]                -= fluxM * areaR_star;
                        resid[cIndexPhStart + kSp + numEqnsCell_] += fluxM * areaR_star;
                    }
                }

                /*
                 *  Change in the total amount of moles of phase jPh
                 */
                old_stuff = concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies] * 4.0 / 3.0 * Pi * (cbR3_init  - cbL3_init) * particleNumberToFollow_;
                new_stuff = concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies] * 4.0 / 3.0 * Pi * (cbR3_final - cbL3_final)* particleNumberToFollow_;
                resid[cIndexPhStart] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;
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
                resid[cIndexPhStart] -= DphMolesSrc_final_[iPh];
                SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart] * DspMoles_final_[iStart];
                for (int kSp = 1; kSp < nSpecies; kSp++) {
                    resid[cIndexPhStart + kSp] -= DspMoles_final_[iStart + kSp];
                    SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart + kSp] *  DspMoles_final_[iStart + kSp];
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
//==================================================================================================================
int  Electrode_SimpleDiff::predictSoln()
{
    /* ---------------------------------------------------------------------------------------
     *                                Update the Internal State of ThermoPhase Objects
     * --------------------------------------------------------------------------------------- */
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
     *    
     */    
    int retn = predictSolnResid();
    return retn;
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
void Electrode_SimpleDiff::unpackNonlinSolnVector(const double* const y)
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
	    concKRSpecies_Cell_final_[iCell * numKRSpecies_ + 0] =  
		concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh];
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
//  Gather the predicted solution values and the predicted integrated source terms
/*
 *  (virtual from Electrode_Integrator)
 *
 *  Both the predicted solution values and the predicted integrated source terms are used
 *  in the time step control
 */
void  Electrode_SimpleDiff::gatherIntegratedSrcPrediction()
{
    extractInfo();
    updateSpeciesMoleChangeFinal();
 
    for (int isp = 0; isp < m_NumTotSpecies; isp++) {
        IntegratedSrc_Predicted[isp] = DspMoles_final_[isp] * deltaTsubcycleCalc_;
    }
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
void Electrode_SimpleDiff::setNLSGlobalSrcTermTolerances(double rtolResid)
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
void  Electrode_SimpleDiff::setResidAtolNLS()
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
int  Electrode_SimpleDiff::evalResidNJ(const doublereal t, const doublereal delta_t,
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
int Electrode_SimpleDiff::integrateResid(const doublereal t, const doublereal delta_t,
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

    /*
     * Calculate the residual
     */
    calcResid(resid, evalType);

  
    int index = 1;
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {

        if (phaseID_TimeDeathMin_ >= 0) {

            printf("\t\t minPH = %3d minCell = %3d delTcalc = %16.7E    Res = %16.7E   \n",
                   phaseID_TimeDeathMin_, cellID_TimeDeathMin_, deltaTsubcycleCalc_,  resid[0]);
        }
        printf("\t\t        PhaseName        Moles_Init    Moles_final    |   Src_Moles  Pred_Moles_Final  |    Resid     |\n");
        for (int iph = 0; iph < m_NumTotPhases; iph++) {
            double src =   DphMolesSrc_final_[iph] * deltaTsubcycleCalc_;
            printf("\t\t %20.20s  %12.4e  %12.4e            | %12.4e %12.4e", PhaseNames_[iph].c_str(), phaseMoles_init_[iph],
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
 */
int Electrode_SimpleDiff::calcResid(double* const resid, const ResidEval_Type_Enum evalType)
{
    // Indexes
    int iCell, iPh, jPh, jCell;
    double I_j;
    double I_jp1;
    double Lleft, L3left, rjCBR3, rjCBL3;
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
    double r0R2_final = 0.0;
    double r0L2_final = 0.0;
    // Reference radius cubed, Right and Left, at final time

    double r0R3_final = 0.0;
    double r0L3_final = 0.0;


    double deltaX, fluxTC, fluxM;
    double dcadx[10];

    // Diffusive flux
    double flux[10];
    // Lattice Velocity on the left side of the cell
    double vLatticeCBL = 0.0;
    // Lattice Velocity on the right side of the cell
    double vLatticeCBR = 0.0;
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
    std::vector<doublereal> cellBoundR_init(numRCells_);
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

    /*
     *  Calculate the cell boundaries in a pre - loop
     */
    for (iCell = 0; iCell < numRCells_; iCell++) {
        rnodeVeloc[iCell] = (rnodePos_final_[iCell] - rnodePos_init_[iCell]) / deltaTsubcycleCalc_;
	if (iCell == numRCells_ - 1) {
	    cellBoundR_init[iCell]   = rnodePos_init_[iCell];
	} else {
	    cellBoundR_init[iCell]   = 0.5 * (rnodePos_init_[iCell] + rnodePos_init_[iCell+1]);
	}
        cellBoundRVeloc[iCell] = (cellBoundR_final_[iCell] - cellBoundR_init[iCell]) / deltaTsubcycleCalc_;

	rLatticeCBR_ref_[iCell] =  cellBoundR_init[iCell];
	if (iCell ==  0) {
	    rLatticeCBR_ref_[iCell] =  m_rbot0_;
	} else {
	    rLatticeCBR_ref_[iCell] =  cellBoundR_init[iCell-1];
	}
    }

    // ---------------------------  Main Loop Over Cells ----------------------------------------------------

    for (int iCell = 0; iCell < numRCells_; iCell++) {

        /*
         *  Copy right side previous to left side current quantities
         */
        vLatticeCBL = vLatticeCBR;
        cbL3_final  = cbR3_final;
        cbL3_init   = cbR3_init;
        r0L2_final  = r0R2_final;
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
        cbR3_init  = cellBoundR_init[iCell]  * cellBoundR_init[iCell]  * cellBoundR_init[iCell];
        cbR3_final = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
        /*
         * Calculate powers of the lattice velocity at the right boundary
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
         */
        double rhs = r0L3_final +  vbarLattice_init / vbarLattice_final * (cbR3_final - cbL3_final);
        resid[rindex] = r0R_final - pow(rhs, 0.333333333333333333);

        /*
         *  Calculate the time derivative of the molar volume
         */
        double vbarDot = (vbarLattice_final - vbarLattice_init) / deltaTsubcycleCalc_;

        /*
         * Find the lattice velocity of the reference radius at the right cell boundary
         */
        double vLatticeCBR = -1.0 / (3. * r0R2_final) *
			     (3.0 *r0L2_final * vLatticeCBL - MolarVolume_Ref_ / (vbarLattice_final * vbarLattice_final) *
			      vbarDot * (cbR3_final - cbL3_final));

	//rLtarget =  rLatticeCBR_ref_[iCell];
	if (rLatticeCBR_final_[iCell] > rLatticeCBR_ref_[iCell]) {
	    // we are here to when the lattice velocity has put the previous lattice CBR into the current cell
	    // Assign the cell id the current cell.
	    jCell = iCell;
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
	    rjCBR3 = cellBoundR_final_[jCell] * cellBoundR_final_[jCell] *  cellBoundR_final_[jCell];
	    rjCBL3 = cellBoundL_final_[jCell] * cellBoundL_final_[jCell] *  cellBoundL_final_[jCell];
	    vbarLattice_final_jcell = 1.0 / concTot_SPhase_Cell_final_[iCell*numSPhases_];
	    vbarLattice_init_jcell = 1.0 /  concTot_SPhase_Cell_init_[iCell*numSPhases_];
	    rnow3 =  L3left +  vbarLattice_init_jcell / vbarLattice_final_jcell * ( rjCBR3 - rjCBL3);
	    rnow = pow(rnow3, 0.333333333333333);
	    vLatticeCBR = (rnow -  rLatticeCBR_ref_[iCell]) / deltaTsubcycleCalc_;
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
	    rjCBR3 = cellBoundR_final_[jCell] * cellBoundR_final_[jCell] *  cellBoundR_final_[jCell];
	    rjCBL3 = cellBoundL_final_[jCell] * cellBoundL_final_[jCell] *  cellBoundL_final_[jCell];
	    vbarLattice_final_jcell = 1.0 / concTot_SPhase_Cell_final_[iCell*numSPhases_];
	    vbarLattice_init_jcell = 1.0 /  concTot_SPhase_Cell_init_[iCell*numSPhases_];
	    rnow3 =  L3left +  vbarLattice_init_jcell / vbarLattice_final_jcell * ( rjCBR3 - rjCBL3);
	    rnow = pow(rnow3, 0.333333333333333);
	    vLatticeCBR = (rnow -  rLatticeCBR_ref_[iCell]) / deltaTsubcycleCalc_;
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
	I_j =  rnodePos_final_[iCell] *  rnodePos_final_[iCell] *  rnodePos_final_[iCell];
	if (iCell == numRCells_-1) {
	} else {
	    if (iCell == 0) {
		rnodePos_final_[iCell] = 0;
	    }
	    I_jp1 =  rnodePos_final_[iCell+1] *  rnodePos_final_[iCell+1] *  rnodePos_final_[iCell+1];
	}
        resid[xindex] = I_jp1 - I_j - fracVolNodePos_[iCell+1];

        /*
         *  Calculate the area of the outer cell which conserves constant functions under mesh movement wrt the Reynolds transport theorum
         */
        double cellBoundR_star2 = (cellBoundR_final_[iCell] * cellBoundR_final_[iCell]
                                   - cellBoundR_final_[iCell] * cellBoundR_init[iCell] + cellBoundR_final_[iCell] * cellBoundR_final_[iCell])/3.0;
        double areaR_star = 4.0 * Pi * cellBoundR_star2 * particleNumberToFollow_;


        int indexMidKRSpecies =  iCell    * numKRSpecies_;
        int indexTopKRSpecies = (iCell+1) * numKRSpecies_;
        int kstart = 0;

        for (jPh = 0; jPh < numSPhases_; jPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jPh];
            ThermoPhase* th = & thermo(iPh);
            int nSpecies = th->nSpecies();
            double old_stuff = concTotalVec_SPhase_init[jPh]  * 4.0 / 3.0 * Pi * (cbR3_init -  cbL3_init)  * particleNumberToFollow_;
            double new_stuff = concTotalVec_SPhase_final[jPh] * 4.0 / 3.0 * Pi * (cbR3_final - cbL3_final) * particleNumberToFollow_;

            /*
             *  Change in the total amount of moles of phase jPh
             */
            resid[cIndexPhStart] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;
            /*
             * Convective flux - mesh movement with NO material movement
             */
            double vtotalR = cellBoundRVeloc[iCell] - vLatticeCBR;
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
                 * Convective flux - mesh movement with NO material movement
                 */
                if (iCell < (numRCells_- 1)) {
                    if (vtotalR >= 0.0) {
                        fluxM = vtotalR * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell+1) + iKRSpecies];
                        resid[cIndexPhStart + kSp]                -= fluxM * areaR_star;
                        resid[cIndexPhStart + kSp + numEqnsCell_] += fluxM * areaR_star;
                    } else {
                        fluxM = vtotalR * concKRSpecies_Cell_final_[numKRSpecies_ * (iCell)   + iKRSpecies];
                        resid[cIndexPhStart + kSp]                -= fluxM * areaR_star;
                        resid[cIndexPhStart + kSp + numEqnsCell_] += fluxM * areaR_star;
                    }
                }

                /*
                 *  Change in the total amount of moles of phase jPh
                 */
                old_stuff = concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies] * 4.0 / 3.0 * Pi * (cbR3_init  - cbL3_init) * particleNumberToFollow_;
                new_stuff = concKRSpecies_Cell_init_[numKRSpecies_ * iCell + iKRSpecies] * 4.0 / 3.0 * Pi * (cbR3_final - cbL3_final)* particleNumberToFollow_;
                resid[cIndexPhStart] += (new_stuff - old_stuff) / deltaTsubcycleCalc_;
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
                resid[cIndexPhStart] -= DphMolesSrc_final_[iPh];
                SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart] * DspMoles_final_[iStart];
                for (int kSp = 1; kSp < nSpecies; kSp++) {
                    resid[cIndexPhStart + kSp] -= DspMoles_final_[iStart + kSp];
                    SolidVolCreationRate += partialMolarVolKRSpecies_Cell_final_[iStart + kSp] *  DspMoles_final_[iStart + kSp];
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
//=======================================================================================================================
void  Electrode_SimpleDiff::showSolution(int indentSpaces)
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
void Electrode_SimpleDiff::showResidual(int indentSpaces, const double * const residual) 
{

    drawline(124, indentSpaces);

    std::string title = "Residual for Lattice Radius CBR (m)";

    const double * const res = residual + 1;
    int iTerm = 0;

    showOneResid(title, indentSpaces, &cellBoundR_final_[0], numRCells_, 1, 0, &rLatticeCBR_init_[0],
		 &rLatticeCBR_final_[0], numEqnsCell_, iTerm, &errorLocalNLS_[0], &atolResidNLS_[0],
		 res);

    title = "Residual for position (m)";
    iTerm = 1;
    showOneResid(title, indentSpaces, &rnodePos_final_[0], numRCells_, 1, 0, &rnodePos_init_[0],
		 &rnodePos_final_[0], numEqnsCell_, iTerm, &errorLocalNLS_[0], &atolResidNLS_[0],
		 res);


    for (int k = 0; k <  numKRSpecies_; k++) {
	title = "Residual for Species Concentration 0f "  + KRsolid_speciesNames_[k] + " (kmol/m3)";
	iTerm = 1;
	iTerm = 2 + k;
	showOneResid(title, indentSpaces, &rnodePos_final_[0], numRCells_, numKRSpecies_, k, &concKRSpecies_Cell_final_[0],
		     &concKRSpecies_Cell_final_[0], numEqnsCell_, iTerm, &errorLocalNLS_[0], &atolResidNLS_[0],
		     res);
    }

    drawline(124, indentSpaces);

}

//====================================================================================================================
void  Electrode_SimpleDiff::showOneField(const std::string &title, int indentSpaces,
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
	    printf("\n%s    %-10.4E ", indent.c_str(), r);
	    int istart = iCell * numFields;
	    for (n = 0; n < 5; n++) {
		v = vals[istart + iBlock * 5 + n];
		printf(" %-10.4E ", v);
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
	    printf("\n%s    %-10.4E ", indent.c_str(), r);
	    int istart = iCell * numFields;
	    for (n = 0; n < nrem; n++) {
		v = vals[istart + numBlockRows * 5 + n];
		printf(" %-10.4E ", v);
	    }
	}
	printf("\n");
    }
}
//====================================================================================================================
void  Electrode_SimpleDiff::showOneFieldInitFinal(const std::string &title, int indentSpaces, 
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
	    printf("\n%s    %-10.4E |", indent.c_str(), r);
	    int istart = iCell * numFields;
	    for (n = 0; n < 4; n++) {
		v_init = vals_init[istart + iBlock * 4 + n];
		v_final = vals_final[istart + iBlock * 4 + n];
		printf(" %-10.4E %-10.4E |", v_final, v_init);
	    }
	}
	printf("\n");
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
	    printf("\n%s    %-10.4E |", indent.c_str(), r);
	    int istart = iCell * numFields;
	    for (n = 0; n < nrem; n++) {
		v_init = vals_init[istart + numBlockRows * 4 + n];
		v_final = vals_init[istart + numBlockRows * 4 + n];
		printf(" %-10.4E %-10.4E |", v_final, v_init);
	    }
	}
	printf("\n");
    }
}

//====================================================================================================================
void  Electrode_SimpleDiff::showOneResid(const std::string &title, int indentSpaces,
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
    printf("%s        z   ", indent.c_str());
    
    printf(" %15s", "Init_Value");
    printf(" %15s", "Final_Value");
    if (resid_error) {
	printf(" %15s", "Error");
    }
    if (solnError_tol) {
	printf(" %15s", "Toler");
    }
    printf(" %15s\n", "Resid");
    

 
    drawline(indentSpaces, 80);

    for (iCell = 0; iCell < numRadialVals; iCell++) {
	doublereal r = radialValues[iCell];
	printf("\n%s    %-10.4E ", indent.c_str(), r);
	int istart = iCell * numFields + iField;
	printf(" %-10.4E ", val_init[istart]);
	printf(" %-10.4E ", val_final[istart]);
	istart = iCell *  numEqnsCell + iEqn;
	if (resid_error) {
	    printf(" %-10.4E ", resid_error[istart]);
	}
	if (solnError_tol) {
	    printf(" %-10.4E ", solnError_tol[istart]);
	}
	printf(" | ");
	printf(" %-10.4E ", residual[istart]);
	printf("\n");
    }
    printf("\n");
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
void  Electrode_SimpleDiff::extractInfo()
{
    /*
     *  this is an initial copy from Electrode_CSTR, since it should be nearly the same
     */
    int maxNumRxns = RSD_List_[0]->nReactions();
    std::vector<double> netROP(maxNumRxns, 0.0);
    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    for (int isk = 0; isk < numSurfaces_; isk++) {
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
                 *  Get the species production rates for the reacting surface, spNetProdPerArea
                 *  This is used to integrate cell model equations.
                 */
                getNetProductionRates(isk, spNetProdPerArea);
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
void Electrode_SimpleDiff::updateSpeciesMoleChangeFinal()
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
void Electrode_SimpleDiff::initialPackSolver_nonlinFunction()
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
	deltaBoundsMagnitudesNLS_[index] = Radius_exterior_init_;
	index++;

	yvalNLS_[index] = rnodePos_final_[iCell];
	ylowNLS_[index] = -1.0E300;
	yhighNLS_[index] = 1.0E300;
	deltaBoundsMagnitudesNLS_[index] = Radius_exterior_init_;
	index++;

	for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
	    int nsp = numSpeciesInKRSolidPhases_[jRPh];
	    yvalNLS_[index] = concTot_SPhase_Cell_final_[iCell * numSPhases_ + jRPh];
	  
	    for (int kSp = 1; kSp < nsp; kSp++) {
		iKRSpecies = kstart + kSp;
		yvalNLS_[index + kSp] = concKRSpecies_Cell_final_[iCell * numKRSpecies_ + iKRSpecies];
		ylowNLS_[index + kSp]  = 0.0;
                yhighNLS_[index + kSp] = 1.0;
                deltaBoundsMagnitudesNLS_[index + kSp] = 1.0E-16;
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
void  Electrode_SimpleDiff::setInitStateFromFinal(bool setInitInit)
{
    /*
     * Call the parent object
     */
    Electrode_Integrator::setInitStateFromFinal_Oin(setInitInit);

    int iCell, i;
 
    int ntotal = numRCells_ * numKRSpecies_;
    for (i = 0; i < ntotal; ++i) {
	spMoles_KRsolid_Cell_init_[i] = spMoles_KRsolid_Cell_final_[i];
	concKRSpecies_Cell_init_[i] =  concKRSpecies_Cell_final_[i];
	spMf_KRSpecies_Cell_init_[i] = spMf_KRSpecies_Cell_final_[i];
	actCoeff_Cell_init_[i] = actCoeff_Cell_final_[i];
    }

    int iTotal =  numSPhases_ * numRCells_;
    for (i = 0; i < iTotal; ++i) {
	concTot_SPhase_Cell_init_[i] = concTot_SPhase_Cell_final_[i];
    }

    for (iCell = 0; iCell < numRCells_; ++iCell) {
	rLatticeCBR_init_[iCell] = rLatticeCBR_final_[iCell];
	rnodePos_init_[iCell] = rnodePos_final_[iCell];
    }

    /*
     *  Now transfer that to other states
     */
    setFinalFinalStateFromFinal();

    if (setInitInit) {
	for (i = 0; i < ntotal; ++i) {
	    spMoles_KRsolid_Cell_init_init_[i] = spMoles_KRsolid_Cell_final_[i];
	    concKRSpecies_Cell_init_init_[i] = concKRSpecies_Cell_final_[i];
	}
	
	for (i = 0; i < iTotal; ++i) {
	    concTot_SPhase_Cell_init_init_[i] = concTot_SPhase_Cell_final_[i];
	}

	for (iCell = 0; iCell < numRCells_; ++iCell) {
	    rLatticeCBR_init_init_[iCell] = rLatticeCBR_final_[iCell];
	    rnodePos_init_init_[iCell] = rnodePos_final_[iCell];
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
void Electrode_SimpleDiff::setInitInitStateFromFinalFinal()
{
    Electrode_Integrator::setInitInitStateFromFinalFinal();

    int iCell, i; 
    int ntotal = numRCells_ * numKRSpecies_;
    for (i = 0; i < ntotal; ++i) {
	spMoles_KRsolid_Cell_init_init_[i] = spMoles_KRsolid_Cell_final_final_[i];
	concKRSpecies_Cell_init_init_[i] =  concKRSpecies_Cell_final_final_[i];
	spMoles_KRsolid_Cell_init_init_[i] = spMoles_KRsolid_Cell_final_final_[i];
	concKRSpecies_Cell_init_init_[i] =  concKRSpecies_Cell_final_final_[i];
    }

    int iTotal =  numSPhases_ * numRCells_;
    for (i = 0; i < iTotal; ++i) {
	concTot_SPhase_Cell_init_init_[i] = concTot_SPhase_Cell_final_final_[i];
	concTot_SPhase_Cell_init_init_[i] = concTot_SPhase_Cell_final_[i];
    }

    for (iCell = 0; iCell < numRCells_; ++iCell) {
	rLatticeCBR_init_init_[iCell] = rLatticeCBR_final_final_[iCell];
	rnodePos_init_init_[iCell] = rnodePos_final_final_[iCell];
	rLatticeCBR_init_init_[iCell] = rLatticeCBR_final_final_[iCell];
	rnodePos_init_init_[iCell] = rnodePos_final_final_[iCell];
    }
}
//====================================================================================================================
void Electrode_SimpleDiff::setFinalStateFromInit()
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
	concKRSpecies_Cell_final_[i] = concKRSpecies_Cell_init_[i];
	actCoeff_Cell_final_[i] = actCoeff_Cell_init_[i];
    }
    
    for (i = 0; i < iTotal; ++i) {
	concTot_SPhase_Cell_final_[i] = concTot_SPhase_Cell_init_[i];
    }
    
    for (iCell = 0; iCell < numRCells_; ++iCell) {
	rLatticeCBR_final_[iCell] = rLatticeCBR_init_[iCell];
	rnodePos_final_[iCell] = rnodePos_init_[iCell];
    }
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
void Electrode_SimpleDiff::setInitStateFromInitInit(bool setFinal)
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
void Electrode_SimpleDiff::setFinalFinalStateFromFinal() 
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
void Electrode_SimpleDiff::printElectrode(int pSrc, bool subTimeStep)
{
    int iph;
    vector<std::string> colNames;
    double* netROP = new double[m_NumTotSpecies];
    double egv = TotalVol();
    printf("   ===============================================================\n");
    if (subTimeStep) {
        printf("      Electrode_SimpleDiff at intermediate-step time final = %12.5E\n", tfinal_);
        printf("                              intermediate-step time init  = %12.5E\n", tinit_);
    } else {
        printf("      Electrode_SimpleDiff at time final = %12.5E\n", t_final_final_);
        printf("                              time init  = %12.5E\n", t_init_init_);
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

    colNames.push_back("LatticeRadius");
    std::string title = "";
    int indentSpaces = 10;
    showOneField(title, indentSpaces, &rnodePos_final_[0], numRCells_, &rLatticeCBR_final_[0], colNames, 1);

    for (iph = 0; iph < m_NumTotPhases; iph++) {
        printElectrodePhase(iph, pSrc);
        printf("     ===============================================================\n");
    }
    delete [] netROP;
}
//===================================================================================================================

void Electrode_SimpleDiff::printElectrodePhase(int iph, int pSrc, bool subTimeStep)
{
    int isph;
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
	for (int iCell = 0; iCell < numRCells_; iCell++) {
	    int istart = iCell * nsp;
	    int jstart = iCell * numKRSpecies_;
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

	string title = "    Species Cell Moles (final and init)";
	showOneFieldInitFinal(title, 14, &rnodePos_final_[0], numRCells_, &spMoles_iph_init[0], &spMoles_iph_final[0],
			      speciesNames, nsp);


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
