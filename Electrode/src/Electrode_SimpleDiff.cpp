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
    phaseIndeciseKRsolidPhases_(0),
    MolarVolume_Ref_(0),
 
    volPP_Cell_final_(0),
    NTflux_final_(0.0),
    DiffCoeff_(1.0E-12),
    DiffCoeff_default_(1.0E-12)
{


}

//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_SimpleDiff::Electrode_SimpleDiff(const Electrode_SimpleDiff& right) :
    Electrode_Integrator(),
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
    phaseIndeciseKRsolidPhases_(0),
    phaseIndeciseNonKRsolidPhases_(0),
    MolarVolume_Ref_(0),
   
    volPP_Cell_final_(0),
    NTflux_final_(0.0),
    DiffCoeff_(1.0E-12),
    DiffCoeff_default_(1.0E-12)
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
    phaseIndeciseKRsolidPhases_         = right.phaseIndeciseKRsolidPhases_;
    concTot_SPhase_Cell_final_final_    = right.concTot_SPhase_Cell_final_final_;
    concTot_SPhase_Cell_final_          = right.concTot_SPhase_Cell_final_;
    concTot_SPhase_Cell_init_           = right.concTot_SPhase_Cell_init_;  
    concTot_SPhase_Cell_init_init_      = right.concTot_SPhase_Cell_init_init_;
    concKRSpecies_Cell_init_            = right.concKRSpecies_Cell_init_;
    concKRSpecies_Cell_final_           = right.concKRSpecies_Cell_final_;
    concKRSpecies_Cell_init_init_       = right.concKRSpecies_Cell_init_init_;
    concKRSpecies_Cell_final_final_     = right.concKRSpecies_Cell_final_final_;
    spMf_KRSpecies_Cell_final_          = right.spMf_KRSpecies_Cell_final_;
    MolarVolume_Ref_                    = right.MolarVolume_Ref_;

    rnodePos_final_final_               = right.rnodePos_final_final_;
    rnodePos_final_                     = right.rnodePos_final_;
    rnodePos_init_                      = right.rnodePos_init_;
    rnodePos_init_init_                 = right.rnodePos_init_init_;
    rRefPos_final_final_                = right.rRefPos_final_final_;
    rRefPos_final_                      = right.rRefPos_final_;
    rRefPos_init_                       = right.rRefPos_init_;
    rRefPos_init_init_                  = right.rRefPos_init_init_;
    cellBoundR_final_                   = right.cellBoundR_final_;
    cellBoundL_final_                   = right.cellBoundL_final_;
    volPP_Cell_final_                   = right.volPP_Cell_final_;

    fracNodePos_                        = right.fracNodePos_;
    fracVolNodePos_                     = right.fracVolNodePos_;
    partialMolarVolKRSpecies_Cell_final_= right.partialMolarVolKRSpecies_Cell_final_;
    DspMoles_final_                     = right.DspMoles_final_;
    m_rbot0_                            = right.m_rbot0_;
    Diff_Coeff_KRSolid_                 = right.Diff_Coeff_KRSolid_;
    DphMolesSrc_final_                  = right.DphMolesSrc_final_;
    surfIndexExteriorSurface_           = right.surfIndexExteriorSurface_;

    NTflux_final_                       = right.NTflux_final_;
    DiffCoeff_                          = right.DiffCoeff_;
    DiffCoeff_default_                  = right.DiffCoeff_default_;
    actCoeff_Cell_final_                = right.actCoeff_Cell_final_;
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
//======================================================================================================================
int
Electrode_SimpleDiff::electrode_model_create(ELECTRODE_KEY_INPUT* eibase)
{
    int iPh;
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
    for (int i = 0; i < (int) phaseIndeciseKRsolidPhases_.size(); i++) {
        iPh =  phaseIndeciseKRsolidPhases_[i];
        ThermoPhase& th = thermo(iPh);
        int nsp = th.nSpecies();
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
    int KRsolid = 0;
    for  (int i = 0; i < (int) phaseIndeciseKRsolidPhases_.size(); i++) {
	iPh =  phaseIndeciseKRsolidPhases_[i];
	ThermoPhase& th = thermo(iPh);
	thermoSPhase_List_[i] = &th;
	int nsp = th.nSpecies();
	int kstart = getGlobalSpeciesIndex(iPh);
	for (int kk = 0; kk < nsp; kk++) {
	    KRsolid_speciesList_[KRsolid] = kstart + kk;
	    KRsolid++;
	}
    }
    
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
     *  Initialize all of the variables on the domain
     *    We take the grid size and the mole numbers from the base Electrode representation
     *    and create a distributed representation.
     *    Any disparity in mole numbers creates an error.
     */
    initializeAsEvenDistribution();

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

    spMf_KRSpecies_Cell_final_.resize(kspCell, 0.0);

    MolarVolume_Ref_ = 1.0;

    rnodePos_final_final_.resize(numRCells_, 0.0);
    rnodePos_final_.resize(numRCells_, 0.0);
    rnodePos_init_.resize(numRCells_, 0.0);
    rnodePos_init_init_.resize(numRCells_, 0.0);

    rRefPos_final_final_.resize(numRCells_, 0.0);
    rRefPos_final_.resize(numRCells_, 0.0);
    rRefPos_init_.resize(numRCells_, 0.0);
    rRefPos_init_init_.resize(numRCells_, 0.0);

    cellBoundR_final_.resize(numRCells_, 0.0);
    cellBoundL_final_.resize(numRCells_, 0.0);

    volPP_Cell_final_.resize(numRCells_, 0.0);
    fracNodePos_.resize(numRCells_, 0.0);
    fracVolNodePos_.resize(numRCells_, 0.0);

    partialMolarVolKRSpecies_Cell_final_.resize(kspCell, 0.0);

    DspMoles_final_.resize(m_NumTotSpecies, 0.0);

    Diff_Coeff_KRSolid_.resize(numKRSpecies_, 0.0);

    DphMolesSrc_final_.resize(m_NumTotPhases, 0.0);


    actCoeff_Cell_final_.resize(kspCell, 1.0);
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
	 rRefPos_final_[iCell] = rnodePos_final_[iCell];
	 rRefPos_final_final_[iCell] = rnodePos_final_[iCell];
	 rRefPos_init_[iCell] = rnodePos_final_[iCell];
	 rRefPos_init_init_[iCell] = rnodePos_final_[iCell];

	 double cbl3 = cellBoundL_final_[iCell] * cellBoundL_final_[iCell] * cellBoundL_final_[iCell];
	 double cbR3 = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
	 volPP_Cell_final_[iCell] = 4. * Pi / 3. * (cbR3 - cbl3);
     }
}
//====================================================================================================================
void
Electrode_SimpleDiff::initializeAsEvenDistribution()
{
    int iCell, i, k, KRSolid,  kspCell, iphCell;
    /*
     *  First, get the mole fractions 
     */
    for (iCell = 0; iCell < numRCells_; ++iCell) {
	for (KRSolid = 0; KRSolid <  numKRSpecies_; KRSolid++) {
	    k = KRsolid_speciesList_[KRSolid];
	    i = KRSolid + numKRSpecies_ * iCell;
	    spMf_KRSpecies_Cell_final_[i] = spMf_final_[k];
	}
    }

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
//    The internal state of the electrode must be kept for the initial and final times of an integration step.
/*
 *  This function advances the initial state to the final state that was calculated
 *  in the last integration step.
 *
 * @param Tinitial   This is the New initial time. This time is compared against the "old"
 *                   final time, to see if there is any problem.
 */
void  Electrode_SimpleDiff::resetStartingCondition(double Tinitial, bool doResetAlways)
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
/*!
 *     This is for distributed phases
 *    We don't calculate a mole fraction vector here. It doesn't make sense to do so.
 *    Instead we take the value of the exterior cell's mole fraction vector.
 *
 *  uupdate phaseMoles_final_[]
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
 * 
 *  update 
 *             spMoles_final_ [] -> sum solid phase species
 *              spMf_final_[]  -> Use exterior cell values
 *
 *
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

    /*
     *  Calculate the voltage field
     */
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
 *       concKRSpecies_Cell_final_[]
 *       rnodePos_final_[iCell]
 *        rRefPos_final_[iCell]
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
 */
void Electrode_SimpleDiff::updateState()
{
    // Indexes
    int iCell, iPh, jRPh, kspStart;
    double tmp;
    // Cubes of the cell boundaries, Right and Left, for the initial and final times
    double cbR3_final = 0.0;
    double cbL3_final = 0.0;

    /*
     *   We now have a vector of cells.
     *   Each of the cells must have their own conditions.
     *   We need to loop over the conditions and have their activivities calculated
     */
    /*
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
        double volCell = 4./3. * (cbR3_final  - cbL3_final);

        double volTotalCell = volCell * particleNumberToFollow_;

        int indexMidKRSpecies =  iCell * numKRSpecies_;
        int kstart = 0;
	int indexCellPhase = iCell * numSPhases_;
        for (jRPh = 0; jRPh < numSPhases_; jRPh++) {
            iPh = phaseIndeciseKRsolidPhases_[jRPh];
            kspStart = m_PhaseSpeciesStartIndex[iPh];
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
                spMoles_KRsolid_Cell_final_[indexMidKRSpecies + iKRSpecies] = volTotalCell * concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies];
            }
            /*
             * Find the mole fractions
             *     from spMoles_KRsolid_Cell_final_;
             */
	    if (concPhase > 1.0E-200) {
		for (int kSp = 0; kSp < nSpecies; kSp++) {
		    int iKRSpecies = kstart + kSp;
		    spMf_KRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] = concKRSpecies_Cell_final_[indexMidKRSpecies + iKRSpecies] / concPhase;
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



     throw CanteraError("Electrode_SimpleDiff::updateState()", "unfinished");
}
//====================================================================================================================
//   Evaluate the residual function
/*
 * (virtual from NonlinearSolver)
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
 *
 * @return
 */
int  Electrode_SimpleDiff::evalResidNJ(const doublereal t, const doublereal delta_t,
                                       const doublereal* const y,
                                       const doublereal* const ydot,
                                       doublereal* const resid,
                                       const ResidEval_Type_Enum evalType,
                                       const int id_x,
                                       const doublereal delta_x)
{

    throw CanteraError("", "");
    return 0;
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

int Electrode_SimpleDiff::calcResid(double* const resid, const ResidEval_Type_Enum evalType)
{
    // Indexes
    int iCell, iPh, jPh;
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
    // Velocity of the reference radius on the left side of the cell
    double vrefL = 0.0;
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
     *  Calculate the cell boundaries in a pre - loop
     */
    for (iCell = 0; iCell < numRCells_; iCell++) {
        rnodeVeloc[iCell] = (rnodePos_final_[iCell] - rnodePos_init_[iCell]) / deltaTsubcycleCalc_;

        cellBoundR_init[iCell]   = 0.5 * (rnodePos_init_[iCell] + rnodePos_init_[iCell+1]);
        cellBoundRVeloc[iCell] = (cellBoundR_final_[iCell] - cellBoundR_init[iCell]) / deltaTsubcycleCalc_;
    }

    // ---------------------------  Main Loop Over Cells ----------------------------------------------------

    for (int iCell = 0; iCell < numRCells_; iCell++) {

        /*
         *  Copy left side to right side
         */
        vrefL       = vrefR;
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
         *  Temporary variables for total concentrations of each of the phases
         */
        concTotalVec_SPhase_init  = &(concTot_SPhase_Cell_init_[iCell*numSPhases_]);
        concTotalVec_SPhase_final = &(concTot_SPhase_Cell_final_[iCell*numSPhases_]);
        /*
         * Calculate the cubes of the cell boundary radii
         */
        cbR3_init  = cellBoundR_init[iCell]  * cellBoundR_init[iCell]  * cellBoundR_init[iCell];
        cbR3_final = cellBoundR_final_[iCell] * cellBoundR_final_[iCell] * cellBoundR_final_[iCell];
        /*
         * Calculate powers of the reference velocity at the right boundary
         */
        double r0R_final = rRefPos_final_[iCell];
        r0R2_final = r0R_final * r0R_final;
        /*
         * Calculate the molar volume of the first phase, final and init
         */
        double vbar_final = 1.0 / concTotalVec_SPhase_final[0];
        double vbar_init  = 1.0 / concTotalVec_SPhase_init[0];
        /*
         * Residual calculation - Value of the reference radius at the right cell boundary
         */
        double rhs = r0L3_final +  MolarVolume_Ref_ / vbar_final * (cbR3_final - cbL3_final);
        resid[rindex] = r0R_final - pow(rhs, 0.333333333333333333);
        /*
         *  Calculate the time derivative of the molar volume
         */
        double vbarDot = (vbar_final - vbar_init) / deltaTsubcycleCalc_;

        /*
         * Find the velocity of the reference radius at the right cell boundary
         */
        vrefR = -1.0 / (3. * r0R2_final) * (3.0 *r0L2_final * vrefL - MolarVolume_Ref_ / (vbar_final * vbar_final) * vbarDot * (cbR3_final - cbL3_final));

        /*
         * Node position residual - spline equations, it all depends on the top node's equation formulation.
         * Everything else get's dragged along with it
         */
        resid[xindex] = rnodePos_final_[iCell] - m_rbot0_
                        - fracNodePos_[iCell] / fracNodePos_[iCell + 1] * (rnodePos_final_[iCell+1] - m_rbot0_);

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
//================================================================================================
/*
 * There is a small dependence on mf_external and mf_internal exhibited by this function
 */
void  Electrode_SimpleDiff::extractInfo(std::vector<int>& justBornMultiSpecies)
{



    updateState();

}


//====================================================================================================================
void Electrode_SimpleDiff::printElectrode(int pSrc, bool subTimeStep)
{
    int iph;
    double* netROP = new double[m_NumTotSpecies];
    double egv = TotalVol();
    printf("   ===============================================================\n");
    if (subTimeStep) {
        printf("      Electrode at intermediate-step time final = %g\n", tfinal_);
        printf("                   intermediate-step time init  = %g\n", tinit_);
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
    delete [] netROP;

}


} // End of namespace Cantera
//======================================================================================================================
