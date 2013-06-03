/*
 * $Id: Electrode_RadialDiffRegions.cpp 298 2012-08-08 20:15:48Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"



#include "cantera/base/mdp_allo.h"

#include "Electrode_RadialDiffRegions.h"
#include "cantera/integrators.h"

#include "BlockEntryGlobal.h"

using namespace Cantera;
using namespace std;
using namespace BEInput;
using namespace TKInput;
using namespace mdpUtil;


#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

namespace Cantera
{

//====================================================================================================================
ELECTRODE_RadialDiffRegions_KEY_INPUT::ELECTRODE_RadialDiffRegions_KEY_INPUT(int printLvl) :
    ELECTRODE_KEY_INPUT(printLvl),
    numRegions_(1),
    solidDiffusionModel_(0),
    rxnPerturbRegions_(0)
{

}
//====================================================================================================================
ELECTRODE_RadialDiffRegions_KEY_INPUT::~ELECTRODE_RadialDiffRegions_KEY_INPUT()
{

}
//====================================================================================================================
ELECTRODE_RadialDiffRegions_KEY_INPUT::ELECTRODE_RadialDiffRegions_KEY_INPUT(const ELECTRODE_RadialDiffRegions_KEY_INPUT& right) :
    ELECTRODE_KEY_INPUT(right),
    numRegions_(right.numRegions_),
    solidDiffusionModel_(right.solidDiffusionModel_)
{
    diffusionCoeffRegions_          = right.diffusionCoeffRegions_;

    rxnPerturbRegions_              = right.rxnPerturbRegions_;
}
//====================================================================================================================
ELECTRODE_RadialDiffRegions_KEY_INPUT&
ELECTRODE_RadialDiffRegions_KEY_INPUT::operator=(const ELECTRODE_RadialDiffRegions_KEY_INPUT& right)
{
    if (this == &right) {
        return *this;
    }

    ELECTRODE_KEY_INPUT::operator=(right);

    numRegions_                     = right.numRegions_;
    solidDiffusionModel_            = right.solidDiffusionModel_;
    diffusionCoeffRegions_          = right.diffusionCoeffRegions_;


    rxnPerturbRegions_              = right.rxnPerturbRegions_;

    return *this;
}
//====================================================================================================================
void ELECTRODE_RadialDiffRegions_KEY_INPUT::setup_input_child1(BEInput::BlockEntry* cf)
{
    /*
     * Obtain the number of regions
     */
    LE_OneInt* s1 = new LE_OneInt("Number of Regions", &(numRegions_), 0, "numRegions");
    s1->set_default(1);
    cf->addLineEntry(s1);
    BaseEntry::set_SkipUnknownEntries(true);

    /*
     * Obtain the number of regions
     */
    LE_OneInt* sdm = new LE_OneInt("Solid Diffusion Model", &(solidDiffusionModel_), 0, "solidDiffusionModel");
    sdm->set_default(0);
    cf->addLineEntry(sdm);
    BaseEntry::set_SkipUnknownEntries(true);


}
//====================================================================================================================
void ELECTRODE_RadialDiffRegions_KEY_INPUT::setup_input_child2(BEInput::BlockEntry* cf)
{

    /*
     * Obtain the number of regions
     */
    int reqd = 0;
    if (solidDiffusionModel_) {
        reqd = 1;
    }


    /*
     * Make the diffusion coefficients optional unless the solid diffusion model is turned on
     */
    reqd = 0;
    if (solidDiffusionModel_) {
        reqd = numRegions_  + 1;
    }



    reqd = 0;

    LE_StdVecDblVarLength* rr1 = new LE_StdVecDblVarLength("Reaction Rate Constant Perturbations for Regions",
            &(rxnPerturbRegions_),
            numRegions_, reqd, "rxnPerturbRegions_");
    rr1->set_default(1.0);
    rr1->set_limits(1.0E30, 1.0E-30);
    cf->addLineEntry(rr1);

    BaseEntry::set_SkipUnknownEntries(false);
}

//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_RadialDiffRegions::Electrode_RadialDiffRegions() :
    Electrode_Integrator(),
    numRadialRegions_(-1),
    radialRegionList_(0),
    SurfaceRegionList_(0)

{


}

//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_RadialDiffRegions::Electrode_RadialDiffRegions(const Electrode_RadialDiffRegions& right) :
    Electrode_Integrator(),
    numRadialRegions_(-1),
    radialRegionList_(0),
    SurfaceRegionList_(0)

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
Electrode_RadialDiffRegions&
Electrode_RadialDiffRegions::operator=(const Electrode_RadialDiffRegions& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode_Integrator::operator=(right);

    numRadialRegions_   = right.numRadialRegions_;
    radialRegionList_   = right.radialRegionList_;
    SurfaceRegionList_  = right.SurfaceRegionList_;

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
Electrode_RadialDiffRegions::~Electrode_RadialDiffRegions()
{


}
//======================================================================================================================
//    Return the type of electrode
/*
 *  Returns the enum type of the electrode. This is used in the factory routine.
 *
 *  @return Returns an enum type, called   Electrode_Types_Enum
 */
Electrode_Types_Enum Electrode_RadialDiffRegions::electrodeType() const
{
    return SIMPLE_DIFF_ET;
}
//======================================================================================================================
int
Electrode_RadialDiffRegions::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{

    Electrode_Integrator::electrode_model_create(ei);




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
Electrode_RadialDiffRegions::setInitialConditions(ELECTRODE_KEY_INPUT* eibase)
{

    /*
     *  Downcast the Key input to make sure we are being fed the correct child object
     */
    ELECTRODE_RadialDiffRegions_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialDiffRegions_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_RadialDiffRegions::electrode_model_create()",
                           " Expecting a child ELECTRODE_KEY_INPUT object and didn't get it");
    }

    Electrode::setInitialConditions(ei);



    setInitStateFromFinal_Oin(true);

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
Electrode_RadialDiffRegions::init_sizes()
{

}
//====================================================================================================================
void
Electrode_RadialDiffRegions::initializeAsEvenDistribution()
{

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
void  Electrode_RadialDiffRegions::resetStartingCondition(double Tinitial, bool doResetAlways)
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
 *
 */
void Electrode_RadialDiffRegions::updateState()
{

}
//================================================================================================
/*
 * There is a small dependence on mf_external and mf_internal exhibited by this function
 */
void  Electrode_RadialDiffRegions::extractInfo(std::vector<int>& justBornMultiSpecies)
{

    updateState();
}
//====================================================================================================================
void Electrode_RadialDiffRegions::printElectrode(int pSrc, bool subTimeStep)
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

void Electrode_RadialDiffRegions::printElectrodePhase(int iph, int pSrc, bool subTimeStep)
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
        mdp::mdp_zero_dbl_1(spNetProdPerArea, m_NumTotSpecies);
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
