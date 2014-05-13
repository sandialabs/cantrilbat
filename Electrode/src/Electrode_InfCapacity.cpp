/*
 * $Id: Electrode_InfCapacity.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */

#include "Electrode_InfCapacity.h"

using namespace std;

namespace Cantera
{
//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_InfCapacity::Electrode_InfCapacity() :
    Electrode()
{

}

//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_InfCapacity::Electrode_InfCapacity(const Electrode_InfCapacity& right) :
    Electrode()
{
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
}
//======================================================================================================================

//! Return the type of electrode
/*!
 *  Returns the enum type of the electrode. This is used in the factory routine.
 *
 *  @return Returns an enum type, called   Electrode_Types_Enum
 */
Electrode_Types_Enum  Electrode_InfCapacity::electrodeType() const
{
    return INF_CAPACITY_ET;
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_InfCapacity& Electrode_InfCapacity::operator=(const Electrode_InfCapacity& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode::operator=(right);


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
Electrode_InfCapacity::~Electrode_InfCapacity()
{
}
//======================================================================================================================
//  Setup the electrode
/*
 * @param ei    ELECTRODE_KEY_INPUT pointer object
 */
int
Electrode_InfCapacity::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{

    Electrode::electrode_model_create(ei);

    setPhaseExistenceForReactingSurfaces(true);
    return 0;
}

//===================================================================================================================
//  Calculate the change in the state of the system when integrating from Tinitial to Tfinal
/*
 *  All information is kept internal within this routine. This may be done continuously
 *  and the solution is not updated.
 *
 *  Note the tolerance parameters refere to the nonlinear solves within the calculation
 *  They do not refer to time step parameters.
 *
 *  @param deltaT        DeltaT for the integration step.
 *  @param GlobalRtolSrcTerm    Relative tolerance allowed for the electron source term over the interval.
 *                       This is a unitless quantity
 *                       Defaults to 1.0E-3
 *  @param fieldInterpolationType Type of interpolation of field variables defaults to 0
 *  @param subIntegrationType     Type of subintegration. defaults to 0
 *
 *  @return Returns the number of subcycle steps it took to complete the full step.
 */
int Electrode_InfCapacity::integrate(double deltaT, double  GlobalRtolSrcTerm,
                                     Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType,
                                     Subgrid_Integration_RunType_Enum subIntegrationType)
{

    double prodTmp;
    double sa_final;
    counterNumberIntegrations_++;
    counterNumberSubIntegrations_++;
    std::vector<doublereal> phaseMoles_tmp(m_NumTotPhases, 0.0);
    std::vector<doublereal> phaseMoles_init_init(m_NumTotPhases, 0.0);
    std::vector<doublereal> spMf_tmp(m_NumTotSpecies, 0.0);
    std::vector<doublereal> spMoles_tmp(m_NumTotSpecies, 0.0);
    std::vector<doublereal> Xf_tmp(m_NumTotSpecies, 0.0);
    std::vector<doublereal> delta(m_NumTotSpecies, 0.0);
    std::vector<int> justBornMultiSpecies(0);


    std::copy(spMoles_init_init_.begin(), spMoles_init_init_.end(), spMoles_init_.begin());
    std::copy(spMoles_init_init_.begin(), spMoles_init_init_.end(), spMoles_final_.begin());
    std::copy(phaseMoles_init_init_.begin(), phaseMoles_init_init_.end(), phaseMoles_init_.begin());
    std::copy(phaseMoles_init_init_.begin(), phaseMoles_init_init_.end(), phaseMoles_final_.begin());
    /*
     *   Set the internal objects to the correct conditions
     *    -> This will be the final conditions.
     *    -> Because we do this within an iterator loop, we are essentially trying to emulate
     *       an implicit algorithm
     */
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];


    double deltaTsubcycle = deltaT;
    t_final_final_ = t_init_init_ + deltaT;
    tinit_ = t_init_init_;
    tfinal_ = tinit_;

    tfinal_ += deltaTsubcycle;
    if (tfinal_ > (t_final_final_)) {
        tfinal_ = t_final_final_;
    }

    justBornPhase_.clear();
    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    setPhaseExistenceForReactingSurfaces(true);


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
            const vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetSurfaceProductionRateDensities();


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
                int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
                for (int k = 0; k < nsp; k++) {
                    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                    if (rsSpeciesProductionRates[kIndexKin] > 0.0) {
                        if ((phaseMoles_init_[jph] <= 0.0) && (jph != metalPhase_)) {
                            bool notFound = true;
                            for (int iiph = 0; iiph < (int)justBornMultiSpecies.size(); iiph++) {
                                if (jph == justBornMultiSpecies[iiph]) {
                                    notFound = false;
                                }
                            }
                            if (notFound) {
                                justBornMultiSpecies.push_back(jph);
                            }
                        }
                    }
                    kIndexKin++;
                }
            }
        }
    }


    /*
     *  Calculate the change in the moles of all of the species
     */
    std::fill(spMoleIntegratedSourceTerm_.begin(), spMoleIntegratedSourceTerm_.end(), 0.);
    for (int k = 0; k < m_NumTotSpecies; k++) {
        spMoles_tmp[k] = spMoles_init_[k];
        for (int isk = 0; isk < numSurfaces_; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                sa_final = surfaceAreaRS_final_[isk];
                double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
                prodTmp = deltaTsubcycle * (sa_final) * spNetProdPerArea[k];
                spMoles_tmp[k] += prodTmp;
                spMoleIntegratedSourceTerm_[k] += prodTmp;
            }
        }
    }


    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase& tp = thermo(iph);
        string pname = tp.id();
        int istart = m_PhaseSpeciesStartIndex[iph];
        int nsp = tp.nSpecies();
        for (int ik = 0; ik < nsp; ik++) {
            int k = istart + ik;
            // HKM -> Note this technique lead to round off error for this case. spMoleIntegratedSourceTerm_[k]
            //        had numerical roundoff error in the 5th or 6th digit.
            //    spMoleIntegratedSourceTerm_[k] = (spMoles_tmp[k] - spMoles_init_[k]);

            if (iph == solnPhase_) {
                continue;
            }
            if (iph == metalPhase_) {
                continue;
            }
            // NOTHING CHANGES HERE -> this is what we mean by infinite capacity
            spMoles_final_[k] = spMoles_init_[k];
        }
    }
    /*
     * indicate that we have a pending integrated step
     */
    pendingIntegratedStep_ = 1;
    return 1;
}
//====================================================================================================================
// The internal state of the electrode must be kept for the initial and
// final times of an integration step.
/*
 *  This function advances the initial state to the final state that was calculated
 *  in the last integration step.
 *
 * @param Tinitial   This is the New initial time. This time is compared against the "old"
 *                   final time, to see if there is any problem.
 */
void  Electrode_InfCapacity::resetStartingCondition(double Tinitial, bool doTestsAlways)
{
    /*
     * If the initial time is input, then the code doesn't advance
     */
    double tbase = std::max(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-9 * tbase) && !doTestsAlways) {
        return;
    }
    Electrode::resetStartingCondition(Tinitial);
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
void Electrode_InfCapacity::setInitStateFromFinal(bool setInitInit)
{
    Electrode::setInitStateFromFinal(setInitInit);
}
//====================================================================================================================
void Electrode_InfCapacity::getIntegratedPhaseMoleTransfer(doublereal* const phaseMolesTransfered)
{

    if (!pendingIntegratedStep_) {
        throw CanteraError(" Electrode_InfCapacity::getIntegratedPhaseMoleTransfer",
                           "no pending integration step");
    }
    double sum = 0.0;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
        phaseMolesTransfered[iph] = 0.0;
        Cantera::ThermoPhase& tp = thermo(iph);
        string pname = tp.id();
        int istart = m_PhaseSpeciesStartIndex[iph];
        int nsp = tp.nSpecies();
        for (int ik = 0; ik < nsp; ik++) {
            int k = istart + ik;
            phaseMolesTransfered[iph] += spMoleIntegratedSourceTerm_[k];
        }
        sum += fabs(phaseMolesTransfered[iph]);
    }
}
//====================================================================================================================
} // End of namespace Cantera
//======================================================================================================================
