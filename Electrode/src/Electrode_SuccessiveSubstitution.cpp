/*
 * $Id: Electrode.cpp 533 2013-02-21 19:50:50Z vebruni $
 */

#include "Electrode_SuccessiveSubstitution.h"

using namespace std;

#ifndef SAFE_DELETE
#define SAFE_DELETE(x)  if ((x)) { delete (x) ; x = 0 ; }
#endif

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//======================================================================================================================
// Constructor
Electrode_SuccessiveSubstitution::Electrode_SuccessiveSubstitution() :
    Electrode()
{
}
//======================================================================================================================
// Destructor
Electrode_SuccessiveSubstitution::~Electrode_SuccessiveSubstitution()
{
}
//======================================================================================================================
//! Copy Constructor
/*!
 * @param right Object to be copied
 */
Electrode_SuccessiveSubstitution::Electrode_SuccessiveSubstitution(const Electrode_SuccessiveSubstitution& right) :
    Electrode()
{
    operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_SuccessiveSubstitution&
Electrode_SuccessiveSubstitution::operator=(const Electrode_SuccessiveSubstitution& right)
{
    if (this == &right) {
        return *this;
    }

    Electrode::operator=(right);

    return *this;
}
//======================================================================================================================
// Duplicator function
/*
 *  Duplicate the current electrode object, returning a base electrode pointer
 *
 *  (Virtual function from Electrode.h)
 */
Electrode* Electrode_SuccessiveSubstitution::duplMyselfAsElectrode() const
{
    Electrode_SuccessiveSubstitution* pp = new Electrode_SuccessiveSubstitution(*this);
    return (Electrode*) pp;
}
//======================================================================================================================
//! Set the electrode ID information
Electrode_Types_Enum Electrode_SuccessiveSubstitution::electrodeType() const
{
    return SUCCESSIVE_SUBSTITUTION_ET;
}
//======================================================================================================================
//===================================================================================================================
//  Calculate the change in the state of the system when integrating from T_global_initial to T_global_final
/*
 *  All information is kept internal within this routine. This may be done continuously
 *  and the solution is not updated.
 *
 *  Note the tolerance parameters refers to the nonlinear solves within the calculation
 *  They do not refer to time step parameters.
 *
 *  @param deltaT        DeltaT for the global integration step.
 *  @param rtolResid     Relative tolerance for nonlinear solves within the calculation
 *                       Defaults to 1.0E-3
 *  @param atolResid     Absolute tolerance for nonlinear solves within the calculation
 *                       Defaults to 1.0E-12
 *
 *  @return Returns the number of subcycle steps it took to complete the full step.
 */
int Electrode_SuccessiveSubstitution::integrate(double deltaT, double rtolResid,
        Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType,
        Subgrid_Integration_RunType_Enum subIntegrationType)
{

    counterNumberIntegrations_++;
    std::vector<double> phaseMoles_tmp(m_NumTotPhases, 0.0);
    std::vector<double> phaseMoles_init_init(m_NumTotPhases, 0.0);
    std::vector<double> spMf_tmp(m_NumTotSpecies, 0.0);
    std::vector<double> spMoles_tmp(m_NumTotSpecies, 0.0);
    std::vector<double> Xf_tmp(m_NumTotSpecies, 0.0);
    std::vector<double> delta(m_NumTotSpecies, 0.0);
    std::vector<size_t> justBornMultiSpecies(0);

    std::copy(spMoles_init_.begin(), spMoles_init_.end(), spMoles_init_init_.begin());
    std::copy(spMoles_final_.begin(), spMoles_final_.end(), spMoles_final_final_.begin());
    std::copy(phaseMoles_init_.begin(), phaseMoles_init_.end(), phaseMoles_init_init.begin());

#ifdef OLD_FOLLOW
    followElectrolyteMoles_ = 1;
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    followElectrolyteMoles_ = 0;
#else
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
#endif

    int iterSubCycle = 0;
    bool notDone = true;


    deltaTsubcycle_ = deltaTsubcycle_init_init_;
    if (deltaTsubcycleMax_ <  deltaTsubcycle_) {
        deltaTsubcycle_  = deltaTsubcycleMax_;
    }
    if (choiceDeltaTsubcycle_init_ == 2) {
        deltaTsubcycle_ = deltaT / 10.;
    }
    if (deltaTsubcycle_ > deltaT) {
        deltaTsubcycle_ = deltaT;
    }


    int numDampPhaseTS = 0;
    tfinal_ = t_init_init_;
    tinit_  = t_init_init_;
    t_final_final_ = t_init_init_ + deltaT;
    std::fill(spMoleIntegratedSourceTerm_.begin(), spMoleIntegratedSourceTerm_.end(), 0.);
    std::fill(spMoleIntegratedSourceTermLast_.begin(), spMoleIntegratedSourceTermLast_.end(), 0.);

    /*
     *  When we call this routine successfully we have an integration for the current step pending
     *  Tempory data now exist for t_final_
     */
    pendingIntegratedStep_ = 1;

    /* ----------------------------------------------------------------------------------------------------
     *   LOOP OVER SUBGRID TIME STEPS
     * ---------------------------------------------------------------------------------------------------- */


    do {
        iterSubCycle++;
        counterNumberSubIntegrations_++;

        /*
         *  Advance the time counter
         */
        tinit_ = tfinal_;
        tfinal_ += deltaTsubcycle_;
        if (tfinal_ > (t_final_final_)) {
            tfinal_ = t_final_final_;
            deltaTsubcycle_ = tfinal_ - tinit_;
        }
        /*
         * Advance the initial state replacing it with the final state from the last time iteration
         */
        spMoles_init_ =  spMoles_final_;
        phaseMoles_init_ =  phaseMoles_final_;



        justBornMultiSpecies.clear();

redoTS:
        /*
         * We start a predictor corrector damping cycle here
         */
        for (int iPredCorr = 0; iPredCorr < 200; iPredCorr++) {
            size_t bornMultiSpecies = npos;

restartStep:

            /*
             *   Set the internal objects to the correct conditions
             *    -> This will be the final conditions.
             *    -> Because we do this within an iterator loop, we are essentially trying to emulate
             *       an implicit algorithm
             */
            updateState();

            /*
             * Loop over surface phases, filling in the phase existence fields within the
             * kinetics operator
             */

            for (size_t isk = 0; isk < numSurfaces_; isk++) {
                /*
                 *  Loop over phases, figuring out which phases have zero moles.
                 *  Volume phases exist if the initial or final mole numbers are greater than zero
                 *  Surface phases exist if the initial or final surface areas are greater than zero.
                 */
                if (ActiveKineticsSurf_[isk]) {
                    ReactingSurDomain* rsd = RSD_List_[isk];
                    size_t nph = rsd->nPhases();
                    for (size_t jph = 0; jph < nph; jph++) {
                        size_t iph = rsd->kinOrder[jph];
                        if (iph == metalPhase_) {
                            continue;
                        }
                        double mm = phaseMoles_init_[iph];
                        double mmf = phaseMoles_final_[iph];
                        ThermoPhase& tp = thermo(iph);
                        size_t nsp = tp.nSpecies();
                        if (iph >=  NumVolPhases_) {
                            // we are in a surface phase
                            size_t isur = iph -  NumVolPhases_;
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
                        if (iph == bornMultiSpecies) {
                            rsd->setPhaseExistence(jph, true);
                        }
                        for (size_t iiph = 0; iiph < justBornMultiSpecies.size(); iiph++) {
                            if (iph == justBornMultiSpecies[iiph]) {
                                rsd->setPhaseExistence(jph, true);
                            }
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

            for (size_t isk = 0; isk < numSurfaces_; isk++) {
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
                    size_t nphRS = RSD_List_[isk]->nPhases();
                    size_t jph, kph;
                    size_t kIndexKin = 0;
                    for (kph = 0; kph < nphRS; kph++) {
                        jph = RSD_List_[isk]->kinOrder[kph];
                        size_t istart = m_PhaseSpeciesStartIndex[jph];
                        size_t nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
                        for (size_t k = 0; k < nsp; k++) {
                            spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                            if (rsSpeciesProductionRates[kIndexKin] > 0.0) {
                                if ((phaseMoles_init_[jph] <= 0.0) && (jph != metalPhase_)) {
                                    bool notFound = true;
                                    for (size_t iiph = 0; iiph < justBornMultiSpecies.size(); iiph++) {
                                        if (jph == justBornMultiSpecies[iiph]) {
                                            notFound = false;
                                        }
                                    }
                                    if (notFound) {
                                        if (nsp > 1) {
                                            bornMultiSpecies = jph;
                                        } else {
                                            justBornMultiSpecies.push_back(jph);
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
            double sa_final = calcSurfaceAreaChange(deltaT);
            surfaceAreaRS_final_[0] = sa_final;

            /*
             *  Calculate the change in the moles of all of the species
             */

            for (size_t k = 0; k < m_NumTotSpecies; k++) {
                spMoles_tmp[k] = spMoles_init_[k];
                for (size_t isk = 0; isk < m_NumSurPhases; isk++) {
                    if (ActiveKineticsSurf_[isk]) {
                        sa_init =  surfaceAreaRS_init_[isk];
                        sa_final = surfaceAreaRS_final_[isk];
                        double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
                        spMoles_tmp[k]  += 0.5 * deltaTsubcycle_ * (sa_init + sa_final) * spNetProdPerArea[k];
                    }
                }
            }

            updatePhaseNumbersTmp(spMoles_tmp, phaseMoles_tmp, spMf_tmp);

            /*
             * For phases which are just born, we need to start with a seed in order for the
             * algorithm to work.
             *  The seed needs to be set at a fraction of the initial mole number
             */
            if (bornMultiSpecies != npos) {
                ThermoPhase* tp = PhaseList_[bornMultiSpecies];
                tp->getMoleFractions(DATA_PTR(Xf_tmp));
                int retn = phasePop(bornMultiSpecies, DATA_PTR(Xf_tmp), deltaTsubcycle_);
                if (retn == 0) {
                    tp->setMoleFractions(DATA_PTR(Xf_tmp));
                    size_t istart = m_PhaseSpeciesStartIndex[bornMultiSpecies];
                    size_t nsp =  m_PhaseSpeciesStartIndex[bornMultiSpecies+1] - istart;
                    for (size_t kp = 0; kp < nsp; kp++) {
                        size_t k = istart + kp;
                        spMf_tmp[k] = Xf_tmp[kp];
                    }
                }
                justBornMultiSpecies.push_back(bornMultiSpecies);
                bornMultiSpecies = npos;
                goto   restartStep;
            }

            int doMoleFraction = false;
            double damp = 1.0;
            double dampTS = 1.0;
            bool doPhase = false;

            for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
                int istart = m_PhaseSpeciesStartIndex[iph];
                doPhase = false;
                // ThermoPhase &tp = thermo(solnPhase_);
                int nsp =  m_PhaseSpeciesStartIndex[iph+1] - istart;
                for (int k = 0; k < nsp; k++) {
                    int kk = istart + k;
                    if (spMoles_tmp[kk] < -1.0E-200) {
                        if (spMoles_final_[kk] > 0.0) {
                            double ratio = fabs(spMoles_tmp[kk] + spMoles_final_[kk]) / spMoles_final_[kk];
                            double ldamp = 0.01 / ratio;
                            damp = std::min(ldamp, damp);
                        }
                        /*
                          if (nsp > 1) {
                          double pmW = phaseMoles_tmp[iph] - spMoles_tmp[kk];
                          if (pmW < -1.0E-16) {
                          doPhase = true;
                          } else {
                          doMoleFraction = true;
                          }
                          } else {

                          doPhase = true;
                          }
                        */
                    } else {
                        double ratio;
                        if (spMf_final_[kk] <= 1.0E-200) {
                            ratio = 1500.;
                        } else {
                            ratio= 2.0 * fabs(spMf_tmp[kk] - spMf_final_[kk]) / (fabs(spMf_final_[kk]));
                        }
                        if (ratio > 0.2) {
                            double ldamp = 0.15 / ratio;
                            damp = std::min(ldamp, damp);
                        }
                        ratio = fabs(spMf_tmp[kk] - spMf_final_[kk]);
                        if (ratio > 0.1) {
                            double ldamp = 0.075 / ratio;
                            damp = std::min(ldamp, damp);
                        }

                    }
                    if (doMoleFraction) {
                        double rate =  spMoles_tmp[k] - spMoles_init_[k];
                        double dampRate = -0.6 * spMoles_init_[k] / rate;
                        damp = std::min(dampRate, damp);
                    }
                    if (doPhase) {
                        double rate = phaseMoles_tmp[iph]- phaseMoles_init_[iph];
                        double dampTStmp = phaseMoles_init_[iph] / -rate;
                        dampTS = std::min(dampTS, dampTStmp);
                    }
                }
            }

            if (dampTS < 0.99) {
                if (numDampPhaseTS < 2) {
                    numDampPhaseTS++;
                    deltaTsubcycle_ *= 0.5;
                    goto redoTS;
                } else {
                    numDampPhaseTS = 0.0;
                }
            }

            if (dampTS < 1.0) {
                tfinal_ -= deltaTsubcycle_;
                deltaTsubcycle_ *=  deltaTsubcycle_;
                tfinal_ += deltaTsubcycle_;
            }


            /*
             *  Now apply the damping.
             */
            spMoles_tmp = spMoles_init_;
            for (size_t isk = 0; isk < numSurfaces_; isk++) {
                if (ActiveKineticsSurf_[isk]) {
                    sa_init =  surfaceAreaRS_init_[isk];
                    sa_final = surfaceAreaRS_final_[isk];
                    double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
                    for (size_t k = 0; k < m_NumTotSpecies; k++) {
                        spMoles_tmp[k] += 0.5 * deltaTsubcycle_ * (sa_init + sa_final) * spNetProdPerArea[k];
                    }
                }
            }
            for (size_t k = 0; k < m_NumTotSpecies; k++) {
                delta[k] = spMoles_tmp[k] - spMoles_final_[k];
            }
            for (size_t k = 0; k < m_NumTotSpecies; k++) {
                spMoles_tmp[k] = damp * delta[k] + spMoles_final_[k];
            }

            /*
             *  Accept the tmp vector and assign it to the final_ vector
             */
            for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
                ThermoPhase& tp = thermo(iph);
                string pname = tp.id();
                int istart = m_PhaseSpeciesStartIndex[iph];
                int nsp = tp.nSpecies();
                for (int ik = 0; ik < nsp; ik++) {
                    int k = istart + ik;

                    if (iph == (size_t) solnPhase_) {
                        continue;
                    }
                    if (iph == (size_t) metalPhase_) {
                        continue;
                    }
                    spMoles_final_[k] = spMoles_tmp[k];
                }
            }
            for (size_t i = 0; i < m_NumTotPhases; i++) {
                updateState_Phase(i);
            }


            /*
             *  Check to see whether any phase has gone negative - adjust the mole numbers accordingly.
             *     TODO
             */

            /*
             *  Some mole numbers are not changes by this object. In particular if variables are defined
             *  to be external variables, these variables are defined to be constant over the time step
             *  interval (first cut, here) and therefore not changed.
             *  Variables which fall into this category are the electrolyte mole numbers, and the voltages.
             */

            /*
             *  Check for any condition which may cause the time step to be subdivided. If there
             *  is such a condition, then insert the code here.
             */

            /*
             * Calculate the change in moles of all of the phases. We take the spMoles_final_ value and propagate
             * it to the other values in the object
             */


            if (damp == 1 && iPredCorr > 0) {
                break;
            }
        } // ------------------------- End of Predictor-Corrector 2-cycle ------------------------------------------------

        /*
         *  We are
         */
        for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
            ThermoPhase& tp = thermo(iph);
            string pname = tp.id();
            int istart = m_PhaseSpeciesStartIndex[iph];
            int nsp = tp.nSpecies();
            for (int ik = 0; ik < nsp; ik++) {
                int k = istart + ik;
                spMoleIntegratedSourceTermLast_[k] = (spMoles_tmp[k] - spMoles_init_[k]);
                spMoleIntegratedSourceTerm_[k] += spMoleIntegratedSourceTermLast_[k];

                if (iph == (size_t) solnPhase_) {
                    continue;
                }
                if (iph == (size_t) metalPhase_) {
                    continue;
                }
                spMoles_final_[k] = spMoles_tmp[k];
            }
        }
        if (tfinal_ >= t_final_final_* (1.0 - 1.0E-8)) {
            notDone = false;
            tfinal_ = t_final_final_;
        }
        //-------------------------------- End of Subcycle ---------------------------------------------------------------
    } while (notDone);

    numIntegrationSubCycles_final_final_ = iterSubCycle;
    return iterSubCycle;
}
//===================================================================================================================
//  Residual calculation for the solution of the Nonlinear integration problem
/*
 *    Given tfinal and delta_t, and given y and ydot which are estimates of the solution
 *    and solution derivative at tfinal, this function calculates the residual equations.
 *    It is the residual function used in the nonlinear solver that relaxes the equations
 *    at each time step.
 *
 *    This is typically called from evalResidNJ(), which is called directly from the
 *    nonlinear solver. However, we expose this routine so that the residual can be queried
 *    given all of the inputs.
 *
 *
 *
 * @param[in] t             Time                    (input)
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
int Electrode_SuccessiveSubstitution::integrateResid(const double tfinal, const double deltaTsubcycle,
        const double* const y, const double* const ydot, double* const resid,
        const ResidEval_Type_Enum evalType, const int id_x, const double delta_x)
{

    /*
    double tinit = tfinal - deltaTsubcycle;
    bool newStep= false;
    if (fabs(tinit - tinit_) > 1.0E-14) {
        newStep = true;
    }
    */

    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        spMoles_final_[k] = y[k];
    }

    std::vector<double> phaseMoles_tmp(m_NumTotPhases, 0.0);
    std::vector<double> spMf_tmp(m_NumTotSpecies, 0.0);
    std::vector<double> spMoles_tmp(m_NumTotSpecies, 0.0);
    std::vector<double> Xf_tmp(m_NumTotSpecies, 0.0);

    std::vector<double> srcTerm(m_NumTotSpecies, 0.0);
    static  std::vector<size_t> justBornMultiSpecies;

#ifdef OLD_FOLLOW
    followElectrolyteMoles_ = 1;
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    followElectrolyteMoles_ = 0;
#else
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
#endif

    tfinal_ = tinit_;
    tfinal_ += deltaTsubcycle;

    /*
     * Advance the initial state so that it matches the final state from the last time iteration
     */
    spMoles_init_ =  spMoles_final_;
    phaseMoles_init_ =  phaseMoles_final_;


    justBornMultiSpecies.clear();

    /*
     * We start a predictor corrector damping cycle here
     */
    size_t bornMultiSpecies = npos;

    /*
     *   Set the internal objects to the correct conditions
     *    -> This will be the final conditions.
     */
    updateState();

    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    for (size_t isk = 0; isk < numSurfaces_; isk++) {
        /*
         *  Loop over phases, figuring out which phases have zero moles.
         *  Volume phases exist if the initial or final mole numbers are greater than zero
         *  Surface phases exist if the initial or final surface areas are greater than zero.
         */
        if (ActiveKineticsSurf_[isk]) {
            ReactingSurDomain* rsd = RSD_List_[isk];
            size_t nph = rsd->nPhases();
            for (size_t jph = 0; jph < nph; jph++) {
                size_t iph = rsd->kinOrder[jph];
                if (iph == metalPhase_) {
                    continue;
                }
                double mm = phaseMoles_init_[iph];
                double mmf = phaseMoles_final_[iph];
                if (iph >=  NumVolPhases_) {
                    // we are in a surface phase
                    size_t isur = iph -  NumVolPhases_;
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
                    } else {
                        rsd->setPhaseExistence(jph, true);
                    }
                }
                if (iph == bornMultiSpecies) {
                    rsd->setPhaseExistence(jph, true);
                }
                for (size_t iiph = 0; iiph <  justBornMultiSpecies.size(); iiph++) {
                    if (iph == justBornMultiSpecies[iiph]) {
                        rsd->setPhaseExistence(jph, true);
                    }
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

    for (size_t isk = 0; isk < numSurfaces_; isk++) {
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
            size_t nphRS = RSD_List_[isk]->nPhases();
            size_t jph, kph;
            int kIndexKin = 0;
            for (kph = 0; kph < nphRS; kph++) {
                jph = RSD_List_[isk]->kinOrder[kph];
                size_t istart = m_PhaseSpeciesStartIndex[jph];
                size_t nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
                for (size_t k = 0; k < nsp; k++) {
                    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                    if (rsSpeciesProductionRates[kIndexKin] > 0.0) {
                        if ((phaseMoles_init_[jph] <= 0.0) && (jph != metalPhase_)) {
                            bool notFound = true;
                            for (size_t iiph = 0; iiph < justBornMultiSpecies.size(); iiph++) {
                                if (jph == justBornMultiSpecies[iiph]) {
                                    notFound = false;
                                }
                            }
                            if (notFound) {
                                if (nsp > 1) {
                                    bornMultiSpecies = jph;
                                } else {
                                    justBornMultiSpecies.push_back(jph);
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
    double sa_final = calcSurfaceAreaChange(deltaTsubcycle);
    surfaceAreaRS_final_[0] = sa_final;

    /*
     *  Calculate the change in the moles of all of the species
     */

    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        spMoles_tmp[k] = spMoles_init_[k];
        for (size_t isk = 0; isk < m_NumSurPhases; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                sa_init =  surfaceAreaRS_init_[isk];
                sa_final = surfaceAreaRS_final_[isk];
                double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
                spMoles_tmp[k]  += 0.5 * deltaTsubcycle * (sa_init + sa_final) * spNetProdPerArea[k];
                srcTerm[k] = 0.5 * (sa_init + sa_final) * spNetProdPerArea[k];
            }
        }
    }

    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        resid[k] = (spMoles_final_[k] - spMoles_init_[k]) / deltaTsubcycle - srcTerm[k];
    }

    return 0;
}
//======================================================================================================================
// Return the number of equations in the Nonlinear equation system used to solve the system
// at every time step
/*
 *  This is also equal to the number of state variables in the problem
 */
int Electrode_SuccessiveSubstitution::nEquations() const
{
    int nsp = nSpecies();
    return nsp;
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

