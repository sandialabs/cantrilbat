/*
 * $Id: Electrode_Equilibrium.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */

#include "Electrode_Equilibrium.h"

#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/equilibrium.h"


using namespace std;

#ifndef SAFE_DELETE
#define SAFE_DELETE(x)  if (x) { delete x;  x = 0;}
#endif
//#ifndef MAX
//#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
//#endif

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{
//==================================================================================================================================
Electrode_Equilibrium::Electrode_Equilibrium(Electrode* elect) :
    ee_(elect),
    m_mp(0),
    LiFixed_(0),
    printLvl_(0)
{
    if (ee_) {
        printLvl_ = std::max(0, printLvl_ - 3);
    }
}
//==================================================================================================================================
Electrode_Equilibrium::Electrode_Equilibrium(const Electrode_Equilibrium& right) :
    ee_(right.ee_),
    m_mp(0),
    LiFixed_(0),
    printLvl_(0)
{
    operator=(right);
}
//==================================================================================================================================
// Destructor
Electrode_Equilibrium::~Electrode_Equilibrium()
{
    delete m_mp;
    m_mp = 0;
    delete LiFixed_;
    LiFixed_ = 0;
}
//==================================================================================================================================
Electrode_Equilibrium& Electrode_Equilibrium::operator=(const Electrode_Equilibrium& right)
{
    if (this == &right) {
        return *this;
    }
    /*
     *  Do a shallow copy of the Electrode pointer. This is all that is necessary
     */
    ee_ = right.ee_;

    if (LiFixed_) {
        delete LiFixed_;
    }
    if (m_mp) {
        delete m_mp;
    }
    setupEquilibriumProblem();

    if (right.LiFixed_) {
        LiFixed_ = (ZZCantera::FixedChemPotSSTP*) (right.LiFixed_)->duplMyselfAsThermoPhase();
    }

    printLvl_ = right.printLvl_;

    return *this;
}
//==================================================================================================================================
int Electrode_Equilibrium::setupEquilibriumProblem()
{
    if (m_mp) {
        delete m_mp;
    }
    m_mp = new MP_EquilStatic();
    PhaseIndex_mp_.resize(ee_->m_NumVolPhases);
    for (size_t iph = 0; iph < ee_->m_NumVolPhases; iph++) {
        ThermoPhase* tp = ee_->VolPhaseList_[iph];
        m_mp->addPhase(tp, ee_->phaseMoles_final_[iph]);
        PhaseIndex_mp_[iph] = iph;
    }

    /*
     *   The multispecies object as it is constituted above is actually useless, because
     *   it only represents a half-cell reaction which doesn't allow for the movement of the ion
     *   between phases. Instead an element potential for the ion must be created.
     */
    if (LiFixed_) {
        delete LiFixed_;
    }

    printLvl_ = std::max(0, printLvl_ - 3);
    return 0;
}
//===========================================================================================================
int Electrode_Equilibrium::addLiIonPlusElectronBathPhase()
{
    if (LiFixed_) {
        delete LiFixed_;
    }
    LiFixed_ = new FixedChemPotSSTP("Li", -2.3E7);

    /*
     * Find the electrochemical potential of the Li+ in the electrolyte
     */
    const ThermoPhase* tElec = &ee_->thermo(ee_->metalPhase_);
    const ThermoPhase* tSalt = &ee_->thermo(ee_->solnPhase_);
    double mu[30];
    tSalt->getElectrochemPotentials(mu);
    //  We assume that the ion in the salt is called Li+
    size_t ilt_Lip = tSalt->speciesIndex("Li+");
    if (ilt_Lip == npos) {
        throw Electrode_Error("Electrode_Equilibrium::addLiIonPlusElectron()",
                              "No species named Li+. Assumptions are violated");
    }
    double muLip = mu[ilt_Lip];
    /*
     * Find the electrochemical potential of the electron in the electrode
     */
    tElec->getElectrochemPotentials(mu);
    double muElec = mu[0];
    /*
     * Define the elemental Li chemical potential as the sum of Li+ and e-
     *     -> the voltage drop will be built into this value.
     */
    double muLi = muLip + muElec;
    LiFixed_->setChemicalPotential(muLi);

    double tMoles = 0.0;
    for (size_t k = 0; k < ee_->m_NumTotSpecies; k++) {
        tMoles += ee_->spMoles_final_[k];
    }
    /*
     *  Add this bath phase into the MultiSpecies
     */
    m_mp->addPhase(LiFixed_, 10. * tMoles);
    PhaseIndex_mp_.push_back(-1);

    return 0;
}
//===========================================================================================================
//    Set the state of the captive ThermoPhase objects and the multiphase object from the current
//    final state of this Electrode.
void Electrode_Equilibrium::update_MP()
{
    int iph;
    int mp_nPhases = m_mp->nPhases();
    m_mp->setState_TP(ee_->temperature_, ee_->pressure_);

    /*
     * Do the volume phases
     */
    for (int mph = 0; mph < mp_nPhases; mph++) {
        iph = PhaseIndex_mp_[mph];
        if (iph >= 0) {
            //ThermoPhase *tp = &(m_mp->phase(mph));
            int eStart = ee_->m_PhaseSpeciesStartIndex[iph];
            m_mp->setPhaseMoles(mph, ee_->phaseMoles_final_[iph]);
            m_mp->setPhaseMoleFractions(mph, &(ee_->spMf_final_[eStart]));
        } else {
            ThermoPhase* tp = &(m_mp->phase(mph));
            string pname = tp->id();
            if (pname == "LiFixed") {
                FixedChemPotSSTP* LiFixed = dynamic_cast<FixedChemPotSSTP*>(tp);
                const ThermoPhase* tElec = &ee_->thermo(ee_->metalPhase_);
                const ThermoPhase* tSalt = &ee_->thermo(ee_->solnPhase_);
                double mu[30];
                tSalt->getElectrochemPotentials(mu);
                //  We assume that the ion in the salt is called Li+
                int ilt_Lip = tSalt->speciesIndex("Li+");
                double muLip = mu[ilt_Lip];
                tElec->getElectrochemPotentials(mu);
                double muElec = mu[0];
                double muLi = muLip + muElec;
                LiFixed->setChemicalPotential(muLi);
            } else {
                throw Electrode_Error("downloadMP()", "unknown event");
            }
        }
    }
}
//====================================================================================================================
// Set the state of the electrode object, from the captive
// MultiPhase object
void Electrode_Equilibrium::uploadMP()
{
    int iph, k;

    ee_->temperature_ = m_mp->temperature();
    ee_->pressure_ = m_mp->pressure();

    int mp_nPhases = m_mp->nPhases();
    for (int mph = 0; mph < mp_nPhases; mph++) {
        iph = PhaseIndex_mp_[mph];
        ThermoPhase* tp = &(m_mp->phase(mph));
        int mStart = m_mp->speciesIndex(0, mph);
        int nSpecies = tp->nSpecies();

        int eStart = ee_->m_PhaseSpeciesStartIndex[iph];

        for (k = 0; k < nSpecies; k++) {
            ee_->spMoles_final_[eStart + k] = m_mp->speciesMoles(mStart + k);
            ee_->spMf_final_[eStart + k] = m_mp->moleFraction(mStart + k);
        }
        ee_->setPhaseMoleNumbers(iph, &(ee_->spMoles_final_[eStart]));

        tp->getPartialMolarVolumes(&(ee_->VolPM_[eStart]));
        tp->getElectrochemPotentials(&(ee_->spElectroChemPot_[eStart]));

        ee_->phaseVoltages_[iph] = tp->electricPotential();
        ee_->phaseMoles_final_[iph] = m_mp->phaseMoles(mph);

        ee_->phaseMolarVolumes_[mph] = tp->molarVolume();
    }

    ee_->deltaVoltage_ = ee_->phaseVoltages_[ee_->metalPhase_] - ee_->phaseVoltages_[ee_->solnPhase_];

    ee_->ElectrodeSolidVolume_ = ee_->SolidVol();

    // Now call the virtual function that takes care of all of the rest of the update
    ee_->updateState();
}
//==============================================================================================================
// Predict the stability of a single phase that is currently zeroed.
/*
 *  This routine fills in an estimate for the solution
 *  Return 1 if the phases are stable and 0 if they are not
 */
int Electrode_Equilibrium::predictStabilitySinglePhase(int iphase, double& funcStab)
{

    double solidMoles = ee_->SolidTotalMoles();

    int pElec = ee_->globalPhaseIndex("metal_Li_LiCl_electrons");
    int pSalt = ee_->globalPhaseIndex("LiKCl_electrolyte");
    ThermoPhase* tElec = &ee_->thermo(pElec);
    ThermoPhase* tSalt = &ee_->thermo(pSalt);

    /*
     *  We  call the general function to calculate phase stability using the MultiPhase object.
     *   In order to do this correctly an element potential for the ion must be implemented.
     */

    /*
     * Access an element potential ThermoPhase obejct
     */
    int mp_nPhases = m_mp->nPhases();
    ThermoPhase* tp = &(m_mp->phase(mp_nPhases - 1));
    FixedChemPotSSTP* LiFixed = dynamic_cast<FixedChemPotSSTP*>(tp);
    if (LiFixed) {
        /*
         * Find the chemical potential of the Lithium
         */
        double mu[10];
        tSalt->getElectrochemPotentials(mu);
        int ilt_Lip = tSalt->speciesIndex("Li+");
        double muLip = mu[ilt_Lip];
        tElec->getElectrochemPotentials(mu);
        double muElec = mu[0];
        double muLi = muLip + muElec;
        LiFixed->setChemicalPotential(muLi);
        m_mp->setPhaseMoles(mp_nPhases - 1, 10. * solidMoles);
        //string pname = ee_->phaseName(iphase);
    }
    int loglevel = 0;
    int iStab = vcs_determine_PhaseStability(*m_mp, iphase, funcStab, printLvl_, loglevel);

    return iStab;
}
//====================================================================================================================
// Returns the mole fractions of a phase just calculated within the routine.
/*
 *
 *   @param iphase      Phase id in the list of phases within the electrode object
 *   @param moleFractions
 */
void Electrode_Equilibrium::getPhaseMoleFractions(int iphase, double* const moleFractions)
{
    int pi = PhaseIndex_mp_[iphase];
    AssertTrace(m_mp != 0);
    ThermoPhase* tp = &(m_mp->phase(pi));
    size_t ns = tp->nSpecies();
    for (size_t k = 0; k < ns; ++k) {
        moleFractions[k] = tp->moleFraction(k);
    }
}
//====================================================================================================================
// Return a complete multiphase object representing equilibrium for the electrode
/*
 * @return returns a multiphase object
 */
MP_EquilStatic* Electrode_Equilibrium::MultiPhase_Obj()
{
    return m_mp;
}
//====================================================================================================================
}// End of namespace
//----------------------------------------------------------------------------------------------------------------------------------
