/*
 * @file Electrode_PhaseStability.cpp 
 */


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"
#include <iomanip>

#include "Electrode_MultiPlateau_NoDiff.h"
#include "zuzax/numerics/NonlinearSolver.h"
#include "zuzax/thermo/FixedChemPotSSTP.h"

#include "Electrode_PhaseStability.h"

using namespace std;

#ifndef SAFE_DELETE
#define SAFE_DELETE(x)  if (x) { delete x;  x = 0;}
#endif

namespace Zuzax
{
//====================================================================================================================
Electrode_PhaseStability::Electrode_PhaseStability(Electrode_MultiPlateau_NoDiff* elect) :
    emp_(elect),
    m_resid(0),
    neq_(0),
    fValue_(0.0),
    nPhasesToPop_(0),
    phasePopNames_(0),
    jac_(0),
    pSolve_(0),
    printLvl_(0),
    enableExtraPrinting_(0),
    detailedResidPrintFlag_(0)
{
    // m_resid = new calcPhaseStabFunc_ResidJacEval(this);
}

//====================================================================================================================
Electrode_PhaseStability::Electrode_PhaseStability(const Electrode_PhaseStability& right) :
    emp_(right.emp_),
    m_resid(0),
    neq_(0),
    fValue_(0.0),
    nPhasesToPop_(0),
    phasePopNames_(0),
    jac_(0),
    pSolve_(0),
    printLvl_(0),
    enableExtraPrinting_(0),
    detailedResidPrintFlag_(0)
{
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
}
//==================================================================================================================================
Electrode_PhaseStability::~Electrode_PhaseStability()
{
    SAFE_DELETE(m_resid);
    SAFE_DELETE(jac_);
    SAFE_DELETE(pSolve_);
}
//==================================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_PhaseStability& Electrode_PhaseStability::operator=(const Electrode_PhaseStability& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    /*
     *  Do a shallow copy of the Electrode pointer. This is all that is necessary
     */
    emp_ = right.emp_;

    SAFE_DELETE(m_resid);
    m_resid = new calcPhaseStabFunc_ResidJacEval(this);

    neq_ = right.neq_;
    fValue_ = right.fValue_;
    nPhasesToPop_ = right.nPhasesToPop_;
    phasePopIndexList_ = right.phasePopIndexList_;
    phasePopSolidIndexList_ = right.phasePopSolidIndexList_;
    phasePopNames_ = right.phasePopNames_;
    mfVector_pl_ = right.mfVector_pl_;
    ZVector_pl_ = right.ZVector_pl_;
    CDotVector_pl_ = right.CDotVector_pl_;
    DDotVector_pl_ = right.DDotVector_pl_;
    DDotLinVector_pl_ = right.DDotLinVector_pl_;
    fValue_pl_ = right.fValue_pl_;
    fValue_pl_lagged_ = right.fValue_pl_lagged_;
    phaseMFBig_ = right.phaseMFBig_;
    residVector_ = right.residVector_;
    solnVector_ = right.solnVector_;
    speciesCreationRatesElectrode_ = right.speciesCreationRatesElectrode_;
    speciesDestructionRatesElectrode_ = right.speciesDestructionRatesElectrode_;
    speciesCreationRatesReactingSurface_ = right.speciesCreationRatesReactingSurface_;
    speciesDestructionRatesReactingSurface_ = right.speciesDestructionRatesReactingSurface_;
    atol_= right.atol_;
    ylow_ = right.ylow_;
    yhigh_ = right.yhigh_;
    yval_ = right.yval_;

    SAFE_DELETE(jac_);
    jac_ = new SquareMatrix(*right.jac_);
    SAFE_DELETE(pSolve_);
    pSolve_ = new NonlinearSolver_JAC(*right.pSolve_);

    phaseMoles_final_Orig_ = right.phaseMoles_final_Orig_;
    printLvl_  = right.printLvl_;
    enableExtraPrinting_ = right.enableExtraPrinting_;
    detailedResidPrintFlag_ = right.detailedResidPrintFlag_;

    /*
     * Return the reference to the current object
     */
    return *this;
}
//==================================================================================================================================
void Electrode_PhaseStability::setup(const std::vector<int>& phasePopIndexList)
{
    bool found = false;
    nPhasesToPop_ = phasePopIndexList.size();
    phasePopSolidIndexList_.clear();
    mfVector_pl_.clear();
    ZVector_pl_.clear();
    phaseMFBig_.clear();
    phasePopNames_.clear();
    phasePopIndexList_ = phasePopIndexList;
    neq_ = 0;
    for (int iii = 0; iii < nPhasesToPop_; iii++) {
        for (int np = 0; np < (int) emp_->phaseIndexSolidPhases_.size(); np++) {
            int iph = emp_->phaseIndexSolidPhases_[np];
            if (phasePopIndexList_[iii] == iph) {
                found = true;
                phasePopSolidIndexList_.push_back(np);
                ThermoPhase* tp = emp_->VolPhaseList[iph];
                int nsp = tp->nSpecies();
                neq_ += nsp - 1;
                std::vector<double> mole_fractions(nsp);
                int iStart = emp_->globalSpeciesIndex(iph, 0);

                std::copy(emp_->spMf_final_[iStart], emp_->spMf_final_[iStart] + nsp, mole_fractions.begin());
                mfVector_pl_.push_back(mole_fractions);
                ZVector_pl_.push_back(mole_fractions);
                CDotVector_pl_.push_back(mole_fractions);
                DDotVector_pl_.push_back(mole_fractions);
                DDotLinVector_pl_.push_back(mole_fractions);
                std::string sss = tp->id();
                phasePopNames_.push_back(sss);

                int bigone = 0;
                double btmp = mole_fractions[0];
                for (int isp = 1; isp < nsp; isp++) {
                    if (mole_fractions[isp] > btmp) {
                        bigone = isp;
                        btmp = mole_fractions[isp];
                    }
                }
                phaseMFBig_.push_back(bigone);

                residVector_.resize(neq_, 0.0);
                solnVector_.resize(neq_, 0.0);
                break;
            }

        }
        if (!found) {
            throw Electrode_Error(" Electrode_PhaseStability", "error");
        }
    }

    fValue_pl_.resize(nPhasesToPop_, 0.5);
    fValue_pl_lagged_.resize(nPhasesToPop_, 0.5);

    /*
     *  Resize arrays based on the number of species in the base Electrode object
     */
    int nsp = emp_->m_NumTotSpecies;
    speciesCreationRatesElectrode_.resize(nsp, 0.0);
    speciesDestructionRatesElectrode_.resize(nsp, 0.0);
    speciesCreationRatesReactingSurface_.resize(nsp, 0.0);
    speciesDestructionRatesReactingSurface_.resize(nsp, 0.0);

    /*
     * Resize arrays based on the number of species in the nonlinear problem
     */
    atol_.resize(neq_, 1.0E-50);
    ylow_.resize(neq_, 0.0);
    yhigh_.resize(neq_, 1.0);
    yval_.resize(neq_, 0.0);

    phaseMoles_final_Orig_.resize(nPhasesToPop_, 0.0);

    jac_ = new SquareMatrix(neq_, 0.0);
    m_resid = new calcPhaseStabFunc_ResidJacEval(this);
    pSolve_ = new NonlinearSolver_JAC(m_resid);
    pSolve_->setAtol(DATA_PTR(atol_));
    pSolve_->setRtol(1.0E-5);
    pSolve_->setBoundsConstraints(&ylow_[0], &yhigh_[0]);
}
//==================================================================================================================================
int Electrode_PhaseStability::determinePhaseStability(double& retnFunc)
{
    printLvl_ = 9;
    /*
     *  In order to get correct rates of progress for phases, it's necessary to put some mass into them
     *  in order to get them to potentially to go backwards. We then have to fix this later.
     */
    for (int iii = 0; iii <  nPhasesToPop_; iii++) {
        int iph = phasePopIndexList_[iii];
        phaseMoles_final_Orig_[iii] = emp_->phaseMoles_final_[iph];
        if (emp_->phaseMoles_final_[iph] <= 0.0) {
            emp_->phaseMoles_final_[iph] = 1.0;
        }
        std::vector<double>& mfVector = mfVector_pl_[iii];
        /*
         *  The default value didn't converge. However, these initial conditions converged robustly.
         *  Therefore, we have a problem with this method. It seems to be very dependent on the initial conditions
         *  of the solution. I tried a bunch of methods but nothing made it robust. These initial conditions
         *  came from the equilibrium solve of the system.
         */
        if (iii == 0) {
            mfVector[0] = 0.34;
            mfVector[1] = 0.66;
        }
        if (iii == 1) {
            mfVector[1] = 6.55E-2;
            mfVector[0] = 1.0 - mfVector[1];
        }
    }
    emp_->setPhaseExistenceForReactingSurfaces(true);

    packNonlinSolnVector(DATA_PTR(yval_));
    double time_curr = 0.0;
    int num_newt_its = 0;
    int num_linear_solves = 0;
    int numBacktracks = 0;
    double* ydot = 0;
    pSolve_->setRtol(1.0E-8);
    int nonlinearFlag = pSolve_->solve_nonlinear_problem(NSOLN_TYPE_STEADY_STATE, &yval_[0], ydot, 0.0,
                        time_curr, *jac_,  num_newt_its, num_linear_solves,
                        numBacktracks, printLvl_);
    if (nonlinearFlag < 0) {
        printf(" Electrode_PhaseStability::determinePhaseStability():  Unsuccessful Nonlinear Solve\n");
        exit(-1);
    }

    for (int iii = 0; iii <  nPhasesToPop_; iii++) {
        int iph = phasePopIndexList_[iii];
        emp_->phaseMoles_final_[iph] = phaseMoles_final_Orig_[iii];
    }
    emp_->setPhaseExistenceForReactingSurfaces(true);

    retnFunc = fValue_;

    if (fValue_ < 1.0) {
        return 1;
    }
    return 0;
}
//==================================================================================================================================
std::vector<std::vector<double> >& Electrode_PhaseStability::moleFractions()
{
    return mfVector_pl_;
}
//==================================================================================================================================
int Electrode_PhaseStability::nResidEquations() const
{
    return neq_;
}
//==================================================================================================================================
// Unpack the soln vector
/*
 *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
 */
int Electrode_PhaseStability::unpackNonlinSolnVector(const double* const y)
{
    int index = 0;
    for (int iii = 0; iii < nPhasesToPop_; iii++) {
        int ph = phasePopSolidIndexList_[iii];
        //int iph = phasePopIndexList_[iii];
        std::vector<double>& mfVector =  mfVector_pl_[iii];
        // int isp = emp_->globalSpeciesIndex(iph, 0);
        if (emp_->numSpecInSolidPhases_[ph] == 1) {
            mfVector[0] = 1.0;
        }  else {
            double BigMF = 1.0;
            for (int sp = 0; sp < emp_->numSpecInSolidPhases_[ph]; sp++) {
                if (sp != phaseMFBig_[iii]) {
                    mfVector[sp] = y[index];
                    BigMF -= mfVector[sp];
                    index++;
                }
            }
            mfVector[phaseMFBig_[iii]] = BigMF;
        }
    }

    if (index != neq_) {
        throw Electrode_Error(" Electrode_PhaseStability::unpackNonlinSolnVector", "nsp error");
    }
    return 0;
}
//==================================================================================================================================
// Pack the soln vector
/*
 *  This function packs the solution vector
 */
int Electrode_PhaseStability::packNonlinSolnVector(double* const y)
{
    int index = 0;
    for (int iii = 0; iii < nPhasesToPop_; iii++) {
        int ph = phasePopSolidIndexList_[iii];
        //int iph = phasePopIndexList_[iii];
        std::vector<double>& mfVector =  mfVector_pl_[iii];
        // int isp = emp_->globalSpeciesIndex(iph, 0);
        if (emp_->numSpecInSolidPhases_[ph] == 1) {
            mfVector[0] = 1.0;
        }  else {
            double BigMF = 1.0;
            for (int sp = 0; sp < emp_->numSpecInSolidPhases_[ph]; sp++) {
                if (sp != phaseMFBig_[iii]) {
                    y[index] =  mfVector[sp];
                    BigMF -= mfVector[sp];
                    index++;
                }
            }
            mfVector[phaseMFBig_[iii]] = BigMF;
        }
    }

    if (index != neq_) {
        throw Electrode_Error(" Electrode_PhaseStability::packNonlinSolnVector", "nsp error");
    }
    return 0;
}
//==================================================================================================================================
// Seed the soln vector
void Electrode_PhaseStability::seedMfVector()
{
    for (int iii = 0; iii < nPhasesToPop_; iii++) {
        int ph = phasePopSolidIndexList_[iii];
        int iph = phasePopIndexList_[iii];
        std::vector<double>& mfVector =  mfVector_pl_[iii];
        int iStart = emp_->globalSpeciesIndex(iph, 0);

        for (int sp = 0; sp < emp_->numSpecInSolidPhases_[ph]; sp++) {
            mfVector[sp] = emp_->spMf_final_[iStart+sp];
        }
    }
}
//==================================================================================================================================
void Electrode_PhaseStability::updatePhaseMoleFractions(int iii, int iph)
{

    int istart = emp_->m_PhaseSpeciesStartIndex[iph];
    ThermoPhase& tp = emp_->thermo(iph);
    int nsp =  emp_->m_PhaseSpeciesStartIndex[iph+1] - istart;

    std::vector<double>& mfVector = mfVector_pl_[iii];
    /*
     * Install a new version of the mole fractions into the electrode objects mole fraction vector
     */
    for (int k = 0; k < nsp; k++) {
        emp_->spMf_final_[istart + k] = mfVector[k];
    }
    tp.setState_TPX(emp_->temperature_, emp_->pressure_, &(emp_->spMf_final_[istart]));
    tp.setElectricPotential(emp_->phaseVoltages_[iph]);
    tp.getPartialMolarVolumes(&(emp_->VolPM_[istart]));
    tp.getElectrochemPotentials(&(emp_->spElectroChemPot_[istart]));

}
//==================================================================================================================================
/*
 *  Calculate CDot, DDot, DDotLin
 */
void Electrode_PhaseStability::extractInfo()
{
    int iph, jph, kph;

#ifdef DEBUG_HKM
    int maxNumRxns = 20;
    std::vector<double> netROP(maxNumRxns, 0.0);
    std::vector<double> fwdROP(maxNumRxns, 0.0);
    std::vector<double> revROP(maxNumRxns, 0.0);
    std::vector<double> deltaGrxn(maxNumRxns, 0.0);
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printf("\t\t  extractInfo() detailed information:\n");
    }
#endif

    /*
     *  First thing we need to do is to set the mole fractions in the thermophases
     */
    for (int iii = 0; iii < nPhasesToPop_; iii++) {
        iph = phasePopIndexList_[iii];
        updatePhaseMoleFractions(iii, iph);
    }



    std::fill(speciesCreationRatesElectrode_.begin(), speciesCreationRatesElectrode_.end(), 0.);
    std::fill(speciesDestructionRatesElectrode_.begin(), speciesDestructionRatesElectrode_.end(), 0.);


    // for (int isk = 0; isk < emp_->m_NumSurPhases; isk++) {
    for (int isk = 0; isk < 2; isk++) {

        double surfaceArea = emp_->surfaceAreaRS_init_[isk];

        bool active = emp_->ActiveKineticsSurf_[isk];
        if (active) {
            ReactingSurDomain* rsd = emp_->RSD_List_[isk];

            rsd->getCreationRates(DATA_PTR(speciesCreationRatesReactingSurface_));
            rsd->getDestructionRates(DATA_PTR(speciesDestructionRatesReactingSurface_));

#ifdef DEBUG_HKM
            if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
                rsd->getNetRatesOfProgress(DATA_PTR(netROP));
                rsd->getFwdRatesOfProgress(DATA_PTR(fwdROP));
                rsd->getRevRatesOfProgress(DATA_PTR(revROP));
                rsd->getDeltaElectrochemPotentials(DATA_PTR(deltaGrxn));
                printf("\t\t   Outer surface ROP   %d  surfaceArea = %10.5g                                "
                       "        ROP      fwdROP        revROP     DeltaElectroChem\n", isk, surfaceArea);
                for (int i = 0; i < emp_->numRxns_[isk]; i++) {
                    std::string ss = rsd->reactionString(i);
                    printf("\t\t    %-80.80s % -12.4e  % -12.4e  % -12.4e  % -12.4e\n", ss.c_str(), netROP[i],
                           fwdROP[i], revROP[i], deltaGrxn[i]);
                }
            }
#endif
            /*
             *  Get the number of phases in the reacting surface and loop over them
             */
            int nphRS = rsd->nPhases();
            int kIndexKin = 0;
            for (kph = 0; kph < nphRS; kph++) {
                /*
                 * Find the phase index in the Electrode Object
                 */
                jph = rsd->kinOrder[kph];
                int istart = emp_->m_PhaseSpeciesStartIndex[jph];
                int nsp = emp_->m_PhaseSpeciesStartIndex[jph+1] - istart;
                for (int k = 0; k < nsp; k++) {
                    speciesCreationRatesElectrode_[istart+k]    += speciesCreationRatesReactingSurface_[kIndexKin] * surfaceArea;
                    speciesDestructionRatesElectrode_[istart+k] += speciesDestructionRatesReactingSurface_[kIndexKin] * surfaceArea;
                    kIndexKin++;
                }
            }


        }

    }

    for (int iii = 0; iii < nPhasesToPop_; iii++) {
        iph = phasePopIndexList_[iii];

        std::vector<double>& CDotVector = CDotVector_pl_[iii];
        std::vector<double>& DDotLinVector = DDotLinVector_pl_[iii];
        std::vector<double>& DDotVector = DDotVector_pl_[iii];
        std::vector<double>& mfVector =     mfVector_pl_[iii];

        int istart = emp_->globalSpeciesIndex(iph, 0);

        int nsp = emp_->m_PhaseSpeciesStartIndex[iph+1] - istart;
        for (int k = 0; k < nsp; k++) {
            CDotVector[k] = speciesCreationRatesElectrode_[istart+k];
            DDotVector[k] = speciesDestructionRatesElectrode_[istart+k];
            if (mfVector[k] > 0.0) {
                DDotLinVector[k] = speciesDestructionRatesElectrode_[istart+k] / mfVector[k];
            } else {
                if (speciesDestructionRatesElectrode_[istart+k] > 0.0) {
                    printf("Unknown situation\n");
                    DDotLinVector[k] = speciesDestructionRatesElectrode_[istart+k];
                } else {
                    DDotLinVector[k] = 0.0;
                }
            }
        }


        if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
            std::string sss = phasePopNames_[iii];
            ThermoPhase tp = emp_->thermo(iph);
            printf("\t\t   Phase_%-15.15s:    "
                   "    X            NetROP       Creation     Destruction  Dest/X       Z\n", sss.c_str());
            for (int k = 0; k < nsp; k++) {
                sss = tp.speciesName(k);
                double val = 0.0;
                if (DDotLinVector[k] > 0.0) {
                    val =  CDotVector[k]/DDotLinVector[k];
                }
                printf("\t\t           %-15.15s     % -12.4e % -12.4e % -12.4e % -12.4e % -12.4e % -12.4e \n",
                       sss.c_str(), mfVector[k], CDotVector[k] - DDotVector[k], CDotVector[k], DDotVector[k], DDotLinVector[k],
                       val);
            }
        }

    }

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
int Electrode_PhaseStability::optResid(const double tdummy, const double delta_t_dummy,
                                       const double* const y,
                                       const double* const ySolnDot,
                                       double* const resid,
                                       const ResidEval_Type_Enum evalType,
                                       const int id_x,
                                       const double delta_x)
{
    int iph, isp;

    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printf("\t\t================================================================================================="
               "==============================\n");
        printf("\t\t  EXTRA PRINTING FROM NONLINEAR RESIDUAL: ");
        if (evalType == Base_ResidEval || evalType == Base_LaggedSolutionComponents) {
            printf(" BASE RESIDUAL");
        } else if (evalType == JacBase_ResidEval) {
            printf(" BASE JAC RESIDUAL");
        } else  if (evalType == JacDelta_ResidEval) {
            printf(" DELTA JAC RESIDUAL");
            printf(" var = %d delta_x = %12.4e Y_del = %12.4e Y_base = %12.4e", id_x, delta_x, y[id_x], y[id_x] - delta_x);
        } else  if (evalType == Base_ShowSolution) {
            printf(" BASE RESIDUAL - SHOW SOLUTION");
        }
        printf("\n");
    }
    /*
     *  Current the solution vector are the nsp -1 mole fractions in multiple phases
     */
    unpackNonlinSolnVector(y);


    extractInfo();

    for (int iii = 0; iii < nPhasesToPop_; iii++) {
        iph = phasePopIndexList_[iii];
        std::vector<double>& ZVector = ZVector_pl_[iii];
        std::vector<double>& CDotVector = CDotVector_pl_[iii];
        std::vector<double>& DDotLinVector = DDotLinVector_pl_[iii];

        /*
         * Calculate ZVector and f
         */
        int ph = phasePopSolidIndexList_[iii];
        int nsp = emp_->numSpecInSolidPhases_[ph];

        fValue_pl_[iii] = 0.0;
        for (isp = 0; isp < nsp; isp++) {
            ZVector[isp] = CDotVector[isp] / DDotLinVector[isp];
            fValue_pl_[iii]  += ZVector[isp];
        }
    }


    if (evalType != JacDelta_ResidEval && (evalType != Base_LaggedSolutionComponents)) {
        std::vector<double>::const_iterator begin = fValue_pl_.begin();
        std::copy(begin, begin+nPhasesToPop_, fValue_pl_lagged_.begin());
    }

    int index = 0;
    for (int iii = 0; iii < nPhasesToPop_; iii++) {
        iph = phasePopIndexList_[iii];
        std::vector<double>& ZVector = ZVector_pl_[iii];
        std::vector<double>& CDotVector = CDotVector_pl_[iii];
        std::vector<double>& DDotVector = DDotVector_pl_[iii];
        //std::vector<double> &DDotLinVector = DDotLinVector_pl_[iii];

        /*
         * Calculate ZVector and f
         */
        int ph = phasePopSolidIndexList_[iii];
        int nsp = emp_->numSpecInSolidPhases_[ph];

        std::vector<double>& mfVector =  mfVector_pl_[iii];
        if (emp_->numSpecInSolidPhases_[ph] == 1) {
            mfVector[0] = 1.0;
        }  else {
            for (isp = 0; isp < nsp; isp++) {
                if (isp != phaseMFBig_[iii]) {

                    // This appears to be the best method. It's not very good but it is better than the ones below.
                    //  - it appears that this has the highest range of convergence
                    resid[index] = mfVector[isp] * fValue_pl_[iii] - ZVector[isp];

                    // This method never converged
                    // resid[index] = mfVector[isp] * fValue_pl_lagged_[iii] - ZVector[isp];

                    // The method below appears to have worse convergence than the method above.
                    // resid[index] = DDotVector[isp] * fValue_pl_[iii] - CDotVector[isp];

                    // The method below appears to have worse convergence than the method above.
                    //resid[index] = DDotVector[isp] * fValue_pl_lagged_[iii] - CDotVector[isp];

                    index++;
                }
            }
        }

        if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
            std::string sss =  phasePopNames_[iii];
            double net = 0.0;
            double fnet = 0.0;
            double rnet = 0.0;
            for (isp = 0; isp < nsp; isp++) {
                net += CDotVector[isp] - DDotVector[isp];
                fnet += CDotVector[isp];
                rnet += DDotVector[isp];
            }
            index -= (nsp - 1);
            ThermoPhase tp = emp_->thermo(iph);
            printf("\t\t       PhaseName                        Net         Creation     Destruction |   f_value - 1.0\n");
            printf("\t\t         %-15.15s              % -12.4e % -12.4e % -12.4e | % -12.4e    \n",
                   sss.c_str(), net, fnet, rnet , fValue_pl_[iii] - 1.0);
            printf("\t\t       SpeciesName              X       Net         Creation     Destruction"
                   " |    Zcalc        ZMeas     |     Resid    |\n");
            for (isp = 0; isp < nsp; isp++) {
                sss = tp.speciesName(isp);
                printf("\t\t         %-15.15s % -12.4e % -12.4e % -12.4e % -12.4e | % -12.4e % -12.4e ",
                       sss.c_str(),  mfVector[isp], CDotVector[isp] -  DDotVector[isp],   CDotVector[isp], DDotVector[isp],
                       mfVector[isp] * fValue_pl_[iii] , ZVector[isp]);
                if (isp != phaseMFBig_[iii]) {
                    printf("| % -12.4e |", resid[index]);
                    index++;
                } else {
                    printf("|              |");
                }
                printf("\n");
            }

        }


    }

    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
        printf("\t\t=============================================================="
               "=================================================================\n");
    }

    return 1;
}
//==================================================================================================================================
// Constructor of the functional
/*
 *  @param elec pointer to the electrode model. This is usually the current Electrode class
 */
Electrode_PhaseStability::calcPhaseStabFunc_ResidJacEval::calcPhaseStabFunc_ResidJacEval(Electrode_PhaseStability* eps) :
    ResidJacEval(),
    eps_(eps)
{
    neq_ = eps_->nResidEquations();
}
//==================================================================================================================================
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
int Electrode_PhaseStability::calcPhaseStabFunc_ResidJacEval::
evalResidNJ(const double tdummy, const double delta_t_dummy,
            const double* const y,
            const double* const ySolnDot,
            double* const resid,
            const ResidEval_Type_Enum evalType,
            const int id_x,
            const double delta_x)
{
    /*
     *  Return control for calculating the residual to the controlling Electrode_PhaseStability object.
     *  This object will calculate the residual at the current conditions
     */
    int retn =  eps_->optResid(tdummy, delta_t_dummy, y, ySolnDot, resid, evalType, id_x, delta_x);
    return retn;
}
//  -----------------------------------------------------------------------------------------------------------------
int Electrode_PhaseStability::calcPhaseStabFunc_ResidJacEval::
getInitialConditions(const double t0, double* const y, double* const ydot)
{
    for (int k = 0; k < neq_; k++) {
        y[k] = 0.0;
    }
    return 1;
}
//  -----------------------------------------------------------------------------------------------------------------
int Electrode_PhaseStability::calcPhaseStabFunc_ResidJacEval::nEquations() const
{
    return neq_;
}
//==================================================================================================================================
//  Apply a filtering process to the step
/*
 *  @param timeCurrent    Current value of the time
 *  @param ybase          current value of the solution
 *  @param step0          Value of the step in the solution vector that will be filtered.
 *                        The filter is applied to the step values.
 *
 *  @return Returns the norm of the value of the amount filtered
 */
double Electrode_PhaseStability::calcPhaseStabFunc_ResidJacEval::
filterNewStep(const double timeCurrent, const double* const ybase,
              double* const step0)
{
    return 0.0;
}
//==================================================================================================================================
//   Determine the big mole fraction in the phase
void  Electrode_PhaseStability::determineBigMoleFractions()
{
    for (int iii = 0; iii < nPhasesToPop_; iii++) {
        int iph = phasePopIndexList_[iii];

        std::vector<double>& mfVector =     mfVector_pl_[iii];

        size_t istart = emp_->globalSpeciesIndex(iph, 0);

        int nsp = emp_->m_PhaseSpeciesStartIndex[iph+1] - istart;
        phaseMFBig_[iii] = 0;

        double xBig = mfVector[0];
        for (int sp = 1; sp < nsp; sp++) {
            if (mfVector[sp] > xBig) {
                phaseMFBig_[iii] = sp;
                xBig = mfVector[sp];
            }
        }
    }
}
//==================================================================================================================================
} 
//==================================================================================================================================

