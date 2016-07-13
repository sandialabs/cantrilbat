/**
 * Predict the stability of a single phase that is currently zereod.
 *
 *  This routine fills in an estimate for the solution
 *  Return 1 if the phases are stable and 0 if they are not
 *
 */
/*
 * $Id: Electrode_PhaseStability.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_PHASESTABILITY_H
#define _ELECTRODE_PHASESTABILITY_H


#include "cantera/equilibrium.h"
#include "cantera/numerics/ResidEval.h"
#include "cantera/numerics/ResidJacEval.h"

#include "tok_input_util.h"

#include "ReactingSurDomain.h"

//#include "ExtraGlobalRxn.h"

#include "BlockEntry.h"

#include "Electrode_input.h"

#include <string>
#include <vector>

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

//! Class which determines the stability of phases due to kinetics
/*!
 *  Note, I believe I can make this class simpler
 *
 */
class Electrode_PhaseStability
{

public:
    //! Constructor
    /*!
     *
     * @param elect  Electrode object pertaining to the stability problem
     *               The electrode object assigns this object as a "friend"
     */
    Electrode_PhaseStability(Electrode_MultiPlateau_NoDiff* elect);

    //! Destructor
    virtual ~Electrode_PhaseStability();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_PhaseStability(const Electrode_PhaseStability& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_PhaseStability& operator=(const Electrode_PhaseStability& right);

    void setup(const std::vector<int>& phasePopIndexList);

    int determinePhaseStability(doublereal& retnFunc);

    std::vector<std::vector<double> >& moleFractions();

    int  nResidEquations() const;

    //! Unpack the soln vector
    /*!
     *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
     */
    int unpackNonlinSolnVector(const double* const y);

    int packNonlinSolnVector(double* const y);

    void  seedMfVector();
    void updatePhaseMoleFractions(int iii, int iph);

    void extractInfo();

    //! Residual Evaluation program
    /*!
     *
     */
    int optResid(const doublereal tdummy, const doublereal delta_t_dummy,
                 const doublereal* const y,
                 const doublereal* const ySolnDot,
                 doublereal* const resid,
                 const ResidEval_Type_Enum evalType,
                 const int id_x,
                 const doublereal delta_x);

    void determineBigMoleFractions();

    //! Nonlinear solver functional
    /*!
     *   This is the interface to the nonlinear solver for advancing the ODE integrals that this
     *   object uses.  This interface is mostly a self-directed passthrough to routines that
     *   are members of this object.
     */
    class calcPhaseStabFunc_ResidJacEval : public ZZCantera::ResidJacEval
    {
    public:

        //! Constructor of the functional
        /*!
         *  @param elec pointer to the electrode model. This is usually the current Electrode class
         */
        calcPhaseStabFunc_ResidJacEval(Electrode_PhaseStability* eps);

        //! Evaluate the residual function
        /*!
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
        virtual int evalResidNJ(const doublereal t, const doublereal delta_t,
                                const doublereal* const y,
                                const doublereal* const ydot,
                                doublereal* const resid,
                                const ResidEval_Type_Enum evalType = Base_ResidEval,
                                const int id_x = -1,
                                const doublereal delta_x = 0.0);

        //! Fill in the initial conditions
        /*!
         * Values for both the solution and the value of ydot may be provided.
         *
         * @param t0            Time                    (input)
         * @param y             Solution vector (output)
         * @param ydot          Rate of change of solution vector. (output)
         */
        virtual int getInitialConditions(const doublereal t0, doublereal* const y,
                                         doublereal* const ydot);


        //! Return the number of equations in the equation system
        virtual  int nEquations() const;


        //! Apply a filtering process to the step
        /*!
         *  @param timeCurrent    Current value of the time
         *  @param ybase          current value of the solution
         *  @param step0          Value of the step in the solution vector that will be filtered.
         *                        The filter is applied to the step values.
         *
         *  @return Returns the norm of the value of the amount filtered
         */
        virtual doublereal filterNewStep(const doublereal timeCurrent, const doublereal* const ybase,
                                         doublereal* const step0);

    protected:

        //! This is a self-reference so that the routine can call member functions of itself
        Electrode_PhaseStability* eps_;
    };

protected:

    //! This is a reference to the friend object, where we will pull most of the data from.
    Electrode_MultiPlateau_NoDiff* emp_;

    //! Residual function that loops back.
    calcPhaseStabFunc_ResidJacEval* m_resid;

    int neq_;

    double fValue_;


    //! Number of phases to check
    int nPhasesToPop_;

    //! Vector of indexes of phases to check in the phase pop problem
    std::vector<int> phasePopIndexList_;

    //! Vector of indexes of phases to check in the phase pop problem
    std::vector<int> phasePopSolidIndexList_;

    std::vector<std::string> phasePopNames_;


    std::vector<std::vector<double> > mfVector_pl_;

    std::vector<std::vector<double> > ZVector_pl_;

    std::vector<std::vector<double> > CDotVector_pl_;
    std::vector<std::vector<double> > DDotVector_pl_;
    std::vector<std::vector<double> > DDotLinVector_pl_;

    std::vector<double > fValue_pl_;

    std::vector<double > fValue_pl_lagged_;

    std::vector<int> phaseMFBig_;

    std::vector<doublereal> residVector_;

    std::vector<doublereal> solnVector_;

    std::vector<doublereal> speciesCreationRatesElectrode_;
    std::vector<doublereal> speciesDestructionRatesElectrode_;


    std::vector<doublereal> speciesCreationRatesReactingSurface_;
    std::vector<doublereal> speciesDestructionRatesReactingSurface_;

    std::vector<doublereal> atol_;
    std::vector<doublereal> ylow_;
    std::vector<doublereal> yhigh_;
    std::vector<doublereal> yval_;

    SquareMatrix* jac_;
    ZZCantera::NonlinearSolver* pSolve_;

    std::vector<doublereal> phaseMoles_final_Orig_;
    int printLvl_;
public:
    int enableExtraPrinting_;
    int detailedResidPrintFlag_;
};

}
#endif
/*****************************************************************************/

