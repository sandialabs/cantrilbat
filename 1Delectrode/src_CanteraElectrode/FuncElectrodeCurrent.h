/**
 *  @file funcElectrodeCurrent.h
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#ifndef _FUNCELECTRODECURRENT_H
#define _FUNCELECTRODECURRENT_H

#include "m1d_SolNonlinear.h"
#include "m1d_SolNonlinear_CurrentSolve.h"
#include "m1d_ProblemResidEval.h"

#include "Electrode.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace Zuzax
{

//==================================================================================================================================
//! Class to solve the battery equations at a constant voltage condition, and then report the resulting current as a function value
/*!
 *   This class is a child of the ResidEval class
 */
class CurrentFunc: public Zuzax::ResidEval
{
public:

    //! Constructor
    /*!
     *
     *
     */
    CurrentFunc(m1d::SolNonlinear_CurrentSolve* solverC, m1d::SolNonlinear* solverV, m1d::ProblemResidEval* func,
                m1d::Epetra_Vector_Ghosted* y_comm, m1d::Epetra_Vector_Ghosted* ydot_comm, double CJ, double time_curr,
                double deltaT);

    //! Destructor
    ~CurrentFunc();

    //! Return the number of equations
    /*!
     *  @return                                  Returns a value of 1
     */
    virtual int nEquations() const override;

    //! Evaluate the steady state residual function at time t
    /*!
     *  (virtual from ResidEval)
     *
     *  @param[in]           t                   Time
     *  @param[in]           y                   Solution vector
     *  @param[out]          resid               Residual vector
     *  @param[in]           residType           Residual type. The default is Base_ResidEval.
     *  @param[in]           solveType           Type of the problem being solved expressed as a Solve_Type_Enum. 
     *                                           Defaults to Solve_Type::SteadyState_Solve
     *
     *  @return                                  Returns a flag to indicate that operation is successful.
     *                                           -  1 = ZZ_RESIDEVAL_SUCCESS  Means a successful operation
     *                                           - -1 = ZZ_RESIDEVAL_FAILURE  or neg value means an unsuccessful operation
     */
    virtual int evalResidSS(const doublevalue t, const doublevalue* const y, doublevalue* const resid,
                            const ResidEval_Type residType = ResidEval_Type::Base_ResidEval,
                            const Solve_Type solveType = Solve_Type::SteadyState_Solve);


    void set_deltaT(double deltaT);

    m1d::SolNonlinear_CurrentSolve* m_solver_constantCurr;
    m1d::SolNonlinear* m_solver_constantVoltage;

    m1d::ProblemResidEval* m_func_constV;

    //! Current value of the solution vector
    m1d::Epetra_Vector_Ghosted* m_y_comm;

    //! Current value of the solution Dot vector
    m1d::Epetra_Vector_Ghosted* m_ydot_comm;

    //! Current value of deltaT
    double m_deltaT;

    //! multiplier before the time derivative terms
    double m_CJ;

    //! Current time
    double m_time_curr;

    //! Print level
    int printLvl_;

    //! used to signal the initialization of ivResultsFile_
    static int fileInit_;

    //! provides output of currents given voltage guesses so you can watch the iterative process.
    FILE* ivResultFile_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
