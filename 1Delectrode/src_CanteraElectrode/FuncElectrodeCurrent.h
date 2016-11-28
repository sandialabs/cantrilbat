/*
 * $Id: FuncElectrodeCurrent.h 534 2013-02-22 21:33:41Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#ifndef _FUNCELECTRODECURRENT_H
#define _FUNCELECTRODECURRENT_H

//#include "cantera/numerics/ResidEval.h"
//#include "cantera/numerics/RootFind.h"

#include "m1d_SolNonlinear.h"
#include "m1d_SolNonlinear_CurrentSolve.h"
#include "m1d_ProblemResidEval.h"

#include "Electrode.h"


#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{


//! Class to solve the battery equations at a constant voltage condition, and then report the resulting
//! current as a function value
/*!
 *   This class is a child of the ResidEval class
 */
class CurrentFunc: public ZZCantera::ResidEval
{
public:

    //! Constructor
    /*!
     *
     *
     */
    CurrentFunc(m1d::SolNonlinear_CurrentSolve * solverC, m1d::SolNonlinear * solverV, m1d::ProblemResidEval *func,
                m1d::Epetra_Vector_Ghosted *y_comm, m1d::Epetra_Vector_Ghosted *ydot_comm, double CJ, double time_curr,
                double deltaT);

    ~CurrentFunc();

    int nEquations() const;

    //!  Function to calculate the current given a voltage
    /*!
     *
     */
    virtual int evalSS(const double t, const double* const x, double* const r);

    void set_deltaT(double deltaT);

    m1d::SolNonlinear_CurrentSolve* m_solver_constantCurr;
    m1d::SolNonlinear* m_solver_constantVoltage;

    m1d::ProblemResidEval *m_func_constV;
    m1d::Epetra_Vector_Ghosted *m_y_comm;
    m1d::Epetra_Vector_Ghosted *m_ydot_comm;

    double m_deltaT;

    double m_CJ;

    double m_time_curr;

    int printLvl_;

    //! used to signal the initialization of ivResultsFile_
    static int fileInit_;

    //! provides output of currents given voltage guesses so you can watch the iterative process.
    FILE *ivResultFile_;
};

}

#endif
