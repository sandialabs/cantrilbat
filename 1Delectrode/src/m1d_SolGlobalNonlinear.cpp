/**
 *  @file m1d_SolNonlinear.cpp
 *  Damped Newton solver for 1D multi-domain problems
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_SolGlobalNonlinear.h"

using namespace std;
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Wrapper error class to denote a parent class that has not been implemented
class errBC : public m1d_Error
{
public:
    //! Constructor
    /*!
     *  @param[in]           procedure           Name of the procedure that is called
     */
    errBC(std::string procedure) :
        m1d_Error(procedure, "Base Class SolGlobalNonlinear called")
    {
    }
};
//==================================================================================================================================
SolGlobalNonlinear::SolGlobalNonlinear() :
    m_print_flag(3)
{
}
//==================================================================================================================================
SolGlobalNonlinear::~SolGlobalNonlinear()
{
}
//==================================================================================================================================
double
SolGlobalNonlinear::soln_error_norm(const Epetra_Vector& delta_y, const bool printLargest,
                                    const char* const title, const int typeYsoln, const double dampFactor) const
{
    throw errBC("soln_error_norm()");
    return 0.0;
}
//==================================================================================================================================
double
SolGlobalNonlinear::res_error_norm(const Epetra_Vector_Owned& resid, const char* title, const int printLargest) const
{

    throw errBC("res_error_norm()");
    return 0.0;
}
//==================================================================================================================================
void
SolGlobalNonlinear::get_res(const double time_curr, const double rdelta_t,
                            const Epetra_Vector_Ghosted* const solnBase_ptr, const Epetra_Vector_Ghosted* const solnDotBase_ptr)
{
    throw errBC("scaleMatrix()");
}

//==================================================================================================================================
void
SolGlobalNonlinear::scaleMatrix(Epetra_Vector_Owned& delta_soln, const Epetra_Vector_Ghosted& y_curr,
                                const Epetra_Vector_Ghosted& ydot_curr, const double time_curr,
                                const double rdelta_t, int loglevel)
{
    throw errBC("scaleMatrix()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::doNewtonSolve(Epetra_Vector_Owned& delta_soln, const Epetra_Vector_Ghosted& y_curr,
                                  const Epetra_Vector_Ghosted& ydot_curr, const double time_curr,
                                  const double rdelta_t, int loglevel)
{
    throw errBC("doNewtonSolve()");
}
//==================================================================================================================================
int
SolGlobalNonlinear::doHardBounds(const Epetra_Vector_Ghosted& y_old, Epetra_Vector_Owned& step, double& fbound)
{
    throw errBC("doHardBounds()");
    return 1;
}
//==================================================================================================================================
double
SolGlobalNonlinear::deltaBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0)
{
    throw errBC("deltaBoundStep()");
    return 0.0;
}
//==================================================================================================================================
double
SolGlobalNonlinear::highLowBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0, int loglevel)
{
    throw errBC("highLowBoundStep()");
    return 1.0;
}
//==================================================================================================================================
void
SolGlobalNonlinear::setDefaultColumnScaleVector()
{
    throw errBC("setDefaultColumnScaleVector()");
}
//==================================================================================================================================
bool
SolGlobalNonlinear::getColumnScaleVector(Epetra_Vector_Owned& colScales) const
{
    throw  errBC("getColumnScaleVector()");
    return true;
}
//==================================================================================================================================
void
SolGlobalNonlinear::setColumnScaleVector(const Epetra_Vector_Owned& colScales)
{
    throw errBC("setColumnScaleVector()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setup_problem(Solve_Type_Enum solveType, const Epetra_Vector_Ghosted* const y_init,
                                  const Epetra_Vector_Ghosted* const ydot_init, double time_curr,
                                  ProblemResidEval& problem, EpetraJac& jac)
{
    throw errBC("setup_problem()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setPredicted_soln(const Epetra_Vector& y_pred)
{
    throw errBC("setPredicted_soln()");
}
//==================================================================================================================================
int
SolGlobalNonlinear::solve_nonlinear_problem(Solve_Type_Enum solveType, Epetra_Vector_Ghosted* const y_comm,
                                            Epetra_Vector_Ghosted* const ydot_comm, double CJ,
                                            double time_curr, int& num_newt_its_comm,
                                            int& num_linear_solves, int& num_backtracks)
{
    throw errBC("solve_nonlinear_problem()");
    return 1;
}
//==================================================================================================================================
void
SolGlobalNonlinear::calc_ydot(int order, const Epetra_Vector& y_curr, Epetra_Vector& ydot_curr)
{
    throw errBC("calc_ydot()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setTolerances(double reltol, int n, const double* const abstol)
{
    throw errBC("setTolerances()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setTolerances_deltaDamping(double reltol_dd, int n, const double* const abstol_dd)
{
    throw errBC("setTolerances_deltaDamping()");
}
//==================================================================================================================================
void SolGlobalNonlinear:: setMaxNewtIts(const int maxNewtIts)
{
    throw errBC("setMaxnewtIts()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setProblemType(int jacFormMethod)
{
    throw errBC("setProblemType()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setDefaultSolnWeights()
{
    throw errBC("setDefaultSolnWeights()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setPrintFlag(int print_flag)
{
    m_print_flag = print_flag;
}
//==================================================================================================================================
void
SolGlobalNonlinear::setPrintSolnOptions(bool dumpJacobians)
{
    throw errBC("setPrintSolnOptions()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setNonLinOptions(int min_newt_its, bool matrixConditioning, bool colScaling, bool rowScaling,
                                     int colScaleUpdateFrequency)
{
    throw errBC("setNonLinOptions()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setPreviousTimeStep(const double timeStep_comm, const Epetra_Vector& y_nm1,
                                        const Epetra_Vector& ydot_nm1)
{
    throw errBC("setPreviousTimeStep()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::computeResidWts()
{
    throw errBC("computeResidWts()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::getResidWts(Epetra_Vector_Owned& residWts)
{
    throw errBC("getResidWts()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setRowScaling(bool onoff)
{
    throw errBC("setRowScaling()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setColScaling(bool onoff, int colScaleUpdateFrequency)
{
    throw errBC("setColScaling()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setDampingToggles(const bool residSolnDamping, const bool deltaDamping, const bool highLowDamping)
{
    throw errBC("setDampingToggles()");
}
//==================================================================================================================================
void
SolGlobalNonlinear::setSolutionBounds(const Epetra_Vector_Owned& lowBounds, const Epetra_Vector_Owned& highBounds)
{
    throw errBC("setSolutionBounds()");
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
