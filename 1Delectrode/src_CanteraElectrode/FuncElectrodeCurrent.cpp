/**
 *  @file FuncElectrodeCurrent.cpp
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "FuncElectrodeCurrent.h"

using namespace std;
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
int CurrentFunc::fileInit_ = 0;
//undef this to stop writing to ivResults.dat
#define WRITE_IV_ITER
//==================================================================================================================================
CurrentFunc::CurrentFunc(m1d::SolNonlinear_CurrentSolve* solverC, m1d::SolNonlinear* solverV, m1d::ProblemResidEval* func,
                         m1d::Epetra_Vector_Ghosted* y_comm, m1d::Epetra_Vector_Ghosted* ydot_comm, double CJ,
                         double time_curr, double deltaT) :
    m_solver_constantCurr(solverC),
    m_solver_constantVoltage(solverV),
    m_func_constV(func),
    m_y_comm(y_comm),
    m_ydot_comm(ydot_comm),
    m_deltaT(deltaT),
    m_CJ(CJ),
    m_time_curr(time_curr)
{
    printLvl_ = 1;
#ifdef WRITE_IV_ITER
    if (printLvl_) {
        if (!fileInit_) {
            fileInit_ = 1;
            ivResultFile_ = fopen("ivResults.dat", "w");
            fprintf(ivResultFile_, "VARIABLES = \"time [s]\" \"Current\" \"Voltage\"\n");
        } else {
            ivResultFile_ = fopen("ivResults.dat", "a");
        }
    }
#endif //WRITE_IV_ITER
}
//==================================================================================================================================
CurrentFunc::~CurrentFunc()
{
    if (ivResultFile_) {
        fclose(ivResultFile_);
    }
}
//==================================================================================================================================
int CurrentFunc::nEquations() const
{
    return 1;
}
//==================================================================================================================================
int CurrentFunc::evalResidSS(const double t, const double* const x, double* const r)
{
    /*
     *   Set the cathode collector voltage in the main function
     */
    double voltsC = x[0];

    int retn = m_func_constV->setSolutionParam("CathodeCollectorVoltage", voltsC);
    if (retn < 0) {
        exit(-1);
    }

    if (printLvl_) {
        printf("CurrentFunc: call integrate with voltage = %20.13g\n", voltsC);
    }

    /*
     *  Set the print level for the inner iteration.
     *   For small print levels we want to suppress the printing in the inner iteration
     *   For larger print levels we want to maintain the printing level that the user expects from the input deck
     */
    if (printLvl_ < 2) {
        m_solver_constantVoltage->setPrintFlag(0);
    } else if (printLvl_ <= 3) {
        m_solver_constantVoltage->setPrintFlag(2);
    } else  {
        m_solver_constantVoltage->setPrintFlag(printLvl_);
    }

    /*
     *  Solve the system of equations at constant voltage
     */
    Zuzax::Solve_Type ss = Zuzax::Solve_Type::TimeDependentAccurate_Solve;
    int num_newt_its, num_linear_solves, numbacktracks;
    m_solver_constantVoltage->solve_nonlinear_problem(ss, m_y_comm, m_ydot_comm, m_CJ, m_time_curr, num_newt_its,
                                                      num_linear_solves, numbacktracks);

    /*
     * Optional print out the intermediate solution to the logfile
     *  event 3 is an intermediate nonlinear step
     */
    if (printLvl_ > 6) {
        m_func_constV->showProblemSolution(3, true, m_time_curr, m_deltaT, *m_y_comm, m_ydot_comm);
    }

    /*
     *  Extract the total amps from the simulation and put it into the residual
     */
    double amps;
    retn =  m_func_constV->getSolutionParam("CathodeCollectorCurrent", &amps);
    r[0] = amps;


#ifdef DEBUG_ELECTRODE_MODE_NOT
    FILE* fp = fopen("iv.txt", "w");
    for (int n = 0; n < 200; n++) {
        voltsC = x[0] - 0.001 + 0.00001 * n;
        retn = m_func_constV->setSolutionParam("CathodeCollectorVoltage", voltsC);
        m_solver_constantVoltage->solve_nonlinear_problem(ss, m_y_comm, m_ydot_comm, m_CJ, m_time_curr, num_newt_its,
                                                          num_linear_solves,
                                                          numbacktracks);
        retn =  m_func_constV->getSolutionParam("CathodeCollectorCurrent", &amps);
        fprintf(fp," %20.13g  %d %20.13g \n", voltsC, amps);
    }
    fclose(fp);
    exit(-1);
#endif
    if (printLvl_) {
        printf("CurrentFunc(time = %10g): Curr(voltage = %20.13g) = %20.13g\n", m_time_curr, x[0], amps);
#ifdef WRITE_IV_ITER
        fprintf(ivResultFile_, "%20.13g \t%20.13g \t%20.13g\n",m_time_curr, amps, x[0]);
#endif //WRITE_IV_ITER
    }
    return 0;
}
//==================================================================================================================================
void  CurrentFunc::set_deltaT(double deltaT)
{
    m_deltaT = deltaT;
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
