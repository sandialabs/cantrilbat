/**
 *  @file BEulerInt_Battery.cpp
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BEulerInt_Battery.h"

#include "m1d_BatteryResidEval.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_CanteraElectrodeGlobals.h"

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

namespace beuler
{

  //=====================================================================================================================
  /*
   *  Constructor. Default settings: dense jacobian, no user-supplied
   *  Jacobian function, Newton iteration.
   */
  BEulerInt_Battery::BEulerInt_Battery() :
    BEulerInt()
  {
  }
  //=====================================================================================================================
  // Copy constructor
  /*
   *
   * @param r  Object to be copied
   */
  BEulerInt_Battery::BEulerInt_Battery(const BEulerInt_Battery &r)
  {
    *this = r;
  }
  //=====================================================================================================================
  BEulerInt_Battery &
  BEulerInt_Battery::operator=(const BEulerInt_Battery &r)
  {
    if (this == &r) {
      return *this;
    }
    ZZCantera::Integrator::operator=(r);

    return *this;
  }
  //=====================================================================================================================
  // Destructor
  BEulerInt_Battery::~BEulerInt_Battery()
  {

  } 
  //=====================================================================================================================
  int BEulerInt_Battery::calcConsistentInitialDerivs() {
    int retn = 0;

    m1d::BatteryResidEval *batResid = dynamic_cast<m1d::BatteryResidEval *>(m_func);
    if (!batResid) {
      throw m1d::m1d_Error("BEulerInt_Battery::calcConsistentInitialDerivs()",
		          " Function doesn't derive from BatteryResidEval");
    }


    // A more general way to get the boundary conditioxsn
    batResid->reportCathodeVoltageBC(time_n, BC_Type_current_, CurrentValueBC_, BC_TimeDep_, TimeDep_);

    // If the bc type is even (specified voltage), we are good to go
    if ( (BC_Type_current_ % 2 ) == 0) {
      retn = calcConsistentInitialDerivs_CV();
      return retn;
    }

   // Let's take the simple case
    currentNeeded_ =  CurrentValueBC_;

    // Go get the current value of the cathode voltage and use that as the initial value of the voltage boundary condition
    //  double cathodeVoltage =  batResid->reportCathodeVoltage();
    //! Specify the voltage at the cathode
    CathodeVoltageBest_ = m1d::PSCinput_ptr->CathodeVoltageSpecified_;

    // Change the problem to a dirichlet condition
    batResid->changeCathodeVoltageBC(0, CathodeVoltageBest_);
    retn = calcConsistentInitialDerivs_CV();

    if (m_print_flag > 0){
      std::string snn = "CV Solution";
      m_func->showSolutionVector(snn, time_n, delta_t_n, *m_y_n);
    }
    if (m_print_flag > 0){
      std::string snn = "CV Solution Time Derivative";
      m_func->showSolutionVector(snn, time_n, delta_t_n, *m_ydot_n);
    }
    
    
    // Change the problem to a flux boundary condition
    batResid->changeCathodeVoltageBC(BC_Type_current_, CurrentValueBC_, BC_TimeDep_, TimeDep_);
    retn = calcConsistentInitialDerivs_CV();

    if (m_print_flag > 0){
      std::string snn = "CC Solution";
      m_func->showSolutionVector(snn, time_n, delta_t_n, *m_y_n);
    }
    if (m_print_flag > 0){
      std::string snn = "CC Solution Time Derivative";
      m_func->showSolutionVector(snn, time_n, delta_t_n, *m_ydot_n);
    }
    
    return retn;

  }

  //====================================================================================================================
  // Solve for the consistent initial conditions and consistent initial time derivatives.
  /*
   *  A special nonlinear problem is solved to find the consistent initial time derivatives and DAE 
   *  values that cause the residual system to be solved at t = time_n. 
   *
   *  This algorithm is described in the notes and utilizes the m_isAlgebraic[] field variables. 
   *
   *  Given values for the non-algabraic unknowns this routine seeks to calculate a consistent value of 
   *  the algebraic unknowns and ydot for the algabraic unknowns to solve the residual equations. 
   *  On a successful solution, the answer is storred in the vector variables m_y_n and m_ydot_n.
   *  Note this answer depends on delta_t_n weakly. The solution depends on a integrated source and the
   *  integrated source depends on delta_t_n.
   *
   *  Note that it can be called at any time in the calculation. However, it should always be called at the
   *  start of the calculation.
   *
   *  @return 0  Always returns 0
   */
  int BEulerInt_Battery::calcConsistentInitialDerivs_CV() {

    /*
     *  Set the max newton iterations high, as this is a do or die calculation
     */
    m_nonlin->setMaxNewtIts(150);

    m_func->setAtolVector_DAEInit(1.0E-3, *m_y_n, *m_ydot_n);

    m_func->setAtolDeltaDamping_DAEInit(1.0, *m_y_n, *m_ydot_n);

    /*
     *  Change the absolute error tolerances to that of the DAE init tolerances. Note, the time derivatives will
     *  probably have a different scaling based on the expected time response of the system.
     *
     *  -> go get atolDAE vector
     */
    const m1d::Epetra_Vector_Owned & atolDAEInitRef = m_func->atolVector_DAEInit();
    /*
     *  Tell nonlinear solver to use atolDAE
     */
    m_nonlin->setTolerances(m_reltol, m_NumLcEqns, & (atolDAEInitRef[0]));

    // not needed ???
    m_nonlin->setPredicted_soln(*m_y_n);

    // not needed ??? -> needed to specify delta_t_n
    m_nonlin->setPreviousTimeStep(delta_t_n, *m_y_nm1, *m_ydot_nm1);

    const m1d::Epetra_Vector_Owned &abs_dd = m_func->atolDeltaDamping();
    m_nonlin->setTolerances_deltaDamping(m_reltol, m_NumLcOwnedEqns, &(abs_dd[0]));
    /*
     * Solve the system of equations at the current time step.
     * Note - x_corr_n and x_dot_n are considered to be updated,
     * on return from this solution.
     */
    int num_linear_solves;
    int numbacktracks;
    int num_newt_its;
    m1d::Solve_Type_Enum ss = m1d::DAESystemInitial_Solve;
    double CJ = 1.0/delta_t_n;

#ifdef DEBUG_INIT_CALCULATION
    int pold = m_nonlin->m_print_flag;
    //m_nonlin->m_print_flag = 10;
    //m1d::SolNonlinear::s_print_NumJac = 10;
#endif


    int ierror = m_nonlin->solve_nonlinear_problem(ss, m_y_n, m_ydot_n, CJ, time_n, num_newt_its, num_linear_solves,
                                                   numbacktracks);

#ifdef DEBUG_INIT_CALCULATION
    m_nonlin->m_print_flag = pold;
    m1d::SolNonlinear::s_print_NumJac = 0;
#endif
    /*
     *  If we have experienced an error in the DAE initial solve calculation, we currently have no recourse but
     *  to end the calculation in failure. We may change this behavior in the future, I don't know.
     */
    if (ierror < 0) {
      throw CanteraError("BEulerInt::calcConsistentInitialDerivs", "Nonlinear solver failed to converge");
    }

    /*
     *   Always call writeSolution to write out the initial conditions
     */
    if (m_print_flag > 0){
      m_func->writeSolution(0, true, time_n, delta_t_n, 0.0, *m_y_n_owned, m_ydot_n_owned);
    }
    if (m_print_flag > 0){
      std::string snn = "Solution Time Derivative";
      m_func->showSolutionVector(snn, time_n, delta_t_n, *m_ydot_n);
    }
    
    /*
     * Call a different user routine at the end of each step,
     * that will probably print to a file.
     */
    m_func->user_out(0, 0.0, 0.0, 0, *m_y_n, m_ydot_n);

    /*
     *  Here we show that we have in fact solved the residual equation for delta_t well. There may be some variation
     *  in the residual due to a change in the equation set. However, barring that, the residual should be near zero.
     */
    if (m_print_flag > 5) {
      double rdelta_t = 1.0 /  delta_t_n;
      m_func->residEval(m_resid, true, m_y_n, m_ydot_n, 0.0, rdelta_t, m1d::Base_ResidEval, m1d::DAESystemInitial_Solve);
      double rnorm = res_error_norm(*m_resid, "DAE_Init_Residual", 10);
      printf("rnorm (DAE) = %g\n", rnorm);
    }

    return 0;
  }
  //====================================================================================================================
  //  Check to see that the predicted solution satisfies proper requirements.
  /*
   *       @return Returns a negative number if the step is inappropriate. Then the stepsize is reduced
   *               and the method is checked again.
   */
  int BEulerInt_Battery::check_predicted_soln(m1d::Epetra_Vector_Ghosted & y_n, m1d::Epetra_Vector_Ghosted & ydot_n,
			                      double CJ, double time_n)
  {
     // took this out as it doesn't do anything
     //m_func->residEval(m_resid, true, &y_n, &ydot_n, time_n,  CJ, m1d::Base_ResidEval, m1d::TimeDependentAccurate_Solve);

     return 0;
  }
  //====================================================================================================================
} // End of beuler namespace
//====================================================================================================================
