/**
 *  @file m1d_ProblemResidEval.cpp
 *
 **/
/*
 * $Author: hkmoffa $
 * $Revision: 564 $
 * $Date: 2013-03-08 16:35:51 -0700 (Fri, 08 Mar 2013) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_BatteryResidEval.h"

#include "m1d_SurDomain_AnodeCollector.h"
#include "m1d_SurDomain_CathodeCollector.h"
#include "m1d_porousLiKCl_dom1D.h"
#include "m1d_porousElectrode_dom1D.h"

#include "m1d_DomainLayout.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_EpetraExtras.h"
#include "m1d_Comm.h"
#include "m1d_GlobalIndices.h"
#include "m1d_globals.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

extern m1d::ProblemStatementCell PSinput;

namespace m1d
{

  //=====================================================================================================================

  // Default constructor
  /*
   *
   * @param atol   Absolute tolerance for calculation
   */
  BatteryResidEval::BatteryResidEval(double atol) :
    ProblemResidEval(atol),
    doHeatSourceTracking_(0),
    maxSubGridTimeSteps_(0),
    QdotPerArea_n_(0.0),
    QdotPerArea_nm1_(0.0),
    QdotAnodePerArea_n_(0.0),
    QdotSeparatorPerArea_n_(0.0),
    QdotCathodePerArea_n_(0.0),
    electrodeCrossSectionalArea_(0.0),
    capacityAnodePA_(0.0),
    capacityCathodePA_(0.0),
    capacityLeftAnodePA_(0.0),
    capacityLeftCathodePA_(0.0),
    capacityDischargedAnodePA_(0.0),
    capacityDischargedCathodePA_(0.0),
    capacityPAEff_(0.0),
    capacityLeftPAEff_(0.0),
    timeLeft_(0.0),
    Crate_current_(0.0)

  {
     doHeatSourceTracking_ = PSinput.doHeatSourceTracking_;
     electrodeCrossSectionalArea_ = PSinput.cathode_input_->electrodeGrossArea;
  }
  //=====================================================================================================================
  // Destructor
  /*
   *
   */
  BatteryResidEval::~BatteryResidEval()
  {
  }
  //=====================================================================================================================
  //! Default copy constructor
  /*!
   *
   * @param r  Object to be copied
   */
  BatteryResidEval::BatteryResidEval(const BatteryResidEval &r) :
    ProblemResidEval(r.m_atol)
  {
    *this = r;
  }
  //=====================================================================================================================
  // Assignment operator
  /*
   *
   * @param r  Object to be copied
   * @return   Returns a copy of the current problem
   */
  BatteryResidEval &
  BatteryResidEval::operator=(const BatteryResidEval &r)
  {
    if (this == &r) {
      return *this;
    }

    ProblemResidEval::operator=(r);

    doHeatSourceTracking_              = r.doHeatSourceTracking_;
    maxSubGridTimeSteps_               = r.maxSubGridTimeSteps_;
    QdotPerArea_n_                     = r.QdotPerArea_n_;
    QdotPerArea_nm1_                   = r.QdotPerArea_nm1_;
    QdotAnodePerArea_n_                = r.QdotAnodePerArea_n_;
    QdotSeparatorPerArea_n_            = r.QdotSeparatorPerArea_n_;
    QdotCathodePerArea_n_              = r.QdotCathodePerArea_n_;

    electrodeCrossSectionalArea_         = r.electrodeCrossSectionalArea_;
    capacityAnodePA_                     = r.capacityAnodePA_;
    capacityCathodePA_                   = r.capacityCathodePA_;
    capacityLeftAnodePA_                 = r.capacityLeftAnodePA_;
    capacityLeftCathodePA_               = r.capacityLeftCathodePA_;
    capacityDischargedAnodePA_           = r.capacityDischargedAnodePA_;
    capacityDischargedCathodePA_         = r.capacityDischargedCathodePA_;
    capacityPAEff_                       = r.capacityPAEff_;
    capacityLeftPAEff_                   = r.capacityLeftPAEff_;
    timeLeft_                            = r.timeLeft_;
    Crate_current_                       = r.Crate_current_;

    return *this;
  }
  //=====================================================================================================================

  //! Calculate the initial conditions
  /*!
   *   This calls the parent class initialConditions method to loop over the volume and surface domains.
   *   Then the method tried to better estimate the electrolyte potential.
   *
   * @param doTimeDependentResid    Boolean indicating whether we should
   *                                formulate the time dependent residual
   * @param soln                    Solution vector. This is the input to
   *                                the residual calculation.
   * @param solnDot                 Solution vector. This is the input to
   *                                the residual calculation.
   * @param t                       Time
   * @param delta_t                 delta_t for the initial time step
   */
  void
  BatteryResidEval::initialConditions(const bool doTimeDependentResid,
				      Epetra_Vector_Ghosted *soln,
				      Epetra_Vector_Ghosted *solnDot,
				      double &t,
				      double &delta_t,
                                      double &delta_t_np1)
  {
      //
      //   We obtain the cross-sectional area from the cathode. 
      //   However, we require the cross-sectional area to be consistent across inputs
      //
      electrodeCrossSectionalArea_ = PSinput.cathode_input_->electrodeGrossArea;
      //
      //   loop over all volume and surface domains providing initial guess
      //
      ProblemResidEval::initialConditions(doTimeDependentResid, soln, solnDot, t, delta_t, delta_t_np1);
      //
      //  improveInitialConditions(soln);
      //
  }
  //====================================================================================================================
  //! Improve upon initial conditions by computing first order potential losses, equilibrium reaction voltages, etc.
  /*!
   * Not currently finished
   */
  void
  BatteryResidEval::improveInitialConditions(Epetra_Vector_Ghosted *soln_ptr)
  {
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    if ( mypid ) return; //using global variables, we only need to do this once. Need to sync.

    Epetra_Vector &soln = *soln_ptr;

    //Domain Layout ptr DL_ptr_
    //Domain layout has the information on the global nodes and the x-positions
    DomainLayout &DL = *DL_ptr_;

    // This algorithm really only works if you have three or more bulk domains
    if ( DL.NumBulkDomains < 3 ) return;

    //   Loop over the Volume Domains
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {

      //positions at boundaries
      double xleft  = DL_ptr_->StartXLoc_Domain[iDom] ;
      double xright = DL_ptr_->EndXLoc_Domain[iDom] ;
      //nodes at boundaries
      int leftGbNode  = DL_ptr_->StartGBNode_Domain[iDom];
      int rightGbNode = DL_ptr_->EndGBNode_Domain[iDom];
      fprintf(stderr, "Boundaries of domain %d are %g and %g\n", iDom, xleft, xright);

      BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
      int iVar_Voltage_BD = d_ptr->BDD_.VariableIndexStart_VarName[Voltage];

      for ( int iGbNode = leftGbNode; iGbNode < rightGbNode; iGbNode++ ) {
	double xNode = (*GI_ptr_->XNodePos_GbNode)[iGbNode];	
	int voltEqn = GI_ptr_->IndexStartGbEqns_GbNode[iGbNode] + iVar_Voltage_BD ;
	fprintf(stderr, "Old var value was %g at x=%g\n", soln[voltEqn], xNode );

      }
    }

    if ( PSinput.cathodeBCType_ == 0 
	 || PSinput.cathodeBCType_ == 2 
	 || PSinput.cathodeBCType_ == 4 
	 || PSinput.cathodeBCType_ == 6 
	 || PSinput.cathodeBCType_ == 8 ) {
      //! CASE FOR SPECIFIED VOLTAGE 
      //! Prior to this, the voltage in the electrolyte should have been specified
      //! as the specified BC voltage minus the open circuit voltage in the edge domains.
      //! Interpolate voltage between guessed voltages in edge domains.
      //
      BulkDomain1D *d_ptr;
      int iVar_Voltage_BD;
      int voltEqn;
      
      //! Left-most bulk domain
      int iDom = 0;
      //The left hand node for interpolation is the N-1st node of the left domain (near it's right edge)
      int iNodeLeft = DL_ptr_->EndGBNode_Domain[iDom] - 1;
      double xleft = (*GI_ptr_->XNodePos_GbNode)[iNodeLeft];	
      //find the voltage at this node
      d_ptr = DL.BulkDomain1D_List[iDom];
      iVar_Voltage_BD = d_ptr->BDD_.VariableIndexStart_VarName[Voltage];
      voltEqn = GI_ptr_->IndexStartGbEqns_GbNode[iNodeLeft] + iVar_Voltage_BD ;
      double voltLeft = soln[voltEqn];
      
      //! Right-most bulk domain
      iDom = DL.NumBulkDomains - 1 ;
      //The right hand node for interpolation is the 2nd node of the right domain (near it's left edge)
      int iNodeRight = DL_ptr_->StartGBNode_Domain[iDom] + 1;
      double xright = (*GI_ptr_->XNodePos_GbNode)[iNodeRight];	
      //find the voltage at this node
      d_ptr = DL.BulkDomain1D_List[iDom];
      iVar_Voltage_BD = d_ptr->BDD_.VariableIndexStart_VarName[Voltage];
      voltEqn = GI_ptr_->IndexStartGbEqns_GbNode[iNodeRight] + iVar_Voltage_BD ;
      double voltRight = soln[voltEqn];
      
      //! Loop over intermediate domains
      //   Loop over the Volume Domains
      for (iDom = 1; iDom < DL.NumBulkDomains - 1; iDom++) {
	
	//nodes at boundaries
	int leftGbNode  = DL_ptr_->StartGBNode_Domain[iDom];
	int rightGbNode = DL_ptr_->EndGBNode_Domain[iDom];
	
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	int iVar_Voltage_BD = d_ptr->BDD_.VariableIndexStart_VarName[Voltage];
	
	for ( int iGbNode = leftGbNode; iGbNode < rightGbNode; iGbNode++ ) {
	  double xNode = (*GI_ptr_->XNodePos_GbNode)[iGbNode];	
	  voltEqn = GI_ptr_->IndexStartGbEqns_GbNode[iGbNode] + iVar_Voltage_BD ;
	  fprintf(stderr, "Old var value was %g", soln[voltEqn] );
	  double voltInterp = voltLeft + ( voltRight - voltLeft ) * ( xNode - xleft ) / ( xright - xleft );
	  soln[voltEqn] = voltInterp;
	  fprintf(stderr, " and new value is %g at x=%g\n", soln[voltEqn], xNode );
	  
	}
      }
      //! First point for Right edge domain
      iDom = DL.NumBulkDomains - 1 ;
      //nodes at boundaries
      int iGbNode  = DL_ptr_->StartGBNode_Domain[iDom];
      d_ptr = DL.BulkDomain1D_List[iDom];
      iVar_Voltage_BD = d_ptr->BDD_.VariableIndexStart_VarName[Voltage];
      
      double xNode = (*GI_ptr_->XNodePos_GbNode)[iGbNode];	
      voltEqn = GI_ptr_->IndexStartGbEqns_GbNode[iGbNode] + iVar_Voltage_BD ;
      fprintf(stderr, "Old var value was %g", soln[voltEqn] );
      double voltInterp = voltLeft + ( voltRight - voltLeft ) * ( xNode - xleft ) / ( xright - xleft );
      soln[voltEqn] = voltInterp;
      fprintf(stderr, " and new value is %g at x=%g\n", soln[voltEqn], xNode );
      
      //exit(1);
    }  else {
      double cathodeCurrent;
      getSolutionParam( "CathodeCollectorCurrent", &cathodeCurrent );

      //! CASE FOR SPECIFIED CURRENT 
      //! Prior to this, the voltage in the left bulk domain electrolyte should have been specified
      //! as the specified BC voltage minus the open circuit voltage in the edge domains.
      //! Take an estimate for the conductivity and the specified current 
      //! to estimated the voltage change across intermediate domains.  
      //! For last domain, use OCV or voltage difference to estimate voltage in electrode (metallic) phase.
      //
    }
  }

//=====================================================================================================================
// Calculate a residual vector
/*
 *   The basic algorithm is to loop over the volume domains.
 *   Then, we loop over the surface domains
 *
 * @param res                     residual output
 * @param doTimeDependentResid    Boolean indicating whether the time dependent residual is requested
 * @param soln                    Pointer to the solution vector. This is the input to the residual calculation.
 * @param solnDot                 Pointer to the solution Dot vector. This is the input to the residual calculation.
 * @param t                       current time
 * @param rdelta_t                delta t inverse
 * @param residType               Residual type
 * @param solveType               Solve type
 */
void
BatteryResidEval::residEval(Epetra_Vector_Owned* const & res,
                            const bool doTimeDependentResid,
                            const Epetra_Vector *soln_ptr,
                            const Epetra_Vector *solnDot_ptr,
                            const double t,
                            const double rdelta_t,
                            const ResidEval_Type_Enum residType,
                            const Solve_Type_Enum solveType)
{
    if (!resInternal_ptr_) {
	resInternal_ptr_ = new Epetra_Vector(*res);
    }

    if (residType == Base_ResidEval) {
	counterResBaseCalcs_++;
    } else if (residType == JacBase_ResidEval) {
	counterJacBaseCalcs_++;
    } else if (residType == JacDelta_ResidEval) {
	counterJacDeltaCalcs_++;
    } else if (residType == Base_ShowSolution) {
	counterResShowSolutionCalcs_++;
    }
    // Get a local copy of the domain layout
    DomainLayout &DL = *DL_ptr_;
    /*
     *   Zero the residual vector
     */
    res->PutScalar(0.0);
  
    /*
     * We calculate solnOld_ptr_ here
     */
    if (doTimeDependentResid) {
	calcSolnOld(*soln_ptr, *solnDot_ptr, rdelta_t);
    }
    /*
     *   Propagate the solution of the system down to the underlying objects where necessary.
     *   It is necessary to do this for mesh unknowns.
     *   Also do a loop over the nodes to carry out any precalculations that are necessary
     */
    setStateFromSolution(doTimeDependentResid, soln_ptr, solnDot_ptr, t, rdelta_t);
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	d_ptr->incrementCounters(residType);
	d_ptr->residEval_PreCalc(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t, residType, solveType);
    }
    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
	SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	d_ptr->incrementCounters(residType);
	d_ptr->residEval_PreCalc(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t, residType, solveType);
    }
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	d_ptr->residEval(*res, doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t, residType, solveType);
    }
    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
	SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	d_ptr->residEval(*res, doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t, residType, solveType);
    }

    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	d_ptr->residEval_PostCalc(*res, doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t, residType, solveType);
    }
    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
	SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	d_ptr->residEval_PostCalc(*res, doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t, residType, solveType);
    }

    maxSubGridTimeSteps_ = 0;
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	porousElectrode_dom1D* e_ptr =  dynamic_cast<porousElectrode_dom1D*>(d_ptr);
	if (e_ptr) {
	  int maxE = e_ptr->getMaxSubGridTimeSteps();
	  maxSubGridTimeSteps_ = MAX(maxSubGridTimeSteps_, maxE);
        }
    }

}
//=====================================================================================================================
// This function may be used to create output at various points in the
// execution of an application.
/*
 *   These functions are not affected by the print controls of the nonlinear solver
 *   and the time stepper.
 *
 *      ievent is a description of the event that caused this
 *      function to be called.
 *
 *      @param ievent  Event that's causing this routine to be called.
 *                     =  0 Initial conditions for a calculation
 *                     =  1 Completion of a successful intermediate time step.
 *                     =  2 Completion of a successful Final time or final calculation.
 *                     =  3 Completion of a successful Intermediate nonlinear step
 *                     = -1 unsuccessful time step that converged, but failed otherwise
 *                     = -2 unsuccessful nonlinear step.
 *
 *      @param time_current      Current time
 *      @param delta_t_n         Current value of delta_t
 *      @param istep             Current step number
 *      @param y_n               Current value of the solution vector
 *      @param ydot_n_ptr        Current value of the time deriv of the solution vector
 */
void
BatteryResidEval::user_out(const int ievent,
			     const double time_current,
			   const double delta_t_n,
			   const int istep,
			   const Epetra_Vector_Ghosted &y_n,
			   const Epetra_Vector_Ghosted * const ydot_n_ptr)
{
    ProblemResidEval::user_out(ievent, time_current, delta_t_n, istep, y_n, ydot_n_ptr);
    
    
}
//=====================================================================================================================
static void
sprint_line(char * buf, const char * const st, const int num)
{
  int n = strlen(buf);
  buf += n;
  for (int k = 0; k < num; k++, buf++) {
    sprintf(buf, "%s", st);
  }
  sprintf(buf, "\n");
}

//=====================================================================================================================
// Write the solution to either the screen or to a log file
/*
 *
 * @param ievent  Type of the event. The following form is used:
 *             0 Initial conditions
 *             1 Completion of a successful intermediate step.
 *             2 Final successful conditions.
 *             3 Intermediate nonlinear step
 *            -1 unsuccessful step
 *
 * @param m_y_n    Current value of the solution vector
 * @param m_ydot_n  Current value of the derivative of the solution vector
 */
void
BatteryResidEval::showProblemSolution(const int ievent,
                                      bool doTimeDependentResid,
                                      const double t,
                                      const double delta_t,
                                      const Epetra_Vector_Ghosted &y_n,
                                      const Epetra_Vector_Ghosted * const ydot_n,
                                      const Solve_Type_Enum solveType,
                                      const double delta_t_np1)
{

  //bool doAllProcs = false;
  time_t aclock;
  ::time(&aclock); /* Get time in seconds */
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  char buf[256];
  bool duplicateOnAllProcs = false;
  // Get a local copy of the domain layout
  DomainLayout &DL = *DL_ptr_;

  Epetra_Vector *solnAll = GI_ptr_->SolnAll;
  m1d::gather_nodeV_OnAll(*solnAll, y_n);

  Epetra_Vector *soln_dot_All = GI_ptr_->SolnDotAll;
  if (ydot_n) {
    m1d::gather_nodeV_OnAll(*soln_dot_All, *ydot_n);
  }
  double rdelta_t = 0.0;
  if (delta_t > 0.0) {
    rdelta_t = 1.0 / delta_t;
  }

  residEval(resInternal_ptr_, doTimeDependentResid, &y_n, ydot_n, t, rdelta_t, Base_ShowSolution, solveType);
 
  gatherCapacityStatistics();


  int indentSpaces = 4;
  string indent = "    ";
  const char *ind = indent.c_str();
  if (!mypid || duplicateOnAllProcs) {
    sprintf(buf, "%s", ind);
    sprint_line(buf, "-", 100);
    Cantera::writelog(buf);
    sprintf(buf, "%s ShowProblemSolution : Time       %-12.3E\n", ind, t);
    Cantera::writelog(buf);
    if (solveType == TimeDependentInitial) {
       sprintf(buf, "%s                       Delta_t    %-12.3E  (initial solution with no previous solution)\n", ind, delta_t);
    } else {
       sprintf(buf, "%s                       Delta_t    %-12.3E\n", ind, delta_t);
    }
    Cantera::writelog(buf);
    sprintf(buf, "%s                       StepNumber %6d\n", ind, m_StepNumber);
    Cantera::writelog(buf);
    sprintf(buf, "%s                       Delta_t_p1 %-12.3E\n", ind, delta_t_np1);
    Cantera::writelog(buf);
    sprintf(buf, "%s                       Capacity_Anode        %-12.3E  Amp hr / m2\n", ind, capacityAnodePA_ / 3600.);
    Cantera::writelog(buf);
    sprintf(buf, "%s                       Capacity_Anode_Left   %-12.3E  Amp hr / m2\n", ind, capacityLeftAnodePA_ / 3600.);
    Cantera::writelog(buf);
    sprintf(buf, "%s                       Capacity_Cathode      %-12.3E  Amp hr / m2\n", ind, capacityCathodePA_ / 3600.);
    Cantera::writelog(buf);
    sprintf(buf, "%s                       Capacity_Cathode_Left %-12.3E  Amp hr / m2\n", ind, capacityLeftCathodePA_ / 3600.);
    Cantera::writelog(buf);
    sprintf(buf, "%s                       Crate_current         %-12.3E  \n", ind, Crate_current_);
    Cantera::writelog(buf);

    sprintf(buf, "%s", ind);
    sprint_line(buf, "-", 100);
    Cantera::writelog(buf);
  }

  Domain1D *d_ptr = DL.SurDomain1D_List[0];
  do {
    d_ptr->showSolution(solnAll, soln_dot_All, &y_n, ydot_n, solnOld_ptr_,
                        resInternal_ptr_, t, rdelta_t, indentSpaces + 2,
                        duplicateOnAllProcs);
    BulkDomain1D *bd_ptr = dynamic_cast<BulkDomain1D *> (d_ptr);
    if (bd_ptr) {
      //BulkDomainDescription &BDD_;
      SurfDomainDescription *sdd = bd_ptr->BDD_.RightSurf;
      if (sdd) {
        int idS = sdd->ID();
        d_ptr = DL.SurDomain1D_List[idS];
      } else {
        d_ptr = 0;
      }
    } else {
      SurDomain1D *sd_ptr = dynamic_cast<SurDomain1D *> (d_ptr);
      DomainDescription *dd = sd_ptr->SDD_.RightDomain;
      if (dd) {
        BulkDomainDescription *bdd = dynamic_cast<BulkDomainDescription *> (dd);
        int idS = bdd->ID();
        d_ptr = DL.BulkDomain1D_List[idS];
      } else {
        d_ptr = 0;
      }
    }
  } while (d_ptr);

  if (!mypid || duplicateOnAllProcs) {
    sprintf(buf, "%s", ind);
    sprint_line(buf, "-", 100);
    Cantera::writelog(buf);
  }

  Epetra_Comm *c = LI_ptr_->Comm_ptr_;
  c->Barrier();

}

  //=====================================================================================================================
  // Write the solution to either the screen or to a log file
  /*
   *
   * @param ievent  Type of the event. The following form is used:
   *             0 Initial conditions
   *             1 Completion of a successful intermediate step.
   *             2 Final successful conditions.
   *             3 Intermediate nonlinear step
   *            -1 unsuccessful step
   *
   * @param m_y_n    Current value of the solution vector
   * @param m_ydot_n  Current value of the derivative of the solution vector
   */
  void
  BatteryResidEval::writeSolution(const int ievent,
				  const bool doTimeDependentResid,
				  const double time_current,
				  const double delta_t_n,
				  int istep,
				  const Epetra_Vector_Ghosted &y_n,
				  const Epetra_Vector_Ghosted * const ydot_n_ptr,
				  const Solve_Type_Enum solveType, 
                                  const double delta_t_np1)
  {
    ProblemResidEval::writeSolution(ievent, doTimeDependentResid, time_current, delta_t_n, istep, y_n,
				    ydot_n_ptr, solveType, delta_t_np1);
    if (ievent == 0 || ievent == 1 || ievent == 2) {
      write_IV(ievent, doTimeDependentResid, time_current, delta_t_n, istep, y_n, ydot_n_ptr);
    }
  }
  //=====================================================================================================================
  void
  BatteryResidEval::write_IV(const int ievent,
			     const bool doTimeDependentResid,
			     const double time_current,
			     const double delta_t_n,
			     int istep,
			     const Epetra_Vector_Ghosted &y_n,
			     const Epetra_Vector_Ghosted * const ydot_n_ptr)
  {
    DomainLayout &DL = *DL_ptr_;

    // we want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    SurDomain_CathodeCollector *c_ptr = dynamic_cast<SurDomain_CathodeCollector *>(d_ptr);

    double c_icurr = c_ptr->icurrCollector_;
    double phiCath = c_ptr->phiCathode_;

    SurDomain1D *ad_ptr = DL.SurDomain1D_List[0];
    SurDomain_AnodeCollector *ac_ptr = dynamic_cast<SurDomain_AnodeCollector *>(ad_ptr);

    double a_icurr = ac_ptr->icurrCollector_;
    //double phiAnode = ac_ptr->phiAnode_;

    int procID = Comm_ptr->MyPID();

    if (!procID) {
    //looking for cathode capacity and depth of discharge
    BulkDomain1D *cd_ptr = DL.BulkDomain1D_List.back();
    //porousLiKCl_FeS2Cathode_dom1D *cc_ptr = dynamic_cast<porousLiKCl_FeS2Cathode_dom1D *>(cd_ptr);
    double capacityZeroDoD, spec_capacityZeroDoD;
    double dischargedCapacity, spec_dischargedCapacity;
    cd_ptr->getSolutionParam( "CapacityZeroDoD", &capacityZeroDoD );
    cd_ptr->getSolutionParam( "DepthOfDischarge", &dischargedCapacity );
    cd_ptr->getSolutionParam( "SpecificCapacityZeroDoD", &spec_capacityZeroDoD );
    cd_ptr->getSolutionParam( "SpecificDepthOfDischarge", &spec_dischargedCapacity );

    FILE *fp = 0;
    if (ievent == 0) {
      fp = fopen("timeDep_IV.dat", "w");
      fprintf(fp, "TITLE = \"Time versus Current or Voltage\"\n");
      fprintf(fp, "VARIABLES = \" T [s]\"\n");
      fprintf(fp, "\"Voltage [volts] \"\n");
      fprintf(fp, "\"CathodeCurrent [mA/cm2]\"\n");
      fprintf(fp, "\"Initial Specific Cathode Capacity [mA-hr/g] \"\n");
      fprintf(fp, "\"Discharged Specific Cathode Capacity [mA-hr/g] \"\n");
      fprintf(fp, "\"Initial Cathode Capacity [A-s/m2] \"\n");
      fprintf(fp, "\"Discharged Cathode Capacity [A-s/m2] \"\n");
      fprintf(fp, "\"AnodeCurrent [mA/cm2]\"\n");
    } else {
      fp = fopen("timeDep_IV.dat", "a");
    }
    fprintf(fp, "   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E \n", 
	    time_current, phiCath, 0.1 * c_icurr, spec_capacityZeroDoD/3.6, spec_dischargedCapacity/3.6, 
	    capacityZeroDoD, dischargedCapacity, 0.1 * a_icurr );
    fflush(fp);
    fclose(fp);
    }
    Comm_ptr->Barrier();
  }
//====================================================================================================================
// Evaluate a supplemental set of equations that are not part of the solution vector, but are considered
// to be time dependent
/*
 *   Equations in this system are evaluated using the time discretization scheme of the nonlinear solver.
 *   It can be used to track supplemental quantities in the calculation, especially if they need to be
 *   integrated in time.
 *
 *   An example of this may be total flux quantites that are dumped into a phase.
 *
 *   This routine is called at the beginning of the time stepping, in order to set up any initializations,
 *   and it is called after every successful time step, once.
 *
 * @param ifunc   0 Initial call to the function, done whenever the time stepper is entered
 *                1 Called after every successful time step.
 *                2 Called 
 * @param t       Current time
 * @param deltaT  Current value of deltaT
 * @param y       Current value of the solution vectors
 * @param ydot    Current value of time derivative of the solution vectors.
 */
void
BatteryResidEval::evalTimeTrackingEqns(const int ifunc,
                                       const double t,
                                       const double deltaT,
                                       const Epetra_Vector_Ghosted & y,
                                       const Epetra_Vector_Ghosted * const solnDot_ptr)
{
    static FILE *fstep = 0;
    FILE *facc = 0;
    static double tinit_calculation = 0.0;
    const Epetra_Vector_Ghosted *soln_ptr = &y;
    if (doHeatSourceTracking_) {

        double rdelta_t = 1.0E8;
        if (deltaT > 0.0) {
	    rdelta_t = 1.0 / deltaT;
        }

        bool doTimeDependentResid = true;
        /*
         * We calculate solnOld_ptr_ here
         */
        if (doTimeDependentResid) {
            calcSolnOld(*soln_ptr, *solnDot_ptr, rdelta_t);
        }

        if (ifunc == 0) {
            fstep = fopen("stepwiseHeatSource.txt", "w");
            fprintf(fstep, "      tinit  ,     tfinal  ,         anode        ,        separator     ,  "
		    "      cathode       ,         total\n");
            fflush(fstep);
	    tinit_calculation = t;
 
            facc = fopen("accumulatedHeatSource.txt", "w");
            fprintf(facc, "      tinit  ,     tfinal  ,         anode        ,        separator     ,  "
		    "      cathode       ,         total\n");
            fflush(facc);
	    fclose(facc);
	    tinit_calculation = t;
        
        }

 
	DomainLayout &DL = *DL_ptr_;
	//
	//   Loop over the Volume Domains
	//
	for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	    BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	    d_ptr->eval_PostSoln(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t);
	}
	//
	//    Loop over the Surface Domains
	//
	for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
	    SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	    d_ptr->eval_PostSoln(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t);
	}

        if (ifunc == 1)  {
            double qtot = 0.0;
	    double tinit = t - deltaT;
            fprintf(fstep, "% 13.5E , % 13.5E ", tinit,  t);
            for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	        BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
                porousFlow_dom1D* p_ptr = dynamic_cast<porousFlow_dom1D*>(d_ptr);
                double qstep = p_ptr->heatSourceLastStep();
                qtot += qstep;
                fprintf(fstep, " , % 20.7E", qstep);
            }
            fprintf(fstep, " ,  % 20.7E\n", qtot); 
            fflush(fstep);
	}

        
        //	
        // 
        if (ifunc == 2) {
            double q = heatSourceAccumulated();          
            // MP issues
            facc = fopen("accumulatedHeatSource.txt", "a");
	    fprintf(facc, "% 13.5E , % 13.5E ", tinit_calculation,  t);
            for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	        BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
                porousFlow_dom1D* p_ptr = dynamic_cast<porousFlow_dom1D*>(d_ptr);
                double qstep = p_ptr->heatSourceAccumulated();
                fprintf(facc, " , % 20.7E", qstep);
            }
            fprintf(facc, " , % 20.7E\n", q);
            fclose(facc);
        }

    }

}
//================================================================================================================
double
BatteryResidEval::heatSourceLastStep() const
{
    double q = 0.0;
    DomainLayout &DL = *DL_ptr_;
    //
    //   Loop over the Volume Domains
    //
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	porousFlow_dom1D* p_ptr = dynamic_cast<porousFlow_dom1D*>(d_ptr);
	if (p_ptr) {
	    q += p_ptr->heatSourceLastStep();
	}
    }
    //
    //    Loop over the Surface Domains (no terms yet)
    //
    // for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
    //	SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
    //	d_ptr->eval_PostSoln(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t);
    //  }
    return q;
}
//================================================================================================================
double
BatteryResidEval::heatSourceAccumulated() const 
{
    double q = 0.0;
    DomainLayout &DL = *DL_ptr_;
    //
    //   Loop over the Volume Domains
    //
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	porousFlow_dom1D* p_ptr = dynamic_cast<porousFlow_dom1D*>(d_ptr);
	if (p_ptr) {
	    q += p_ptr->heatSourceAccumulated();
	}
    }
    return q;
}
//================================================================================================================
void
BatteryResidEval::heatSourceZeroAccumulated() const
{
    DomainLayout &DL = *DL_ptr_;
    //
    //   Loop over the Volume Domains
    //
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	porousFlow_dom1D* p_ptr = dynamic_cast<porousFlow_dom1D*>(d_ptr);
	if (p_ptr) {
	   p_ptr->heatSourceZeroAccumulated();
	}
    }
}
//====================================================================================================================
int
BatteryResidEval::setSolutionParam(std::string paramName, double paramVal) {
    if (paramName != "CathodeCollectorVoltage") {
	return -1;
    } 
    // Go find the Cathode Collector object
    DomainLayout &DL = *DL_ptr_;
    // we want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    SurDomain_CathodeCollector *c_ptr = dynamic_cast<SurDomain_CathodeCollector *>(d_ptr);
    VarType v1(Voltage, 2, "CathodeVoltage");
    int num = c_ptr->changeDirichletConditionValue(v1, paramVal);
    if (num != 1) {
      throw m1d_Error("", "");
    }
    return num;
  }
  //====================================================================================================================
  int
  BatteryResidEval::getSolutionParam(std::string paramName, double * const paramVal) {
    if (paramName == "CathodeCollectorCurrent") {
      // Go find the Cathode Collector object -- the last domain of the DomainLayout
      DomainLayout &DL = *DL_ptr_;
      SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
      SurDomain_CathodeCollector *c_ptr = dynamic_cast<SurDomain_CathodeCollector *>(d_ptr);
      double icurr = c_ptr->icurrCollector_;
      paramVal[0] = icurr; 
      return 1;
    }
    else if (paramName == "SeparatorConductivity") {
      // Go find the Separator object -- the second domain of the DomainLayout
      DomainLayout &DL = *DL_ptr_;
      BulkDomain1D *d_ptr = DL.BulkDomain1D_List[1];
      porousLiKCl_dom1D *c_ptr = dynamic_cast<porousLiKCl_dom1D *>(d_ptr);
      double conductivity = c_ptr->getConductivity();
      paramVal[0] = conductivity; 
    }
    else {
      throw m1d_Error("BatteryResidEval::getSolutionParam", "unknown parameter");
    } 
    return 1;
  }
  //====================================================================================================================
  //   Report on the boundary condition applied to the cathode voltage equation
  /*
   *
   *   @param[in] time     Current time for evaluating time dependent BC
   *   @param[out] BC_Type  Type of the boundary condition
   *   @param[out] value    Value of the dirichlet condition or flux - default 0.0
   *   @param[out] BC_TimeDep BoundaryCondition Pointers for time dependent BC for BC_Tppe = 3,4,5
   *                   (default 0)
   *   @param[out] TimeDep  Function pointer to a function that returns a double given a single parameter (the time).
   *                   Defaults to a NULL pointer.
   *
   */
  void
  BatteryResidEval::reportCathodeVoltageBC(double time, int &BC_Type, double &value, BoundaryCondition * &BC_TimeDep,
					   TimeDepFunctionPtr &TimeDep) const {
    VarType v1(Voltage, 2, "CathodeVoltage");

    DomainLayout &DL = *DL_ptr_;
    // we want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    SurDomain_CathodeCollector *c_ptr = dynamic_cast<SurDomain_CathodeCollector *>(d_ptr);


    int num = c_ptr->reportBoundaryCondition(time, v1, BC_Type, value, BC_TimeDep, TimeDep);
    if (num != 1) {
      throw m1d_Error("", "");
    }
  }

  //====================================================================================================================
  void
  BatteryResidEval::changeCathodeVoltageBC(int BC_Type, double value, BoundaryCondition * BC_TimeDep,
					   TimeDepFunctionPtr TimeDep) {
    VarType v1(Voltage, 2, "CathodeVoltage");

    DomainLayout &DL = *DL_ptr_;
    // we want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    SurDomain_CathodeCollector *c_ptr = dynamic_cast<SurDomain_CathodeCollector *>(d_ptr);


    int num = c_ptr->changeBoundaryCondition(v1, BC_Type, value, BC_TimeDep, TimeDep);
    if (num != 1) {
      throw m1d_Error("", "");
    }
  }
  //====================================================================================================================
  double
  BatteryResidEval::reportCathodeVoltage() const {
    DomainLayout &DL = *DL_ptr_;
    // we want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    SurDomain_CathodeCollector *c_ptr = dynamic_cast<SurDomain_CathodeCollector *>(d_ptr);
    // might have to update the SurDomain.
    double phi =  c_ptr->phiCathode_;

    return phi;
  } 
  //====================================================================================================================
  double
  BatteryResidEval::reportCathodeCurrent() const {
    DomainLayout &DL = *DL_ptr_;
    // we want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    SurDomain_CathodeCollector *c_ptr = dynamic_cast<SurDomain_CathodeCollector *>(d_ptr);
    double icurr = c_ptr->icurrCollector_;
    return icurr;
  }
//====================================================================================================================
// Get the max value of the sub grid time step number from the last residual calculation
/*
 *   @return   Returns the max subgrid time step number from all of the electrodes. Special steps
 *             aren't counted in this number.
 */
int   BatteryResidEval::getMaxSubGridTimeSteps() const	
{
    return maxSubGridTimeSteps_;
}

//====================================================================================================================
void BatteryResidEval::gatherCapacityStatistics() 
{
    DomainLayout &DL = *DL_ptr_;
    BulkDomain1D *d_ptr = DL.BulkDomain1D_List[0];
    porousElectrode_dom1D*  d_anode_ptr= dynamic_cast<porousElectrode_dom1D*>(d_ptr);
    d_ptr = DL.BulkDomain1D_List[2];
    porousElectrode_dom1D* d_cathode_ptr = dynamic_cast<porousElectrode_dom1D*>(d_ptr);

    capacityAnodePA_   =   d_anode_ptr->capacityPA();
    capacityCathodePA_ = d_cathode_ptr->capacityPA();

    capacityLeftAnodePA_   =   d_anode_ptr->capacityLeftPA();
    capacityLeftCathodePA_ = d_cathode_ptr->capacityLeftPA();

    capacityDischargedAnodePA_   =   d_anode_ptr->capacityDischargedPA();
    capacityDischargedCathodePA_ = d_cathode_ptr->capacityDischargedPA();

    capacityPAEff_ = std::min(capacityAnodePA_, capacityCathodePA_);
    capacityLeftPAEff_ = std::min(capacityLeftAnodePA_, capacityLeftCathodePA_);

  
    double icurr = reportCathodeCurrent();

    double timeToDischargeWhole = 1.0E10;

    if (icurr  > 1.0E-5) {
	timeToDischargeWhole = capacityPAEff_ / icurr; 
    } else if (icurr < 1.0E-5) {
	timeToDischargeWhole = - capacityPAEff_ / icurr; 
    }
   
    timeLeft_ =  capacityLeftPAEff_ / icurr;
 
    Crate_current_ = 3600. /  timeToDischargeWhole;
}

//=====================================================================================================================
} // end of m1d namespace
//=====================================================================================================================
