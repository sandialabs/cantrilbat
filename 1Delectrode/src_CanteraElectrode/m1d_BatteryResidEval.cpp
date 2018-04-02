/**
 *  @file m1d_BatteryResidEval.cpp 
 *     definitions for the particular functions involved with calculating and analysing 1D battery problems
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
// @deprecated -> shouldn't be needed or called here porousElectrode is the only one
#include "m1d_porousLiKCl_dom1D.h"
#include "m1d_porousElectrode_dom1D.h"

#include "m1d_DomainLayout.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_EpetraExtras.h"
#include "m1d_GlobalIndices.h"
#include "m1d_globals.h"
#include "m1d_CanteraElectrodeGlobals.h"
#include "m1d_SurfDomainDescription.h"

#include "m1d_NodalVars.h"

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! boolean indicating that the polarization index calculations are current
bool polIndecisesCurrent = false;

//==================================================================================================================================
//! Draw a line on the Zuzax log file
/*!
 *  No flushing is carried out.
 *
 *  @param[in]              sp                   Number of indent characters
 *  @param[in]              ll                   Number "-" characters in the line
 */
static void drawline(int sp, int ll)
{
    for (int i = 0; i < sp; i++) {
        ZZCantera::writelog(" ");
    }
    for (int i = 0; i < ll; i++) {
        ZZCantera::writelog("-");
    }
    ZZCantera::writelog("\n");
}
//==================================================================================================================================
// Default constructor
/*
 *
 * @param atol   Absolute tolerance for calculation
 */
BatteryResidEval::BatteryResidEval(double atol) :
    ProblemResidEval(atol),
    doHeatSourceTracking_(0),
    doPolarizationAnalysis_(0),
    doResistanceTracking_(0),
    anodeType_(0),
    cathodeType_(0),
    maxSubGridTimeSteps_(0),
    QdotPerArea_n_(0.0),
    QdotPerArea_nm1_(0.0),
    QdotAnodePerArea_n_(0.0),
    QdotSeparatorPerArea_n_(0.0),
    QdotCathodePerArea_n_(0.0),
    capacityAnodePA_(0.0),
    capacityCathodePA_(0.0),
    capacityLeftAnodePA_(0.0),
    capacityLeftCathodePA_(0.0),
    capacityDischargedAnodePA_(0.0),
    capacityDischargedCathodePA_(0.0),
    dodAnodePA_(0.0),
    dodCathodePA_(0.0),
    capacityPAEff_(0.0),
    capacityLeftPAEff_(0.0),
    timeLeft_(0.0),
    Crate_current_(0.0),
    ocvAnode_(0.0),
    ocvCathode_(0.0),
    hasPressureEquation_(false),
    hasGasResevoir_(false),
    Anode_BDI_(0),
    Separator_BDI_(1),
    Cathode_BDI_(2)
{
     doHeatSourceTracking_ = PSCinput_ptr->doHeatSourceTracking_;
     doPolarizationAnalysis_ = PSCinput_ptr->doPolarizationAnalysis_;
     doResistanceTracking_ = PSCinput_ptr->doResistanceTracking_;
     crossSectionalArea_ = PSCinput_ptr->cathode_input_->electrodeGrossArea;

     // Determine whether there is a darcy formulation
     if (PSCinput_ptr->Pressure_formulation_prob_type_ >= 1) {
	 hasPressureEquation_ = true;
     }

     // Determine whether there is a gas resevoir
     if (PSCinput_ptr->Pressure_formulation_prob_type_ >= 2) {
	 hasGasResevoir_ = true;
     }
}
//================================================================================================================================
BatteryResidEval::~BatteryResidEval()
{
}
//================================================================================================================================
BatteryResidEval::BatteryResidEval(const BatteryResidEval &r) :
    ProblemResidEval(r.m_atolDefault)
{
    *this = r;
}
//==================================================================================================================================
BatteryResidEval &
BatteryResidEval::operator=(const BatteryResidEval &r)
{
    if (this == &r) {
	return *this;
    }

    ProblemResidEval::operator=(r);

    doHeatSourceTracking_              = r.doHeatSourceTracking_;
    doPolarizationAnalysis_            = r.doPolarizationAnalysis_;
    doResistanceTracking_              = r.doResistanceTracking_;
    anodeType_                         = r.anodeType_;
    cathodeType_                       = r.cathodeType_;
    maxSubGridTimeSteps_               = r.maxSubGridTimeSteps_;
    QdotPerArea_n_                     = r.QdotPerArea_n_;
    QdotPerArea_nm1_                   = r.QdotPerArea_nm1_;
    QdotAnodePerArea_n_                = r.QdotAnodePerArea_n_;
    QdotSeparatorPerArea_n_            = r.QdotSeparatorPerArea_n_;
    QdotCathodePerArea_n_              = r.QdotCathodePerArea_n_;

    capacityAnodePA_                     = r.capacityAnodePA_;
    capacityCathodePA_                   = r.capacityCathodePA_;
    capacityLeftAnodePA_                 = r.capacityLeftAnodePA_;
    capacityLeftCathodePA_               = r.capacityLeftCathodePA_;
    capacityDischargedAnodePA_           = r.capacityDischargedAnodePA_;
    capacityDischargedCathodePA_         = r.capacityDischargedCathodePA_;
    dodAnodePA_                          = r.dodAnodePA_;
    dodCathodePA_                        = r.dodCathodePA_;
    capacityPAEff_                       = r.capacityPAEff_;
    capacityLeftPAEff_                   = r.capacityLeftPAEff_;
    timeLeft_                            = r.timeLeft_;
    Crate_current_                       = r.Crate_current_;
    ocvAnode_                            = r.ocvAnode_;
    ocvCathode_                          = r.ocvCathode_;
    hasPressureEquation_                 = r.hasPressureEquation_;
    hasGasResevoir_                      = r.hasGasResevoir_;
    Anode_BDI_                           = r.Anode_BDI_;
    Separator_BDI_                       = r.Separator_BDI_;
    Cathode_BDI_                         = r.Cathode_BDI_;

    return *this;
}
//==================================================================================================================================
void
BatteryResidEval::residSetupTmps()
{
    //
    //  This routine needs a lot more sophistication 
    //    Need to run it whenever the tmp info needs to be recalculated.
    //
    DomainLayout &DL = *DL_ptr_;
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	porousFlow_dom1D* pd_ptr = dynamic_cast<porousFlow_dom1D*>(d_ptr);
	if (pd_ptr != 0) {
	    pd_ptr->residSetupTmps();
	}
    }
}
//==================================================================================================================================
void
BatteryResidEval::setStateFromSolution(const bool doTimeDependentResid, const Epetra_Vector_Ghosted *soln, 
				       const Epetra_Vector_Ghosted *solnDot,
				       const double t, const double delta_t, const double t_old)
{

    if (doTimeDependentResid) {
	if (delta_t > 0.0) {
	    if (fabs(t - (t_old + delta_t)) > 0.001 * delta_t) {
		throw m1d_Error("BatteryResidEval::setStateFromSolution", "case of t != (t_old + delta_t) not handled yet");
	    }
	}
    }
    //
    // Call the base class
    //
    ProblemResidEval::setStateFromSolution(doTimeDependentResid, soln, solnDot, t, delta_t, t_old);


    //Domain Layout ptr DL_ptr_
    //Domain layout has the information on the global nodes and the x-positions
    DomainLayout &DL = *DL_ptr_;
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	d_ptr->setStateFromSolution(doTimeDependentResid, soln, solnDot, t, delta_t, t_old);
    }
    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
	SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	d_ptr->setStateFromSolution(doTimeDependentResid, soln, solnDot, t, delta_t, t_old);
    }

}
//==================================================================================================================================
void
BatteryResidEval::initialConditions(const bool doTimeDependentResid, Epetra_Vector_Ghosted *soln,
				    Epetra_Vector_Ghosted *solnDot, double &t, double &delta_t, double &delta_t_np1)
{
    //
    //   We obtain the cross-sectional area from the cathode. 
    //   However, we require the cross-sectional area to be consistent across inputs
    //
    crossSectionalArea_ = PSCinput_ptr->cathode_input_->electrodeGrossArea;
    //
    //   loop over all volume and surface domains providing initial guess
    //
    ProblemResidEval::initialConditions(doTimeDependentResid, soln, solnDot, t, delta_t, delta_t_np1);
    //
    //  improveInitialConditions(soln);
    //
}
//==================================================================================================================================
// Improve upon initial conditions by computing first order potential losses, equilibrium reaction voltages, etc.
/* 
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
	int iVar_Voltage_BD = d_ptr->BDD_ptr_->VariableIndexStart_VarName[Voltage];

	for ( int iGbNode = leftGbNode; iGbNode < rightGbNode; iGbNode++ ) {
	    double xNode = (*GI_ptr_->XNodePos_GbNode)[iGbNode];	
	    int voltEqn = GI_ptr_->IndexStartGbEqns_GbNode[iGbNode] + iVar_Voltage_BD ;
	    fprintf(stderr, "Old var value was %g at x=%g\n", soln[voltEqn], xNode );

	}
    }

    if (    PSCinput_ptr->cathodeBCType_ == 0 
	    || PSCinput_ptr->cathodeBCType_ == 2 
	    || PSCinput_ptr->cathodeBCType_ == 4 
	    || PSCinput_ptr->cathodeBCType_ == 6 
	    || PSCinput_ptr->cathodeBCType_ == 8 ) {
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
	iVar_Voltage_BD = d_ptr->BDD_ptr_->VariableIndexStart_VarName[Voltage];
	voltEqn = GI_ptr_->IndexStartGbEqns_GbNode[iNodeLeft] + iVar_Voltage_BD ;
	double voltLeft = soln[voltEqn];
      
	//! Right-most bulk domain
	iDom = DL.NumBulkDomains - 1 ;
	//The right hand node for interpolation is the 2nd node of the right domain (near it's left edge)
	int iNodeRight = DL_ptr_->StartGBNode_Domain[iDom] + 1;
	double xright = (*GI_ptr_->XNodePos_GbNode)[iNodeRight];	
	//find the voltage at this node
	d_ptr = DL.BulkDomain1D_List[iDom];
	iVar_Voltage_BD = d_ptr->BDD_ptr_->VariableIndexStart_VarName[Voltage];
	voltEqn = GI_ptr_->IndexStartGbEqns_GbNode[iNodeRight] + iVar_Voltage_BD ;
	double voltRight = soln[voltEqn];
      
	//! Loop over intermediate domains
	//   Loop over the Volume Domains
	for (iDom = 1; iDom < DL.NumBulkDomains - 1; iDom++) {
	
	    //nodes at boundaries
	    int leftGbNode  = DL_ptr_->StartGBNode_Domain[iDom];
	    int rightGbNode = DL_ptr_->EndGBNode_Domain[iDom];
	
	    BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	    int iVar_Voltage_BD = d_ptr->BDD_ptr_->VariableIndexStart_VarName[Voltage];
	
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
	iVar_Voltage_BD = d_ptr->BDD_ptr_->VariableIndexStart_VarName[Voltage];
      
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

//==================================================================================================================================
void
BatteryResidEval::residEval(Epetra_Vector_Owned* const & res,
                            const bool doTimeDependentResid,
                            const Epetra_Vector *soln_ptr,
                            const Epetra_Vector *solnDot_ptr,
                            const double t,
                            const double rdelta_t,
                            const Zuzax::ResidEval_Type residType,
                            const Zuzax::Solve_Type solveType)
{
    if (!resInternal_ptr_) {
	resInternal_ptr_ = new Epetra_Vector(*res);
    }

    if (residType == Zuzax::ResidEval_Type::Base_ResidEval) {
	counterResBaseCalcs_++;
    } else if (residType == Zuzax::ResidEval_Type::JacBase_ResidEval) {
	counterJacBaseCalcs_++;
    } else if (residType == Zuzax::ResidEval_Type::JacDelta_ResidEval) {
	counterJacDeltaCalcs_++;
    } else if (residType == Zuzax::ResidEval_Type::Base_ShowSolution) {
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
    double delta_t = 0.0;
    double t_old = t;
    if (rdelta_t != 0.0) {
        delta_t = 1.0 / rdelta_t;
        t_old = t - delta_t;
    }
    if (doTimeDependentResid) {
	calcSolnOld(*soln_ptr, *solnDot_ptr, rdelta_t);
    }
    /*
     *   Propagate the solution of the system down to the underlying objects where necessary.
     *   It is necessary to do this for mesh unknowns.
     *   Also do a loop over the nodes to carry out any precalculations that are necessary
     */
    setStateFromSolution(doTimeDependentResid, soln_ptr, solnDot_ptr, t, delta_t, t_old);
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
	    maxSubGridTimeSteps_ = std::max(maxSubGridTimeSteps_, maxE);
        }
    }

}
//==================================================================================================================================
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
    doPolarizationAnalysis_ = true; 
    if (doPolarizationAnalysis_) {

        doPolarizationAnalysis(ievent, time_current, delta_t_n, y_n, ydot_n_ptr);

    }
    // When heat source tracking is turned on and the enthalpy equation is used, we do a full tracking of species and enthalpy 
    if (energyEquationProbType_ == 3 && doHeatSourceTracking_) {

        doSpeciesAnalysis(ievent, time_current, delta_t_n, y_n, ydot_n_ptr);

	doHeatAnalysis(ievent, time_current, delta_t_n, y_n, ydot_n_ptr);
    }
}
//==================================================================================================================================
//! print an appended line into a buffer
/*!
 *  @param[in, out]          buf                 Buffer to append the line to
 *  @param[in]               st                  string to write the line
 *  @param[in]               num                 Number of times to write
 */
static void sprint_line(char * buf, const char * const st, const int num)
{
    int n = strlen(buf);
    buf += n;
    for (int k = 0; k < num; k++, buf++) {
	sprintf(buf, "%s", st);
    }
    sprintf(buf, "\n");
}
//==================================================================================================================================
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
                                      const Zuzax::Solve_Type solveType,
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

    residEval(resInternal_ptr_, doTimeDependentResid, &y_n, ydot_n, t, rdelta_t,
              Zuzax::ResidEval_Type::Base_ShowSolution, solveType);

    gatherCapacityStatistics();


    int indentSpaces = 4;
    string indent = "    ";
    const char *ind = indent.c_str();
    if (!mypid || duplicateOnAllProcs) {
	sprintf(buf, "%s", ind);
	sprint_line(buf, "-", 100);
	ZZCantera::writelog(buf);
	sprintf(buf, "%s ShowProblemSolution : Time       %-12.3E\n", ind, t);
	ZZCantera::writelog(buf);
	if (solveType == Zuzax::Solve_Type::TimeDependentInitial) {
	    sprintf(buf, "%s                       Delta_t    %-12.3E  (initial solution with no previous solution)\n", ind, delta_t);
	} else {
	    sprintf(buf, "%s                       Delta_t    %-12.3E\n", ind, delta_t);
	}
	ZZCantera::writelog(buf);
	sprintf(buf, "%s                       StepNumber %6d\n", ind, m_StepNumber);
	ZZCantera::writelog(buf);
	sprintf(buf, "%s                       Delta_t_p1 %-12.3E\n", ind, delta_t_np1);
	ZZCantera::writelog(buf);
	sprintf(buf, "%s                       Capacity_Anode        %-12.3E  Amp hr / m2\n", ind, capacityAnodePA_ / 3600.);
	ZZCantera::writelog(buf);
	sprintf(buf, "%s                       Capacity_Anode_Left   %-12.3E  Amp hr / m2\n", ind, capacityLeftAnodePA_ / 3600.);
	ZZCantera::writelog(buf);
	sprintf(buf, "%s                       DepthDischarge_Anode  %-12.3E  Amp hr / m2\n", ind, dodAnodePA_ / 3600.);
	ZZCantera::writelog(buf);
	sprintf(buf, "%s                       Capacity_Cathode      %-12.3E  Amp hr / m2\n", ind, capacityCathodePA_ / 3600.);
	ZZCantera::writelog(buf);
	sprintf(buf, "%s                       Capacity_Cathode_Left %-12.3E  Amp hr / m2\n", ind, capacityLeftCathodePA_ / 3600.);
	ZZCantera::writelog(buf);
	sprintf(buf, "%s                       DepthDischarge_Cathode%-12.3E  Amp hr / m2\n", ind, dodCathodePA_ / 3600.);
	ZZCantera::writelog(buf);
	sprintf(buf, "%s                       Crate_current         %-12.3E  \n", ind, Crate_current_);
	ZZCantera::writelog(buf);

	sprintf(buf, "%s", ind);
	sprint_line(buf, "-", 100);
	ZZCantera::writelog(buf);
    }

    Domain1D *d_ptr = DL.SurDomain1D_List[0];
    do {
	d_ptr->showSolution(solnAll, soln_dot_All, &y_n, ydot_n, solnOld_ptr_,
			    resInternal_ptr_, t, rdelta_t, indentSpaces + 2,
			    duplicateOnAllProcs);
	BulkDomain1D *bd_ptr = dynamic_cast<BulkDomain1D *> (d_ptr);
	if (bd_ptr) {
	    //BulkDomainDescription &BDD_;
	    SurfDomainDescription *sdd = bd_ptr->BDD_ptr_->RightSurf;
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
	ZZCantera::writelog(buf);
    }

    Epetra_Comm *c = LI_ptr_->Comm_ptr_;
    c->Barrier();
}
//==================================================================================================================================
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
				const Zuzax::Solve_Type solveType, 
				const double delta_t_np1)
{
    ProblemResidEval::writeSolution(ievent, doTimeDependentResid, time_current, delta_t_n, istep, y_n,
				    ydot_n_ptr, solveType, delta_t_np1);
    if (ievent == 0 || ievent == 1 || ievent == 2) {
	write_IV(ievent, doTimeDependentResid, time_current, delta_t_n, istep, y_n, ydot_n_ptr);
    }

    writeGlobalTecplot(ievent, doTimeDependentResid, time_current, delta_t_n, istep, y_n,
				    ydot_n_ptr, solveType, delta_t_np1);
}
//==================================================================================================================================
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
    double TempCathodeCollector = c_ptr->TempCollector;

    //
    // The first domain is the surface anode collector
    //  get the current -> its a member data
    //
    SurDomain1D *ad_ptr = DL.SurDomain1D_List[0];
    SurDomain_AnodeCollector *ac_ptr = dynamic_cast<SurDomain_AnodeCollector *>(ad_ptr);
    double a_icurr = ac_ptr->icurrCollector_;

    int procID = Comm_ptr->MyPID();

    if (!procID) {
	//
	//  looking for cathode capacity and depth of discharge
	//
	//BulkDomain1D *cd_ptr = DL.BulkDomain1D_List.back();

	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[0];
	porousElectrode_dom1D*  d_anode_ptr= dynamic_cast<porousElectrode_dom1D*>(d_ptr);
	d_ptr = DL.BulkDomain1D_List.back();
	porousElectrode_dom1D* d_cathode_ptr = dynamic_cast<porousElectrode_dom1D*>(d_ptr);

	double capacityZeroDoD, spec_capacityZeroDoD;
	double dischargedCapacity, spec_dischargedCapacity;
	d_cathode_ptr->reportSolutionParam( "CapacityZeroDoD", &capacityZeroDoD );
	d_cathode_ptr->reportSolutionParam( "DepthOfDischarge", &dischargedCapacity );
	d_cathode_ptr->reportSolutionParam( "SpecificCapacityZeroDoD", &spec_capacityZeroDoD );
	d_cathode_ptr->reportSolutionParam( "SpecificDepthOfDischarge", &spec_dischargedCapacity );

	ocvAnode_ = d_anode_ptr->openCircuitPotentialQuick();
	ocvCathode_ = d_cathode_ptr->openCircuitPotentialQuick();
	double ocvQuick = ocvCathode_ - ocvAnode_;   

	FILE *fp = 0;
	bool doOldFormat = false;
	if (doOldFormat) {
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
		    capacityZeroDoD,  dischargedCapacity, 0.1 * a_icurr);
	} else {
	    if (ievent == 0) {
		fp = fopen("timeDep_IV.dat", "w");
		fprintf(fp, "TITLE = \"Time versus Current or Voltage\"\n");
		fprintf(fp, "VARIABLES = \" T [s]\"\n");
		fprintf(fp, "\"Voltage [volts] \"\n");
		fprintf(fp, "\"CathodeCurrent [A/m2]\"\n");
		fprintf(fp, "\"Cathode Capacity [A-s/m2] \"\n");
		fprintf(fp, "\"Cathode DepthOfDischarge [A-s/m2] \"\n");
		fprintf(fp, "\"Discharged Cathode Capacity [A-s/m2] \"\n");
		fprintf(fp, "\"Cathode Capacity Left [A-s/m2] \"\n");
		fprintf(fp, "\"AnodeCurrent [A/m2]\"\n");
		fprintf(fp, "\"Anode Capacity [A-s/m2] \"\n");
		fprintf(fp, "\"Anode DepthOfDischarge [A-s/m2] \"\n");
		fprintf(fp, "\"Discharged Anode Capacity [A-s/m2] \"\n");
		fprintf(fp, "\"Anode Capacity Left [A-s/m2] \"\n");
		fprintf(fp, "\"OCV_Quick [volts] \"\n");
                fprintf(fp, "\"TemperatureCathode [K]\"\n");

	    } else {
		fp = fopen("timeDep_IV.dat", "a");
	    }
	    fprintf(fp, "   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E ,  %15.5E ,"
		    "  %15.5E,  %15.5E, %15.5E ,  %15.5E,   %15.5E\n", 
		    time_current, phiCath, c_icurr, capacityCathodePA_, dodCathodePA_,
		    capacityDischargedCathodePA_, capacityLeftCathodePA_, a_icurr,
		    capacityAnodePA_, dodAnodePA_, capacityDischargedAnodePA_, capacityLeftAnodePA_, ocvQuick, TempCathodeCollector);
	}
	fclose(fp);
	//
	if (ievent == 0) { 
	    fp = fopen("Capacity_Starting.dat", "w"); 
	    fprintf(fp, " starting time = %g\n", time_current);
	    fprintf(fp, "  Quantity             Units           Anode                Separator           Cathode\n");
	    fprintf(fp, "  InitialCapacity      [A-s/m2]    %20.10E %20.10E %20.10E\n", capacityAnodePA_, 0.0, capacityCathodePA_);
	    fprintf(fp, "  InitialDepthDischar  [A-s/m2]    %20.10E %20.10E %20.10E\n", dodAnodePA_,      0.0, dodCathodePA_);

	    fclose(fp);
      
	}
    }
    Comm_ptr->Barrier();
}
//==================================================================================================================================
void
BatteryResidEval::writeGlobalTecplot(const int ievent, const bool doTimeDependentResid, const double time_current,
				     const double delta_t_n, int istep, const Epetra_Vector_Ghosted &y_n, 
                                     const Epetra_Vector_Ghosted * const ydot_n_ptr, const Zuzax::Solve_Type solveType, 
				     const double delta_t_np1)
{
    static int headerWritten = false;
    int numRtn;
    // Create a communications vector
    std::vector<double> volInfoVector;
    std::vector<double> varsVector;
    std::string requestID;
    std::string name;
    VarType vt;
    int requestType;
    if (!headerWritten) {
	headerWritten = true;
	writeGlobalTecplotHeader(ievent, doTimeDependentResid, time_current, delta_t_n, istep, y_n,
				 ydot_n_ptr, solveType, delta_t_np1);
    }
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid ; //only proc 0 should write out the header
    DomainLayout* dl = DL_ptr_;
    BulkDomainDescription* bdd0 = dl->BulkDomainDesc_global[0];
    if (doWrite) {
	// Open Tecplot file
	FILE* ofp;
	string sss = getBaseFileName();
	char filename[30];
	sprintf(filename,"%s%s",sss.c_str(),"_multiBulk.dat");
	ofp = fopen(filename, "a");
	// Number of equations per node
	int numVar = bdd0->NumEquationsPerNode;
	// Vector containing the variable names as they appear in the solution vector
	std::vector<VarType> &variableNameList = bdd0->VariableNameList;
	//! First global node of this bulk domain

	double numNodes = 0;
	for (size_t iDom = 0; iDom < (size_t) dl->NumBulkDomains; ++iDom) {
	    // BulkDomain1D *d_ptr = dl->BulkDomain1D_List[iDom];
	    BulkDomainDescription* bdd_ptr = dl->BulkDomainDesc_global[iDom];
	    int numNodesDD = bdd_ptr->LastGbNode - bdd_ptr->FirstGbNode + 1;
	    numNodes += numNodesDD;
	}
	//
	//  Write the Zone header information
	//
        fprintf(ofp, "ZONE T = \"All, %10.4E\", I = %d, SOLUTIONTIME = %12.6E\n", time_current,
		(int)  numNodes, time_current);
        fprintf(ofp, "ZONETYPE = ORDERED\n");
        fprintf(ofp, "DATAPACKING = BLOCK\n");
        fprintf(ofp, "STRANDID = 1\n");

       

	for (int k = -1; k < numVar; k++) {
	    if (k == -1) {
		name = "x [m]";
		requestType = 1;
	    } else {
		vt = variableNameList[k];
	        name = vt.VariableName();
		requestType = 0;
	    }
	    requestID = name;
	  
	    /*
	     *   Loop over the Volume Domains
	     */
	    varsVector.clear();
	    for (int iDom = 0; iDom < dl->NumBulkDomains; iDom++) {
		BulkDomain1D *d_ptr = dl->BulkDomain1D_List[iDom];
		BulkDomainDescription* bdd_ptr = dl->BulkDomainDesc_global[iDom];
		int numNodesDD = bdd_ptr->LastGbNode - bdd_ptr->FirstGbNode + 1;
		requestID = name;
		if (name ==  "Volt(AnodeVoltage)") {
		    if (iDom == 2) {
			requestID = "Volt(CathodeVoltage)";
		    }
		}
		if ((requestID ==  "Volt(AnodeVoltage)") && iDom == 1) {
		    numRtn =  numNodesDD;
		    volInfoVector.clear();
		    volInfoVector.resize(numNodesDD, 0.0);
		} else {
		    numRtn = d_ptr->reportSolutionVector(requestID, requestType, &y_n, volInfoVector);
		}
		/// CBL This is a hack. Fix it for real \todo
		if(numRtn != numNodesDD && ( name == "Displacement_Axial" || name == "Solid Stress Axial"))
		  numRtn = numNodesDD;
		if (numRtn != numNodesDD) {
		    throw m1d_Error("BatteryResidEval::writeGlobalTecplot", 
				    "Can't process requestID " + requestID + " from domain " + int2str(iDom));
		}
		varsVector.insert(varsVector.end(), volInfoVector.begin(), volInfoVector.end());	
	    }
	    fwriteTecplotVector(ofp, varsVector, 13, 10);
	}
	fprintf(ofp, "\n" );
	fclose(ofp);
    }
}
//==================================================================================================================================
void
BatteryResidEval::writeGlobalTecplotHeader(const int ievent,
					   const bool doTimeDependentResid,
					   const double time_current,
					   const double delta_t_n,
					   int istep,
					   const Epetra_Vector_Ghosted &y_n,
					   const Epetra_Vector_Ghosted * const ydot_n_ptr,
					   const Zuzax::Solve_Type solveType, 
					   const double delta_t_np1)
{
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid ; //only proc 0 should write out the header
    DomainLayout* dl = DL_ptr_;
    BulkDomainDescription* bdd0 = dl->BulkDomainDesc_global[0];
    if (doWrite) {
	//open tecplot file
	FILE* ofp;
	string sss = getBaseFileName();
	char filename[30];
	sprintf(filename,"%s%s",sss.c_str(),"_multiBulk.dat");
	ofp = fopen(filename, "w");

	//write title and variable list
	fprintf( ofp, "TITLE = \"Tecplot Solution on multiple bulk domains for %s\"\n",sss.c_str());

	// Number of equations per node
	int numVar = bdd0->NumEquationsPerNode;
	// Vector containing the variable names as they appear in the solution vector
	std::vector<VarType> &variableNameList = bdd0->VariableNameList;
	//! First global node of this bulk domain

	fprintf( ofp, "VARIABLES = ");
	fprintf( ofp, "\"x [m]\"  \n" );

	for (int k = 0; k < numVar; k++) {
	    VarType &vt = variableNameList[k];
	    string name = vt.VariableName(128);
	    if (name == "Volt(AnodeVoltage)") {
		name = "Volt(Electrode)";
	    }
	    fprintf( ofp, "\"%s\" \n", name.c_str() );
	}
	fprintf(ofp, "\n" );
	fclose(ofp);
    }
}
//==================================================================================================================================
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
    static FILE *fstepE = 0;
    FILE *facc = 0;
    FILE *faccE = 0;
    static double tinit_calculation = 0.0;
    const Epetra_Vector_Ghosted *soln_ptr = &y;
    DomainLayout &DL = *DL_ptr_;
    bool doTimeDependentResid = true;
    double rdelta_t = 0.0;
    if (deltaT > 0.0) {
	rdelta_t = 1.0 / deltaT;
    }
    /*
     * We calculate solnOld_ptr_ here
     */
    if (doTimeDependentResid) {
        calcSolnOld(*soln_ptr, *solnDot_ptr, rdelta_t);
    }

    bool need_Post =true;
    if (doHeatSourceTracking_ || doResistanceTracking_) {
	need_Post = true;
    }
    if (need_Post) {

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
    }

    if (doHeatSourceTracking_) {
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
        
    if (doResistanceTracking_) {

	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[0];
	porousElectrode_dom1D*  d_anode_ptr= dynamic_cast<porousElectrode_dom1D*>(d_ptr);
	d_ptr = DL.BulkDomain1D_List.back();
	porousElectrode_dom1D* d_cathode_ptr = dynamic_cast<porousElectrode_dom1D*>(d_ptr);
	ocvAnode_ = d_anode_ptr->openCircuitPotentialQuick();
	ocvCathode_ = d_cathode_ptr->openCircuitPotentialQuick();
	d_ptr = DL.BulkDomain1D_List[1];
	porousFlow_dom1D* d_separator_ptr= dynamic_cast<porousFlow_dom1D*>(d_ptr);

	double tinit = t - deltaT;
        if (ifunc == 0) {
            fstepE = fopen("stepwiseElectricalOutput.txt", "w");
          
            fprintf(fstepE, "       tinit  ,       tfinal  ,  Volts_Anode ,    Volts_Sep , Volts_Cathode,"
		    "  Volts_Total ,    Volts_OCV ,   OCV_Anode  ,  OCV_Cathode , "
		    "     Current , Resist_anode ,  Resist_sep  , Resist_cathod,  Resist_Total\n");
            fflush(fstepE);
            fclose(fstepE);
        
	    faccE = fopen("accumulatedElectricalOutput.txt", "w");
            fprintf(faccE, "       tinit  ,       tfinal  ,  Volts_Anode ,    Volts_Sep , Volts_Cathode,"
		    "  Volts_Total ,    Volts_OCV ,   OCV_Anode  ,  OCV_Cathode , "
		    "     Current , Resist_anode ,  Resist_sep  , Resist_cathod,  Resist_Total\n");
	    fclose(faccE);
	}

        if (ifunc == 1) {
	    double pot1, pot2, current;
	    double  volts_OCV_Anode,  volts_OCV_Cathode;
            double resistAnode = d_anode_ptr->effResistanceLayer(pot1, pot2, volts_OCV_Anode, current);
	    double pot_AnodeColl = pot1;
	    double volts_Anode = pot1 - pot2;
	    double resistCathode = d_cathode_ptr->effResistanceLayer(pot1, pot2, volts_OCV_Cathode, current);
	    double volts_Cathode = pot2 - pot1;
	    double pot_CathodeColl = pot2;
	    double volts_OCV_Separator;
            double volts_Total = pot_CathodeColl - pot_AnodeColl;
	    double resistSeparator = d_separator_ptr->effResistanceLayer(pot1, pot2, volts_OCV_Separator, current);
	    double volts_Separator = pot1 - pot2;

	    double volts_OCV = volts_OCV_Cathode - volts_OCV_Anode;
	    double resistTotal = 0.0;
	    if (fabs(current) > 1.0E-100) {
		resistTotal = (volts_OCV - (pot_CathodeColl - pot_AnodeColl)) / current;
	    } 

            fstepE = fopen("stepwiseElectricalOutput.txt", "a");
	    fprintf(fstepE, "% 13.5E , % 13.5E ,", tinit,  t);

	    fprintf(fstepE, "% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,",
		    volts_Anode, volts_Separator, volts_Cathode, volts_Total, volts_OCV,  volts_OCV_Anode,  volts_OCV_Cathode);
	    fprintf(fstepE, "% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E\n",
		    current, resistAnode, resistSeparator,  resistCathode, resistTotal);

            fclose(fstepE);
        }

	if (ifunc == 2) {   
            // MP issues
            faccE = fopen("accumulatedElectricalOutput.txt", "a");
	    fprintf(faccE, "% 13.5E , % 13.5E ,", tinit_calculation,  t);
            double pot1, pot2, current;
	    double  volts_OCV_Anode,  volts_OCV_Cathode;
            double resistAnode = d_anode_ptr->effResistanceLayer(pot1, pot2, volts_OCV_Anode, current);
	    double pot_AnodeColl = pot1;
	    double volts_Anode = pot1 - pot2;
	    double resistCathode = d_cathode_ptr->effResistanceLayer(pot1, pot2, volts_OCV_Cathode, current);
	    double volts_Cathode = pot2 - pot1;
	    double pot_CathodeColl = pot2;
	    double volts_OCV_Separator;
            double volts_Total = pot_CathodeColl - pot_AnodeColl;
	    double resistSeparator = d_separator_ptr->effResistanceLayer(pot1, pot2, volts_OCV_Separator, current);
	    double volts_Separator = pot1 - pot2;

	    double volts_OCV = volts_OCV_Cathode - volts_OCV_Anode;
	    double resistTotal = 0.0;
	    if (fabs(current) > 1.0E-100) {
		resistTotal = (volts_OCV - (pot_CathodeColl - pot_AnodeColl)) / current;
	    } 
	    fprintf(faccE, "% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,",
		    volts_Anode, volts_Separator, volts_Cathode, volts_Total, volts_OCV,  volts_OCV_Anode,  volts_OCV_Cathode);
	    fprintf(faccE, "% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E ,% 13.5E\n",
		    current, resistAnode, resistSeparator,  resistCathode, resistTotal);
            fclose(faccE);
        }
    }
}
//==================================================================================================================================
void BatteryResidEval::doHeatAnalysis(const int ifunc, const double t, const double deltaT, const Epetra_Vector_Ghosted &y,
	                              const Epetra_Vector_Ghosted * const solnDot_ptr)
{
    class globalHeatBalValsBat dVals;
    DomainLayout &DL = *DL_ptr_;
    double HeatFluxRight = 0.0;
    double JHelecRight = 0.0;
    double currRight = 0.0;
    double phiCath = 0.0;
    double HeatFluxLeft = 0.0;
    double JHelecLeft = 0.0;
    double currLeft = 0.0;
    double phiAnode = 0.0;
    double nEnthOldB[10];
    double nEnthNewB[10];
    double nEnthOldTotal = 0.0;
    double nEnthNewTotal = 0.0;
    //double enthalpyIVfluxRight  = 0.0;
    //double enthalpyIVfluxLeft  = 0.0;
    double enthalpyFlowOut = 0.0;
    double sourceTermExtra = 0.0;
    double sourceTermExtraCath = 0.0;

    //
    //   Loop over the Volume Domains
    //
    double totalHeatCapacity = 0.0;
    for (size_t iDom = 0; iDom < (size_t) DL.NumBulkDomains; iDom++) {
        dVals.zero();
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	d_ptr->eval_HeatBalance(ifunc, t, deltaT, &y, solnDot_ptr, solnOld_ptr_, dVals);
	totalHeatCapacity += dVals.totalHeatCapacity;
        nEnthNewB[iDom] = dVals.newNEnthalpy;
        nEnthOldB[iDom] = dVals.oldNEnthalpy;
        nEnthNewTotal +=  dVals.newNEnthalpy;
        nEnthOldTotal +=  dVals.oldNEnthalpy;
	sourceTermExtra += dVals.sourceTermExtra;
        if (iDom == 0) {
            JHelecLeft = dVals.JHelecLeft;
        }
	if (iDom == 2) {
            JHelecRight = dVals.JHelecRight;
	    //enthalpyIVfluxRight = dVals.enthalpyIVfluxRight;
	    HeatFluxRight = dVals.HeatFluxRight;
	    enthalpyFlowOut = dVals.enthFluxOut;
	}
    }
    //
    //    Loop over the Surface Domains
    //
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
        dVals.zero();
	SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	d_ptr->eval_HeatBalance(ifunc, t, deltaT, &y, solnDot_ptr, solnOld_ptr_, dVals);

        totalHeatCapacity += dVals.totalHeatCapacity;
        //nEnthOldS[iDom] = dVals.oldNEnthalpy;
        //nEnthNewS[iDom] = dVals.newNEnthalpy;
        nEnthNewTotal +=  dVals.newNEnthalpy;
        nEnthOldTotal +=  dVals.oldNEnthalpy;
	sourceTermExtra += dVals.sourceTermExtra;

        if (iDom == 0) {
            HeatFluxLeft = dVals.HeatFluxLeft;
            currLeft = dVals.currentLeft;
            phiAnode = dVals.phiSolid;
        }

        if (iDom == 3) {
	    HeatFluxRight = dVals.HeatFluxRight;
	    currRight = dVals.currentRight;
	    phiCath = dVals.phiSolid;
            sourceTermExtraCath = dVals.sourceTermExtra;
        }
    }

    printf(" Total Heat Capacity = %g\n", 	totalHeatCapacity);

    printf(" Heat flux Cathode = %g\n", HeatFluxRight);

    printf(" JHelecRight = %g\n", JHelecRight);

    printf(" current Volts = %g\n", currRight * phiCath); 

    printf(" Heat flux Anode = %g\n", HeatFluxLeft);

    printf(" JHelecLeft = %g\n", JHelecLeft);

    printf(" current Volts anode = %g\n", currLeft * phiAnode);

    printf(" Old Enthalpy = %g \n", nEnthOldTotal); 

    printf(" New Enthalpy = %g \n", nEnthNewTotal);

    // extend this to surfaces when they start to have enthalpy

    printf("\n");
    printf("   Region      NewTotalEnthalpy  OldTotalEnthalpy  DeltaTotalEnth    IVWork_Out*delT"
	   "Heat_Out*delT       enthFlowOut*delT      SourceExtra\n");
    printf("  ---------------------------------------------------------------------------------"
	   "----------------------------------------------------\n");
    for (size_t iDom = 0; iDom < (size_t) DL.NumBulkDomains; iDom++) {
        if (iDom == 0) printf("  Anode        ");
        if (iDom == 1) printf("  Separator    ");
        if (iDom == 2) printf("  Cathode      ");
        printf(" % 12.6E    % 12.6E       % 12.6E" , nEnthNewB[iDom], nEnthOldB[iDom],   nEnthNewB[iDom] - nEnthOldB[iDom]);
	if (iDom == 0) {
	    printf(" % 12.6E    % 12.6E       % 12.6E    % 12.6E" ,   - deltaT * JHelecLeft, deltaT * HeatFluxLeft, 0.0, 0.0);
	    //printf(" % 12.6E    % 12.6E       % 12.6E" ,   deltaT * enthalpyIVfluxLeft, deltaT * HeatFluxLeft, 0.0);
	}
	if (iDom == 1) {
	    printf(" % 12.6E    % 12.6E       % 12.6E    % 12.6E" , 0.0, 0.0, 0.0, 0.0);
	}
	if (iDom == 2) {
	    printf(" % 12.6E    % 12.6E       % 12.6E    % 12.6E" ,  
		   deltaT * JHelecRight, deltaT * HeatFluxRight, deltaT * enthalpyFlowOut, deltaT * sourceTermExtraCath);
	    // printf(" % 12.6E    % 12.6E       % 12.6E" ,  
	    //        deltaT * enthalpyIVfluxRight, deltaT * HeatFluxRight, deltaT * enthalpyFlowOut);
	}
	printf("\n");
    }
    printf("-----------------------------------------------------------------------------------------"
	   "----------------------------------------------\n");
    //double IVWorkdT =  deltaT * (enthalpyIVfluxRight +  enthalpyIVfluxLeft);
    double IVWorkdT =  deltaT * (JHelecRight - JHelecLeft);
    double HeatFluxdT =  deltaT * (HeatFluxRight +  HeatFluxLeft);
    double enthFlowdT = deltaT *  enthalpyFlowOut;
    double sourceExtradT = deltaT * sourceTermExtraCath;
    //if (sourceExtradT != 0.0) {
    //	printf(" sourceExtradT = %g\n", sourceExtradT);
    //}
    double enthLossdT =  IVWorkdT + HeatFluxdT + enthFlowdT;
   
    double delEnthtotal = nEnthNewTotal - nEnthOldTotal;

    double delHtotal = delEnthtotal + enthLossdT - sourceExtradT;
      
    printf("  Total     ");
    printf("    % 12.6E    % 12.6E       % 12.6E" , nEnthNewTotal, nEnthOldTotal,   nEnthNewTotal - nEnthOldTotal);
    printf(" % 12.6E    % 12.6E       % 12.6E    % 12.6E", IVWorkdT , HeatFluxdT, enthFlowdT, sourceExtradT);
    printf("\n");
    printf("                                 % 12.6E \n",  delHtotal);

    printf("\n\n");
    
}
//==================================================================================================================================
void BatteryResidEval::doPolarizationAnalysis(const int ifunc, const double t, const double deltaT, const Epetra_Vector_Ghosted &soln,
	                                      const Epetra_Vector_Ghosted * const solnDot_ptr)
{
    // Global index for the electrode voltage at the anode - anode collector interface.
    static size_t gindex_VoltageSolid_ACA = npos;
    static size_t gindex_Voltage_AS = npos;
    static size_t gindex_Voltage_SC = npos;
    static size_t gindex_VoltageSolid_CCC = npos;
    /*
     *  In this analysis we will assume the following about the domains:
     *          domain 0 is  anode
     *          domain 1 is separator
     *          domain 2 is cathode
     *  We may also have to analyse the current collectors to get more polarization resistance
     */
    //class globalHeatBalValsBat pVals;
    //DomainLayout &DL = *DL_ptr_;

    //const size_t domAnode = 0;
    //const size_t domSeparator = 1;
    //const size_t domCathode = 2;

    /*
     *  at any time step we know the current.`At any time step we know the current through each electrode object
     *  Use that ratio to create an averaging operation to yield the polarization data.
     *
     *  At any single electrode object we know all of the polarization voltage losses from the separate start
     *  to the current collector. We know the open circuit voltage of the surface where the current came across.
     *  If there are two surfaces over which the current came across we can effectively create two electrode
     *  objects and average further. 
     *
     *  At any time step we can obtain the average extent of reaction, and find the open circuit voltage based
     *  on that number. 
     */

    /*
     *  Identify 6 points in the domain. 
     *          ACA -> anode-collector to anode point.
     *                 We will obtain the voltage of the solid-phase here
     *
     *          AS  -> anode-separator interface. We will obtain the electric potential of the electrolyte at this point.
     *                 We obtain the instantaneous mole fractions here
     *
     *          SC  -> separator-cathode interface.  We will obtain the electric potential of the electrolyte at this point.
     *                 We obtain the instantaneous mole fractions here
     *
     *          CCC -> cathode-collector to cathode point.
     *                 We will obtain the voltage of the solid-phase here
     *
     *          Anode_Voltage -> this is the voltage at the anode terminal , It may be different than the ACA point, due to losses in 
     *                           the anode collector itself.
     *
     *          Cathode_Voltage -> this is the voltage at the cathode terminal, It may be different than the CCC point, due to losses in 
     *                           the anode collector itself.
     *
     *          S  -> Middle of the separator with average electrolyte conditions. v_S = (AS + SC) / 2
     */

    // Calculation of the OCV

    /*
     *  Calculate an average value of the electrolyte concentration homogenized over all parts of the battery.
     *   (or assume a value based on input conditions).
     *   -> Call this condition "S".
     */

    /*
     *  Calculation of the average value of the anode concentrations
     *     -> Take an electrode object that is representative of the entire anode and calculate a homogeneous half-cell OCV
     *        We can do this based on having an eletrode object and having a definition of "S".
     */
     
    /*
     *   Loop over the Anode Electrodes.
     *        All loops over the anode electrodes will sum up to the current i through the battery.
     *        For each Electrode:
     *             We have V_EE_solid and V_EE_lyte.
     *             Identify deltaV_solid = V_EE_solid - V_ACA
     *                      deltaV_lyte_aa = V_AS - V_EE_lyte.
     *             Fill out polarization losses from these terms.
     * 
     *
     *   Loop over the Cathode Electrodes.
     *        All loops over the anode electrodes will sum up to the current i through the battery.
     *        For each Electrode:
     *             We have V_EE_solid and V_EE_lyte.
     *             Identify deltaV_solid = V_CCC - V_EE_solid
     *                      deltaV_lyte_aa = V_EE_lyte - V_CS.
     *             Fill out polarization losses from these terms.
     *
     *   Calculate polarization losses through the separator.
     *
     *   Do an average over all anodes to get an average result -> probabably separate by plateau!
     *
     *   Do an average over all cathodes to get an average result -> probably separatete by plateau for special processing
     *
     *   Create output results as a 0-D global result.
     *
     */

     // Identify the ACA point

     // Assume it is located at the left of the first domain
     if (polIndecisesCurrent == false) {
         polIndecisesCurrent = true;
         DomainLayout &DL = *DL_ptr_;
         BulkDomain1D *d_anode = DL.BulkDomain1D_List[0];
         //porousFlow_dom1D* p_anode = dynamic_cast<porousFlow_dom1D*>(d_anode);

         // Get the pointer to the Bulk domain descriptor
         m1d::BulkDomainDescription* bdd_anode_ptr = d_anode->BDD_ptr_;

         // Pointer to the LocalNodeIndices object for the domain
         LocalNodeIndices* lni = d_anode->LI_ptr_;

         // Pointer to the global 
         GlobalIndices* gi_ptr = lni->GI_ptr_;

         // Get the index of the  first global node in the anode domain
         int fgn = bdd_anode_ptr->FirstGbNode;

         // get the NodalVars structure pertaining to that global node
         NodalVars* node = gi_ptr->NodalVars_GbNode[fgn];

         // Loop up the stating equation index for that global node
         int index_EqnStart = node->EqnStart_GbEqnIndex;

         // Look up the offset for the voltage unknown
         size_t vs = node->Offset_VarType[Voltage];
         if (vs == npos) {
             throw m1d_Error("", "error");
         }

         // Check to see that there are two voltages at that node and select the second one for the solid voltage
         size_t num = node->Number_VarType[Voltage];
         if (num != 2) {
             throw m1d_Error("", "error");
         }
         gindex_VoltageSolid_ACA = index_EqnStart + vs + 1;

 
         BulkDomain1D *d_separator = nullptr; 
         m1d::BulkDomainDescription* bdd_separator_ptr = nullptr;
         int index_EqnStart_AS = -1;
         if (Separator_BDI_ != npos) {
            d_separator = DL.BulkDomain1D_List[Separator_BDI_];
            bdd_separator_ptr = d_separator->BDD_ptr_;
         }
         //porousFlow_dom1D* p_separator = dynamic_cast<porousFlow_dom1D*>(d_separator);
         //LocalNodeIndices* lni_separator = d_separator->LI_ptr_;
         int gn_AS = bdd_separator_ptr->FirstGbNode;
         NodalVars* node_AS = gi_ptr->NodalVars_GbNode[gn_AS];
         index_EqnStart_AS = node_AS->EqnStart_GbEqnIndex;
         vs = node_AS->Offset_VarType[Voltage];
         gindex_Voltage_AS = index_EqnStart_AS + vs;

         // calculate SC
         int gn_SC = bdd_separator_ptr->LastGbNode;
         NodalVars* node_SC = gi_ptr->NodalVars_GbNode[gn_SC];
         int index_EqnStart_SC = node_SC->EqnStart_GbEqnIndex;
         vs = node_SC->Offset_VarType[Voltage];
         gindex_Voltage_SC = index_EqnStart_SC + vs;

          BulkDomain1D *d_cathode = DL.BulkDomain1D_List[2];
          m1d::BulkDomainDescription* bdd_cathode_ptr = d_cathode->BDD_ptr_;
          int gn_CCC = bdd_cathode_ptr->LastGbNode;
          NodalVars* node_CCC = gi_ptr->NodalVars_GbNode[gn_CCC];
          int index_EqnStart_CCC = node_CCC->EqnStart_GbEqnIndex;
          vs = node_CCC->Offset_VarType[Voltage];
          gindex_VoltageSolid_CCC = index_EqnStart_CCC + vs + 1;

     }

     //double vSolid_AC =  soln[gindex_VoltageSolid_ACA];
     //double vLyte_AS =  soln[gindex_Voltage_AS];
     //double vLyte_SC =  soln[gindex_Voltage_SC];
     //double vSolid_CCC =  soln[gindex_VoltageSolid_CCC];
   

     //double vCathode = reportCathodeVoltage();
     reportCathodeVoltage();

     //double vAnode = 0.0;

     // Loop over the domains gathering the information

     



     // Loop over the anode, filling in the missing pieces that are part of the anode
     
     
      // for each electrode in the anode



       //       void addSolidPol(double phiCurrentCollector, int region);



}
//==================================================================================================================================
void
BatteryResidEval::doSpeciesAnalysis(const int ifunc, const double t, const double deltaT, const Epetra_Vector_Ghosted &y,
				    const Epetra_Vector_Ghosted* const solnDot_ptr)
{
    int iLip = 1;
    static std::vector<double> elem_Start_Total(30, 0.0);
    //static double elemLi_Start_Dom[3];
    static double elemLi_Solid_Start_Dom[3];
    static std::vector<double> elem_Lyte_Start_Total(10, 0.0);

    double elemLi_Solid_New_Dom[3];
    double elemLi_Solid_Old_Dom[3];
    static int iset = 0;
    class globalHeatBalValsBat dVals;
    DomainLayout &DL = *DL_ptr_;
    dVals.sizeLyte(3); 
    std::vector<double> elem_Lyte_New_Total(10, 0.0);
    std::vector<double> elem_Lyte_Old_Total(10, 0.0);
    std::vector<double> elem_Solid_New_Total(30, 0.0);
    std::vector<double> elem_Solid_Old_Total(30, 0.0);
    

    //
    //   Loop over the Volume Domains
    //
    for (size_t iDom = 0; iDom < (size_t) DL.NumBulkDomains; iDom++) {
        dVals.zero();
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	d_ptr->eval_SpeciesElemBalance(ifunc, t, deltaT, &y, solnDot_ptr, solnOld_ptr_, dVals);
	elem_Lyte_New_Total[0] += dVals.elem_Lyte_New[0];
	elem_Lyte_Old_Total[0] += dVals.elem_Lyte_Old[0];
	elem_Lyte_New_Total[1] += dVals.elem_Lyte_New[1];
	elem_Lyte_Old_Total[1] += dVals.elem_Lyte_Old[1];
	elem_Lyte_New_Total[2] += dVals.elem_Lyte_New[2];
	elem_Lyte_Old_Total[2] += dVals.elem_Lyte_Old[2];
        elem_Solid_New_Total[iLip] += dVals.elem_Solid_New[iLip];
        elem_Solid_Old_Total[iLip] += dVals.elem_Solid_Old[iLip];
 
	elemLi_Solid_New_Dom[iDom] = dVals.elem_Solid_New[iLip];
	elemLi_Solid_Old_Dom[iDom] = dVals.elem_Solid_Old[iLip];

	if (!iset) {
	    //elemLi_Start_Dom[iDom] = dVals.elem_Lyte_New[iLip] + dVals.elem_Solid_New[iLip];
	    elemLi_Solid_Start_Dom[iDom] = dVals.elem_Solid_New[iLip];
	}
      
    }
    //
    //    Loop over the Surface Domains
    //
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
        dVals.zero();
	SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	d_ptr->eval_SpeciesElemBalance(ifunc, t, deltaT, &y, solnDot_ptr, solnOld_ptr_, dVals);

    }
    double elem_li_old_total = elem_Lyte_Old_Total[iLip] + elem_Solid_Old_Total[iLip];
    double elem_li_new_total = elem_Lyte_New_Total[iLip] + elem_Solid_New_Total[iLip];
    if (!iset) {
        iset = 1;
        elem_Start_Total[iLip] =  elem_li_old_total;
	elem_Lyte_Start_Total[0] = elem_Lyte_New_Total[0];
	elem_Lyte_Start_Total[1] = elem_Lyte_New_Total[1];
	elem_Lyte_Start_Total[2] = elem_Lyte_New_Total[2];
    }
 
    drawline(0, 132);
    printf("                       ELEMENT BALANCES \n");
    drawline(0, 132);
    //
    // Solvent Species Evolution with Time
    //
    double delta = elem_Lyte_New_Total[0] - elem_Lyte_Start_Total[0];
    printf("Solvent Lyte:  Curr = % 12.5E   Start =  % 12.5E  Delta = % 12.5E\n",
	   elem_Lyte_New_Total[0], elem_Lyte_Start_Total[0], delta);

    delta = elem_Lyte_New_Total[0] - elem_Lyte_Old_Total[0];
    printf("               Curr = % 12.5E   Old   =  % 12.5E  Delta = % 12.5E\n",
	   elem_Lyte_New_Total[0], elem_Lyte_Old_Total[0], delta);
    printf("\n");
    //
    // PF6- Species Evolution with Time
    //
    delta = elem_Lyte_New_Total[2] - elem_Lyte_Start_Total[2];
    printf("PF6- Lyte:     Curr = % 12.5E   Start =  % 12.5E  Delta = % 12.5E\n",
	   elem_Lyte_New_Total[2], elem_Lyte_Start_Total[2], delta);

    delta = elem_Lyte_New_Total[2] - elem_Lyte_Old_Total[2];
    printf("               Curr = % 12.5E   Old   =  % 12.5E  Delta = % 12.5E\n",
	   elem_Lyte_New_Total[2], elem_Lyte_Old_Total[2], delta);
    printf("\n");
    //
    // Lip Species Evolution with Time
    //
    delta = elem_Lyte_New_Total[1] - elem_Lyte_Start_Total[1];
    printf("Lip Lyte:      Curr = % 12.5E   Start =  % 12.5E  Delta = % 12.5E\n",
	   elem_Lyte_New_Total[1], elem_Lyte_Start_Total[1], delta);

    delta = elem_Lyte_New_Total[1] - elem_Lyte_Old_Total[1];
    printf("               Curr = % 12.5E   Old   =  % 12.5E  Delta = % 12.5E\n",
	   elem_Lyte_New_Total[1], elem_Lyte_Old_Total[1], delta);
    printf("\n\n");
    //
    // Anode Solid Phase Evolution of Lithium
    //
    delta = elemLi_Solid_New_Dom[0] - elemLi_Solid_Start_Dom[0]; 
    printf("Anode Solid:   Curr = % 12.5E   Start =  % 12.5E  Delta = % 12.5E\n",
	   elemLi_Solid_New_Dom[0],  elemLi_Solid_Start_Dom[0], delta);
    delta = elemLi_Solid_New_Dom[0] - elemLi_Solid_Old_Dom[0]; 
    printf("               Curr = % 12.5E   Old   =  % 12.5E  Delta = % 12.5E\n",
	   elemLi_Solid_New_Dom[0], elemLi_Solid_Old_Dom[0], delta);
    //
    // Cathode Solid Phase Evolution of Lithium
    //
    delta = elemLi_Solid_New_Dom[2] - elemLi_Solid_Start_Dom[2]; 
    printf("Cath Solid:    Curr = % 12.5E   Start =  % 12.5E  Delta = % 12.5E\n",
	   elemLi_Solid_New_Dom[2],  elemLi_Solid_Start_Dom[2], delta);
    delta = elemLi_Solid_New_Dom[2] - elemLi_Solid_Old_Dom[2]; 
    printf("               Curr = % 12.5E   Old   =  % 12.5E  Delta = % 12.5E\n",
	   elemLi_Solid_New_Dom[2], elemLi_Solid_Old_Dom[2], delta);
    printf("\n");
    //
    // Solid Phase Balance of Lithium
    //
    double newS =  elemLi_Solid_New_Dom[0] + elemLi_Solid_New_Dom[2];
    double oldS =  elemLi_Solid_Old_Dom[0] + elemLi_Solid_Old_Dom[2];
    double startS= elemLi_Solid_Start_Dom[0] + elemLi_Solid_Start_Dom[2];
    delta = newS - startS;
    printf("Total Solid:   Curr = % 12.5E   Start =  % 12.5E  Delta = % 12.5E\n", newS, startS,  delta);
    delta =  newS - oldS;
    printf("               Curr = % 12.5E   Old   =  % 12.5E  Delta = % 12.5E\n\n", newS, oldS, delta);
    //
    // Total Balance of Lithium
    //
    delta = elem_li_new_total - elem_Start_Total[iLip];
    printf("Total Li:      Curr = % 12.5E   Start =  % 12.5E  Delta = % 12.5E\n", elem_li_new_total, elem_Start_Total[iLip], delta);
    delta = elem_li_new_total - elem_li_old_total;
    printf("               Curr = % 12.5E   Old   =  % 12.5E  Delta = % 12.5E\n\n", elem_li_new_total, elem_li_old_total, delta);
    printf("\n");

    drawline(0, 132);
}
//==================================================================================================================================
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
//==================================================================================================================================
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
//==================================================================================================================================
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
//==================================================================================================================================
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
//==================================================================================================================================
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
//==================================================================================================================================
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
//================================================================================================================================
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
//================================================================================================================================
double
BatteryResidEval::reportCathodeVoltage() const
{
    DomainLayout &DL = *DL_ptr_;
    // we want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    SurDomain_CathodeCollector *c_ptr = dynamic_cast<SurDomain_CathodeCollector *>(d_ptr);
#ifdef DEBUG_MODE
    if (!c_ptr) {
        throw m1d_Error("BatteryResidEval::reportCathodeVoltage()", "SurDomain_CathodeCollector failed");
    }
#endif
    // might have to update the SurDomain.
    double phi = c_ptr->phiCathodeCC_;
    return phi;
} 
//==================================================================================================================================
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

    if (anodeType_ == 0) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[0];
	porousElectrode_dom1D*  d_anode_ptr= dynamic_cast<porousElectrode_dom1D*>(d_ptr);
	capacityAnodePA_   =   d_anode_ptr->capacityPA();
	capacityLeftAnodePA_   =   d_anode_ptr->capacityLeftPA();
	capacityDischargedAnodePA_   =   d_anode_ptr->capacityDischargedPA();
	dodAnodePA_   = d_anode_ptr->depthOfDischargePA();
    }
 
    if (cathodeType_ == 0) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[2];
	porousElectrode_dom1D* d_cathode_ptr = dynamic_cast<porousElectrode_dom1D*>(d_ptr);
	capacityCathodePA_ = d_cathode_ptr->capacityPA();
	capacityLeftCathodePA_ = d_cathode_ptr->capacityLeftPA();
	capacityDischargedCathodePA_ = d_cathode_ptr->capacityDischargedPA();
	dodCathodePA_ = d_cathode_ptr->depthOfDischargePA();
    }

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

//==================================================================================================================================
} // end of m1d namespace
//==================================================================================================================================
