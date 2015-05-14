/**
 *  @file m1d_ProblemResidEval.cpp
 *
 **/
/*
 * $Author: hkmoffa $
 * $Revision: 592 $
 * $Date: 2013-05-13 10:57:58 -0600 (Mon, 13 May 2013) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_ProblemResidEval.h"

#include "m1d_DomainLayout.h"
#include "m1d_GlobalIndices.h"
#include "m1d_globals.h"
#include "m1d_EpetraJac.h"
#include "m1d_EpetraExtras.h"
#include "m1d_NodalVars.h"

#include "m1d_ProblemStatement.h"
#include "m1d_SurfDomainDescription.h"
#include "m1d_BulkDomain1D.h"
#include "m1d_SurDomain1D.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include <fstream>

using namespace std;


namespace m1d
{

//=====================================================================================================================
  /*
   *  Print flag that is used to turn on extra printing 
   */ 
  int ProblemResidEval::s_printFlagEnv = 0; 
//=====================================================================================================================

// Default constructor
/*
 *
 * @param atol   Absolute tolerance for calculation
 */
ProblemResidEval::ProblemResidEval(double atol) :
    m_atol(atol), 
    m_neq(0), DL_ptr_(0), GI_ptr_(0), LI_ptr_(0), m_jac(0), m_baseFileName("solution"), m_StepNumber(0),
    solnOld_ptr_(0), resInternal_ptr_(0),
  m_atolVector(0), 
  m_atolVector_DAEInit(0),
  m_atolDeltaDamping(0),
  m_atolDeltaDamping_DAEInit(0),
  psInput_ptr_(0),
  m_numTimeRegions(1),
  m_currentTimeRegion(0),
  m_writeStartEndFile(0),
  counterResBaseCalcs_(0), counterJacBaseCalcs_(0), 
  counterJacDeltaCalcs_(0), counterResShowSolutionCalcs_(0),
  SolutionBehavior_printLvl_(0),
  Residual_printLvl_(0),
  coordinateSystemType_(Rectilinear_Coordinates),
  energyEquationProbType_(0)
{
  /*
   *  Read an environmental variable to set the static variable s_printFlagEnv
   */
  char *resp_str = getenv("PRE_printFlagEnv");
  if (resp_str) {
    if (resp_str[0] == 'y') {
      s_printFlagEnv = 1; 
    } else {
     try {
       string rr(resp_str) ;
       double ff = fpValueCheck(rr);
       s_printFlagEnv  = ff;
     } catch (Cantera::CanteraError &cE) {
       Cantera::showErrors();
       Cantera::popError();
     }
    }
  }
  if (PSinput_ptr == 0) {
     printf("ProblemResidEval constructor: Error expect the input file to have been read and processed\n");
     exit(-1);
  }
  m_writeStartEndFile = PSinput_ptr->writeStartEndFile_;

  if (PSinput_ptr->Energy_equation_prob_type_ == 3) {
       energyEquationProbType_ = PSinput_ptr->Energy_equation_prob_type_;
  } else if (PSinput_ptr->Energy_equation_prob_type_ != 0) {
       printf("unimplemetned\n");
       exit(-1);
  }

}
//=====================================================================================================================
// Destructor
/*
 *
 */
ProblemResidEval::~ProblemResidEval()
{
  // We are responsible now for deleting the global data about the problem
  safeDelete(DL_ptr_);
  safeDelete(LI_ptr_);
  safeDelete(GI_ptr_);

  //safeDelete(m_jac);

  safeDelete(solnOld_ptr_);
  safeDelete(resInternal_ptr_);
  safeDelete(m_atolVector);
  safeDelete(m_atolVector_DAEInit);
  safeDelete(m_atolDeltaDamping);
  safeDelete(m_atolDeltaDamping_DAEInit);
}
//=====================================================================================================================
//! Default copy constructor
/*!
 *
 * @param r  Object to be copied
 */
ProblemResidEval::ProblemResidEval(const ProblemResidEval &r) :
  m_atol(r.m_atol), m_neq(0), DL_ptr_(0), GI_ptr_(0), LI_ptr_(0), m_jac(0), m_baseFileName("solution"),
  m_StepNumber(0), solnOld_ptr_(0), resInternal_ptr_(0),
  m_atolVector(0),
  m_atolVector_DAEInit(0),
  m_atolDeltaDamping(0),
  psInput_ptr_(0),
  m_numTimeRegions(1),
  m_currentTimeRegion(0),
  counterResBaseCalcs_(0), counterJacBaseCalcs_(0),
  counterJacDeltaCalcs_(0), counterResShowSolutionCalcs_(0),
  coordinateSystemType_(r.coordinateSystemType_),
  energyEquationProbType_(r.energyEquationProbType_)
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
ProblemResidEval &
ProblemResidEval::operator=(const ProblemResidEval &r)
{
  if (this == &r) {
    return *this;
  }

  m_atol = r.m_atol;

  safeDelete(DL_ptr_);
  DL_ptr_ = new DomainLayout(*(r.DL_ptr_));

  safeDelete(GI_ptr_);
  GI_ptr_ = new GlobalIndices(*(r.GI_ptr_));

  safeDelete(LI_ptr_);
  LI_ptr_ = new LocalNodeIndices(*(r.LI_ptr_));

  safeDelete(solnOld_ptr_);
  solnOld_ptr_ = new Epetra_Vector(*(r.solnOld_ptr_));

  safeDelete(resInternal_ptr_);
  resInternal_ptr_ = new Epetra_Vector(*(r.resInternal_ptr_));

  safeDelete(m_atolVector);
  m_atolVector = new Epetra_Vector(*(r.m_atolVector));

  safeDelete(m_atolVector_DAEInit);
  m_atolVector_DAEInit = new Epetra_Vector(*(r.m_atolVector_DAEInit));

  safeDelete(m_atolDeltaDamping);
  m_atolDeltaDamping = new Epetra_Vector(*(r.m_atolDeltaDamping));

  safeDelete(m_atolDeltaDamping_DAEInit);
  m_atolDeltaDamping_DAEInit = new Epetra_Vector(*(r.m_atolDeltaDamping_DAEInit));

  safeDelete(psInput_ptr_);
  psInput_ptr_  = new ProblemStatement(*(r.psInput_ptr_));


  counterResBaseCalcs_ = r.counterResBaseCalcs_;
  counterJacBaseCalcs_ = r.counterJacBaseCalcs_;
  counterJacDeltaCalcs_ = r.counterJacDeltaCalcs_;
  counterResShowSolutionCalcs_ = r.counterResShowSolutionCalcs_;

  SolutionBehavior_printLvl_ = r.SolutionBehavior_printLvl_;
  Residual_printLvl_ = r.Residual_printLvl_;
  coordinateSystemType_ = r.coordinateSystemType_;

  return *this;
}
//=====================================================================================================================
/*
 *
 *  Initialize the domain structure for the problem
 *
 *   This routine will grow as we add more types.
 */
void
ProblemResidEval::specifyProblem(int problemType, ProblemStatement *ps_ptr)
{
  psInput_ptr_ = ps_ptr;
  if (problemType == 1) {
    DL_ptr_ = new SimpleDiffusionLayout(1, ps_ptr);
  } else if (problemType == 2) {
    DL_ptr_ = new SimpleTimeDependentDiffusionLayout(1, ps_ptr);
  } else {
    throw m1d_Error("specifyProblem", "Unknown problem type");
  }
  DL_ptr_->setProblemResid(this);
  //
  // Specify the type of the coordinate system
  //
  coordinateSystemType_ = psInput_ptr_->coordinateSystemType_;

  SolutionBehavior_printLvl_ = psInput_ptr_->SolutionBehavior_printLvl_;
}
//=====================================================================================================================
/*
 *
 *  Initialize the domain structure for the problem
 *
 *   This routine will grow as we add more types.
 */
void
ProblemResidEval::specifyProblem(DomainLayout *dl, ProblemStatement *ps_ptr)
{
    psInput_ptr_ = ps_ptr;
    DL_ptr_ = dl;
    DL_ptr_->setProblemResid(this);

    //
    // Specify the type of the coordinate system
    //
    coordinateSystemType_ = psInput_ptr_->coordinateSystemType_;
    SolutionBehavior_printLvl_ = psInput_ptr_->SolutionBehavior_printLvl_;
}
//=====================================================================================================================
/*
 *
 *  Initialize the domain structure for the problem
 *   Initialize the global indices
 *   This routine will grow as we add more types.
 */
void
ProblemResidEval::generateGlobalIndices()
{

  GI_ptr_ = new m1d::GlobalIndices(Comm_ptr);

  m1d::GlobalIndices &GI = *GI_ptr_;

  // (Now that MPI is initialized, we can access/use MPI_COMM_WORLD.)
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &(GI.NumProc));
  MPI_Comm_rank(MPI_COMM_WORLD, &(GI.MyProcID));
#else
  GI.NumProc = 1;
  GI.MyProcID = 0;
#endif

  /*
   *   Initialize the number of global nodes
   * Size all arrays based on the total number of global nodes
   */
  GI.init(DL_ptr_);

  /*
   *  Initialize the mesh positions
   */
  GI.InitMesh();

  /*
   *  Initialize the global NodalVars object
   *     We figure out what domains are located at what nodes.
   *    We determine the number of equations that are located at each
   *    node. This is done on a global basis.
   *     we figure out what equations are located at what nodes. Then, we
   *     fill up the nodal data structure with this information
   */
  GI.discoverNumEqnsPerNode();

  /*
   * Figure out the division of the nodes on the processor
   *   we determine what nodes and equations are on each processor
   */
  GI.procDivide();
}
//=====================================================================================================================
// Generate and fill up the local node vectors on this processor
/*
 *  These indices have to do with accessing the local nodes on the processor
 *  At the end of this routine, the object LocalNodeIncides is fully formed
 *   on this processor.
 *  Within this object, all numbers are in local Row Node Format.
 *  Local Row node format is
 *  defined as the following. All owned nodes come first. They are
 *  ordered in terms of increasing global node number.
 *  Then the right ghost node is listed.
 *  Then the left ghost node is listed
 *  Then, the "globally-all-connected node is listed, if available.
 */
void
ProblemResidEval::generateLocalIndices()
{
  m1d::GlobalIndices &GI = *GI_ptr_;
  /*
   *  Once we have the basic vectors, we can generate the node maps
   *  necessary to formulate the problem using Epetra.
   *  Form a Map to describe the distribution of our owned equations
   *  and our owned nodes on the processors.
   */
  GI.initNodeMaps();
  /*
   *  Create the object to hold the local node data
   */
  LI_ptr_ = new m1d::LocalNodeIndices(Comm_ptr, GI_ptr_);
  /*
   *  Determine the external nodes on each of the processors. With this step,
   *  we have a determination of all of the nodes that each processor needs.
   *  This is also an essential part of determining the matrix stencil for
   *  eventual solution of the problem.
   */
  LI_ptr_->determineLcNodeMaps(DL_ptr_);
  /*
   * Generate the nodal vars data at the local node level by coping the pointers
   */
  LI_ptr_->GenerateNodalVars();
  /*
   *  Update the equation count and the vector of number of equations per node
   *  and the value of the  local equation start index at each local local node.
   */
  LI_ptr_->UpdateEqnCount();
  LI_ptr_->generateEqnMapping();
  /*
   * Once we have found the external nodes, then we can generate the block
   * maps that are part of establishing the matrix
   */
  GI.initBlockNodeMaps(DATA_PTR(LI_ptr_->NumEqns_LcNode));

  LI_ptr_->determineLcEqnMaps();

  /*
   * Determine total number of unknowns in the global equation system
   */
  m_neq = GI.NumGbEqns;

  /*
   * determine number of unknowns on each processor
   */

  //Number of equations on the processor including ghost equations
  m_NumLcEqns = LI_ptr_->NumLcEqns;

  // Number of equations on the processor not including ghost equations
  m_NumLcOwnedEqns = LI_ptr_->NumLcOwnedEqns;

}
//=====================================================================================================================
  void ProblemResidEval::fillIsAlgebraic(Epetra_IntVector & isAlgebraic)
  { 
    DomainLayout &DL = *DL_ptr_;
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
      BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
      d_ptr->fillIsAlgebraic(isAlgebraic);
    }

    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
      SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
      d_ptr->fillIsAlgebraic(isAlgebraic);
    }
  }
  //=====================================================================================================================
  void ProblemResidEval::fillIsArithmeticScaled(Epetra_IntVector & isArithmeticScaled)
  { 
    DomainLayout &DL = *DL_ptr_;
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
      BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
      d_ptr->fillIsArithmeticScaled(isArithmeticScaled);
    }

    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
      SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
      d_ptr->fillIsArithmeticScaled(isArithmeticScaled);
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
ProblemResidEval::residEval(Epetra_Vector_Owned* const & res,	   
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
  double delta_t = 0.0;
  double t_old = t;
  if (rdelta_t > 0.0) {
    delta_t = 1.0/rdelta_t;
    t_old = t - delta_t;
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
  setStateFromSolution(doTimeDependentResid, soln_ptr, solnDot_ptr, t, delta_t, t_old);

  /*
   *   Loop over the Volume Domains
   */
  for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
    BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
    d_ptr->incrementCounters(residType);
    d_ptr->residEval_PreCalc(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t, residType);
  }
  /*
   *    Loop over the Surface Domains
   */
  for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
    SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
    d_ptr->incrementCounters(residType);
    d_ptr->residEval_PreCalc(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr_, t, rdelta_t, residType);
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

}
//=====================================================================================================================
void
ProblemResidEval::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector *soln_ptr,
				      const Epetra_Vector *solnDot_ptr, const Epetra_Vector *solnOld_ptr,
				      const double t, const double t_old)
{
  // Get a local copy of the domain layout
  DomainLayout &DL = *DL_ptr_;

  double rdelta_t = 1.0E200;
  double delta_t = 0.0;
  if (t > t_old) {
    delta_t = t - t_old;
    rdelta_t = 1.0 / delta_t;
  }

  /*
   * We calculate solnOld_ptr_ here
   */
  if (doTimeDependentResid) {
    calcSolnOld(*soln_ptr, *solnDot_ptr, rdelta_t);
  }
  /*
   *   Propagate the solution of the system down to the underlying objects where necessary.
   *   It is necessary to do this for mesh unknowns.
   *   Also do a loop over nodes, calculating quantities.
   */
  setStateFromSolution(doTimeDependentResid, soln_ptr, solnDot_ptr, t, delta_t, t_old);
  /*
   *   Loop over the Volume Domains
   */
  for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
    //BulkDomainDescription *bdd_ptr = DomainDesc_global[iDom];
    BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
    d_ptr->advanceTimeBaseline(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr, t, t_old);
  }
  /*
   *    Loop over the Surface Domains
   */
  for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
    SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
    d_ptr->advanceTimeBaseline(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr, t, t_old);
  }
}
//=====================================================================================================================
// Revert the Residual object's conditions to the conditions at the start of the global time step
/*
 *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init_init.
 *  We get rid of the pendingIntegratedFlags_ flag here as well.
 */
  void ProblemResidEval::revertToInitialGlobalTime()
{ 
    // Get a local copy of the domain layout
    DomainLayout &DL = *DL_ptr_;
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	d_ptr->revertToInitialGlobalTime();
    }
    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
	SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	d_ptr->revertToInitialGlobalTime();
    }
}
//=====================================================================================================================
// Set the underlying state of the system from the solution vector
/*
 *   Note this is an important routine for the speed of the solution.
 *   It would be great if we could supply just exactly what is changing here.
 *   This routine is always called at the beginning of the residual evaluation process.
 *   It is also called during the advanceTimeBaseline() process. This is appropriate as we make sure that
 *   the calculation of nodal quantities is done using the end solution values, and therefore the _old values
 *   will be correct,  before advancing to the next time step. 
 *
 *   This is a natural place to put any precalculations of nodal quantities that
 *   may be needed by the residual before its calculation.
 *
 * @param doTimeDependentResid
 * @param soln
 * @param solnDot
 * @param t
 * @param delta_t
 * @param t_old
 */
void
ProblemResidEval::setStateFromSolution(const bool doTimeDependentResid,
                                       const Epetra_Vector *soln,
                                       const Epetra_Vector *solnDot,
                                       const double t,
                                       const double rdelta_t, const double t_old)
{
  // Get a local copy of the domain layout
  // DomainLayout &DL = *DL_ptr_;

  LI_ptr_->ExtractPositionsFromSolution(soln);

  GI_ptr_->updateGlobalPositions(LI_ptr_->Xpos_LcOwnedNode_p);

  /*
   *  Update the nodal values in the LocalNodeIndices structure
   */
  //  void
  // LocalNodeIndices::UpdateNodalVarsPositions()
}
//=====================================================================================================================
// Calculate the initial conditions
/*
 *   The basic algorithm is to loop over the volume domains.
 *   Then, we loop over the surface domains
 *
 * @param doTimeDependentResid    INPUT Boolean indicating whether we should formulate the time dependent problem
 * @param soln                    OUTPUT Solution vector.  On return this contains the initial solution vector
 * @param solnDot                 OUTPUT Solution Dot Vector. On return this contains the initial time derivative
 *                                      of the solution vector.
 * @param t                       INPUT Time
 * @param delta_t                 INPUT/OUTPUT delta_t for the initial time step
 * @param delta_t                 OUTPUT       delta_t_np1 for the initial time step
 */
void
ProblemResidEval::initialConditions(const bool doTimeDependentResid, Epetra_Vector_Ghosted *soln,
                                    Epetra_Vector_Ghosted *solnDot, double& t,
                                    double& delta_t, double& delta_t_np1)
{
  if (!solnOld_ptr_) {
    solnOld_ptr_ = new Epetra_Vector(*soln);
  }

  // Get a local copy of the domain layout
  DomainLayout &DL = *DL_ptr_;
  /*
   *   Zero the solution vector as a start
   */
  soln->PutScalar(0.0);
  if (doTimeDependentResid) {
    solnDot->PutScalar(0.0);
  }

  /*
   * We set initial conditions here that make sense to do by looping over nodes
   * instead of cells.
   */
  LI_ptr_->setInitialConditions(doTimeDependentResid, soln, solnDot, t, delta_t);
  /*
   *   Loop over the Volume Domains
   */
  for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
    //BulkDomainDescription *bdd_ptr = DomainDesc_global[iDom];
    BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
    d_ptr->initialConditions(doTimeDependentResid, soln, solnDot, t, delta_t);

    //JCH adding tecplot output call here.
    //probably want to change this later
    //d_ptr->writeSolutionTecplotHeader();
  }
  /*
   *    Loop over ths Surface Domains
   */
  for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
    SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
    d_ptr->initialConditions(doTimeDependentResid, soln, solnDot, t, delta_t);

    //JCH adding tecplot output call here.
    //probably want to change this later
    //d_ptr->writeSolutionTecplotHeader();
  }

  // Find the restart record number. If this is greater than zero, then the user has
  // signaled that a restart file should be used
  //
  int rrn = psInput_ptr_->restartRecordNumber_; 
  double t_read;
  double delta_t_read = delta_t;
  double delta_t_next_read = delta_t;
  delta_t_np1 = delta_t;
  if (rrn >= 0) {
      string& rfn = psInput_ptr_->restartFileName_;
      size_t locdot = rfn.find('.');
      string baseFileName = rfn.substr(0, locdot);

      readSolution(rrn, baseFileName, *soln, solnDot, t_read, delta_t_read, delta_t_next_read);
      t = t_read;
      delta_t = delta_t_read;
      delta_t_np1 = delta_t_next_read;
      double t_nm1 = t - delta_t;

      setStateFromSolution(doTimeDependentResid, soln, solnDot, t, delta_t, t_nm1);

  } 

  /*
   *  Create an initial value of the atol vector used in the convergence of the time stepping
   *  algorithm and in the nonlinear solver.
   */
  setAtolVector(psInput_ptr_->absTol_, *soln);
  /*
   *  Create an initial vector of delta damping for the nonlinear solver
   */
  setAtolDeltaDamping(1.0, *soln);
}
//=====================================================================================================================
void
ProblemResidEval::createMatrix(RecordTree_base *linearSolver_db)
{
  m_jac = new EpetraJac(*this);
  //m_jac->solverType_=Iterative;
  m_jac->process_BEinput(linearSolver_db);
  
  m_jac->allocateMatrix();
#ifdef DEBUG_MATRIX_STRUCTURE
  m1d::stream0 w0;
  const Epetra_Comm *cc = m_jac->Comm_ptr_;
  print0_sync_start(true, w0, *cc);
  ostringstream ssSave;
  m_jac->queryMatrixStructure(ssSave);
  w0 << ssSave.str();
  w0 << endl;
  ssprint0(w0);
  print0_sync_end(true, w0, *cc);
#endif
  if (!resInternal_ptr_) {
    resInternal_ptr_ = new Epetra_Vector((m_jac->A_)->RangeMap());
  }
}
//=====================================================================================================================
void
ProblemResidEval::filterSolnPrediction(double t, Epetra_Vector_Ghosted & y)
{

}
//=====================================================================================================================
double
ProblemResidEval::delta_t_constraint(const double time_n,
                                     const Epetra_Vector_Ghosted &y_n,
                                     const Epetra_Vector_Ghosted &ydot_n)
{
  return 0.0;
}
//=====================================================================================================================
void
ProblemResidEval::applyFilter(const double timeCurrent,
                              const double delta_t_n,
                              const Epetra_Vector_Ghosted &y_current,
                              const Epetra_Vector_Ghosted &ydot_current,
                              Epetra_Vector_Ghosted &delta_y)
{
  delta_y.PutScalar(0.0);
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
ProblemResidEval::user_out(const int ievent,
                           const double time_current,
                           const double delta_t_n,
                           const int istep,
                           const Epetra_Vector_Ghosted &y_n,
                           const Epetra_Vector_Ghosted * const ydot_n_ptr)
{
  switch (ievent)
    {
    case -1:
      break;
    case -2:
      break;
    case 0:
      break;
    case 1:
      break;
    case 2:
      break;
    case 3:
      break;
    default:
      throw m1d_Error("ProblemResidEval::user_out","print command not defined for ievent of this type");
    }
}
//=====================================================================================================================
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
 * @param t       Current time
 * @param deltaT  Current value of deltaT
 * @param y       Current value of the solution vectors
 * @param ydot    Current value of time derivative of the solution vectors.
 */
void
ProblemResidEval::evalTimeTrackingEqns(const int ifunc,
                                       const double t,
                                       const double deltaT,
                                       const Epetra_Vector_Ghosted & y,
                                       const Epetra_Vector_Ghosted * const ydot)
{
}
  //====================================================================================================================
  // Set a solution parameter 
  /*
   *  @param paramName   String identifying the parameter to be set
   *  @param paramVal    Single double value of the parameter to be set
   *
   *  @return returns a 0 if the parameter makes sense
   *          Returns a negative number if the parameter is unknown
   */
  int
  ProblemResidEval::setSolutionParam(std::string paramName, double paramVal) 
  {
    return -1;
  }
  //====================================================================================================================
  // Get a solution parameter 
  /*
   *  @param paramName   String identifying the parameter to be set
   *  @param paramVal    Vector of parameters returned
   *
   *  @return returns the number of parameters returned.
   */
  int
  ProblemResidEval::getSolutionParam(std::string paramName, double * const paramVal) 
  {
    return -1;
  }
//=====================================================================================================================
/*
 * Return a vector of delta y's for calculation of the
 * numerical Jacobian
 */
void
ProblemResidEval::calcDeltaSolnVariables(const double t, const double * const ySoln,
                                         const double * const ySolnDot, double * const deltaYSoln,
                                         const double * const solnWeights)
{
  int nEq = m_NumLcEqns;
  if (!solnWeights) {
    for (int i = 0; i < nEq; i++) {
      deltaYSoln[i] = m_atol + fabs(1.0E-6 * ySoln[i]);
    }
  } else {
    for (int i = 0; i < nEq; i++) {
      deltaYSoln[i] = m_atol + fmax(1.0E-2 * solnWeights[i], 1.0E-6 * fabs(ySoln[i]));
    }
  }
}
//=====================================================================================================================
// Save the solution to the end of an XML file using
// XML solution format
/*
 *  We write out the solution to a file. 
 *
 *  
 *
  <ctml>
   <simulation id="0">
    <string title="timestamp">
      Mon Aug  8 10:43:58 2011
    </string>
    <time type="time" units="s"> 0.000000000E+00 </time>
    <integer title="StepNumber" type="time"> 0 </integer>
    <domain id="SurDomain1D_0" numVariables="1" points="1" type="surface">
      <X0> 0.000000000E+00 </X0>
      <X> 0.000000000E+00 </X>
      <Conc_sp(Species0)> 1.000000000E+00 </Conc_sp(Species0)>
    </domain>
    <domain id="BulkDomain1D_0" numVariables="1" points="10" type="bulk">
      <grid_data>
        <floatArray size="10" title="X0" type="length" units="m">
          0.000000000E+00,   1.111111111E-01,   2.222222222E-01,
          3.333333333E-01,   4.444444444E-01,   5.555555556E-01,
          6.666666667E-01,   7.777777778E-01,   8.888888889E-01,
          1.000000000E+00
        </floatArray>
        <floatArray size="10" title="Conc_sp(Species0)" type="concentration" units="kmol/m3">
          1.000000000E+00,   0.000000000E+00,   0.000000000E+00,
          0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
          0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
          0.000000000E+00
        </floatArray>
      </grid_data>
    </domain>
    <domain id="SurDomain1D_1" numVariables="1" points="1" type="surface">
      <X0> 1.000000000E+00 </X0>
      <X> 1.000000000E+00 </X>
      <Conc_sp(Species0)> 0.000000000E+00 </Conc_sp(Species0)>
    </domain>
  </simulation>
 </ctml>
 *
 * NOTES ABOUT THE FORMAT
 * -----------------------------------------
 *   RecordNumber = position of the simulation within the xml file.
 *   SimulationID = String name of the simulation xml node
 *
 *   Initially the RecordNumber and SimulationID are the same. the string name is a
 *   defined as a string, while the RecordNumber is defined as an integer starting
 *   with zero.   If the file is  manipulated this relation may not hold. There are
 *   two ways to get a hold of the record.
 *
 *   The RecordNumber is kept as a static int in this file, as the name solNum. This
 *   may change in the future. The basename of the file is also kept as the static string savedBase.
 *
 *   If the file is called with a different name, then a new RecordNumber saved in the
 *   static var solNumAlt is used to keep track of the record number, with the static
 *   string savedAltBase as its basename.
 *
 * @param ievent  Type of the event. The following form is used:
 *             0 Initial conditions
 *             1 Completion of a successful intermediate step.
 *             2 Final successful conditions.
 *             3 Intermediate nonlinear step
 *            -1 unsuccessful step
 *
 * @param baseFileName to be used. .xml is appended onto the filename
 *          processors other than 0 have the pid appended to the name
 * @param m_y_n    Current value of the solution vector
 * @param m_ydot_n  Current value of the derivative of the solution vector
 */
void
ProblemResidEval::saveSolutionEnd(const int itype,
                                  std::string baseFileName,
                                  const Epetra_Vector_Ghosted &y_n_ghosted,
                                  const Epetra_Vector_Ghosted *ydot_n_ghosted,
                                  const double t,
                                  const double delta_t,
                                  const double delta_t_np1)
{
    //bool doAllProcs = false;
    struct tm *newtime;
    time_t aclock;
    static int solNum = 0;
    static int solNumAlt = 0;
    static string savedBase = baseFileName;
    static string savedAltBase;
    int appendAtEnd = 1;
    ::time(&aclock); /* Get time in seconds */
    newtime = localtime(&aclock); /* Convert time to struct tm form */
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    if (mypid != 0) {
	baseFileName += ("_" + Cantera::int2str(mypid));
    }
    std::string fname = baseFileName + ".xml";
    if (solNum == 0) {
	savedBase = fname;
	appendAtEnd = 0;
    }
    if (savedBase != fname) {
	if ((fname != savedAltBase) && (savedAltBase == "")) {
	    solNumAlt = 0;
	    savedAltBase = fname;
	    appendAtEnd = 0;
	}
    }
    if ((fname != savedBase) && (fname != savedAltBase)) {
	appendAtEnd = 0;
    }

    Epetra_BlockMap *nmap = ownedMap();
    Epetra_Vector *y_n_owned_ptr =  new_EpetraVectorView(y_n_ghosted, *nmap);
    Epetra_Vector *ydot_n_owned_ptr = 0; 
    if (ydot_n_ghosted) {
	ydot_n_owned_ptr =  new_EpetraVectorView(*ydot_n_ghosted, *nmap);
    }

    Cantera::XML_Node root("--");

    Cantera::XML_Node& ct = root.addChild("ctml");

    Cantera::XML_Node& sim = ct.addChild("simulation");
    //
    // Initially define the 
    std::string simulationID = "0";
    if (fname == savedBase) {
	simulationID = Cantera::int2str(solNum);
    } else if (fname == savedAltBase) {
	simulationID = Cantera::int2str(solNumAlt);
    }
    sim.addAttribute("id", simulationID);
    ctml::addString(sim, "timestamp", asctime(newtime));
    // if (desc != "")
    //  addString(sim, "description", desc);

    ctml::addFloat(sim, "time", t, "s", "time");
    if (delta_t > 0.0) {
	ctml::addFloat(sim, "delta_t", delta_t, "s", "time");
    } else {
	ctml::addFloat(sim, "delta_t", 0.0, "s", "time");
    }
    ctml::addFloat(sim, "delta_t_np1", delta_t_np1, "s", "time");
    ctml::addInteger(sim, "StepNumber", m_StepNumber, "", "time");

    // Get a local copy of the domain layout
    DomainLayout &DL = *DL_ptr_;

    Epetra_Vector *solnAll = GI_ptr_->SolnAll;
    m1d::gather_nodeV_OnAll(*solnAll, *y_n_owned_ptr);
    //solnAll = m1d::gatherOnAll(y_n);
    Epetra_Vector *soln_dot_All = GI_ptr_->SolnDotAll;
    if (ydot_n_ghosted) {
	m1d::gather_nodeV_OnAll(*soln_dot_All, *ydot_n_owned_ptr);
    }

    Domain1D *d_ptr = DL.SurDomain1D_List[0];
    do {
	d_ptr->saveDomain(sim, solnAll, soln_dot_All, t, false);

	//JCH adding tecplot output call here.
	//probably want to change this later
	//d_ptr->writeSolutionTecplot(solnAll, soln_dot_All, t);

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

    if (mypid == 0) {
	if (!appendAtEnd) {
	    fstream s(fname.c_str(), fstream::in | fstream::out | ios::trunc);
	    root.write(s);
	    s.close();
	} else {

	    int length = 0;
	    ifstream sr(fname.c_str());
	    if (sr) {
		sr.seekg(0, ios::end);
		length = sr.tellg();
		sr.close();
	    }

	    if (length > 0) {
		fstream s(fname.c_str(), fstream::in | fstream::out | fstream::ate);
		//int pos = s.tellp();
		s.seekp(-8, ios_base::end);
		sim.write(s, 2);
		s.write("</ctml>\n", 8);
		s.close();
	    } else {
		fstream s(fname.c_str(), fstream::in | fstream::out | fstream::app);
		root.write(s);
		s.close();
	    }
	}
    }

    if (mypid == 0) {
	Cantera::writelog("Solution saved to file " + fname + " as solution " + simulationID + ".\n");
    }
    Epetra_Comm *c = LI_ptr_->Comm_ptr_;
    c->Barrier();
    if (fname == savedBase) {
	solNum++;
    } else if (fname == savedAltBase) {
	solNumAlt++;
    }
    delete y_n_owned_ptr;
    delete ydot_n_owned_ptr;
}
//=====================================================================================================================
void
ProblemResidEval::readSolutionRecordNumber(const int iNumber,
					   std::string baseFileName,
					   Epetra_Vector_Ghosted &y_n_ghosted,
					   Epetra_Vector_Ghosted * const ydot_n_ghosted,
					   double &t_read,
					   double &delta_t_read,
					   double &delta_t_next_read)
{
    //
    // Flesh out the name of the file
    //
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    if (mypid != 0) {
	baseFileName += ("_" + Cantera::int2str(mypid));
    }
    std::string fname = baseFileName + ".xml";

    //
    //    Read in the XML file 
    //
    XML_Node* xSavedSoln = get_XML_File(fname);
    if (!xSavedSoln) {
	throw m1d_Error("ProblemResidEval::readSolutionRecordNumber()",
			"Error could not read the file " + fname);
    }
 
    XML_Node* simulRecord = selectSolutionRecordNumber(xSavedSoln, iNumber);
    if (!simulRecord) {
	throw m1d_Error("ProblemResidEval::readSolutionRecordNumber()",
			"Error could not find the requested record " +  int2str(iNumber));	
    }
    //
    //  Now call the underlying routine that reads the record
    //
    readSolutionXML(simulRecord, y_n_ghosted,ydot_n_ghosted, t_read, delta_t_read, delta_t_next_read);

    if (mypid == 0) {
	Cantera::writelog("Read saved solution from  file " + fname + " as solution " + int2str(iNumber) + ".\n");
    }
    Epetra_Comm *c = LI_ptr_->Comm_ptr_;
    c->Barrier();

}
//=====================================================================================================================
void
ProblemResidEval::readSolution(const int iNumber,
			       std::string baseFileName,
			       Epetra_Vector_Ghosted &y_n_ghosted,
			       Epetra_Vector_Ghosted * const ydot_n_ghosted,
			       double &t_read,
			       double &delta_t_read,
                               double &delta_t_next_read)
{
    readSolutionRecordNumber(iNumber, baseFileName, y_n_ghosted, ydot_n_ghosted,
			     t_read, delta_t_read, delta_t_next_read);
}
//=====================================================================================================================
// Read the solution from a saved file.
/*
 *  We read only successful final steps
 *
 * @param iNumber      Solution number to read from. Solutions are numbered consequetively from 
 *                     zero to the number of global time steps
 *
 * @param baseFileName to be used. .xml is appended onto the filename
 *          processors other than 0 have the pid appended to the name
 * @param m_y_n        Current value of the solution vector
 * @param m_ydot_n     Current value of the derivative of the solution vector
 * @param t_read       time that is read out
 * @param delta_t_read delta time step for the last time step.
 * @param delta_t_next_read delta time step for the next time step if available
 */
void
ProblemResidEval::readSolutionXML(XML_Node* simulRecord, Epetra_Vector_Ghosted &y_n_ghosted,
				  Epetra_Vector_Ghosted * const ydot_n_ghosted, double &t_read,
				  double &delta_t_read, double &delta_t_next_read)
{    
    // Get a local copy of the domain layout
    DomainLayout &DL = *DL_ptr_;

    //
    //  Go get a vector for the global solution vector -> this has already been malloced  (i hope)
    //
    // Epetra_Vector *solnAll = GI_ptr_->SolnAll;
  
    // 
    //  Go get a vector for the global solution vector dot -> this has already been malloced  (i hope)
    //
    // Epetra_Vector *soln_dot_All = GI_ptr_->SolnDotAll;


    t_read = ctml::getFloat(*simulRecord, "time");
    delta_t_read = ctml::getFloat(*simulRecord, "delta_t");
    delta_t_next_read = ctml::getFloat(*simulRecord, "delta_t_np1");

    //
    //   Loop over the domains reading in the current solution vector
    //
    Domain1D *d_ptr = DL.SurDomain1D_List[0];
    do {
	d_ptr->readDomain(*simulRecord, &y_n_ghosted, ydot_n_ghosted);


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
}
//====================================================================================================================
//  Select the global time step increment record by the consequuatively numbered record index number
/*
 *    @param   xSoln               Solution file for the simulation object
 *    @param   globalRecordNumbe   Time step record number to select
 *
 */
 XML_Node* ProblemResidEval::selectSolutionRecordNumber(XML_Node* xmlTop, int globalRecordNumber)
 {   
     /*
      *  Search for a particular global step number
      */
     XML_Node* xctml = xmlTop;
     if (xctml->name() != "ctml") {
	 xctml = xmlTop->findByName("ctml");
	 if (!xctml) {
	     throw CanteraError("selectSolutionRecordNumber()","Can't find the top ctml node");
	 }
     }
     //
     //  Get a vector of simulation children, and then pick the record number from that
     //
     std::vector<XML_Node*> ccc;
     xctml->getChildren("simulation", ccc);
     int sz = ccc.size();
     if (globalRecordNumber < 0 || globalRecordNumber >= sz) {
        throw CanteraError("selectSolutionRecordNumber()", 
                           "Can't find the global record number " + int2str(globalRecordNumber));
     }
     XML_Node *xs = ccc[globalRecordNumber];

     return xs;
}
//====================================================================================================================
//  Select the global time step record by its string id.
/*
 *    @param   xSoln               Solution file for the simulation object
 *    @param   timeStepID          Time step number to select
 *
 */
 XML_Node* ProblemResidEval::selectSolutionTimeStepID(XML_Node* xSoln, std::string timeStepID)
 {   
     /*
      *  Search for a particular global step number
      */
     XML_Node* eRecord = xSoln->findNameID("simulation", timeStepID);
     if (!eRecord) {
	 throw m1d_Error("selectSolutionTimeStepID()", "Solution record for following id not found : " + timeStepID);
     }
     /*
      * Return a pointer to the record
      */
     return eRecord;
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
ProblemResidEval::showProblemSolution(const int ievent, 
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

  // residEval(resInternal_ptr_, doTimeDependentResid, &y_n, ydot_n, t, rdelta_t, Base_ShowSolution, solveType);

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
void
ProblemResidEval::writeTecplot(const int ievent,
			       std::string baseFileName,
			       bool doTimeDependentResid, 
			       const double t,
			       const double delta_t,
			       const Epetra_Vector_Ghosted &y_n,
			       const Epetra_Vector_Ghosted * const ydot_n,
			       const Solve_Type_Enum solveType,
			       const double delta_t_np1)
{ 

    // Get a local copy of the domain layout
    DomainLayout &DL = *DL_ptr_;

    static int firstTime = 1;

    Epetra_Vector *solnAll = GI_ptr_->SolnAll;
    m1d::gather_nodeV_OnAll(*solnAll, y_n);
    
    Epetra_Vector *soln_dot_All = GI_ptr_->SolnDotAll;
    if (ydot_n) {
	m1d::gather_nodeV_OnAll(*soln_dot_All, *ydot_n);
    }
    
    //
    // Write individual Tecplot files for each domain
    //
    if (firstTime) {
	firstTime = 0;
	for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	    BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	    d_ptr->writeSolutionTecplotHeader();
	}
	for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
	    SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	    d_ptr->writeSolutionTecplotHeader();
	}
	

	//writeGlobalTecplotHeader();
    }
    
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
	BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
	d_ptr->writeSolutionTecplot(solnAll, soln_dot_All, t);
    }
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
	SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
	d_ptr->writeSolutionTecplot(solnAll, soln_dot_All, t);
    }
  
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
ProblemResidEval::showSolutionVector(std::string& solnVecName,
				     const double t,
				     const double delta_t,
				     const Epetra_Vector_Owned &solnVector,
				     FILE *outFile)
{

  //bool doAllProcs = false;
  time_t aclock;
  ::time(&aclock); /* Get time in seconds */
  int mypid = LI_ptr_->Comm_ptr_->MyPID();

  bool duplicateOnAllProcs = false;
  // Get a local copy of the domain layout
  DomainLayout &DL = *DL_ptr_;

  Epetra_Vector *solnVectorAll = GI_ptr_->SolnAll;
  m1d::gather_nodeV_OnAll(*solnVectorAll, solnVector);


  double rdelta_t = 0.0;
  if (delta_t > 0.0) {
    rdelta_t = 1.0 / delta_t;
  }


  int indentSpaces = 4;
  string indent = "    ";
  const char *ind = indent.c_str();
  if (!mypid || duplicateOnAllProcs) {
    fprintf(outFile, "%s", ind);
    for (int i = 0; i < 100; i++) fprintf(outFile, "-");
    fprintf(outFile, "\n");
    fprintf(outFile, "%s showSolutionVector : %s\n", ind, solnVecName.c_str());
    fprintf(outFile, "%s                       Time       %-12.3E\n", ind, t);
    fprintf(outFile, "%s                       Delta_t    %-12.3E\n", ind, delta_t);
    fprintf(outFile, "%s                       StepNumber %6d\n", ind, m_StepNumber);
    fprintf(outFile, "%s", ind);
    for (int i = 0; i < 100; i++) fprintf(outFile, "-");
    fprintf(outFile, "\n");
  }

  Domain1D *d_ptr = DL.SurDomain1D_List[0];
  do {
    d_ptr->showSolutionVector(solnVecName, solnVectorAll,  &solnVector, t, rdelta_t, indentSpaces + 2, 
			      duplicateOnAllProcs, outFile);
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
    fprintf(outFile, "%s", ind);
    for (int i = 0; i < 100; i++) fprintf(outFile, "-");
    fprintf(outFile, "\n");
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
ProblemResidEval::showSolutionIntVector(std::string& solnVecName,
				     const double t,
				     const double delta_t,
				     const Epetra_IntVector &solnIntVector,
				     FILE *outFile)
{

  //bool doAllProcs = false;
  time_t aclock;
  ::time(&aclock); /* Get time in seconds */
  int mypid = LI_ptr_->Comm_ptr_->MyPID();

  bool duplicateOnAllProcs = false;
  // Get a local copy of the domain layout
  DomainLayout &DL = *DL_ptr_;

  Epetra_IntVector *solnIntVectorAll = GI_ptr_->SolnIntAll;
  m1d::gather_nodeIntV_OnAll(*solnIntVectorAll, solnIntVector);


  double rdelta_t = 0.0;
  if (delta_t > 0.0) {
    rdelta_t = 1.0 / delta_t;
  }


  int indentSpaces = 4;
  string indent = "    ";
  const char *ind = indent.c_str();
  if (!mypid || duplicateOnAllProcs) {
    fprintf(outFile, "%s", ind);
    for (int i = 0; i < 100; i++) fprintf(outFile, "-");
    fprintf(outFile, "\n");
    fprintf(outFile, "%s showSolutionVector : %s\n", ind, solnVecName.c_str());
    fprintf(outFile, "%s                       Time       %-12.3E\n", ind, t);
    fprintf(outFile, "%s                       Delta_t    %-12.3E\n", ind, delta_t);
    fprintf(outFile, "%s                       StepNumber %6d\n", ind, m_StepNumber);
    fprintf(outFile, "%s", ind);
    for (int i = 0; i < 100; i++) fprintf(outFile, "-");
    fprintf(outFile, "\n");
  }

  Domain1D *d_ptr = DL.SurDomain1D_List[0];
  do {
    d_ptr->showSolutionIntVector(solnVecName, solnIntVectorAll,  &solnIntVector, t, rdelta_t, indentSpaces + 2, 
				 duplicateOnAllProcs, outFile);
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
    fprintf(outFile, "%s", ind);
    for (int i = 0; i < 100; i++) fprintf(outFile, "-");
    fprintf(outFile, "\n");
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
ProblemResidEval::writeSolution(const int ievent,
                                const bool doTimeDependentResid,
                                const double time_current,
                                const double delta_t_n,
                                int istep,
                                const Epetra_Vector_Ghosted &y_n,
                                const Epetra_Vector_Ghosted * const ydot_n_ptr,
				const Solve_Type_Enum solveType,
                                const double delta_t_np1)
{
    //
    // HKM -> Could gather solnAll, i.e., the entire solution on proc 0,  here
    //

    //
    //  Do a Base_ShowSolution calculation, so that all items are evaluated no matter what the print level.
    //
    double rdelta_t = 0.0;
    if (delta_t_n > 0.0) {
       rdelta_t = 1.0 / delta_t_n;
    }
    residEval(resInternal_ptr_, doTimeDependentResid, &y_n, ydot_n_ptr, time_current, rdelta_t, Base_ShowSolution, solveType);


    static double lastTime = -10000000.;
    if (ievent == 0 || ievent == 1 || ievent == 2) {
        if (ievent == 2) {
           if (fabs(lastTime - time_current) > 1.0E-3 * delta_t_n) {
	       saveSolutionEnd(ievent, m_baseFileName, y_n, ydot_n_ptr, time_current, delta_t_n, delta_t_np1);
            }
        } else {
	    saveSolutionEnd(ievent, m_baseFileName, y_n, ydot_n_ptr, time_current, delta_t_n, delta_t_np1);
        }
    }
    if (m_writeStartEndFile) {
        if (ievent == 0 || ievent == 2) {
  	    saveSolutionEnd(ievent, "solutionStartEnd", y_n, ydot_n_ptr, time_current, delta_t_n, delta_t_np1);
        }
    }
   
    if (ievent == 0 || ievent == 1 || ievent == 2) {
 	if (PSinput_ptr->SolutionBehavior_printLvl_ > 3) {
	    showProblemSolution(ievent,  doTimeDependentResid, time_current, delta_t_n, y_n, ydot_n_ptr, solveType, 
                                delta_t_np1);
	}
    }

    
    writeTecplot(ievent,  m_baseFileName, doTimeDependentResid, time_current, delta_t_n, y_n, ydot_n_ptr, solveType, 
		 delta_t_np1);
    
    lastTime = time_current;
}
//=====================================================================================================================
/*
 *  evalStoppingCriteria()
 *
 * If there is a stopping critera other than time set it here.
 *
 */
bool
ProblemResidEval::evalStoppingCritera(double &time_current,
                                      double &delta_t_n,
                                      const Epetra_Vector_Ghosted &y_n,
                                      const Epetra_Vector_Ghosted &ydot_n)
{
  return false;
}



//=====================================================================================================================
/*
 * matrixConditioning()
 *
 * Multiply the matrix by the inverse of a matrix which lead to a
 * better conditioned system. The default, specified here, is to
 * do nothing.
 */
void
ProblemResidEval::matrixConditioning(double * const matrix, int nrows, double * const rhs)
{
}
//=====================================================================================================================
//  Prepare all of the pointers for a fast, efficient residual calculation
/*
 *
 */
void
ProblemResidEval::domain_prep()
{

  // Get a local copy of the domain layout
  DomainLayout &DL = *DL_ptr_;

  /*
   *   Loop over the Volume Domains
   *
   */
  for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
    //BulkDomainDescription *bdd_ptr = DomainDesc_global[iDom];
    BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];

    d_ptr->domain_prep(LI_ptr_);
  }

  /*
   *    Loop over the Surface Domains
   */
  for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
    SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];

    d_ptr->domain_prep(LI_ptr_);
  }
}
//=====================================================================================================================
  // Link the ProblemStatement class with the ProblemResidEval class and other classes
  // that need input from the user
  /*
   *  @param psInput_ptr   Pointer to the ProblemStatement class
   */
  void ProblemResidEval::link_input(ProblemStatement *psInput_ptr) {
    psInput_ptr_ = psInput_ptr;
  }
//=====================================================================================================================
EpetraJac&
ProblemResidEval::jacobian()
{
  return (*m_jac);
}
//=====================================================================================================================
// Returns the ghosted map of the problem
Epetra_BlockMap *
ProblemResidEval::ghostedMap() const
{
  return LI_ptr_->GbBlockNodeEqnstoLcBlockNodeEqnsColMap;
}
//=====================================================================================================================
// Returns the unique unknown map of the problem
Epetra_BlockMap *
ProblemResidEval::ownedMap() const
{
  return LI_ptr_->GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap;
}
//=====================================================================================================================
int
ProblemResidEval::NumGbEqns() const
{
  if (!GI_ptr_) {
    throw m1d_Error("ProblemResidEval::NumGbEqns()", "pointer is null");
  }
  return GI_ptr_->NumGbEqns;
}
//=====================================================================================================================
int
ProblemResidEval::NumGbNodes() const
{
  if (!GI_ptr_) {
    throw m1d_Error("ProblemResidEval::NumGbNodes()", "pointer is null");
  }
  return GI_ptr_->NumGbNodes;
}
//=====================================================================================================================
  std::string ProblemResidEval::equationID(int ilocalEqn, int& iLcNode, int & iGbNode,  int& iNodeEqnNum, EqnType &eqn,
					   EQ_TYPE_SUBNUM &etsub)
{
  std::string vstring;
  int ilfound = LI_ptr_->NumLcNodes - 1;
  int istart;
  // Find the local node number
  for (int i = 0; i < LI_ptr_->NumLcNodes; i++) {
    istart = LI_ptr_->IndexLcEqns_LcNode[i];
    if (ilocalEqn < istart) {
      ilfound = i - 1;
      break;
    }
  }
  iLcNode = ilfound;
  // Find the global node number
  iGbNode = LI_ptr_->IndexGbNode_LcNode[ilfound];
  // Get the nodal vars structure for the global node
  NodalVars *nv =  GI_ptr_->NodalVars_GbNode[iGbNode];
  // Find the local equation number on that node
  iNodeEqnNum = ilocalEqn - LI_ptr_->IndexLcEqns_LcNode[ilfound];
  // Now query the NodalVars structure for the variable name and subtype
  eqn   = nv->EquationNameList_EqnNum[iNodeEqnNum];
  etsub = nv->EquationSubType_EqnNum[iNodeEqnNum];
  // Get the string which identifies the variable and return
  vstring = eqn.EquationName(200);
  return vstring;
}
//=====================================================================================================================
std::string ProblemResidEval::variableID(int ilocalVar, int& iLcNode, int & iGbNode, int& iNodeEqnNum, VarType &var,
					 VAR_TYPE_SUBNUM &vtsub)
{
  std::string vstring;
  int ilfound = LI_ptr_->NumLcNodes - 1;
  int istart;
  // Find the local node number
  for (int i = 0; i < LI_ptr_->NumLcNodes; i++) {
    istart = LI_ptr_->IndexLcEqns_LcNode[i];
    if (ilocalVar < istart) {
      ilfound = i - 1;
      break;
    }
  }
  iLcNode = ilfound;
  // Find the global node number
  iGbNode = LI_ptr_->IndexGbNode_LcNode[ilfound];
  // Get the nodal vars structure for the global node
  NodalVars *nv =  GI_ptr_->NodalVars_GbNode[iGbNode];
  // Find the local equation number on that node
  iNodeEqnNum = ilocalVar - LI_ptr_->IndexLcEqns_LcNode[ilfound];
  // Now query the NodalVars structure for the variable name and subtype
  var = nv->VariableNameList_EqnNum[iNodeEqnNum];
  vtsub = nv->VariableSubType_EqnNum[iNodeEqnNum];
  // Get the string which identifies the variable and return
  vstring = var.VariableName(200);
  return vstring;
}
//=====================================================================================================================
// Return the global equation number given the global block number
/*
 *
 * @param gbBlockNum
 * @param localRowNumInBlock
 * @return  the global equation number
 */
int
ProblemResidEval::GbEqnNum_From_GbBlock(const int gbBlockNum, const int localRowNumInBlock)
{
#ifdef DEBUG_MODE
  if (gbBlockNum < 0 || gbBlockNum >= GI_ptr_->NumGbNodes) {
    throw m1d_Error("ProblemResidEval::GbEqnNum_From_GbBlock", "gbBlockNum out of range: "
        + Cantera::int2str(gbBlockNum) + "\n");
  }
  int indexStart = GI_ptr_->IndexStartGbEqns_Proc[gbBlockNum];
  if ((localRowNumInBlock < 0) || (localRowNumInBlock >= GI_ptr_->NumEqns_GbNode[gbBlockNum])) {
    throw m1d_Error("ProblemResidEval::GbEqnNum_From_GbBlock", "localRowNumInBlock out of range: "
        + Cantera::int2str(localRowNumInBlock) + "\n");
  }
  return indexStart + localRowNumInBlock;
#else
  return GI_ptr_->IndexStartGbEqns_Proc[gbBlockNum] + localRowNumInBlock;
#endif
}
//=====================================================================================================================
// global node to local node mapping
/*
 * Given a global node, this function will return the local node value.
 * If the global node is not on this processor, then this function returns -1.
 *
 * This is a pass through routine to Local Node Class which already has implemented this
 *
 * @param gbNode  global node
 * @return returns the local node value
 */
int
ProblemResidEval::GbNodeToLcNode(const int gbNode) const
{
  return (LI_ptr_->GbNodeToLcNode(gbNode));
}
//=====================================================================================================================
int
ProblemResidEval::GbEqnToLcEqn(const int gbEqn) const
{
  int rowEqnNum;
  int gbNode = GI_ptr_->GbEqnToGbNode(gbEqn, rowEqnNum);

  int lcnode = LI_ptr_->GbNodeToLcNode(gbNode);

  if (lcnode >= 0) {
    int lceqn = LI_ptr_->IndexLcEqns_LcNode[lcnode];
    return (rowEqnNum + lceqn);
  }
  return -1;
}
//======================================================================================================================
//  Find a delta of a solution component for use in the numerical jacobian
/*
 *    @param soln      Reference to the complete solution vector
 *    @param ieqn      local equation number of the solution vector
 */
double
ProblemResidEval::deltaSolnComp(const Epetra_Vector & soln, const int ieqn)
{
  /*
   *  We are screwing around here with the algorithm. Later we will put in
   *  a serious algorithm.
   */
  double base = soln[ieqn];
  double delta = 1.0E-6 * fabs(base) + 1.0E-9;
  //double delta = 1.0E-3 * fabs(base) + 1.0E-3;
  return (base + delta);
}
//======================================================================================================================
// Evaluates the atol vector used in the time stepper routine and in the nonlinear solver.
/*
 * @param atolDefault  Double containing the default value of atol
 * @param soln         current solution vector.
 */
void
ProblemResidEval::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
                                const Epetra_Vector_Ghosted * const atolV)
{  
  DomainLayout &DL = *DL_ptr_;
  m_atol =  atolDefault;
  if (!m_atolVector) {
    m_atolVector = new Epetra_Vector(soln);
  }
  if (atolV == 0) {
    for (int i = 0; i < m_NumLcEqns; i++) {
      (*m_atolVector)[i] = atolDefault;
    }
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
      BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
      d_ptr->setAtolVector(atolDefault, soln, *m_atolVector);
    }

    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
      SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
      d_ptr->setAtolVector(atolDefault, soln, *m_atolVector);
    }

  } else {
    for (int i = 0; i < m_NumLcEqns; i++) {
      (*m_atolVector)[i] = (*atolV)[i];
    }
  }


}
//======================================================================================================================
// Evaluates the atol vector used in the time stepper routine and in the nonlinear solver.
/*
 * @param atolDefault  Double containing the default value of atol
 * @param soln         current solution vector.
 */
void
ProblemResidEval::setAtolVector_DAEInit(double atolDAEInitDefault, const Epetra_Vector_Ghosted & soln, const Epetra_Vector_Ghosted & solnDot,
				        const Epetra_Vector_Ghosted * const atolVector_DAEInit) const
{
  DomainLayout &DL = *DL_ptr_;
  // m_atol =  atolDefault;
  if (!m_atolVector_DAEInit) {
    m_atolVector_DAEInit = new Epetra_Vector(soln);
  }
  if (atolVector_DAEInit == 0) {

    for (int i = 0; i < m_NumLcEqns; i++) {
      (*m_atolVector_DAEInit)[i] = atolDAEInitDefault;
    }
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
      BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
      d_ptr->setAtolVector_DAEInit(atolDAEInitDefault, soln, solnDot, *m_atolVector_DAEInit);
    }

    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
      SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
      d_ptr->setAtolVector_DAEInit(atolDAEInitDefault, soln, solnDot, *m_atolVector_DAEInit);
    }

  } else {
    for (int i = 0; i < m_NumLcEqns; i++) {
      (*m_atolVector_DAEInit)[i] = (*atolVector_DAEInit)[i];
    }
  }
}
//======================================================================================================================
void
ProblemResidEval::setAtolDeltaDamping(double relCoeff, const Epetra_Vector_Ghosted & soln, 
				      const Epetra_Vector_Ghosted * const atolDeltaDamping) 
{
  DomainLayout &DL = *DL_ptr_;
  double atolDAEInitDefault = m_atol * 1.0E4;
  if (!m_atolDeltaDamping) {
    m_atolDeltaDamping = new Epetra_Vector(soln);
  }
  if (atolDeltaDamping == 0) {

    for (int i = 0; i < m_NumLcEqns; i++) {
      (*m_atolDeltaDamping)[i] = m_atol * relCoeff;
    }
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
      BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
      d_ptr->setAtolDeltaDamping(atolDAEInitDefault, relCoeff, soln, *m_atolDeltaDamping);
    }

    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
      SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
      d_ptr->setAtolDeltaDamping(atolDAEInitDefault, relCoeff, soln, *m_atolDeltaDamping);
    }

  } else {
    for (int i = 0; i < m_NumLcEqns; i++) {
      (*m_atolDeltaDamping)[i] = (*atolDeltaDamping)[i];
    }
  }
}
//======================================================================================================================
void
ProblemResidEval::setAtolDeltaDamping_DAEInit(double relCoeff, const Epetra_Vector_Ghosted & soln,
					      const Epetra_Vector_Ghosted & solnDot,
					      const Epetra_Vector_Ghosted * const atolDeltaDamping) 
{
  DomainLayout &DL = *DL_ptr_;
  double atolDeltaDampingDefault = m_atol * 1.0E4;
  if (!m_atolDeltaDamping_DAEInit) {
    m_atolDeltaDamping_DAEInit = new Epetra_Vector(soln);
  }
  if (atolDeltaDamping == 0) {

    for (int i = 0; i < m_NumLcEqns; i++) {
      (*m_atolDeltaDamping_DAEInit)[i] = m_atol * relCoeff;
    }
    /*
     *   Loop over the Volume Domains
     */
    for (int iDom = 0; iDom < DL.NumBulkDomains; iDom++) {
      BulkDomain1D *d_ptr = DL.BulkDomain1D_List[iDom];
      d_ptr->setAtolDeltaDamping_DAEInit(atolDeltaDampingDefault , relCoeff, soln, solnDot, *m_atolDeltaDamping_DAEInit);
    }

    /*
     *    Loop over the Surface Domains
     */
    for (int iDom = 0; iDom < DL.NumSurfDomains; iDom++) {
      SurDomain1D *d_ptr = DL.SurDomain1D_List[iDom];
      d_ptr->setAtolDeltaDamping_DAEInit(atolDeltaDampingDefault, relCoeff, soln, solnDot, *m_atolDeltaDamping_DAEInit);
    }

  } else {
    for (int i = 0; i < m_NumLcEqns; i++) {
      (*m_atolDeltaDamping_DAEInit)[i] = (*atolDeltaDamping)[i];
    }
  }
}
//======================================================================================================================
// Return a const reference to the atol vector
const Epetra_Vector_Ghosted &ProblemResidEval::atolVector() const
{
  if (!m_atolVector) {
    throw CanteraError("ProblemResidEval:: atolVector()",
                       "m_atolVector vector hasn't been malloced yet. This is done when an initial guess "
                       "to the solution is supplied");
  }
  return *m_atolVector;
}
//======================================================================================================================
// Return a const reference to the atol DAE init vector
const Epetra_Vector_Owned & ProblemResidEval::atolVector_DAEInit() const
{
  if (!m_atolVector_DAEInit) {
     throw CanteraError("ProblemResidEval:: atolVector_DAEInit()",
                       "m_atolVector_DAEInit vector hasn't been malloced yet. This is done when an initial guess "
                       "to the solution is supplied");
  }
 
  
  return *m_atolVector_DAEInit;
}
//=====================================================================================================================
// Return a const reference to the atol vector
const Epetra_Vector_Ghosted &ProblemResidEval::atolDeltaDamping() const
{
  if (!m_atolDeltaDamping) { 
    throw CanteraError("ProblemResidEval:: atol_deltaDamping()",
                       "m_atolDeltaDamping vector hasn't been malloced yet");
  }
  return *m_atolDeltaDamping;
}
//=====================================================================================================================
// Update the ghost unknowns on the processor.
/*
 *
 * @param soln             Ghosted vector
 * @param v      Vector that may be ghosted or not ghosted
 */
void
ProblemResidEval::updateGhostEqns(Epetra_Vector_Ghosted * const soln, const Epetra_Vector * const v)
{
  int myelems = soln->MyLength();
  if (myelems != m_NumLcEqns) {
    throw m1d_Error("ProblemResidEval::updateGhostEqns", "Problem with myelems");
  }
  LI_ptr_->updateGhostEqns(soln, v);
}

//=====================================================================================================================
std::string
ProblemResidEval::getBaseFileName() const
{
  return m_baseFileName;
}
//=====================================================================================================================
void
ProblemResidEval::setBaseFileName(std::string baseFileName)
{
  m_baseFileName = baseFileName;
}
//=====================================================================================================================
void
ProblemResidEval::calcSolnOld(const Epetra_Vector_Ghosted &soln, const Epetra_Vector_Ghosted &solnDot, double rdelta_t)
{
  if (rdelta_t < 1.0E-200) {
    for (int i = 0; i < m_NumLcEqns; i++) {
      (*solnOld_ptr_)[i] = soln[i];
    }
  } else {
    for (int i = 0; i < m_NumLcEqns; i++) {
      (*solnOld_ptr_)[i] = soln[i] - solnDot[i] / rdelta_t;
    }
  }
}
//=====================================================================================================================
} // end of m1d namespace
//=====================================================================================================================
