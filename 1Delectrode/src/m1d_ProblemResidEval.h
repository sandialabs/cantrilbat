/**
 *  @file m1d_ProblemResidEval.h
 *  File that contains the description of a single m1d problem that
 *  will be solved.
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2013-05-13 10:57:58 -0600 (Mon, 13 May 2013) $
 *  $Revision: 592 $
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_PROBLEMRESIDEVAL_H
#define M1D_PROBLEMRESIDEVAL_H

#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "m1d_VBRIndices.h"
#include "m1d_LocalNodeIndices.h"
#include "m1d_EqnVarTypes.h"
#include "cantera/base/xml.h"


namespace m1d
{
// Forward declarations of m1d Classes
class DomainLayout;
class GlobalIndices;
class LocalNodeIndices;
class LocalRowNodeVBRIndices;
class EpetraJac;
class ProblemStatement;
class RecordTree_base;

/*!
 *  Below are the three distributed maps that are used throughout the program. In order
 *  to document the program we typedef the different types of Epetra_Vectors according
 *  to their underlying maps.
 */

//!  This is an epetra vector consisting of just owned nodes/equations
typedef Epetra_Vector Epetra_Vector_Owned;

//! This is an Epetra_Vector consisting of owned and ghost nodes/equations
typedef Epetra_Vector Epetra_Vector_Ghosted;

//! This is an Epetra_Vector consisting of a global vector of all nodes/equations on all processors.
typedef Epetra_Vector Epetra_Vector_GlAll;

//! The types of solution problems that are solved.
enum Solve_Type_Enum {
    SteadyState_Solve = 0,
    TimeDependentAccurate_Solve,
    //! Specify this when we don't have a previous initially solved solution
    //!   -> in other words the electrode objects don't have a history for a previous time step
    TimeDependentInitial,
    TimeDependentRelax_Solve,
    //! Initial conditions for a DAE system
    DAESystemInitial_Solve
};
  

//! Differentiates the type of residual evaluations according to functionality
enum ResidEval_Type_Enum
{
    //! Base residual calculation for the time-stepping function
    Base_ResidEval = 0,
    //! Base residual calculation for the Jacobian calculation
    JacBase_ResidEval,
    //! Delta residual calculation for the Jacbobian calculation
    JacDelta_ResidEval,
    //! Base residual calculation for the showSolution routine
    /*!
     *    We calculate this when we want to display a solution
     */
    Base_ShowSolution
};

//! Differentiates the type of coordinate system
enum CoordinateSystem_Type_Enum
{
    //! Rectilinear coordinate system
    Rectilinear_Coordinates = 0,

    //! Cylindrical coordinate system
    Cylindrical_Coordinates,

    //! Spherical coordinate system
    Spherical_Coordinates
};

//! Class to hold Porosity equation status
struct Porosity_EqnType_Status
{
    //! Enum indicating the porosity equation properties
    enum Porosity_EqnType_Enum
    {
	//! There is no porosity
	None     =                   0x00,
	
	//! There is a constant porosity
	Constant =                   0x01,
	
	//! It is not constant, but calculated outside of equation system
	CalculatedOutOfEqnSystem =   0x02,
	
	//! Calculated as a function of the equation system
	CalculatedInEqnSystem  =     0x04,
	
	//! There are mechanics calculations which change the porosity
	PartOfMechanics =            0x08,
	
	//! Added phases are in the equation system. These may be nucleation phases
	AddedPhasesInEqnSystem =     0x16
    }; 
};

/**
 *  A class for the description of 1D problems that
 *  encompass multiple regions in 1D and multiple time regions
 */
class ProblemResidEval
{

public:

    //! Default constructor
    /*!
     *
     * @param atol   Absolute tolerance for calculation
     */
    ProblemResidEval(double atol = 1.0e-13);

  //! Destructor
  /*!
   *
   */
  virtual
  ~ProblemResidEval();

  //! Default copy constructor
  /*!
   *
   * @param r  Object to be copied
   */
  ProblemResidEval(const ProblemResidEval &r);

  //! Assignment operator
  /*!
   *
   * @param r  Object to be copied
   * @return   Returns a copy of the current problem
   */
  ProblemResidEval &
  operator=(const ProblemResidEval &r);

  //! Specify the problem
  /*!
   *  Initialize the domain structure for the problem.
   *
   *  @param problemType Problem type
   */
  void
  specifyProblem(int problemType, ProblemStatement *ps_ptr);

  //! Specify the problem
  /*!
   *  Initialize the domain structure for the problem.
   *
   *  @param dl DomainLayout object that specifies most of the problem
   */
  void
  specifyProblem(DomainLayout *dl,  ProblemStatement *ps_ptr);

  //! Create the global indices for the problem
  /*!
   *
   */
  void
  generateGlobalIndices();

  //!  Prepare all of the pointers for a fast, efficient residual calculation
  /*!
   *
   */
  void
  domain_prep();

  
  //! Link the ProblemStatement class with the ProblemResidEval class and other classes
  //! that need input from the user
  /*!
   *  @param psInput_ptr   Pointer to the ProblemStatement class
   */
  virtual void link_input(ProblemStatement *psInput_ptr);
  

  //! Generate and fill up the local node vectors on this processor
  /*!
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
  generateLocalIndices();

  //! Fill the vector isAlgebraic with flags indicating whether the degree of freedom
  //! is a DAE or not
  /*!
   *     0    Regular dof with time derivative
   *     1    DAE variable
   *     2    DAE variable in the regular time dependent residual.
   *          However, the variable is a time derivative in the special DAE initial solve system
   *          (e.g., mole fraction sum to 1 variable).
   *
   *  @param isAlgebraic This is the vector to be filled up
   *   
   */
  void fillIsAlgebraic(Epetra_IntVector & isAlgebraic);


  //! Fill the vector isAlgegraicScaled with flags indicating whether the degree of freedom
  //! is on an arithmetic scale and should be handled differently
  /*!
   *     0    Regular dof that can't go below the origin.
   *     1    arithmetic scaled variable
   *
   *  @param isAlgebraicScaled. This is the vector to be filled up
   *   
   */
  void fillIsArithmeticScaled(Epetra_IntVector & isAlgebraicScaled);



  //! Calculate a residual vector
  /*!
   *   The basic algorithm is to loop over the volume domains.
   *   Then, we loop over the surface domains
   *
   * @param res                     residual output
   * @param doTimeDependentResid    Boolean indicating whether the time
   *                                dependent residual is requested
   * @param doTimeDependentResid    Boolean indicating whether we should
   *                                formulate the time dependent residual
   * @param soln                    Pointer to the solution vector. This is the input to the residual calculation.
   * @param solnDot                 Pointer to the solution Dot vector. This is the input to the residual calculation.
   * @param t                       current time
   * @param rdelta_t                delta t inverse
   * @param residType               Residual type
   * @param solveType               Solve type
   */
  virtual void
  residEval(Epetra_Vector* const &  res,  
	    const bool doTimeDependentResid,
            const Epetra_Vector_Ghosted *soln,
            const Epetra_Vector_Ghosted *solnDot,
            const double t,
            const double rdelta_t,
            const ResidEval_Type_Enum residType = Base_ResidEval,
	    const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);

  //! Function that gets called at end the start of every time step
  /*!
   *  This function provides a hook for a residual that gets called whenever a
   *  time step has been accepted and we are about to move on to the next time step.
   *  The call is made with the current time as the time
   *  that is accepted. The old time may be obtained from t and rdelta_t_accepted.
   *
   *  After this call interrogation of the previous time step's results will not be
   *  valid.
   *
   *  This call also calculates all of the "old" cell information for the residual calculation.
   *  The "old" values are storred from calculation of the "current" values.
   *
   *  Note, when t is equal to t_old, soln_ptr should equal solnOld_ptr values. However,
   *  solnDot_ptr values may not be zero.
   *
   *   @param  doTimeDependentResid  This is true if we are solving a time dependent
   *                                 problem.
   *   @param  soln_ptr              Solution value at the current time
   *   @param  solnDot_ptr           derivative of the solution at the current time.
   *   @param  solnOld_ptr           Solution value at the old time step, n-1
   *   @param  t                     current time to be accepted, n
   *   @param  t_old                 previous time step value, t_old may be equal to t, 
   *                                 When we are calculating the initial conditions we
   *                                 require that we have values of "old" cell information.
   *                                 The call to this routine calculates the "old" information.
   */ 
  virtual void
  advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector *soln_ptr,
		      const Epetra_Vector *solnDot_ptr, const Epetra_Vector *solnOld_ptr,
		      const double t, const double t_old);

   //! Revert the Residual object's conditions to the conditions at the start of the global time step
   /*!
    *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init_init.
    *  We get rid of the pendingIntegratedFlags_ flag here as well.
    */
   virtual void 
   revertToInitialGlobalTime();

  //! Set the underlying state of the system from the solution vector
  /*!
   *   Note this is an important routine for the speed of the solution.
   *   It would be great if we could supply just exactly what is changing here.
   *   This routine is always called at the beginning of the residual evaluation process.
   *
   *   This is a natural place to put any precalculations of nodal quantities that
   *   may be needed by the residual before its calculation.
   *
   * @param doTimeDependentResid   If true then we are doing a time step. If false then we are not
   *                               doing a time step. 
   * @param soln                   Solution values
   * @param solnDot                Solution derivative values
   * @param t         Time to apply the changes. Time can be between delta_t+t_old and t_old
   * @param delta_t   Global time step number
   * param  t_old     Old time step value
   */
  virtual void
  setStateFromSolution(const bool doTimeDependentResid,
                       const Epetra_Vector_Ghosted *soln,
                       const Epetra_Vector_Ghosted *solnDot,
                       const double t,
                       const double delta_t,
                       const double t_old);

  //! Calculate the initial conditions
  /*!
   *   The basic algorithm is to loop over the volume domains.
   *   Then, we loop over the surface domains
   *
   * @param doTimeDependentResid    Boolean indicating whether we should
   *                                formulate the time dependent residual
   * @param soln                    Solution vector. This is the input to
   *                                the residual calculation.
   * @param solnDot                 Solution vector. This is the input to
   *                                the residual calculation.
   * @param t                       Time
   * @param delta_t                 delta_t for the current time step
   * @param delta_t_np1             delta_t_np1 for the next step. This 
   *                                may come from a saved solution.
   */
  void
  initialConditions(const bool doTimeDependentResid,
                    Epetra_Vector_Ghosted *soln,
                    Epetra_Vector_Ghosted *solnDot,
                    double& t,
                    double& delta_t,
                    double& delta_t_np1);

  //! Allocate and create storage for the Jacobian matrix for this problem
  void
  createMatrix(RecordTree_base *linearSolver_db);

  //! Return a reference to the Jacobian
  EpetraJac&
  jacobian();

  //! Returns the ghosted map of the problem
  Epetra_BlockMap *
  ghostedMap() const;

  //! Returns the unique unknown map of the problem
  Epetra_BlockMap * ownedMap() const;

  //! Returns the total number of global equations in the equation set
  int
  NumGbEqns() const;

  //! Returns the total number of global nodes in the problem
  int
  NumGbNodes() const;

  //! Returns information about the local equation number
  /*!
   *  @param ilocalEqn     INPUT      Equation number
   *  @param iLcNode       OUTPUT     Local node number
   *  @param iGbNode       OUTPUT     Global node number
   *  @param iNodeEqnNum   OUTPUT     Local node equation offset
   *  @param eqn           OUTPUT     Equation type 
   *  @param vtsub         OUTPUT     Equation subtype
   *
   *  @return returns the string id of the Equation variable
   */
  std::string equationID(int ilocalEqn, int& iLcNode, int &iGbNode, int& iNodeEqnNum, EqnType &var,
			 EQ_TYPE_SUBNUM &vtsub);

  //! Returns information about the local variable number
  /*!
   *  @param ilocalVar    INPUT        Variable number
   *  @param iLcNode      OUTPUT       Local node number
   *  @param iGbNode       OUTPUT     Global node number
   *  @param iNodeEqnNum   OUTPUT    Local node equation offset
   *  @param var type  OUTPUT   
   *  @param vtsub   OUTPUT     Variable subtype
   *
   *  @return returns the string id of the variable.
   */
  std::string variableID(int ilocalVar, int& iLcNode, int &iGbNode, int& iNodeEqnNum, VarType &var,
                         VAR_TYPE_SUBNUM &vtsub);


  //! Update the ghost unknowns on the processor.
  /*!
   *
   * @param soln             Ghosted vector
   * @param v      Vector that may be ghosted or not ghosted
   */
  void
  updateGhostEqns(Epetra_Vector_Ghosted * const soln, const Epetra_Vector * const v);

  /****************************************************************************/
  /*                   Sets of indexing functions                             */
  /****************************************************************************/

  //! Return the global equation number given the global block number
  /*!
   *
   * @param gbBlockNum
   * @param localRowNumInBlock
   * @return  the global equation number
   */
  int
  GbEqnNum_From_GbBlock(const int gbBlockNum, const int localRowNumInBlock);

  //! Global node to Local node mapping
  /*!
   * Given a global node, this function will return the local node value.
   * If the global node is not on this processor, then this function returns -1.
   *
   * This is a pass through routine to Local Node Class which already has implemented this
   *
   * @param gbNode   Global node
   * @return returns the local node value. Returns -1 if the local node is not on this
   *         processor.
   */
  int
  GbNodeToLcNode(const int gbNode) const;

  //! Global equation to local equation mapping
  /*!
   *    Returns the local equation number of the global equation
   *    Note, if the local equation number isn't on the processor, this routine
   *    will return -1.
   *
   * @param gbEqn  Input global equation number
   * @return     Returns the local equation number of the global equation
   *             Note, if the local equation number isn't on the processor, this routine
   *             will return -1.
   */
  int
  GbEqnToLcEqn(const int gbEqn) const;


  //!  Find a delta of a solution component for use in the numerical jacobian
  /*!
   *    @param soln      Reference to the complete solution vector
   *    @param ieqn      local equation number of the solution vector
   *
   *  @return returns the delta for the component
   */
  double
  deltaSolnCompJac(const Epetra_Vector_Ghosted & soln, const int ieqn);

  //! Evaluates the atol vector used in the time stepper routine and in the nonlinear solver.
  /*!
   *  @param atolDefault  Double containing the default value of atol
   *  @param soln         Current solution vector.
   *  @param atolVector   Optional vector containing the individual entries for the
   *                               absolute tolerance. Use this if you want to
   *                              override the default calculation of atol.
   *                              Defaults to zero.
   */
  virtual void
  setAtolVector(double atolDefault, 
		const Epetra_Vector_Ghosted & soln,
		const Epetra_Vector_Ghosted * const atolVector = 0);


  //! Evaluates the atol DAESystemInitial vector used in the time stepper routine and in the nonlinear solver.
  /*!
   *    The default behavior is to set all of the atol values to the constant atolDAEInitDefault.
   *
   *   @param atolDefault         Double containing the default value of atol used in the
   *                              DAESystemInitial solve
   *   @param soln                Current solution vector.
   *   @param solnDot             Current solution vector.
   *   @param atolDAEInitVector   Optional vector containing the individual entries for the time
   *                              derivative absolute tolerance. Use this if you want to
   *                              override the default calculation of atol.
   *                              Defaults to zero.
   */
  virtual void
  setAtolVector_DAEInit(double atolDAEInitDefault, const Epetra_Vector_Ghosted & soln,
			const Epetra_Vector_Ghosted & solnDot,
			const Epetra_Vector_Ghosted * const atolVector_DAEInit = 0) const;

  //! Evaluates the atol vector used in the delta damping process.
  /*!
   *   @param relcoeff     Relative constant to multiply all terms by
   *   @param soln         current solution vector.
   *   @param atolDeltaDamping      If non-zero, this copies the vector into the object as input
   *                      The default is zero.
   */
  virtual void
  setAtolDeltaDamping(double relcoeff, 
		      const Epetra_Vector_Ghosted & soln,
		      const Epetra_Vector_Ghosted * const atolDeltaDamping = 0);

  //! Evaluates the atol vector used in the delta damping process for the DAE problem
  /*!
   *   @param relcoeff     Relative constant to multiply all terms by
   *   @param soln         current solution vector.
   *   @param solnDot      Current solutionDot vector.
   *   @param atolDeltaDamping       If non-zero, this copies the vector into the object as input
   *                       The default is zero.
   */
  virtual void
  setAtolDeltaDamping_DAEInit(double relcoeff, 
		      const Epetra_Vector_Ghosted & soln,
		      const Epetra_Vector_Ghosted & solnDot,
		      const Epetra_Vector_Ghosted * const atolDeltaDamping = 0);

  //! Evaluate the delta damping vector from the abs tol vector
  /*!
   *  The default behavior is to use the atolVector values to fill this vector.
   *
   *  @param relCoeff sets a multiplicative constant to the atol vector
   */
  //virtual void
  //setDeltaDamping(double relCoeff = 1.0);

  //! Return a const reference to the atol vector
  /*!
   *  @return const reference to the atol vector
   */
  const Epetra_Vector_Owned & atolVector() const;

  //! Return a const reference to the atolDAEInit vector
  /*!
   *  @return const reference to the atolDAEInit vector
   */
  const Epetra_Vector_Owned & atolVector_DAEInit() const;

  //! Return a const reference to the atol_deltaDamping vector
  /*!
   *  @return const reference to the atol delta damping vector
   */
  const Epetra_Vector_Owned & atolDeltaDamping() const;

  virtual void
  filterSolnPrediction(double t, Epetra_Vector_Ghosted &y);

  void
  applyFilter(const double timeCurrent,
              const double delta_t_n,
              const Epetra_Vector_Ghosted &y_current,
              const Epetra_Vector_Ghosted &ydot_current,
              Epetra_Vector_Ghosted &delta_y);

  double
  delta_t_constraint(const double time_n, const Epetra_Vector_Ghosted &y_n, const Epetra_Vector_Ghosted &ydot_n);

  //! Evaluate a supplemental set of equations that are not part of the solution vector, but are considered
  //! to be time dependent
  /*!
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
  virtual void
  evalTimeTrackingEqns(const int ifunc,
                       const double t,
                       const double deltaT,
                       const Epetra_Vector_Ghosted & y,
                       const Epetra_Vector_Ghosted * const ydot);

  //! Set a solution parameter 
  /*!
   *  @param paramName   String identifying the parameter to be set
   *  @param paramVal    Single double value of the parameter to be set
   *
   *  @return returns a 0 if the number makes sense
   *          Returns a negative number if the parameter is unknown
   */
  virtual int
  setSolutionParam(std::string paramName, double paramVal);


  //! Get a solution parameter 
  /*!
   *  @param paramName   String identifying the parameter to be set
   *  @param paramVal    Vector of parameters returned
   *
   *  @return returns the number of parameters returned.
   */
  virtual int
  getSolutionParam(std::string paramName, double * const paramVal);


  virtual bool
  evalStoppingCritera(double &time_current,
                      double &delta_t_n,
                      const Epetra_Vector_Ghosted &y_n,
                      const Epetra_Vector_Ghosted &ydot_n);

  /**
   * Return a vector of delta y's for calculation of the
   *  numerical Jacobian
   */
  virtual void
  calcDeltaSolnVariables(const double t, const Epetra_Vector& soln,
			 const Epetra_Vector* solnDot_ptr, Epetra_Vector& deltaSoln,
                         const Solve_Type_Enum solveType = TimeDependentAccurate_Solve,
                         const  Epetra_Vector* solnWeights=0);

  //! Save the solution to the end of an XML file using XML solution format
  /*!
   *  We write out the solution to a file.
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
   *
   *      @param delta_t_np1       Suggested next delta t value (defaults to 0.0 if there isn't a
   *                               good guess for the next delta_t).
   */
  virtual void
  saveSolutionEnd(const int itype,
                  std::string baseFileName,
                  const Epetra_Vector_Ghosted &m_y_n,
                  const Epetra_Vector_Ghosted *m_ydot_n,
                  const double t,
                  const double delta_t,
                  const double delta_t_np1 = 0.0);

  //! Read the solution from a saved file using the record number to identify the solution 
  /*!
   *  We read only successful final steps.
   *
   * @param iNumber            Solution Record number to read from. Solutions are numbered consequetively from 
   *                           zero to the number of global time steps. We will use the last
   *                           time step if there aren't as many time steps in the solution file
   *                           as requested.
   *
   * @param baseFileName       File name to be used. .xml is appended onto the filename.
   *                           Processors other than 0 have the pid appended to the name as well.
   *
   * @param y_n_ghosted        Current value of the solution vector
   * @param ydot_n_ghosted     Current value of the derivative of the solution vector
   * @param t_read             time that is read in
   * @param delta_t_read       delta time step for the last time step.
   * @param delta_t_next_read  delta time step for the next time step if available
   */
  void
  readSolutionRecordNumber (const int itype,
	       std::string baseFileName,
	       Epetra_Vector_Ghosted &y_n_ghosted,
	       Epetra_Vector_Ghosted * const ydot_n_ghosted,
	       double &t_read,
	       double &delta_t_read,
	       double &delta_t_next_read);

  //! Read the solution from a saved file. (deprecated)
  /*!
   *  We read only successful final steps
   *
   * @param iNumber            Solution number to read from. Solutions are numbered consequetively from 
   *                           zero to the number of global time steps. We will use the last
   *                           time step if there aren't as many time steps in the solution file
   *                           as requested.
   *
   * @param baseFileName       File Name to be used. .xml is appended onto the filename
   *                           processors other than 0 have the pid appended to the name
   * @param y_n_ghosted        Current value of the solution vector
   * @param ydot_n_ghosted     Current value of the derivative of the solution vector
   * @param t_read             time that is read in
   * @param delta_t_read       delta time step for the last time step.
   * @param delta_t_next_read  delta time step for the next time step if available
   */
  void
  readSolution(const int itype,
               std::string baseFileName,
               Epetra_Vector_Ghosted &y_n_ghosted,
               Epetra_Vector_Ghosted * const ydot_n_ghosted,
               double &t_read,
               double &delta_t_read,
               double &delta_t_next_read);

  //! Read the solution from a saved solution XML record.
  /*!
   * This is the underlying program that reads the record
   *
   * @param simulRecord        XML element to read from. Name must be a simulation XML_Node.
   * @param y_n_ghosted        Current value of the solution vector
   * @param ydot_n_ghosted     Current value of the derivative of the solution vector
   * @param t_read             time that is read in
   * @param delta_t_read       delta time step for the last time step.
   * @param delta_t_next_read  delta time step for the next time step if available
   */
  void
  readSolutionXML(Cantera::XML_Node* simulRecord, Epetra_Vector_Ghosted &y_n_ghosted,
		  Epetra_Vector_Ghosted * const ydot_n_ghosted, double &t_read,
		  double &delta_t_read, double &delta_t_next_read);

  //!  Select the global time step increment record by the consequuatively numbered record index number
  /*
   *    @param   xSoln               Solution file for the simulation object
   *    @param   globalTimeStepNum   Time step number to select
   *
   *    @return Returns a pointer to the selected record.
   */
  Cantera::XML_Node* selectSolutionRecordNumber(Cantera::XML_Node* xSoln, int globalTimeStepNum);
 
  //!  Select the global time step increment record by the consequuatively numbered record index number
  /*
   *    @param   xSoln               Solution file for the simulation object
   *    @param   solnTimeStepID      string value of the time step ID to select
   *
   *    @return Returns a pointer to the selected record.
   */
  Cantera::XML_Node* selectSolutionTimeStepID(Cantera::XML_Node* xSoln, std::string solnTimeStepID);

  //! Write the solution to either the screen or to a log file
  /*!
   *  This is a general output utility to Cantera's logfile.
   *  It's not hooked into the IO algorithm at all. It should be
   *  conditionally called depending on the whims of the user.
   *
   * @param ievent  Type of the event. The following form is used:
   *             0 Initial conditions
   *             1 Completion of a successful intermediate step.
   *             2 Final successful conditions.
   *             3 Intermediate nonlinear step
   *            -1 unsuccessful step
   * @param doTimeDependentResid   Do the time dependent residual calculation
   * @param t                      Current time
   * @param delta_t                delta t
   * @param y_n    Current value of the solution vector
   * @param ydot_n  Current value of the derivative of the solution vector
   */
  virtual void
  showProblemSolution(const int ievent,
		      bool doTimeDependentResid, 
                      const double t,
                      const double delta_t,
                      const Epetra_Vector_Owned &y_n,
                      const Epetra_Vector_Owned * const ydot_n,
		      const Solve_Type_Enum solveType = TimeDependentAccurate_Solve,
                      const double delta_t_np1 = 0.0);

  virtual void
  writeTecplot(const int ievent,
               std::string m_baseFileName,
	       bool doTimeDependentResid,
	       const double t,
	       const double delta_t,
	       const Epetra_Vector_Ghosted &y_n,
	       const Epetra_Vector_Ghosted * const ydot_n,
	       const Solve_Type_Enum solveType,
	       const double delta_t_np1);


  //! Write a solution vector type to either the screen or a log file
  /*!
   *   @param solnVecName    Name of the solution vector -> will appear in the output
   *   @param t              Current time
   *   @param delta_t        delta t
   *   @param solnVector     Vector of information to be printed
   */
  virtual void
  showSolutionVector(std::string& solnVecName,
                     const double t,
                     const double delta_t,
                     const Epetra_Vector_Owned &solnVector,
		     FILE *outFile = stdout);

 //! Write an int solution vector type to either the screen or a log file
  /*!
   *   @param solnVecName    Name of the solution vector -> will appear in the output
   *   @param t              Current time
   *   @param delta_t        delta t
   *   @param solnVector     Integer Vector of information to be printed
   */
  virtual void
  showSolutionIntVector(std::string& solnVecName,
                     const double t,
                     const double delta_t,
                     const Epetra_IntVector &solnVector,
                     FILE *outFile = stdout);

  //! Write out to a file or to standard output the current solution
  /*!
   *   These functions are affected by the print controls of the nonlinear solver
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
   *      @param doTimeDependentResid  true if solving a time dependent problem
   *      @param time_current      Current time
   *      @param delta_t_n         Current value of delta_t
   *      @param istep             Current step number
   *      @param soln_n               Current value of the solution vector
   *      @param solnDot_n_ptr            Current value of the time deriv of the solution vector
   *      @param delta_t_np1       Suggested next delta t value (defaults to 0.0 if there isn't a
   *                               good guess for the next delta_t).
   */
  virtual void
  writeSolution(const int ievent, 
		const bool doTimeDependentResid,
                const double time_current,
                const double delta_t_n,
                const int istep,
                const Epetra_Vector_Ghosted &soln_n,
                const Epetra_Vector_Ghosted * const solnDot_n_ptr,
		const Solve_Type_Enum solveType = TimeDependentAccurate_Solve,
                const double delta_t_np1 = 0.0);

  //! This function may be used to create output at various points in the
  //! execution of an application.
  /*!
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
  virtual void
  user_out(const int ievent,
           const double time_current,
           const double delta_t_n,
           const int istep,
           const Epetra_Vector_Ghosted &y_n,
           const Epetra_Vector_Ghosted * const ydot_n_ptr);

  virtual void
  matrixConditioning(double * const matrix, int nrows, double * const rhs);

  std::string
  getBaseFileName() const;

  void
  setBaseFileName(std::string m_baseFileName);

  void
  calcSolnOld(const Epetra_Vector_Ghosted &soln, const Epetra_Vector_Ghosted &solnDot, double rdelta_t);

public:

  double m_atol;

  //! Global total of the equations
  int m_neq;

  //! Number of equations on the processor including ghost equations
  int m_NumLcEqns;

  //! Number of equations on the processor not including ghost equations
  int m_NumLcOwnedEqns;

  //! DomainLayout is a light class that describes the overall
  //! disposition of the domains in the problem
  DomainLayout *DL_ptr_;

  //! Global indices class is the same on all processors
  /*!
   *  This is a global structure, containing global information.
   *  All information in this structure, except for my_procID, numLnNodes,
   *  and numLnEqns, is the same on all processors.
   *
   */
  GlobalIndices *GI_ptr_;

  //! Class setting up the local node indices on the current processor
  /*!
   *   These indices have to do with accessing the local nodes on the processor
   *
   *   All numbers are in local Row Node Format. Local Row node format is
   *   defined as the following. All owned nodes come first. They are
   *   ordered in terms of increasing global node number.
   *   Then the right ghost node is listed.
   *   Then the left ghost node is listed
   *   Then, the "globally-all-connected node is listed, if available.
   */
  LocalNodeIndices *LI_ptr_;

  //! Pointer to the jacobian calculator
  /*!
   *   We don't own this
   */
  EpetraJac *m_jac;

  std::string m_baseFileName;

  //! Step Number
  /*!
   *   Step number is the time step number for time dependent calculations.
   *   It's the number of the current solution, for time independent calculations,
   *   which have continuation numbers.
   */
  int m_StepNumber;

  //! Pointer to the old solution vector
  Epetra_Vector_Ghosted *solnOld_ptr_;

  //! Residual used for internal calculation
  Epetra_Vector_Owned *resInternal_ptr_;

  //! Absolute tolerances used for nonlinear solver convergences
  Epetra_Vector_Ghosted *m_atolVector;

  //! Absolute tolerances used for nonlinear solver convergences
  mutable Epetra_Vector_Ghosted *m_atolVector_DAEInit;

  //! Tolerances used for delta damping in the nonlinear solver
  Epetra_Vector_Ghosted *m_atolDeltaDamping;

  //! Tolerances used for delta damping in the nonlinear solver
  Epetra_Vector_Ghosted *m_atolDeltaDamping_DAEInit;


  //! Pointer to the ProblemStatement class
  ProblemStatement *psInput_ptr_;

public:

  //!  Number of time regions defined in the problem
  /*!
   *  A time region is a region of time where the boundary conditions are specified. Between
   *  time regions there can be step discontinuities in the boundary conditions.
   */
  int m_numTimeRegions;

  //! Current time region defined in the problem
  /*!
   *  This number varies from zero to m_numTimeRegions
   */
  int m_currentTimeRegion;

  //! Write a special solution consisting of the starting and ending condition
  /*!
   *  this is used as a restart file for ORNL CAEBAT program
   */
  int m_writeStartEndFile;

public:

  //! counter containing the number of Residual base case evaluations
  int counterResBaseCalcs_;

  //! Counter containing the number of Residual evaluations taken in doing the base residual needed for the jacobian calculation
  //! This is frequently equal to the number of jacobian calculations
  int counterJacBaseCalcs_;

  int counterJacDeltaCalcs_;
  int counterResShowSolutionCalcs_;

  //! Print flag that is used in conjunction with the printlvl. 
  /*!
   *  This variable is set by reading the environmental variable PRE_printFlagEnv at the start of the calculation.
   *  If this variable is set to 'y' or to an integer, the the int below is set to that variable. 
   */
  static int s_printFlagEnv;

  //! Print level that is set through the input file.
  /*!
   *   0 -> Don't print anything
   *   1 -> Print only about significant issues going on
   *   2 -> Print status information at regular intervals.
   *   3 -> Print ShowSolution at regular intervals
   *   4 -> Print ShowSolution at all successful time steps
   *   5 -> Print additional information at each time step
   *   6 -> Print some information about each electrode object at each time step
   *   7 -> Print a lot of information about each electrode object at each time step
   */
  int SolutionBehavior_printLvl_;

  //! Level of residual information printing done to stdout
  /*!
   *   0 -> Don't print anything
   *   1 -> Print only about significant issues going on
   *   2 -> Print status information at regular intervals.
   *   3 -> Print ShowResidual at regular intervals
   *   4 -> Print ShowResidual at all successful time steps
   *   5 -> Print additional information when ShowSolution is called.
   *   6 -> Print additional information when any residual is called.
   *   7 -> Print a lot of information about each when ShowSolution is called.
   *   8 -> Print a lot of information when any base or show residual is called
   *   9 -> Print a lot of information when any residual is called
   */
   int Residual_printLvl_;

   //! The type of coordinate system that is used
   /*!
    *  There are two that are envisioned: Rectinear_Coordinates and Cylindrical_Coordinates
    */
   CoordinateSystem_Type_Enum coordinateSystemType_;

   //! Cross sectional area, if in cartesian coordinates
   /*!
    *  The overwhelming output from the program is on a per-crosssectional area basis
    *  However, there are some times when the cross-section is needed. This is the place
    *  where it is supplied.
    *
    *    units m2
    */
   double crossSectionalArea_;

   //! Cylinder Length, if in cylindrical coordinates
   /*!
    *  The overwhelming output from the program is on a per-crosssectional area basis
    *  However, there are some times when the cross-section is needed. This is the place
    *  where it is supplied.  We assume 2 pi radians always, i.e., a full radius
    *
    *    units m
    */
   double cylinderLength_;

   //!  Int indicating the energy equation treatment
   /*!
    *  0 -> isothermal               Don't solve an energy equation (default)
    *  1 -> Fixed                    Don't solve an energy equation
    *  2 -> Dirichlet Equation       Solve a Dirichlet equation for temperature. 
    *                                This is a way to do the fixed system while keeping the
    *                                matrix structure the same.
    *  3 -> Enthalpy Equation        Solve a full enthalpy equation for the temperature
    *  4 -> Temperature Equation     Solve a Cp dT/dt formulation for the temperature
    */
   int energyEquationProbType_;

   //! Int indicationg the solid mechanics treatment
   /*!
    *   0 -> none                    Don't solve an equation for mesh motion due to stress effects
    *   1 -> Linear Elastic          Solve a simple treatment for treating the stress-strain effects
    */
   int solidMechanicsProbType_;

};

// *****************************************************************
//  end of m1d namespace
}
// ******************************************************************

#endif

