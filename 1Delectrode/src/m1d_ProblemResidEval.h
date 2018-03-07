/**
 *  @file m1d_ProblemResidEval.h
 *  File that contains the description of a single m1d problem that
 *  will be solved.
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_PROBLEMRESIDEVAL_H
#define M1D_PROBLEMRESIDEVAL_H

#include "m1d_defs.h"
#include "m1d_EqnVarTypes.h"


#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_VbrMatrix.h"

#include "cantera/numerics/ResidJacEval.h"
#include "cantera/base/xml.h"

#include "m1d_LocalNodeIndices.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

// Forward declarations of m1d Classes
class GlobalIndices;
class LocalNodeIndices;
class EpetraJac;
class ProblemStatement;
class DomainLayout;
class RecordTree_base;

//==================================================================================================================================
/*!
 *  Below are the three distributed maps that are used throughout the program. In order
 *  to document the program we typedef the different types of Epetra_Vectors according
 *  to their underlying maps.
 */

//!  This is an epetra vector consisting of just owned nodes/equations
//typedef Epetra_Vector Epetra_Vector_Owned;

//! This is an Epetra_Vector consisting of owned and ghost nodes/equations
//typedef Epetra_Vector Epetra_Vector_Ghosted;

//! This is an Epetra_Vector consisting of a global vector of all nodes/equations on all processors.
//typedef Epetra_Vector Epetra_Vector_GlAll;

//==================================================================================================================================
//! The types of solution problems that are solved.
/*
enum class Solve_Type {
    //! Steady state solve
    SteadyState_Solve = 0,
    //! Time dependent solve that is time accurate
    TimeDependentAccurate_Solve,
    //! Specify this when we don't have a previous initially solved solution
    //!   -> in other words the electrode objects don't have a history for a previous time step
    TimeDependentInitial,
    //! We solve a time dependent problem only to relax the system to steady state
    TimeDependentRelax_Solve,
    //! Initial conditions for a DAE system
    DAESystemInitial_Solve
};
*/
// take this out eventually
//using Zuzax::Solve_Type;
//==================================================================================================================================
//! Differentiates the type of coordinate system
enum CoordinateSystem_Type_Enum {
    //! Rectilinear coordinate system
    Rectilinear_Coordinates = 0,

    //! Cylindrical coordinate system
    Cylindrical_Coordinates,

    //! Spherical coordinate system
    Spherical_Coordinates
};
//==================================================================================================================================
//! Class to hold Porosity equation status
struct Porosity_EqnType_Status {
    //! Enum indicating the porosity equation properties
    enum Porosity_EqnType_Enum {
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
//==================================================================================================================================
//!  A class for the description of 1D problems that  encompass multiple regions in 1D and multiple time regions
/*!
 *   This is the base class for 
 *  todo make it inherit from Zuzax::ResidJacEval.  Rename Epetra routines that override
 */
class ProblemResidEval
//class ProblemResidEval : public Zuzax::ResidJacEval  (under construction)
{
public:
    //! Default constructor
    /*!
     *  @param[in]           atol                Absolute tolerance for calculation
     */
    ProblemResidEval(double atol = 1.0e-13);

    //! Virtual Destructor
    /*!
     *
     */
    virtual ~ProblemResidEval();

    //! Default copy constructor
    /*!
     *  @param[in]           r                   Object to be copied
     */
    ProblemResidEval(const ProblemResidEval& r);

    //! Assignment operator
    /*!
     *  @param[in]           r                   Object to be copied
     *  @return                                  Returns a copy of the current problem
     */
    ProblemResidEval& operator=(const ProblemResidEval& r);

    //! Specify the problem
    /*!
     *  Initialize the domain structure for the problem and initialize the parameters for the problem by sending
     *  in the ProblemStatement structure
     *
     *  @param[in]           problemType         Problem type
     *  @param[in]           ps_ptr              Pointer to the ProblemStatement object that contains the parameters for the problem
     */
    void specifyProblem(int problemType, ProblemStatement* ps_ptr);

    //! Specify the problem
    /*!
     *  Initialize the domain structure for the problem.
     *
     *  @param[in]             dl                  DomainLayout object that specifies most of the problem
     *  @param[in]           ps_ptr              Pointer to the ProblemStatement object that contains the parameters for the problem
     */
    void specifyProblem(DomainLayout* dl, ProblemStatement* ps_ptr);

    //! Create the global indices for the problem
    /*!
     *  Initialize the global indices for the problem
     */
    void generateGlobalIndices();

    //!  Prepare all of the pointers for a fast, efficient residual calculation
    /*!
     *
     */
    void domain_prep();


    //! Link the ProblemStatement class with the ProblemResidEval class and other classes
    //! that need input from the user
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  @param[in]             psInput_ptr         Pointer to the ProblemStatement class
     */
    virtual void link_input(ProblemStatement* psInput_ptr);

    //! Generate and fill up the local node vectors on this processor
    /*!
     *  These indices have to do with accessing the local nodes on the processor. At the end of this routine,
     *  the object LocalNodeIncides is fully formed on this processor. Within this object, all numbers are in local Row Node Format.
     *  Local Row node format is defined as the following. All owned nodes come first. They are ordered 
     *  in terms of increasing global node number. Then the right ghost node is listed.
     *  Then the left ghost node is listed. Then, the "globally-all-connected node is listed, if available.
     */
    void generateLocalIndices();

    //! Fill the vector isAlgebraic with flags indicating whether the degree of freedom is a DAE or not
    /*!
     *  The following entries are used:
     *
     *     0    Regular dof with time derivative
     *     1    DAE variable
     *     2    DAE variable in the regular time dependent residual.
     *          However, the variable is a time derivative in the special DAE initial solve system
     *          (e.g., mole fraction sum to 1 variable).
     *
     *  @param[out]          isAlgebraic         This is the vector to be filled up
     */
    void fillIsAlgebraic(Epetra_IntVector& isAlgebraic);

    //! Fill the vector isAlgegraicScaled with flags indicating whether the degree of freedom
    //! is on an arithmetic scale and should be handled differently
    /*!
     *  The following entries are used
     *
     *     0    Regular dof that can't go below the origin.
     *     1    arithmetic scaled variable
     *
     *  @param[out]          isAlgebraicScaled  This is the vector to be filled up
     */
    void fillIsArithmeticScaled(Epetra_IntVector& isAlgebraicScaled);

    //! Calculate a residual vector
    /*!
     *  (virtual from ProblemResidEval)
     *
     *   The basic algorithm is to loop over the volume domains. Then, we loop over the surface domains.
     *
     *  @param[out]          res                 Residual output
     *  @param[in]          doTimeDependentResid Boolean indicating whether the time dependent residual is requested
     *  @param[in]           soln                Pointer to the solution vector. This is the input to the residual calculation.
     *  @param[in]           solnDot             Pointer to the solution Dot vector. This is the input to the residual calculation.
     *  @param[in]           t                   Current time
     *  @param[in]           rdelta_t            Delta t inverse
     *  @param[in]           residType           Residual type. The default is Base_ResidEval.
     *  @param[in]           solveType           Type of the problem being solved expressed as a  Solve_Type_Enum. 
     *                                           Defaults to TimeDependentAccurate_Solve
     */
    virtual void residEval(Epetra_Vector* const&  res, const bool doTimeDependentResid,
                           const Epetra_Vector_Ghosted* const soln, const Epetra_Vector_Ghosted* const solnDot,
                           const double t, const double rdelta_t,
                           const Zuzax::ResidEval_Type residType = Zuzax::ResidEval_Type::Base_ResidEval,
                           const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve);

    //! Function that gets called at end the start of every time step
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  This function provides a hook for a residual that gets called whenever a
     *  time step has been accepted and we are about to move on to the next time step.
     *  The call is made with the current time as the time
     *  that is accepted. The old time may be obtained from t and rdelta_t_accepted.
     *
     *  After this call interrogation of the previous time step's results will not be valid.
     *
     *  This call also calculates all of the "old" cell information for the residual calculation.
     *  The "old" values are storred from calculation of the "current" values.
     *
     *  Note, when t is equal to t_old, soln_ptr should equal solnOld_ptr values. However, solnDot_ptr values may not be zero.
     *
     *   @param[in]          doTimeDependentResid  This is true if we are solving a time dependent problem.
     *   @param[in]          soln_ptr              Solution value at the current time
     *   @param[in]          solnDot_ptr           derivative of the solution at the current time.
     *   @param[in]          solnOld_ptr           Solution value at the old time step, n-1
     *   @param[in]          t                     current time to be accepted, n
     *   @param[in]          t_old                 Previous time step value, t_old may be equal to t,
     *                                             When we are calculating the initial conditions we
     *                                             require that we have values of "old" cell information.
     *                                             The call to this routine calculates the "old" information.
     */
    virtual void advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                                     const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr,
                                     const double t, const double t_old);

    //! Revert the Residual object's conditions to the conditions at the start of the global time step
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init_init.
     *  We get rid of the pendingIntegratedFlags_ flag here as well.
     */
    virtual void revertToInitialGlobalTime();

    //! Set the underlying state of the system from the solution vector
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  Note this is an important routine for the speed of the solution.
     *  It would be great if we could supply just exactly what is changing here.
     *  This routine is always called at the beginning of the residual evaluation process.
     *
     *  This is a natural place to put any precalculations of nodal quantities that
     *  may be needed by the residual before its calculation.
     *
     *  @param[in]           doTimeDependentResid  If true then we are doing a time step. If false then we are not doing a time step
     *  @param[in]           soln                  Solution values
     *  @param[in]           solnDot               Solution derivative values
     *  @param[in]           t                     Time to apply the changes. Time can be between delta_t+t_old and t_old
     *  @param[in]           delta_t               Global time step number
     *  @param[in]           t_old                 Old time step value
     */
    virtual void setStateFromSolution(const bool doTimeDependentResid, const Epetra_Vector_Ghosted* const soln,
                                     const Epetra_Vector_Ghosted* const solnDot, const double t,
                                     const double delta_t, const double t_old);

    //! Calculate the initial conditions
    /*!
     *  (virtual from ProblemResidEval)
     *
     *   Note this is an important routine for the speed of the solution.
     *   The basic algorithm is to loop over the volume domains.
     *   Then, we loop over the surface domains
     *
     *  @param[in]           doTimeDependentResid  Boolean indicating whether we should formulate the time dependent residual
     *  @param[out]          soln                  Solution vector. This is the input to the residual calculation.
     *  @param[out]          solnDot               Solution vector time derivative. This is the input to the residual calculation.
     *  @param[in]           t                     Current time
     *  @param[in]           delta_t               Delta_t for the current time step
     *  @param[in]           delta_t_np1           Delta_t_np1 for the next step. This may come from a saved solution.
     */
    virtual void initialConditions(const bool doTimeDependentResid, Epetra_Vector_Ghosted* const soln,
                                   Epetra_Vector_Ghosted* const solnDot, double& t, double& delta_t,
                                   double& delta_t_np1);

    //! Allocate and create storage for the Jacobian matrix for this problem
    /*!
     *  @param[in]           linearSolver_db     Linear solver database
     */
    void createMatrix(RecordTree_base* linearSolver_db);

    //! Return a reference to the Jacobian
    /*!
     *  @return                                  Returns a reference to the jacobian
     */
    EpetraJac& jacobian();

    //! Returns the ghosted map of the problem
    /*!
     *  @return                                  Returns a pointer to a BlockMap of the ghosted Solution vector
     */
    Epetra_BlockMap* ghostedMap() const;

    //! Returns the unique unknown map of the problem
    /*!
     *  @return                                  Returns a pointer to a BlockMap of the unique unknown Solution vector
     */
    Epetra_BlockMap* ownedMap() const;

    //! Returns the total number of global equations in the equation set
    /*!
     *  @return                                  Returns the number of global equations
     */
    int NumGbEqns() const;

    //! Returns the total number of global nodes in the problem
    /*!
     *  @return                                  Returns the number of global unknowns
     */
    int NumGbNodes() const;

    //! Returns information about the local equation number
    /*!
     *  @param[in]           ilocalEqn           Local equation number
     *
     *  @param[out]          iLcNode             Local node number
     *  @param[out]          iGbNode             Global node number
     *  @param[out]          iNodeEqnNum         Local node equation offset
     *  @param[out]          eqn                 Equation type
     *  @param[out]          vtsub               Equation subtype
     *
     *  @return                                  returns the string id of the Equation variable
     */
    std::string equationID(int ilocalEqn, int& iLcNode, int& iGbNode, int& iNodeEqnNum, EqnType& eqn,
                           EQ_TYPE_SUBNUM& vtsub);

    //! Returns information about the local variable number
    /*!
     *  @param[in]           ilocalVar           Variable number
     *  @param[out]          iLcNode             Local node number
     *  @param[out]          iGbNode             Global node number
     *  @param[out]          iNodeEqnNum         Local node equation offset
     *  @param[out]          var                 Variable type
     *  @param[out]          vtsub               Variable subtype
     *
     *  @return                                  returns the string id of the variable.
     */
    std::string variableID(int ilocalVar, int& iLcNode, int& iGbNode, int& iNodeEqnNum, VarType& var,
                           VAR_TYPE_SUBNUM& vtsub);

    //! Update the ghost unknowns on the processor.
    /*!
     *
     * @param[out]           soln                Ghosted vector
     * @param[in]            v                   Vector that may be ghosted or not ghosted
     */
    void updateGhostEqns(Epetra_Vector_Ghosted* const soln, const Epetra_Vector* const v);

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
    int GbEqnNum_From_GbBlock(const int gbBlockNum, const int localRowNumInBlock);

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
    int GbNodeToLcNode(const int gbNode) const;

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
    int GbEqnToLcEqn(const int gbEqn) const;


    //!  Find a delta of a solution component for use in the numerical jacobian
    /*!
     *    @param soln      Reference to the complete solution vector
     *    @param ieqn      local equation number of the solution vector
     *
     *  @return returns the delta for the component
     */
    double deltaSolnCompJac(const Epetra_Vector_Ghosted& soln, const int ieqn);

    //! Evaluates the atol vector used in the time stepper routine and in the nonlinear solver.
    /*!
     *  (virtual from ProblemResidEval)
     *
     *   Note this is an important routine for the speed of the solution.
     *  @param atolDefault  Double containing the default value of atol
     *  @param soln         Current solution vector.
     *  @param atolVector   Optional vector containing the individual entries for the
     *                               absolute tolerance. Use this if you want to
     *                              override the default calculation of atol.
     *                              Defaults to zero.
     */
    virtual void setAtolVector(double atolDefault, const Epetra_Vector_Ghosted& soln,
                               const Epetra_Vector_Ghosted* const atolVector = 0);


    //! Evaluates the atol DAESystemInitial vector used in the time stepper routine and in the nonlinear solver.
    /*!
     *  (virtual from ProblemResidEval)
     *
     *   Note this is an important routine for the speed of the solution.
     *    The default behavior is to set all of the atol values to the constant atolDAEInitDefault.
     *
     *  @param[in]           atolDAEInitDefault  Double containing the default value of atol used in the DAESystemInitial solve
     *  @param[in]           soln                Current solution vector.
     *  @param[in]           solnDot             Current solution vector.
     *  @param[in]           atolVector_DAEInit  Optional vector containing the individual entries for the time derivative
     *                                           absolute tolerance. Use this if you want to override the default calculation of atol.
     *                                           Defaults to nullptr.
     */
    virtual void setAtolVector_DAEInit(double atolDAEInitDefault, const Epetra_Vector_Ghosted& soln,
                                       const Epetra_Vector_Ghosted& solnDot,
                                       const Epetra_Vector_Ghosted* const atolVector_DAEInit = nullptr) const;

    //! Evaluates the atol vector used in the delta damping process.
    /*!
     *  (virtual from ProblemResidEval)
     *
     *   Note this is an important routine for the speed of the solution.
     *  @param[in]           relcoeff            Relative constant to multiply all terms by
     *  @param[in]           soln                current solution vector.
     *  @param[in]           atolDeltaDamping    If non-zero, this copies the vector into the object as input.
     *                                           The default is nullptr.
     */
    virtual void setAtolDeltaDamping(double relcoeff, const Epetra_Vector_Ghosted& soln,
                                     const Epetra_Vector_Ghosted* const atolDeltaDamping = nullptr);

    //! Evaluates the atol vector used in the delta damping process for the DAE problem
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  Note this is an important routine for the speed of the solution.
     *  @param[in]           relcoeff            Relative constant to multiply all terms by
     *  @param[in]           soln                current solution vector.
     *  @param[in]           solnDot             Current solutionDot vector.
     *  @param[in]           atolDeltaDamping    If non-zero, this copies the vector into the object as input
     *                                           The default is nullptr.
     */
    virtual void setAtolDeltaDamping_DAEInit(double relcoeff, const Epetra_Vector_Ghosted& soln,
                                            const Epetra_Vector_Ghosted& solnDot,
                                            const Epetra_Vector_Ghosted* const atolDeltaDamping = nullptr);

    //! Evaluate the delta damping vector from the abs tol vector
    /*!
     *  The default behavior is to use the atolVector values to fill this vector.
     *
     *  param relCoeff sets a multiplicative constant to the atol vector
     */
    //virtual void
    //setDeltaDamping(double relCoeff = 1.0);

    //! Return a const reference to the atol Epetra_Vector
    /*!
     *  This is an owned only vector.
     *
     *  @return                                  const reference to the atol vector
     */
    const Epetra_Vector_Owned& atolVector() const;

    //! Return a const reference to the atolDAEInit vector
    /*!
     *  @return const reference to the atolDAEInit vector
     */
    const Epetra_Vector_Owned& atolVector_DAEInit() const;

    //! Return a const reference to the atol_deltaDamping vector
    /*!
     *  @return const reference to the atol delta damping vector
     */
    const Epetra_Vector_Owned& atolDeltaDamping() const;

    //! Filter the solution predictions
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  @param[in]           time_current        Current time
     *  @param[in]           y_n                 Current value of the solution vector
     */
    virtual void filterSolnPrediction(double time_current, Epetra_Vector_Ghosted& y_n);

    //! Filter the new solution estimate, getting rid of forbidden values that may have been inserted into the solution during
    //! the nonlinear iteration process
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  This is where we impose the requirements of a positive or zero value for species concentrations and other bounds conditions
     *  required by the solution process
     *
     *  This base class implementation just sets delta_y_filter to zero.
     *
     *  @param[in]           t_n                 Current time
     *  @param[in]           delta_t_n           Current value of delta_t used on this global time step
     *  @param[in]           yCurrent_n          Current proposed value of the solution vector at t_n
     *  @param[in]           ydotCurrent_n       Current proposed value of the solution dot vector at t_n
     *  @param[out]          delta_y_filter      Change in yCurrent_n dictated by the filter to get the solution back into
     *                                           compliance with the bounds.
     */
    virtual void applyFilter(const double t_n, const double delta_t_n, const Epetra_Vector_Ghosted& yCurrent_n,
                             const Epetra_Vector_Ghosted& ydotCurrent_n, Epetra_Vector_Ghosted& delta_y_filter);

    //! Apply any constraints on time region bounds to provide a constraint on the next delta_t 
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  If there are any bounds on going beyound a time period, we apply it here.
     *  @param[in]           t_n                 Current time
     *  @param[in]           y_n                 Value of the solution vector at t_n
     *  @param[in]           ydot_n              Value of the solution dot vector at t_n
     *
     *  @return                                  Returns the constraint on the value of delta_t for the next time step.
     *                                           Returns 0.0, if there is no time constraint.
     */
    virtual double delta_t_constraint(const double t_n, const Epetra_Vector_Ghosted& y_n, const Epetra_Vector_Ghosted& ydot_n);

    //! Evaluate a supplemental set of equations that are not part of the solution vector, but are considered to be time dependent
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  Equations in this system are evaluated using the time discretization scheme of the nonlinear solver.
     *  It can be used to track supplemental quantities in the calculation, especially if they need to be integrated in time.
     *
     *  An example of this may be total flux quantites that are dumped into a phase.
     *
     *  This routine is called at the beginning of the time stepping, in order to set up any initializations,
     *  and it is called after every successful time step, once.
     *
     *  @param[in]           ifunc               Function call value. Why the routine was called.
     *                                             0 Initial call to the function, done whenever the time stepper is entered
     *                                             1 Called after every successful time step.
     *  @param[in]           t                   Current time
     *  @param[in]           deltaT              Current value of deltaT
     *  @param[in]           y                   Current value of the solution vectors
     *  @param[in]           ydot                Current value of time derivative of the solution vectors.
     */
    virtual void evalTimeTrackingEqns(const int ifunc, const double t, const double deltaT, const Epetra_Vector_Ghosted& y,
                                      const Epetra_Vector_Ghosted* const ydot);

    //! Set a solution parameter
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  @param[in]           paramName           String identifying the parameter to be set
     *  @param[in]           paramVal            Single double value of the parameter to be set
     *
     *  @return                                  Returns a 0 if the number makes sense.
     *                                           Returns a negative number if the parameter is unknown.
     */
    virtual int setSolutionParam(std::string paramName, double paramVal);

    //! Get a solution parameter
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  @param[in]           paramName           String identifying the parameter to be set
     *  @param[out]          paramVal            Vector of parameters returned
     *
     *  @return                                  Returns the number of parameters returned.
     */
    virtual int getSolutionParam(std::string paramName, double* const paramVal);

    //! Evaluate the stopping criteria
    /*!
     *  (virtual from ProblemResidEval)
     *  If there is a stopping criteria for the simulation other than time, set it here
     *
     *  @param[out]          time_current        Current time
     *  @param[out]          delta_t_n           Current delta t
     *  @param[in]           y_n                 Const reference to the current solution vector
     *  @param[in]           ydot_n              Const reference to the current solution dot vector
     *
     *  @return                                  Returns true if the simulation should be stopped after the current global step
     */
    virtual bool evalStoppingCritera(double& time_current, double& delta_t_n, const Epetra_Vector_Ghosted& y_n,
                                     const Epetra_Vector_Ghosted& ydot_n);

    //! Return a vector of delta y's for calculation of the numerical Jacobian
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  Routine to find the delta's for all the unknowns at the same time. 
     *  Note that all input vectors to this routine are ghosted. THIS MUST BE CHECKED.
     *  The idea is to calculate all deltas on the ghosted unknowns directly, and have them be exactly the same as on their
     *  owned processor.
     *
     *  @param[in]           t                   time
     *  @param[in]           soln                Reference to the current solution vector.
     *  @param[in]           solnDot_ptr         Pointer to the Current solutionDot vector. Can be 0 for steady state problems
     *  @param[out]          deltaSoln           Refererence to the changeable delta for the solution vector
     *                                           that will be used.
     *  @param[in]           solveType           Type of the problem being solved expressed as a  Solve_Type_Enum. 
     *                                           Defaults to TimeDependentAccurate_Solve
     *  @param[in]           solnWeights         Vector of solution weights used in creating normalized error values.
     *                                           Defaults to 0, indicating none is available.
     */
    virtual void calcDeltaSolnVariables(double t, const Epetra_Vector_Ghosted& soln,
                                        const Epetra_Vector_Ghosted* const solnDot_ptr, Epetra_Vector_Ghosted& deltaSoln,
                                        const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve,
                                        const Epetra_Vector_Ghosted* const solnWeights = nullptr);

    //! Save the solution to the end of an XML file using XML solution format
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  We write out the solution to a file.
     *
     *  @param[in]           ievent              Type of the event. The following form is used:
     *                                              0 Initial conditions
     *                                              1 Completion of a successful intermediate step.
     *                                              2 Final successful conditions.
     *                                              3 Intermediate nonlinear step
     *                                             -1 unsuccessful step
     *
     *  @param[in]           baseFileName        Filename to be used. .xml is appended onto the filename
     *                                           processors other than 0 have the pid appended to the name
     *  @param[in]           m_y_n               Current value of the solution vector
     *  @param[in]           m_ydot_n            Current value of the derivative of the solution vector
     *  @param[in]           t                   Current time
     *  @param[in]           delta_t             delta t
     *  @param[in]           delta_t_np1         Suggested next delta t value (defaults to 0.0 if there isn't a
     *                                           good guess for the next delta_t).
     */
    virtual void
    saveSolutionEnd(const int ievent, std::string baseFileName, const Epetra_Vector_Ghosted& m_y_n,
                    const Epetra_Vector_Ghosted* const m_ydot_n, const double t,
                    const double delta_t, const double delta_t_np1 = 0.0);

    //! Read a solution from a saved file, using the record number to identify the solution
    /*!
     *  Note, we read only successful final steps.
     *
     *  @param[in]           iNumber             Solution Record number to read from. Solutions are numbered consequetively from
     *                                           zero to the number of global time steps. We will use the last
     *                                           time step if there aren't as many time steps in the solution file as requested.
     *
     *  @param[in]           baseFileName        File name to be used. .xml is appended onto the filename.
     *                                           Processors other than 0 have the pid appended to the name as well.
     *
     *  @param[out]          y_n_ghosted         Current value of the solution vector
     *  @param[out]          ydot_n_ghosted      Current value of the derivative of the solution vector
     *  @param[out]          t_read              time that is read in
     *  @param[out]          delta_t_read        delta time step for the last time step.
     *  @param[out]          delta_t_next_read   delta time step for the next time step if available
     */
    void readSolutionRecordNumber(const int iNumber, std::string baseFileName, Epetra_Vector_Ghosted& y_n_ghosted,
                                  Epetra_Vector_Ghosted* const ydot_n_ghosted, double& t_read, double& delta_t_read, 
                                  double& delta_t_next_read);

    //! Read the solution from a saved file. (deprecated)
    /*!
     *  We read only successful final steps
     *
     *  @param[in]           iNumber             Solution Record number to read from. Solutions are numbered consequetively from
     *                                           zero to the number of global time steps. We will use the last
     *                                           time step if there aren't as many time steps in the solution file as requested.
     *
     *  @param[in]           baseFileName        File name to be used. .xml is appended onto the filename.
     *                                           Processors other than 0 have the pid appended to the name as well.
     *
     *  @param[out]          y_n_ghosted         Current value of the solution vector
     *  @param[out]          ydot_n_ghosted      Current value of the derivative of the solution vector
     *  @param[out]          t_read              time that is read in
     *  @param[out]          delta_t_read        delta time step for the last time step.
     *  @param[out]          delta_t_next_read   delta time step for the next time step if available
     */
    void readSolution(const int iNumber, std::string baseFileName, Epetra_Vector_Ghosted& y_n_ghosted, 
                      Epetra_Vector_Ghosted* const ydot_n_ghosted, double& t_read, double& delta_t_read,
                      double& delta_t_next_read);

    //! Read the solution from a saved solution XML record.
    /*!
     *  This is the underlying program that reads the record
     *
     *  @param[in]           simulRecord         XML element to read from. Name must be a simulation XML_Node.
     *  @param[out]          y_n_ghosted         Current value of the solution vector
     *  @param[out]          ydot_n_ghosted      Current value of the derivative of the solution vector
     *  @param[out]          t_read              time that is read in
     *  @param[out]          delta_t_read        delta time step for the last time step.
     *  @param[out]          delta_t_next_read   delta time step for the next time step if available
     */
    void readSolutionXML(ZZCantera::XML_Node* simulRecord, Epetra_Vector_Ghosted& y_n_ghosted,
                         Epetra_Vector_Ghosted* const ydot_n_ghosted, double& t_read,
                         double& delta_t_read, double& delta_t_next_read);

    //!  Select the global time step increment record by the consequatively numbered record index number
    /*!
     *   @param[in]          xSoln               Solution file for the simulation object
     *   @param[in]          globalTimeStepNum   Time step number to select
     *
     *   @return                                 Returns a pointer to the selected record.
     */
    ZZCantera::XML_Node* selectSolutionRecordNumber(ZZCantera::XML_Node* xSoln, int globalTimeStepNum);

    //! Select the global time step increment record by the consequatively numbered record index number
    /*!
     *  @param[in]           xSoln               Solution file for the simulation object
     *  @param[in]           solnTimeStepID      String value of the time step ID to select
     *
     *  @return                                  Returns a pointer to the selected record.
     */
    ZZCantera::XML_Node* selectSolutionTimeStepID(ZZCantera::XML_Node* xSoln, std::string solnTimeStepID);

    //! Write the solution to either the screen or to a log file
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  This is a general output utility to Zuzax's logfile.  It's not hooked into the IO algorithm at all. It should be
     *  conditionally called depending on the whims of the user.
     *
     *  @param[in]           ievent              Type of the event. The following form is used:
     *                                              0  Initial conditions
     *                                              1  Completion of a successful intermediate step.
     *                                              2  Final successful conditions.
     *                                              3  Intermediate nonlinear step
     *                                             -1  unsuccessful step
     *  @param[in]         doTimeDependentResid  Do the time dependent residual calculation
     *  @param[in]           t                   Current time
     *  @param[in]           delta_t             delta t
     *  @param[in]           y_n                 Current value of the solution vector
     *  @param[in]           ydot_n              Current value of the derivative of the solution vector
     *  @param[in]           solveType           Type of the problem being solved expressed as a  Solve_Type_Enum. 
     *                                           Defaults to TimeDependentAccurate_Solve
     *  @param[in]           delta_t_np1         Value of the next time step. Defaults to 0.0, if not available.
     *
     */
    virtual void showProblemSolution(const int ievent, bool doTimeDependentResid, const double t,
                                     const double delta_t, const Epetra_Vector_Owned& y_n,
                                     const Epetra_Vector_Owned* const ydot_n,
                                     const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve,
                                     const double delta_t_np1 = 0.0);

    //! Write out Tecplot solution files
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  @param[in]           ievent              Type of the event. The following form is used:
     *                                              0  Initial conditions
     *                                              1  Completion of a successful intermediate step.
     *                                              2  Final successful conditions.
     *                                              3  Intermediate nonlinear step
     *                                             -1  unsuccessful step
     *  @param[in]           baseFileName        Base file name for writing out the file.
     *  @param[in]         doTimeDependentResid  Do the time dependent residual calculation
     *  @param[in]           t                   Current time
     *  @param[in]           delta_t             delta t
     *  @param[in]           y_n                 Current value of the solution vector
     *  @param[in]           ydot_n              Current value of the derivative of the solution vector
     *  @param[in]           solveType           Type of the problem being solved expressed as a  Solve_Type_Enum. 
     *                                           Defaults to TimeDependentAccurate_Solve
     *  @param[in]           delta_t_np1         Value of the next time step. Defaults to 0.0, if not available.
     *
     *
     */
    virtual void writeTecplot(const int ievent, std::string baseFileName, bool doTimeDependentResid,
                              const double t, const double delta_t, const Epetra_Vector_Ghosted& y_n,
                              const Epetra_Vector_Ghosted* const ydot_n, const Zuzax::Solve_Type solveType,
                              const double delta_t_np1);

    //! Write a solution vector type to either the screen or a log file
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  @param[in]           solnVecName         Name of the solution vector -> will appear in the output
     *  @param[in]           t                   Current time
     *  @param[in]           delta_t             delta t
     *  @param[in]           solnVector          Vector of information to be printed
     *  @param[in]           outFile             FILE pointer. Defaults to stdout if not specified.
     */
    virtual void showSolutionVector(std::string& solnVecName, const double t, const double delta_t,
                                    const Epetra_Vector_Owned& solnVector, FILE* outFile = stdout);

    //! Write an int solution vector type to either the screen or a log file
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  @param[in]           solnVecName         Name of the solution vector -> will appear in the output
     *  @param[in]           t                   Current time
     *  @param[in]           delta_t             delta t
     *  @param[in]           solnVector          Integer Vector of information to be printed
     *  @param[in]           outFile             FILE pointer. Defaults to stdout if not specified.
     */
    virtual void showSolutionIntVector(std::string& solnVecName, const double t, const double delta_t,
                                       const Epetra_IntVector& solnVector, FILE* outFile = stdout);

    //! Write out to a file or to standard output the current solution
    /*!
     *  (virtual from ProblemResidEval)
     *
     *   These functions are affected by the print controls of the nonlinear solver
     *   and the time stepper.
     *
     *   ievent is a description of the event that caused this function to be called.
     *
     *  @param[in]           ievent              Event that's causing this routine to be called.
     *                                            =  0 Initial conditions for a calculation
     *                                            =  1 Completion of a successful intermediate time step.
     *                                            =  2 Completion of a successful Final time or final calculation.
     *                                            =  3 Completion of a successful Intermediate nonlinear step
     *                                            = -1 unsuccessful time step that converged, but failed otherwise
     *                                            = -2 unsuccessful nonlinear step.
     *  @param[in]    doTimeDependentResid       true if solving a time dependent problem
     *  @param[in]           time_current        Current time
     *  @param[in]           delta_t_n           Current value of delta_t
     *  @param[in]           istep               Current step number
     *  @param[in]           soln_n              Current value of the solution vector
     *  @param[in]           solnDot_n_ptr       Current value of the time deriv of the solution vector
     *  @param[in]           solveType           Type of the problem being solved expressed as a  Solve_Type_Enum. 
     *                                           Defaults to TimeDependentAccurate_Solve
     *  @param[in]           delta_t_np1         Suggested next delta t value (defaults to 0.0 if there isn't a
     *                                           good guess for the next delta_t).
     */
    virtual void writeSolution(const int ievent, const bool doTimeDependentResid, const double time_current,
                               const double delta_t_n, const int istep, const Epetra_Vector_Ghosted& soln_n,
                               const Epetra_Vector_Ghosted* const solnDot_n_ptr,
                               const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve,
                               const double delta_t_np1 = 0.0);

    //! This function may be used to create output at various points in the execution of an application.
    /*!
     *  (virtual from ProblemResidEval)
     *
     *   These functions are not affected by the print controls of the nonlinear solver
     *   and the time stepper.
     *
     *  ievent is a description of the event that caused this   function to be called.
     *
     *  @param[in]           ievent              Event that's causing this routine to be called.
     *                                            =  0 Initial conditions for a calculation
     *                                            =  1 Completion of a successful intermediate time step.
     *                                            =  2 Completion of a successful Final time or final calculation.
     *                                            =  3 Completion of a successful Intermediate nonlinear step
     *                                            = -1 unsuccessful time step that converged, but failed otherwise
     *                                            = -2 unsuccessful nonlinear step.
     *
     *  @param[in]           time_current        Current time
     *  @param[in]           delta_t_n           Current value of delta_t
     *  @param[in]           istep               Current step number
     *  @param[in]           y_n                 Current value of the solution vector
     *  @param[in]           ydot_n_ptr          Current value of the time deriv of the solution vector
     */
    virtual void user_out(const int ievent, const double time_current, const double delta_t_n,
                          const int istep, const Epetra_Vector_Ghosted& y_n,
                          const Epetra_Vector_Ghosted* const ydot_n_ptr);

    //! Handle matrix conditioning if necessary here
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  In general, multiply the matrix by the inverse of a matrix which lead to a better conditioned system. 
     *  The default, specified here, is to do nothing.
     *
     *  @param[in,out]       matrix_ptr          Pointer to the matrix 
     *  @param[in]           nrows               Number of rows in the matrix
     *  @param[in,out]       rhs                 Right-hand side of the linear system, for which we are trying to find
     *                                           the solution for
     *  @return                                  Returns a flag to indicate that the operation is successful.
     *                                             - 1  Means a successful operation
     *                                             - 0  or neg value means an unsuccessful operation
     */
    virtual int matrixConditioning(double* const matrix_ptr, int nrows, double* const rhs);

    //! Returns the base file name of the solution file
    /*!
     *  The default name is "solution".
     *  @return                                  Returns the base file name
     */
    std::string getBaseFileName() const;

    //! Sets the base file name of the solution file
    /*!
     *  The default name is "solution". This is used unless it is overwritten here.
     *
     *  @param[in]           baseFileName        Sets the value of the base file name.
     */
    void setBaseFileName(const std::string& baseFileName);

    //! Calculate the solution at the last time step, t_nm1
    /*!
     *  Calculates the old solution using the current solution, and stores it internally in the vector *solnOld_ptr_
     *  Assumes backwards euler time stepping.
     *
     *  @param[in]          soln                 Solution at the current time step
     *  @param[in]          solnDot              SolutionDot at the current time step
     *  @param[in]          rdelta_t             Inverse of the time step
     */
    void calcSolnOld(const Epetra_Vector_Ghosted& soln, const Epetra_Vector_Ghosted& solnDot, double rdelta_t);

public:

    // --------------------------------------------- D A T A ----------------------------------------------------------------

    //! Default value of the atol parameter
    double m_atolDefault;

    //! Global total of the equations
    int m_neq;

    //! Number of equations on the processor including ghost equations
    int m_NumLcEqns;

    //! Number of equations on the processor not including ghost equations
    int m_NumLcOwnedEqns;

    //! DomainLayout is a light class that describes the overall
    //! disposition of the domains in the problem
    DomainLayout* DL_ptr_;

    //! Global indices class is the same on all processors
    /*!
     *  This is a global structure, containing global information.
     *  All information in this structure, except for my_procID, numLnNodes,
     *  and numLnEqns, is the same on all processors.
     *
     */
    GlobalIndices* GI_ptr_;

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
    LocalNodeIndices* LI_ptr_;

    //! Pointer to the jacobian calculator
    /*!
     *   We don't own this
     */
    EpetraJac* m_jac;

    //! Base file name for the output file
    /*!
     *  The default is to use the name, "solution"
     */
    std::string m_baseFileName;

    //! Step Number
    /*!
     *   Step number is the time step number for time dependent calculations.
     *   It's the number of the current solution, for time independent calculations,
     *   which have continuation numbers.
     */
    int m_StepNumber;

    //! Pointer to the old solution vector
    Epetra_Vector_Ghosted* solnOld_ptr_;

    //! Residual used for internal calculation
    Epetra_Vector_Owned* resInternal_ptr_;

    //! Absolute tolerances used for nonlinear solver convergences
    Epetra_Vector_Ghosted* m_atolVector;

    //! Absolute tolerances used for nonlinear solver convergences
    mutable Epetra_Vector_Ghosted* m_atolVector_DAEInit;

    //! Tolerances used for delta damping in the nonlinear solver
    Epetra_Vector_Ghosted* m_atolDeltaDamping;

    //! Tolerances used for delta damping in the nonlinear solver
    Epetra_Vector_Ghosted* m_atolDeltaDamping_DAEInit;

    //! Pointer to the ProblemStatement class
    ProblemStatement* psInput_ptr_;

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

    //! Number of of residual calculations used for delta calculations to create the jacobian.
    int counterJacDeltaCalcs_;

    //! Counter containing the number of residual calculation used for ShowSolution calcs
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
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

#endif

