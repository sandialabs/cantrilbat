/*
 * @file m1d_Domain1D.h
 *        Base class for calculation of residuals from a single domain
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_DOMAIN1D_H_
#define M1D_DOMAIN1D_H_
//! This is a heavyweight base class that provides the function
//! evaluation for a single domain whether a surface or a bulk.

#include "m1d_ProblemResidEval.h"
#include "m1d_LocalNodeIndices.h"

//#include "Epetra_Vector.h"

#include "zuzax/base/xml.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

//==================================================================================================================================
//! Class to accept global heat balance numbers from internal calculations
/*!
 *  We use this class to calculate a macroscopic heat balance. This is a virtual class so that it may be inherited
 *  and specialized for situations
 *
 *  All elements of this class are on a per unit area basis for cartesian coordinate systems. In other words the units
 *  are Joules / m2.
 */
class globalHeatBalVals
{
public:

    //! Default constructor
    globalHeatBalVals() :
        totalHeatCapacity(0.0),
        HeatFluxRight(0.0),
        HeatFluxLeft(0.0),
        oldNEnthalpy(0.0),
        newNEnthalpy(0.0),
        moleFluxOut(0.0),
        enthFluxOut(0.0)
    {
    }

    //! Empty destructor
    virtual ~globalHeatBalVals()
    {
    }

    //! Member function to zero all
    virtual void zero()
    {
        totalHeatCapacity = 0.0;
        HeatFluxRight = 0.0;
        HeatFluxLeft = 0.0;
        oldNEnthalpy = 0.0;
        newNEnthalpy = 0.0;
    }

    //! Total Extrinsic heat capacity per unit area of the
    /*!
     *   Units: Joules/m2/K
     */
    double totalHeatCapacity;

    //! Heat flox out of the domain to the right
    /*!
     *   Units: Joules/m2/s
     */
    double HeatFluxRight;

    //! Heat flox out of the domain to the left
    /*!
     *   Units: Joules/m2/s
     */
    double HeatFluxLeft;

    //! Old extrinsic enthalpy
    /*!
     *  Units: Joules/m2
     */
    double oldNEnthalpy;

    //! New extrinsic enthalpy
    /*!
     *  Units: Joules/m2
     */
    double newNEnthalpy;

    //! Mole flux out of the domain
    /*!
     *  Units: kmol/m2/s
     */
    double moleFluxOut;

    //! Vector of species flux out of the domain
    /*!
     *  Units: kmol/m2/s
     */
    std::vector<double> speciesMoleFluxOut;

    //! Enthalpy flux out of the domain
    /*!
     *  Units: Joules/m2/s
     */
    double enthFluxOut;
};
//==================================================================================================================================

//! Base class for solving residuals for bulk and surface domains
/*!
 * There is a 1 to 1 mapping between the local control volume indexing and
 * the Local Consecutive Ordering indexing
 *
 * There is a 1 to 1 mapping between the global control volume indexing
 * and the Global node number indexing that is given by a single offset.
 *
 */
class Domain1D
{

public:

    //! Default Constructor
    /*!
     */
    Domain1D();

    //! Copy Constructor
    /*!
     *  @param[in]             r                   Object to be copied
     */
    Domain1D(const Domain1D& r);

    //!  Destructor
    virtual ~Domain1D();

    //! Assignment operator
    /*!
     *  @param[in]             r                   object to be copied
     *
     *  @return                                    Returns a reference to an object
     */
    Domain1D& operator=(const Domain1D& r);

    //! Specify an identifying tag for this domain.
    /*!
     *  (virtual from Domain1D)
     *  The identifying tag for the domain is the name of the domain. It will appear on all output files.
     *
     *  @param[in]             s                   Name of the domain
     */
    virtual void setID(const std::string& s);

    //! Returns the identifying tag for the domain
    /*!
     *  (virtual from Domain1D)
     *  @return                                    Returns the string name for the domain
     */
    virtual std::string id() const;

    //! Prepare all of the indices for fast calculation of the residual
    /*!
     *  (virtual from Domain1D)
     *  Ok, at this point, we will have figured out the number of equations
     *  to be calculated at each node point. The object NodalVars will have been fully formed.
     *
     *  We use domain_prep() to figure out what the local node number is and what equations correspond to what unknown.
     *
     *  Child objects of domain_prep() will normally call parent classes in a recursive fashion.
     *
     *  @param[in]             li_ptr              Pointer to the LocalNodeIndices Structure that contains information
     *                                             about how the mesh is layed out within this domain and other domains in the problem
     */
    virtual void domain_prep(LocalNodeIndices* const li_ptr);

    //! Basic function to calculate the residual for the current domain.
    /*!
     *  (virtual from Domain1D)
     *  This base class is used just for volumetric domains.
     *
     *  All residual terms are written with the following sign convention
     *  based on keeping the time derivative term positive.
     *
     *       res = dcdt - dc2 /dx2 - src = 0
     *
     *  @param[out]            res                 Output vector containing the residual
     *  @param[in]             doTimeDependentResid  boolean indicating whether the time dependent residual is requested
     *  @param[in]             soln_ptr            Solution vector at which the residual should be evaluated
     *  @param[in]             solnDot_ptr         Solution dot vector at which the residual should be evaluated.
     *  @param[in]             solnOld_ptr         Pointer to the solution vector at the old time step
     *  @param[in]             t                   time
     *  @param[in]             rdelta_t            inverse of delta_t
     *  @param[in]             residType           Type of evaluation of the residual. Uses the ResidEval_Type_Enum type.
     *                                             Defaults to Base_ResidEval
     *  @param[in]             solveType           Type of solution Type. Uses the Solve_Type_Enum  type.
     *                                             Defualts to  TimeDependentAccurate_Solve
     */
    virtual void
    residEval(Epetra_Vector& res, const bool doTimeDependentResid,
              const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr,
              const Epetra_Vector* const solnOld_ptr, const double t, const double rdelta_t,
              const Zuzax::ResidEval_Type residType = Zuzax::ResidEval_Type::Base_ResidEval,
              const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve);


    //! Evaluate Post-solution quantities of interest after each converged time step
    /*!
     *  (virtual from Domain1D)
     *  @param[in]             doTimeDependentResid  Boolean indicating whether the time dependent residual is requested
     *  @param[in]             soln_ptr            Solution vector at which the residual should be evaluated
     *  @param[in]             solnDot_ptr         Solution dot vector at which the residual should be evaluated.
     *  @param[in]             solnOld_ptr         Pointer to the solution vector at the old time step
     *  @param[in]             t                   time
     *  @param[in]             rdelta_t            inverse of delta_t
     */
    virtual void
    eval_PostSoln(const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                  const Epetra_Vector* const solnDot_ptr,
                  const Epetra_Vector* const solnOld_ptr, const double t, const double rdelta_t);

    //! Evaluate the macroscopic heat balance on the domain given a converged solution of the problem
    /*!
     *  (virtual from Domain1D)
     *  Calculations are done on a per m2 basis. So, the basic units are Joules / m2.
     *
     *  @param[in]             ifunc               situation function parameter, input from doHeatAnalysis()
     *  @param[in]             t                   time
     *  @param[in]             deltaT              deltaT for the just-finished time step
     *  @param[in]             soln_ptr            Solution vector at which the residual should be evaluated
     *  @param[in]             solnDot_ptr         Solution dot vector at which the residual should be evaluated.
     *  @param[in]             solnOld_ptr         Pointer to the solution vector at the old time step
     *  @param[in,out]         dVals               Reference to the globalHeatBalVals structure that will contain the results
     *                                             of the macroscopic heat balance calculation
     */
    virtual void
    eval_HeatBalance(const int ifunc, const double t, const double deltaT, const Epetra_Vector* const soln_ptr,
                     const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr,
                     struct globalHeatBalVals& dVals);

    //! Evalualte the macroscopic Species balance equations on a domain given a converged solution of the problem
    /*!
     *  (virtual from Domain1D)
     *  Calculations are done on a per m2 basis. So, the basic units are kmol/ m2.
     *
     *  @param[in]             ifunc               situation function parameter, input from doHeatAnalysis()
     *  @param[in]             t                   time
     *  @param[in]             deltaT              deltaT for the just-finished time step
     *  @param[in]             soln_ptr            Solution vector at which the residual should be evaluated
     *  @param[in]             solnDot_ptr         Solution dot vector at which the residual should be evaluated.
     *  @param[in]             solnOld_ptr         Pointer to the solution vector at the old time step
     *  @param[in,out]         dVals               Reference to the globalHeatBalVals structure that will contain the results
     *                                             of the macroscopic heat balance and species balance calculations
     */
    virtual void
    eval_SpeciesElemBalance(const int ifunc, const double t, const double deltaT,
                            const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr,
                            const Epetra_Vector* const solnOld_ptr, class globalHeatBalVals& dVals);

    //! Utility function to calculate quantities before the main residual routine.
    /*!
     *  (virtual from Domain1D)
     *  This is used for a loop over nodes in the domain. All calculated quantities must be internally storred within the
     *  domain structure. Currently this is called during the residual evalultion of the problem.
     *
     *  @param[in]             doTimeDependentResid  Boolean indicating whether the time dependent residual is requested
     *  @param[in]             soln_ptr            Solution vector at which the residual should be evaluated
     *  @param[in]             solnDot_ptr         Solution dot vector at which the residual should be evaluated.
     *  @param[in]             solnOld_ptr         Pointer to the solution vector at the old time step
     *  @param[in]             t                   time
     *  @param[in]             rdelta_t            inverse of delta_t
     *  @param[in]             residType           Type of evaluation of the residual. Uses the ResidEval_Type_Enum type.
     *                                             Defaults to Base_ResidEval
     *  @param[in]             solveType           Type of solution Type. Uses the Solve_Type_Enum  type.
     *                                             Defaults to  TimeDependentAccurate_Solve
     */
    virtual void
    residEval_PreCalc(const bool doTimeDependentResid, const Epetra_Vector* soln_ptr, const Epetra_Vector* solnDot_ptr,
                      const Epetra_Vector* solnOld_ptr, const double t, const double rdelta_t,
                      const Zuzax::ResidEval_Type residType, const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve);

    //! Auxiliary function to calculate the residual for the current domain.
    /*!
     *  (virtual from Domain1D)
     *  By default this function does nothing.
     *
     *  @param[out]          res                 Output vector containing the residual
     *  @param[in]           doTimeDependentResid boolean indicating whether the time dependent residual is requested
     *  @param[in]           soln_ptr            solution vector at which the residual should be evaluated
     *  @param[in]           solnDot_ptr         solution dot vector at which the residual should be evaluated.
     *  @param[in]           solnOld_ptr         Pointer to the solution vector at the old time step
     *  @param[in]           t                   time
     *  @param[in]           rdelta_t            inverse of delta_t
     *  @param[in]           residType           Type of evaluation of the residual. Uses the ResidEval_Type_Enum type.
     *                                             Defaults to Base_ResidEval
     *  @param[in]           solveType           Type of solution Type. Uses the Solve_Type  type.
     *                                             Defaults to  TimeDependentAccurate_Solve
     */
    virtual void
    residEval_PostCalc(Epetra_Vector& res, const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                       const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr, const double t,
                       const double rdelta_t, const Zuzax::ResidEval_Type residType = Zuzax::ResidEval_Type::Base_ResidEval,
                       const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve);

    //! Function that gets called at end the start of every time step
    /*!
     *  (virtual from Domain1D)
     *
     *  This function provides a hook for a residual that gets called whenever a
     *  time step has been accepted and we are about to move on to the next time step.
     *  The call is made with the current time as the time
     *  that is accepted. The old time may be obtained from t and rdelta_t_accepted.
     *
     *  After this call interrogation, of the previous time step's results will not be valid.
     *
     *  This call also calculates all of the "old" cell information for the residual calculation.
     *  The "old" values are storred from calculation of the "current" values.
     *
     *  Note, when t is equal to t_old, soln_ptr should equal solnOld_ptr values. However,
     *  solnDot_ptr values may not be zero.
     *
     *  @param[in]           doTimeDependentResid  This is true if we are solving a time dependent problem.
     *  @param[in]           soln_ptr            Solution value at the current time
     *  @param[in]           solnDot_ptr         derivative of the solution at the current time.
     *  @param[in]           solnOld_ptr         Solution value at the old time step, n-1
     *  @param[in]           t                   Current time to be accepted, n
     *  @param[in]           t_old               Previous time step value, t_old may be equal to t,
     *                                           When we are calculating the initial conditions we
     *                                           require that we have values of "old" cell information.
     *                                           The call to this routine calculates the "old" information.
     */
    virtual void
    advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                        const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr,
                        const double t, const double t_old);

    //! Revert the Residual object's conditions to the conditions at the start of the global time step
    /*!
     *  (virtual from Domain1D)
     *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init_init.
     *  We get rid of the pendingIntegratedFlags_ flag here as well.
     */
    virtual void revertToInitialGlobalTime();

    //! Base class for saving the solution on the domain in an xml node.
    /*!
     *  (virtual from Domain1D)
     *
     *  @param[in,out]         oNode               Reference to the XML_Node
     *  @param[in]             soln_GlAll_ptr      Pointer to the Global-All solution vector
     *  @param[in]             solnDot_GlAll_ptr   Pointer to the time derivative of the Global-All solution vector
     *  @param[in]             t                   time
     *  @param[in]             duplicateOnAllProcs If this is true, all processors will include the same XML_Node information
     *                                             as proc 0. If false, the xml_node info will only exist on proc 0.
     *                                             Defaults to false.
     */
    virtual void saveDomain(Zuzax::XML_Node& oNode, const Epetra_Vector* const soln_GlAll_ptr,
                            const Epetra_Vector* const solnDot_GlAll_ptr, const double t, bool duplicateOnAllProcs = false);

    //! Base Class for reading the solution from the saved file
    /*!
     *  (virtual from Domain1D)
     *
     *  @param[in]             simulationNode      Reference to the XML_Node named simulation
     *  @param[out]            soln_GlAll_ptr      Pointer to the Global-All solution vector that will accept the solution
     *  @param[out]            solnDot_GlAll_ptr   Pointer to the time derivative of the Global-All solution vector that will accept
     *                                             the solution
     */
    virtual void readSimulation(const Zuzax::XML_Node& simulationNode, Epetra_Vector* const soln_GlAll_ptr,
                                Epetra_Vector* const solnDot_GlAll_ptr);

    //! Base Class for reading the solution from the saved file
    /*!
     *  (virtual from Domain1D)
     *
     * @param[in]            domainNode          Reference to the XML_Node to read the solution from
     * @param[in]            soln_GlAll_ptr      Pointer to the Global-All solution vector
     * @param[in]            solnDot_GlAll_ptr   Pointer to the time derivative of the Global-All solution vector
     * @param[in]            globalTimeRead      Value of the global time that is read in. This is used for
     *                                           comparison and quality control purposes
     */
    virtual void
    readDomain(const Zuzax::XML_Node& domainNode, Epetra_Vector* const soln_GlAll_ptr,
               Epetra_Vector* const solnDot_GlAll_ptr, double globalTimeRead);

    //! Method for writing the header for the surface domain to a tecplot file.
    /*!
     *  (virtual from Domain1D)
     *
     *  Only proc0 will write tecplot files.
     */
    virtual void writeSolutionTecplotHeader();

    //! Method for writing the solution on the surface domain to a tecplot file.
    /*!
     *  (virtual from Domain1D)
     *
     *  Only proc0 will write tecplot files.
     *
     *  @param[in]           soln_GlAll_ptr      Pointer to the Global-All solution vector
     *  @param[in]           solnDot_GlAll_ptr   Pointer to the time derivative of the Global-All solution vector
     *  @param[in]           t                   time
     */
    virtual void writeSolutionTecplot(const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
                                      const double t);

    //! Base class for writing the solution on the domain to a logfile.
    /*!
     *  (virtual from Domain1D)
     *
     *  @param[in]           soln_GlAll_ptr      Pointer to the Global-All solution vector
     *  @param[in]           solnDot_GlAll_ptr   Pointer to the Global-All solution dot vector
     *  @param[in]           soln_ptr            Pointer to the solution vector
     *  @param[in]           solnDot_ptr         Pointer to the time-derivative of the solution vector
     *  @param[in]           solnOld_ptr         Pointer to the solution vector at the old time step
     *  @param[in]           residInternal_ptr   Pointer to the current value of the residual just calculated
     *                                           by a special call to the residEval()
     *  @param[in]           t                   time
     *  @param[in]           rdelta_t            The inverse of the value of delta_t
     *  @param[in]           indentSpaces        Indentation that all output should have as a starter
     *  @param[in]           duplicateOnAllProcs If this is true, all processors will include the same log information as proc 0. If
     *                                           false, the loginfo will only exist on proc 0.
     */
    virtual void
    showSolution(const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
                 const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr,
                 const Epetra_Vector* const solnOld_ptr, const Epetra_Vector_Owned* const residInternal_ptr,
                 const double t, const double rdelta_t, int indentSpaces, bool duplicateOnAllProcs = false);

    //! Base class for writing a solution vector, not the solution, on the domain to a logfile.
    /*!
     *  (virtual from Domain1D)
     *  @param[in]           solnVecName         String name of the solution vector
     *  @param[in]           solnVector_GlAll_ptr Pointer to the Global-All solution vector
     *  @param[in]           solnVector_ptr       Pointer to the solution vector
     *  @param[in]           t                   time
     *  @param[in]           rdelta_t            The inverse of the value of delta_t
     *  @param[in]           indentSpaces        Indentation that all output should have as a starter
     *  @param[in]           duplicateOnAllProcs If this is true, all processors will include the same log information as proc 0. If
     *                                           false, the loginfo will only exist on proc 0.
     *  @param[in]           of                  FILE pointer to write the contents to. Defaults to stdout.
     */
    virtual void
    showSolutionVector(std::string& solnVecName, const Epetra_Vector* const solnVector_GlAll_ptr,
                       const Epetra_Vector* const solnVector_ptr, const double t, const double rdelta_t,
                       int indentSpaces, bool duplicateOnAllProcs = false, FILE* of = stdout);

    //! Base class for writing an int solution vector, not the solution, on the domain to a logfile.
    /*!
     *  (virtual from Domain1D)
     *
     *  @param[in]           solnVecName         String name of the solution vector
     *  @param[in]           solnIntVector_GlAll_ptr Pointer to the Global-All solution vector
     *  @param[in]           solnIntVector_ptr   Pointer to the solution vector
     *  @param[in]           t                   time
     *  @param[in]           rdelta_t            The inverse of the value of delta_t
     *  @param[in]           indentSpaces        Indentation that all output should have as a starter
     *  @param[in]           duplicateOnAllProcs If this is true, all processors will include the same log information as proc 0. If
     *                                           false, the loginfo will only exist on proc 0.
     *  @param[in]           of                  FILE pointer to write the contents to. Defaults to stdout.
     */
    virtual void
    showSolutionIntVector(std::string& solnVecName, const Epetra_IntVector* const solnIntVector_GlAll_ptr,
                          const Epetra_IntVector* const solnIntVector_ptr, const double t, const double rdelta_t,
                          int indentSpaces, bool duplicateOnAllProcs = false, FILE* of = stdout);

    //! Get solution parameters specified by text strings
    /*!
     *  (virtual from Domain1D)
     *
     *  @param[in]           paramID             String name for the item to be requested
     *  @param[out]          paramVal            Vector of information returned.
     *
     *  @return                                  Returns the number of items returned. A value of -1 signifies a failure.
     */
    virtual int
    reportSolutionParam(const std::string& paramID, double* const paramVal) const;

    //! Get vectors of solution quantities requested by text strings
    /*!
     *  (virtual from Domain1D)
     *
     *  @param[in]           requestID           String name for the item to be requested
     *  @param[in]           requestType         Type of the request
     *                                             0    solution variable
     *                                             1    other variable
     *  @param[in]           soln_ptr            Current solution vector (if null, not available)
     *  @param[out]          vecInfo             Vector of information returned.
     *
     *  @return                                  Returns the number of items returned. A value of -1 signifies a failure.
     */
    virtual int
    reportSolutionVector(const std::string& requestID, const int requestType, const Epetra_Vector* const soln_ptr,
                         std::vector<double>& vecInfo) const;

    //! Set the underlying state of the system from the solution vector
    /*!
     *  (virtual from Domain1D)
     *   Note this is an important routine for the speed of the solution.
     *   It would be great if we could supply just exactly what is changing here.
     *   This routine is always called at the beginning of the residual evaluation process.
     *
     *   This is a natural place to put any precalculations of nodal quantities that
     *   may be needed by the residual before its calculation.
     *
     *   Also, this routine is called with delta_t = 0. This implies that a step isn't being taken. However, the
     *   the initial conditions must be propagated.
     *
     *   Note, in general t may not be equal to t_old + delta_t. If this is the case, then the solution is
     *   interpolated across the time interval and then the solution applied.
     *
     *   If doTimeDependentResid then delta_t > 0.
     *   If !doTimeDependentResid then usually delta_t = 0 but not necessarily
     *
     *  @param[in]           doTimeDependentResid Boolean indicating whether we should formulate the time dependent residual
     *  @param[in]           soln              Solution vector. This is the input to the residual calculation
     *  @param[in]           solnDot           Solution vector. This is the input to the residual calculation.
     *  @param[in]           t                 Current time
     *  @param[in]           delta_t           delta t. If zero then delta_t equals 0.
     *  @param[in]           t_old             Old value of the time
     */
    virtual void
    setStateFromSolution(const bool doTimeDependentResid, const Epetra_Vector_Ghosted* const soln,
                         const Epetra_Vector_Ghosted* const solnDot, const double t, const double delta_t, const double t_old);

    //! Generate the initial conditions for the problem. This routine finds the solution vector and solution dot vector given
    //! the specification of the problem.
    /*!
     *  (virtual from Domain1D)
     *  The basic algorithm is to loop over the volume domains. Then, we loop over the surface domains
     *
     *  @param[in]            doTimeDependentResid Boolean indicating whether we should formulate the time dependent residual
     *  @param[out]           soln                Solution vector. This is the input to the residual calculation.
     *  @param[out]           solnDot             Solution vector. This is the input to the residual calculation.
     *  @param[in]            t                   Time
     *  @param[in]            delta_t             delta_t for the initial time step
     */
    virtual void
    initialConditions(const bool doTimeDependentResid, Epetra_Vector* const soln, Epetra_Vector* const solnDot,
                      const double t, const double delta_t);

    //! Evaluate a vector of delta quantities to use when evaluating the Jacobian by numerical differencing
    /*!
     *  (virtual from Domain1D)
     *  @param[in]           t                   Time
     *  @param[in]           soln                Solution vector. This is the input to the algorithm for picking a delta value
     *  @param[in]           solnDot_ptr         Pointer to the time-derivative of the solution vector
     *  @param[out]          deltaSoln           Reference to the Epetra_Vector that will contain the solution deltas.
     *  @param[in]           atolVector_ptr      Pointer to  the atol vector for the solution unknowns
     *  @param[in]           solveType           Type of solution Type. Uses the Solve_Type_Enum  type.
     *                                           Defaults to  TimeDependentAccurate_Solve
     *  @param[in]           solnWeights         Pointer to the solution weights vector. Defaults to nullptr.
     */
    virtual void
    calcDeltaSolnVariables(const double t, const Epetra_Vector& soln, const Epetra_Vector* const solnDot_ptr,
                           Epetra_Vector& deltaSoln, const Epetra_Vector* const atolVector_ptr,
                           const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve,
                           const  Epetra_Vector* const solnWeights = nullptr);

    //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
    /*!
     *  (virtual from Domain1D)
     *
     *  @param[in]           atolDefault         Default atol value
     *  @param[in]           soln                Solution vector. This is a constant the residual calculation.
     *  @param[out]          atolVector          Reference for the atol vector to fill up
     *  @param[in]           atolV               A previously defined atol vector from a previous context. Defaults to nullptr.
     */
    virtual void setAtolVector(double atolDefault, const Epetra_Vector_Ghosted& soln,
                               Epetra_Vector_Ghosted& atolVector, const Epetra_Vector_Ghosted* const atolV = nullptr);

    //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
    /*!
     *  (virtual from Domain1D)
     *  @param[in]           atolDefault         Default atol value
     *  @param[in]           soln                Solution vector. This is a constant the residual calculation.
     *  @param[in]           solnDot             Current solutionDot vector.
     *  @param[out]          atolVector_DAEInit  Reference for the atol vector to fill up
     *  @param[in]           atolV               A previously defined atol vector from a previous context. Defaults to nullptr.
     */
    virtual void 
    setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted& soln, const Epetra_Vector_Ghosted& solnDot,
                          Epetra_Vector_Ghosted& atolVector_DAEInit, const Epetra_Vector_Ghosted* const atolV = nullptr);

    //! Evaluates the atol vector used in the delta damping process.
    /*!
     *  (virtual from Domain1D)
     *
     *  @param[in]           atolDefault         Default atol value
     *  @param[in]           relcoeff            Relative constant to multiply all terms by
     *  @param[in]           soln                Current solution vector.
     *  @param[out]          atolDeltaDamping    If non-zero, this copies the vector into the object as input
     *                                           The default is zero.
     *  @param[in]           atolV               A previously defined atol vector from a previous context. Defaults to nullptr.
     */
    virtual void
    setAtolDeltaDamping(double atolDefault, double relcoeff, const Epetra_Vector_Ghosted& soln,
                        Epetra_Vector_Ghosted& atolDeltaDamping, const Epetra_Vector_Ghosted* const atolV = nullptr);

    //! Evaluates the atol vector used in the delta damping process for the DAE problem
    /*!
     *  (virtual from Domain1D)
     *
     *  @param[in]           atolDefault         Default atol value
     *  @param[in]           relcoeff            Relative constant to multiply all terms by
     *  @param[in]           soln                Current solution vector.
     *  @param[in]           solnDot             Current solutionDot vector.
     *  @param[out]          atolDeltaDamping    Vector to be filled up.
     *  @param[in]           atolV               A previously defined atol vector from a previous context. Defaults to nullptr.
     */
    virtual void
    setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, const Epetra_Vector_Ghosted& soln,
                                const Epetra_Vector_Ghosted& solnDot, Epetra_Vector_Ghosted& atolDeltaDamping,
                                const Epetra_Vector_Ghosted* const atolV = nullptr);

    //! Increment the counters
    /*!
     *  @param[in]           residType           Type of the residual call
     */
    void incrementCounters(const Zuzax::ResidEval_Type residType);

    // ------------------------------------------------  DATA ---------------------------------------------------------------

    //! Number of equations associated with this domain, whether it be a bulk or surface domain
    int NumDomainEqns;

    //! Identifying tag for the domain
    std::string m_id;

    //! The type of coordinate system that is used
    /*!
     *  There are two that are envisioned: Rectinear_Coordinates and Cylindrical_Coordinates
     */
    CoordinateSystem_Type_Enum coordinateSystemType_;

    //! CrossSectional Area of the domain
    /*!
     *  This get's copied to the object from the problem statement.
     *  The default is 1 m**2
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

    //! Reference Temperature (Kelvin)
    /*!
     *  For each domain, we have a reference temperature. This temperature will be used for property
     *  evaluation as the default temperature within the domain whenever there isn't another source for
     *  the value of the temperature.
     *     The default is 298.15 Kelvin
     */
    double TemperatureReference_;

    //! The initial state of the battery has some initial stress due to the mfg process. This produces no strain/displacement.
    /*!
     *  While this is being imput as a scalar, it in principal could be an initial condition of each node, making it a function
     * of position in the battery layers.
     */
    double SolidStressAxialRef_;

    //! Reference Pressure  (pascal)
    /*!
     *  For each domain, we have a reference thermodynamic pressure. This pressure will be used for property
     *  evaluation as the default pressure within the domain whenever there isn't another source for
     *  the value of the thermodynamic pressure.
     *    The default is 1.01325E5 Pascal = OneAtm
     */
    double PressureReference_;

    //! Integer representing the energy equation problem type
    /*!
     *  0 -> isothermal               Don't solve an energy equation (default)
     *  1 -> Fixed Temperature Profile Don't solve an energy equation
     *  2 -> Dirichlet Equation       Solve a Dirichlet equation for temperature.
     *                                This is a way to do the fixed system while keeping the
     *                                matrix structure the same.
     *  3 -> Enthalpy Equation        Solve a full enthalpy equation for the temperature
     *  4 -> Temperature Equation     Solve a Cp dT/dt formulation for the temperature
     */
    int energyEquationProbType_;

    //! Integer representing the solid mechanics problem type
    /*!
     *  0 -> none                     Don't solve an stress-strain relationship for mesh motion
     *  1 -> LinearElastic            Solve for mesh motion using a global simple stress-strain relationship
     */
    int solidMechanicsProbType_;

    //! Porosity equation type
    /*!
     *  This turns on the calculation of the porosity volume fraction in the equation system.
     *  This also
     *
     *  None     =                   0x00,
     *  Constant =                   0x01,
     *  CalculatedOutOfEqnSystem =   0x02,
     *  CalculatedInEqnSystem  =     0x04,
     *  PartOfMechanics =            0x08,
     *  AddedPhasesInEqnSystem =     0x16
     */
    int porosityEquationProbType_;

protected:

    //! Current value of the residual type
    Zuzax::ResidEval_Type residType_Curr_;

private:
    //! Local error routine
    /*!
     *
     * @param msg Meesage indicating where the error message is originating from
     */
    void err(const char* msg) const;

public:
    //!  Solid mechanics equation type enum
    enum SolidMechEqn {
        None =        0x00,
        All  =        0x01,
        ChemEx =      0x02,
        TempEx =      0x04,
        FluidPr =     0x08
    };

    //! Counter for the total number of base residual calculations undertaken
    int counterResBaseCalcs_;

    //! Counter for the total number of base Jacobian residual calculations undertaken
    int counterJacBaseCalcs_;

    //! Counter for the total number of Jacobian delta residual calculations undertaken
    int counterJacDeltaCalcs_;

    //! Counter for the total number of show residual calculations undertaken
    int counterResShowSolutionCalcs_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
