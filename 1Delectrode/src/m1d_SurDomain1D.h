/**
 * @file m1d_SurDomain1D.h
 *  Basic object to calculate the surface residuals for surface domains.
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#ifndef M1D_SURDOMAIN1D_H_
#define M1D_SURDOMAIN1D_H_
//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_Domain1D.h"
#include "m1d_SurfDomainDescription.h"
#include "m1d_BoundaryCondition.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
// Forward declarations
class NodalVars;

//! Definition of a function pointer that takes the time and returns a double
/*!
 *  This typdef is used to create a hook for making boundary condition functions
 */
typedef double (*TimeDepFunctionPtr)(double time);

//==================================================================================================================================
//! Basic object to calculate the surface residuals for surface domains.
/*!
 *  This is a heavyweight base class that provides the function evaluation for a single surface domain.
 *  Surface domains are identified to exist at one and only one node in the domain. There may be multiple bulk domains
 *  attached to this one node.
 *
 *  It's paired up with a surface domain description class that gets set at input, and provides the editorial control.
 *
 *  Notes on the MP implementation:
 *      All surface domains objects will exist on all processors. However,
 *      the local equation variables corresponding to those surface domains
 *      may not exist on the current processor. Therefore, we need to always check
 *      whether a local node and local equation index exists on the current processor  before proceeding.
 */
class SurDomain1D : public Domain1D
{
public:

    //! Constructor
    /*!
     *  @param[in]           sdd                 Contains the surface domain description, which is a required parameter
     */
    SurDomain1D(SurfDomainDescription& sdd);

    //! Copy constructor
    /*!
     *  @param[in]           r                   Object to be copied
     */
    SurDomain1D(const SurDomain1D& r);

    //! Destructor
    virtual ~SurDomain1D();

    //! Assignment operator
    /*!
     *  @param[in]           r                   Object to be copied into the current object
     *  @return                                  Returns a changeable reference to the current object
     */
    SurDomain1D& operator=(const SurDomain1D& r);

    //! Returns the identifying tag for the domain
    /*!
     *  (virtual from Domain1D)
     *  @return                                    Returns the string name for the domain
     */
    virtual std::string id() const;

    //! Prepare all of the indices for fast calculation of the residual
    /*!
     *  Ok, at this point, we will have figured out the number of equations
     *  to be calculated at each node point. The object NodalVars will have  been fully formed.
     *
     *  We use domain_prep() to figure out what the local node number is and what equations correspond to what unknown.
     *
     *  Child objects of domain_prep() will normally call parent classes in a recursive fashion.
     *
     *  @param[in]           li_ptr              Pointer to the LocalNodeIndices structure
     */
    virtual void
    domain_prep(LocalNodeIndices* const li_ptr);

    //! Basic function to calculate the residual for the current domain.
    /*!
     *  (virtual from Domain1D)
     *
     *  All residual terms are written with the following sign convention
     *  based on keeping the time derivative term positive.
     *
     *       res = dcdt - dc2 /dx2 - src = 0
     *
     *  For SurDomains, we calculate the additions and/or replacement of the residual here for the equations 
     *  that this dirichlet condition is responsible for.
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
     *  @param[in]             solveType           Type of solution Type. Uses the Solve_Type  type.
     *                                             Defualts to  TimeDependentAccurate_Solve
     */
    virtual void
    residEval(Epetra_Vector& res, const bool doTimeDependentResid,
              const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr,
              const Epetra_Vector* const solnOld_ptr, const double t, const double rdelta_t, 
              const Zuzax::ResidEval_Type residType = Zuzax::ResidEval_Type::Base_ResidEval,
              const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve) override;

    //! Generate the initial conditions for the problem. This routine finds the solution vector and solution dot vector given
    //! the specification of the problem.
    /*!
     *  (virtual from Domain1D)
     *  For surface dirichlet conditions, we impose the t = 0- condition.
     *
     *  @param[in]            doTimeDependentResid Boolean indicating whether we should formulate the time dependent residual
     *  @param[out]           soln                Solution vector. This is the input to the residual calculation.
     *  @param[out]           solnDot             Solution vector. This is the input to the residual calculation.
     *  @param[in]            t                   Time
     *  @param[in]            delta_t             delta_t for the initial time step
     */
    virtual void
    initialConditions(const bool doTimeDependentResid, Epetra_Vector* const soln, Epetra_Vector* const solnDot,
                      const double t, const double delta_t) override;

    //! Base class for saving the solution on the domain in an xml node.
    /*!
     *  (virtual from Domain1D)
     *  @param[in,out]         oNode               Reference to the XML_Node
     *  @param[in]             soln_GlAll_ptr      Pointer to the Global-All solution vector
     *  @param[in]             solnDot_GlAll_ptr   Pointer to the time derivative of the Global-All solution vector
     *  @param[in]             t                   time
     *  @param[in]             duplicateOnAllProcs If this is true, all processors will include the same XML_Node information
     *                                             as proc 0. If false, the xml_node info will only exist on proc 0.
     *                                             Defaults to false.
     */
    virtual void
    saveDomain(ZZCantera::XML_Node& oNode, const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
               const double t, bool duplicateOnAllProcs = false) override;

    //! Base Class for reading the solution from the saved file
    /*!
     *  (virtual from Domain1D)
     *  This class assumes that the XML_Node is the domain node in the example below.
     *
     * @verbatim
     *  <simulation id="0">
     *    <time type="time" units="s" vtype="float"> 0.000000000000000E+00 </time>
     *    <delta_t type="time" units="s" vtype="float"> 1.000000000000000E-08 </delta_t>
     *    <StepNumber type="time" vtype="integer"> 0 </StepNumber>
     *    <domain id="SurDomain1D_0" numVariables="6" points="1" type="surface">
     *      <X0 vtype="float"> 0.000000000000000E+00 </X0>
     *      <X vtype="float"> 0.000000000000000E+00 </X>
     *      <Vel_axial(0) vtype="float"> 0.000000000000000E+00 </Vel_axial(0)>
     *      <MF_sp(ECDMC) vtype="float"> 9.200000000000000E-01 </MF_sp(ECDMC)>
     *      <MF_sp(Li+) vtype="float"> 4.000000000000000E-02 </MF_sp(Li+)>
     *      <MF_sp(PF6-) vtype="float"> 4.000000000000000E-02 </MF_sp(PF6-)>
     *      <Volt(Electrolyte) vtype="float"> -7.000000000000001E-02 </Volt(Electrolyte)>
     *      <Volt(AnodeVoltage) vtype="float"> 0.000000000000000E+00 </Volt(AnodeVoltage)>
     *    </domain>
     *  </simulation>
     * @endverbatim
     *
     * @param[in]            domainNode          Reference to the XML_Node to read the solution from
     * @param[in]            soln_GlAll_ptr      Pointer to the Global-All solution vector
     * @param[in]            solnDot_GlAll_ptr   Pointer to the time derivative of the Global-All solution vector
     * @param[in]            globalTimeRead      Value of the global time that is read in. This is used for
     *                                           comparison and quality control purposes
     */
    virtual void
    readDomain(const ZZCantera::XML_Node& domainNode, Epetra_Vector* const soln_GlAll_ptr,
               Epetra_Vector* const solnDot_GlAll_ptr, double globalTimeRead) override;


    //! Fill the vector isAlgebraic with the values from the DomainDescription
    /*!
     *  (virtual from SurDomain1D)
     *  @param[in]           isAlgebraic         Epetra_IntVector to be filled with the IsAlgebraic values
     */
    virtual void fillIsAlgebraic(Epetra_IntVector& isAlgebraic);

    //!  Fill the vector isArithmeticScaled with the values from the DomainDescription
    /*!
     *  (virtual from SurDomain1D)
     *  @param[in]           isArithmeticScaled  Epetra_IntVector to be filled with the IsArithmeticScaled values
     */
    virtual void fillIsArithmeticScaled(Epetra_IntVector& isArithmeticScaled);

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
                           const  Epetra_Vector* const solnWeights=nullptr) override;

    //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
    /*!
     *  (virtual from Domain1D)
     *  @param[in]           atolDefault         Default atol value
     *  @param[in]           soln                Solution vector. This is a constant the residual calculation.
     *  @param[out]          atolVector          Reference for the atol vector to fill up
     *  @param[in]           atolV               A previously defined atol vector from a previous context. Defaults to nullptr.
     */
    virtual void setAtolVector(double atolDefault, const Epetra_Vector_Ghosted& soln,
                               Epetra_Vector_Ghosted& atolVector, const Epetra_Vector_Ghosted* const atolV = nullptr) override;

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
                          Epetra_Vector_Ghosted& atolVector_DAEInit, 
                          const Epetra_Vector_Ghosted* const atolV = nullptr) override;

    //! Evaluates the atol vector used in the delta damping process.
    /*!
     *  (virtual from Domain1D)
     *  @param[in]           atolDefault         Default atol value
     *  @param[in]           relcoeff            Relative constant to multiply all terms by
     *  @param[in]           soln                Current solution vector.
     *  @param[out]          atolDeltaDamping    If non-zero, this copies the vector into the object as input
     *                                           The default is zero.
     *  @param[in]           atolV               A previously defined atol vector from a previous context. Defaults to nullptr.
     */
    virtual void
    setAtolDeltaDamping(double atolDefault, double relcoeff, const Epetra_Vector_Ghosted& soln,
                        Epetra_Vector_Ghosted& atolDeltaDamping, const Epetra_Vector_Ghosted* const atolV = nullptr) override;

    //! Evaluates the atol vector used in the delta damping process for the DAE problem
    /*!
     *  (virtual from Domain1D)
     *  @param[in]           atolDefault         Default atol value
     *  @param[in]           relcoeff            Relative constant to multiply all terms by
     *  @param[in]           soln                Current solution vector.
     *  @param[in]           solnDot             Current time-derivative of the solution vector
     *  @param[out]          atolDeltaDamping    If non-zero, this copies the vector into the object as input
     *                                           The default is zero.
     *  @param[in]           atolV               A previously defined atol vector from a previous context. Defaults to nullptr.
     */
    virtual void
    setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, const Epetra_Vector_Ghosted& soln,
                                const Epetra_Vector_Ghosted& solnDot, Epetra_Vector_Ghosted& atolDeltaDamping,
                                const Epetra_Vector_Ghosted* const atolV = nullptr) override;

    //! Method for writing the header for the surface domain to a tecplot file.
    /*!
     *  (virtual from Domain1D)
     * Only proc0 will write tecplot files.
     */
    virtual void writeSolutionTecplotHeader() override;

    //! Method for writing the solution on the surface domain to a tecplot file.
    /*!
     *  (virtual from Domain1D)
     *  Only proc0 will write tecplot files.Therefore, we must be sure to always use a Epetra_Vector_GlAll solution type.
     *
     *  @param[in]           soln_GlAll_ptr      Pointer to the Global-All solution vector
     *  @param[in]           solnDot_GlAll_ptr   Pointer to the time derivative of the Global-All solution vector
     *  @param[in]           t                   time
     */
    virtual void 
    writeSolutionTecplot(const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
                         const double t) override;

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
                 const double t, const double rdelta_t, int indentSpaces, bool duplicateOnAllProcs = false) override;
 
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
                       int indentSpaces, bool duplicateOnAllProcs = false, FILE* of = stdout) override;

    //! Base class for writing an int solution vector, not the solution, on the domain to a logfile.
    /*!
     *  (virtual from Domain1D)
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
                          int indentSpaces, bool duplicateOnAllProcs = false, FILE* of = stdout) override;

    //! Class for writing the solution on the domain to a logfile.
    /*!
     *  @param[in]           soln_GlAll_ptr      Pointer to the Global-All solution vector
     *  @param[in]           solnDot_GlAll_ptr   Pointer to the Global-All solution dot vector
     *  @param[in]           soln_ptr            Pointer to the solution vector
     *  @param[in]           solnDot_ptr         Pointer to the time-derivative of the solution vector
     *  @param[in]           solnOld_ptr         Pointer to the solution vector at the old time step
     *  @param[in]           residInternal_ptr   Pointer to the current value of the residual just calculated
     *                                             by a special call to the residEval()
     *  @param[in]           t                   time
     *  @param[in]           rdelta_t            The inverse of the value of delta_t
     *  @param[in]           indentSpaces        Indentation that all output should have as a starter
     *  @param[in]           duplicateOnAllProcs If this is true, all processors will include the same log information as proc 0. If
     *                                           false, the loginfo will only exist on proc 0.
     */
    void
    showSolution0All(const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
                     const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr,
                     const Epetra_Vector* const solnOld_ptr, const Epetra_Vector_Owned* const residInternal_ptr,
                     const double t, const double rdelta_t, int indentSpaces, bool duplicateOnAllProcs = false);

    //! Extract the double value out of a solution vector for a particular
    //! variable defined at a node corresponding to the surface domain
    /*!
     *  @param[in]           soln_ptr            Pointer to the ghosted solution vector
     *  @param[in]           v1                  VarType of the variable to be returned
     *
     *  @return                                  Returns the value of the variable. If the variable doesn't exist
     *                                           on the processor this routine returns the value of -1.0E300.
     */
    double
    extractSolnValue(const Epetra_Vector_Ghosted* const soln_ptr, VarType v1);

    //! Transfer the bulk flux vectors to the surface flux vectors
    /*!
     *  (virtual from SurDomain1D)
     *  This routine update the following right flux vectors from the neighboring right bulk domain.
     *
     *    DiffFluxRightBulkDomain_LastResid_NE[i]
     *	  TotalFluxRightBulkDomain_LastResid_NE[i]
     *	  VarVectorRightBulkDomain_LastResid_NE[i]
     *
     *  It then updates the corresponding left flux vectors from the left bulk domain.
     *
     *  This routine is typically called from the residEval routine. However, it's modular
     *  enough to be carved out as its own routine. We make it virtual because it may be overridden in other routines.
     */
    virtual void updateBulkFluxVectors();

    // --------------------------------------- D A T A ----------------------------------------------------------------

    //! Light description of what this domain is about
    /*!
     *   This is currently a 1 to 1 mapping between SurfDomainDescription structures
     *   and SurfDomain1D structures. This may be relaxed in the future.
     */
    SurfDomainDescription& SDD_;

    //! Number of owned nodes in this domain
    /*!
     *  For surfaces, this is either 1 or 0
     */
    int NumOwnedNodes;

    //! First global node that is owned by this processor.
    int FirstGbNode;

    //! Last Global node that is owned by this processor
    int LastOwnedGbNode;

    //! Number of equations this boundary condition is to be
    //! applied on.
    /*!
     *   This number is used to loop over the number of equations.
     *   Subtypes are counted as different equations in this loop.
     */
    int NumBCs;

    //! True if this processor owns the left-most node  of this domain
    bool IOwnLeft;

    //! True if this processor owns the right-most node of this domain
    bool IOwnRight;

    //! Index of the local node that will have this boundary condition applied
    /*!
     *  This is the index of the local node number. If the global node
     *  doesn't exist on this processor, this value will be -1.
     *  Always check for whether the node exists on the processor by checking to see
     *  if this value is greater or equal to zero.
     */
    int Index_LcNode;

    //! Pointer to the NodalVars object for the node that contains this surface domain
    /*!
     *  This is set even if the processor doesn't own the node.
     */
    NodalVars* NodalVarPtr;

    //! This is the index of this surface domain in the list of surface domains
    //! that are located at the current node.
    int Index_NodalSD_;

    //! Pointer to the Class that contains the Local Node Indices for this processor
    LocalNodeIndices* LI_ptr_;

    //! Number of equations defined at the current node on which the SurDomain1D is located
    /*!
     *  This will be equal to or greater than the number of unknowns assigned to the surface domain.
     *  This is set even if this processor doesn't own the node.
     *  Frequently, this will be used a loop variable.
     */
    size_t NumNodeEqns;

    //! Vector of flags to indicate whether the variable has a time derivative
    /*!
     *   -1 No information
     *    0 DAE
     *    1 Has a time derivative
     *
     *  length of the number of equations defined at this node
     */
    std::vector<int> IsAlgebraic_Node;

    //! Vector of flags to indicate whether the variable has a time derivative
    /*!
     *   -1 No information
     *    0 regular
     *    1 arithmetic scaled
     *
     *  length of the number of equations defined at this node
     */
    std::vector<int> IsArithmeticScaled_Node;

    //! Number of Equations per node in the BulkDomain1d to the left
    int NumDomainEqnsLeft_;

    //! Number of Equations per node in the BulkDomain1d to the Right
    int NumDomainEqnsRight_;

    //! Diffusive fluxes into the left bulk domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the diffusive flux calculated
     *  at the left Boundary during the last residual calculation. Flux is dotted with normal
     *  oriented away from the surface.
     *
     *  This quantities are useful (I think) for the specification of boundary
     *  conditions for multi domains.
     *
     *  They are also useful for the specification of global balances, when Dirichlet conditions are set
     *  on some variables. Then, we seek the specification of fluxes which preserve global conservation laws.
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> DiffFluxLeftBulkDomain_LastResid_NE;

    //! Diffusive Fluxes into the right bulk domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the diffusive flux calculated
     *  at the right boundary during the last residual calculation. Flux is dotted with normal
     *  oriented away from the surface.
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> DiffFluxRightBulkDomain_LastResid_NE;

    //! Total Fluxes into the left bulk domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the diffusive flux calculated
     *  at the left boundary during the last residual calculation. Flux is dotted with normal
     *  oriented away from the surface.
     *
     *  This quantities are useful (I think) for the specification of boundary
     *  conditions for multi domains.
     *
     *  They are also useful for the specification of global balances, when Dirichlet conditions are set
     *  on some variables. Then, we seek the specification of fluxes which preserve global conservation laws.
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> TotalFluxLeftBulkDomain_LastResid_NE;

    //! Total Fluxes into the right bulk domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the total flux calculated
     *  at the right boundary during the last residual calculation. Flux is dotted with normal
     *  oriented away from the surface.
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> TotalFluxRightBulkDomain_LastResid_NE;

    //! Variable Vector from the left bulk domain
    /*!
     *  This is the value of the solution vector from the left bulk domain. We store a temporary copy here
     *  Length: NumDomainEqnsLeft_
     */
    //std::vector<double> VarVectorLeftBulkDomain_LastResid_NE;

    //! Variable Vector from the right bulk domain
    /*!
     *  This is the value of the solution vector from the right bulk domain. We store a temporary copy here
     *   Length: NumDomainEqnsRight_
     */
    //std::vector<double> VarVectorRightBulkDomain_LastResid_NE;

    //!  Residual fed into the surface domain before application of BCs
    /*!
     *   This may be used to close balances for cases where Dirichlet conditions are applied at this surface.
     *   Length: Number of equations defined at this surface node = NumNodeEqns
     */
    std::vector<double> Resid_BeforeSurDomain_NE;

    //!  Vector that stores the SurDomain1D's contributions to the residual equations for all node equations defined on the node
    /*!
     *   This is only used for printout purposes when residType == Base_ShowSolution
     *   Length:  NumNodeEqns
     */
    std::vector<double> DomainResidVector_LastResid_NE;

private:
    //! Call an error exit condition
    /*!
     *  @param[in]           msg                 Message to be printed out
     */
    void err(const char* const msg);
};
//==================================================================================================================================
//! Specification of a set of simple Dirichlet conditions on a surface domain
/*!
 *  This boundary condition specifies one or more Dirichlet conditions on variables at an interface.
 *  It sets up and maintains all of the general structures necessary for handling Dirichlet conditions at interfaces.
 */
class SurBC_Dirichlet : public SurDomain1D
{
public:

    //! Constructor
    /*!
     *
     *  @param[in]           sdd                 Contains the surface domain description.
     */
    SurBC_Dirichlet(SurfDomainDescription& sdd);

    //! Copy Constructor
    /*!
     *  @param[in]             r                   Object to be copied
     */
    SurBC_Dirichlet(const SurBC_Dirichlet& r);

    //! Destructor
    virtual ~SurBC_Dirichlet();

    //! Assignment operator
    /*!
     *  @param[in]           r                   Object to be copied into the current object
     *  @return                                  Returns a changeable reference to the current object
     */
    SurBC_Dirichlet& operator=(const SurBC_Dirichlet& r);

    //! Prepare all of the indices for fast calculation of the residual
    /*!
     *  (virtual from Domain1D)
     *  Here we collect all of the information necessary to speedily implement SpecFlag_NE and Value_NE within the
     *  residual calculation.  We transfer the information from SDT_Dirichlet structure to this structure for quick processing.
     *
     *  @param[in]             li_ptr              Pointer to the LocalNodeIndices Structure that contains information
     *                                             about how the mesh is layed out within this domain and other domains in the problem
     */
    virtual void
    domain_prep(LocalNodeIndices* const li_ptr) override;

    //! Change the value of an existing dirichlet condition
    /*!
     *  (virtual from SurBC_Dirichlet)
     *  @param[in]           vtDir               Variable type that is to be changed
     *  @param[in]           newVal              new value of dirichlet condition
     *
     *  @return                                  Number of boundary conditions that are changed.
     */
    virtual int
    changeDirichletConditionValue(VarType vtDir, double newVal);

    //! Change the boundary condition applied to a variable
    /*!
     *  @param[in]           vtDir               Variable type class. Note, general matches are allowed with this parameter
     *  @param[in]           BC_Type             Type of the boundary condition
     *  @param[in]           value               Value of the dirichlet condition or flux - default 0.0
     *  @param[in]           BC_TimeDep          BoundaryCondition Pointers for time dependent BC for BC_Tppe = 3,4,5
     *                                             (default 0)
     *  @param[in]           TimeDep             Function pointer to a function that returns a double given a single
     *                                             parameter (the time). Defaults to a NULL pointer.
     *
     *  @return                                  Returns the number of boundary conditions matched.
     *                                             A negative number means that an error has been encountered
     */
    int
    changeBoundaryCondition(VarType vtDir, int BC_Type, double value = 0.0, BoundaryCondition* BC_TimeDep = nullptr,
                            TimeDepFunctionPtr TimeDep = nullptr);

    //! Report on the boundary condition applied on the first match to VarType
    /*!
     *  @param[in]           time                Current time for evaluating time dependent BC
     *  @param[in]           vtDir               Variable type class. Note, general matches are allowed with this parameter
     *  @param[out]          BC_Type             Type of the boundary condition
     *  @param[out]          value               Value of the Dirichlet condition or flux - default 0.0
     *  @param[out]          BC_TimeDep          BoundaryCondition Pointers for time dependent BC for BC_Type = 3,4,5
     *                                             (default 0)
     *  @param[out]          TimeDep             Function pointer to a function that returns a double given a single parameter (the time).
     *                                           Defaults to a NULL pointer.
     *
     *  @return                                  Returns the number of boundary conditions matched.
     *                                           A negative number means that an error has been encountered
     */
    int
    reportBoundaryCondition(double time, const VarType vtDir, int& BC_Type, double& value, BoundaryCondition*& BC_TimeDep,
                            TimeDepFunctionPtr& TimeDep) const;

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
     *  We calculate the additions and/or replacement of the residual here for the equations that this dirichlet condition
     *  is responsible for.
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
    residEval(Epetra_Vector& res, const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
              const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr,
              const double t, const double rdelta_t, const Zuzax::ResidEval_Type residType = Zuzax::ResidEval_Type::Base_ResidEval,
              const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve) override;

    //! Base class for saving the solution on the domain in an xml node.
    /*!
     *  (virtual from Domain1D)
     *  @param[in,out]         oNode               Reference to the XML_Node
     *  @param[in]             soln_GlAll_ptr      Pointer to the Global-All solution vector
     *  @param[in]             solnDot_GlAll_ptr   Pointer to the time derivative of the Global-All solution vector
     *  @param[in]             t                   time
     *  @param[in]             duplicateOnAllProcs If this is true, all processors will include the same XML_Node information
     *                                             as proc 0. If false, the xml_node info will only exist on proc 0.
     *                                             Defaults to false.
     */
    virtual void
    saveDomain(ZZCantera::XML_Node& oNode, const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
               const double t, bool duplicateOnAllProcs = false) override;

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
                 const double t, const double rdelta_t, int indentSpaces, bool duplicateOnAllProcs = false) override;

    //! Generate the initial conditions for the problem. This routine finds the solution vector and solution dot vector given
    //! the specification of the problem.
    /*!
     *  (virtual from Domain1D)
     *   For surface dirichlet conditions, we impose the t = 0- condition.
     *
     *  @param[in]            doTimeDependentResid Boolean indicating whether we should formulate the time dependent residual
     *  @param[out]           soln                Solution vector. This is the input to the residual calculation.
     *  @param[out]           solnDot             Solution vector. This is the input to the residual calculation.
     *  @param[in]            t                   Time
     *  @param[in]            delta_t             delta_t for the initial time step
     */
    virtual void
    initialConditions(const bool doTimeDependentResid, Epetra_Vector* const soln, Epetra_Vector* const solnDot,
                      const double t, const double delta_t) override;

    // ------------------------------------------- D A T A -------------------------------------------------------------------------

    //! Boolean flag indicating which variables at the node are being specified
    //! with Dirichlet conditions or with other boundary conditions such as flux conditions.
    /*!
     *   Vector has length equal to the number of equations defined at the node
     */
    std::vector<int> SpecFlag_NE;

    //! Value of the variable
    /*!
     *   Vector has length equal to the number of equations defined at the node
     */
    std::vector<double> Value_NE;

    //! Function Pointers for time dependent BC for BC_Type_NE = 2
    /*!
     *   Vector has length equal to the number of equations defined at the node
     *   The default is the null pointer
     */
    std::vector<TimeDepFunctionPtr> TimeDep_NE;

    //! BoundaryCondition Pointers for time dependent BC for BC_Type_NE = 3,4,5
    /*!
     *   Vector has length equal to the number of equations defined at the node
     */
    std::vector<BoundaryCondition*> BC_TimeDep_NE;

    //! Type of the Boundary condition to be applied to that variable
    /*!
     *  -1 Boundary condition is not used
     *   0 Dirichlet
     *   1 Flux condition
     *   2 Time Dependent Dirichlet. A function provides the dependence
     *   3
     *   4
     *   5
     *   6
     *   7
     *   8
     *   9
     *  10  Robin boundary condition
     *
     *   Vector has length equal to the number of equations defined at the node
     */
    std::vector<int> BC_Type_NE;
};

//==================================================================================================================================
} /* End of namespace */
//----------------------------------------------------------------------------------------------------------------------------------
#endif /* M1D_SURDOMAIN1D_H_ */
