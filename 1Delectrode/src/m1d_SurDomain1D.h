/**
 * @file m1d_SurDomain1D.h
 *  Basic object to calculate the surface residuals for surface domains.
 */
/*
 *    $Id: m1d_SurDomain1D.h 363 2012-08-22 03:37:42Z hkmoffa $
 */

#ifndef M1D_SURDOMAIN1D_H_
#define M1D_SURDOMAIN1D_H_
//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_DomainDescription.h"
#include "m1d_Domain1D.h"
#include "m1d_ProblemResidEval.h"
#include "m1d_BoundaryCondition.h"

#include "Epetra_Vector.h"

namespace m1d
{
// Forward declarations
class LocalNodeIndices;


typedef double (*TimeDepFunctionPtr)(double time);


//! Basic object to calculate the surface residuals for surface domains.
/*!
 *  This is a heavyweight base class that provides the function
 *  evaluation for a single surface domain.
 *
 *   It's paired up with a surface domain description class that gets set
 *   at input, and provides the editorial control.
 *
 *  Notes on the MP implementation:
 *      All surface domains objects will exist on all processors. However,
 *      the local equation variables corresponding to those surface domains
 *      may not exist on the current processor. Therefore, we need to always check
 *      whether a local node and local equation index exists on the current processor
 *      before proceeding.
 */
class SurDomain1D : public Domain1D
{

public:

  //! Constructor
  /*!
   * @param sdd   Contains the surface domain description, which is a required parameter
   */
  SurDomain1D(m1d::SurfDomainDescription &sdd);

  //! Copy constructor
  /*!
   *
   * @param r  Object to be copied
   */
  SurDomain1D(const SurDomain1D &r);

  //! Destructor
  virtual
  ~SurDomain1D();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  SurDomain1D &
  operator=(const SurDomain1D &r);

  //! Returns the identifying tag for the domain
  virtual std::string
  id() const;

  //! Prepare all of the indices for fast calculation of the residual
  /*!
   *  Ok, at this point, we will have figured out the number of equations
   *  to be calculated at each node point. The object NodalVars will have
   *  been fully formed.
   *
   *  We use this to figure out what the local node number is and what
   *  equations correspond to what unknown.
   *
   *  Child objects of this one will normally call this routine in a
   *  recursive fashion.
   */
  virtual void
  domain_prep(LocalNodeIndices *li_ptr);

  //! Basic function to calculate the residual for the domain.
  /*!
   *  We calculate the additions and/or replacement of the
   *  residual here for the equations that this dirichlet condition
   *  is responsible for.
   *
   * @param res           Output vector containing the residual
   * @param doTimeDependentResid  boolean indicating whether the time
   *                         dependent residual is requested
   * @param soln         Solution vector at which the residual should be
   *                     evaluated
   */
  virtual void
  residEval(Epetra_Vector &res,
            const bool doTimeDependentResid,
            const Epetra_Vector *soln_ptr,
            const Epetra_Vector *solnDot_ptr,
            const Epetra_Vector *solnOld_ptr,
            const double t,
            const double rdelta_t,
            const ResidEval_Type_Enum residType = Base_ResidEval,
	    const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);

  //! Generate the initial conditions
  /*!
   *   For surface dirichlet conditions, we impose the t = 0- condition.
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
  virtual void
  initialConditions(const bool doTimeDependentResid,
                    Epetra_Vector *soln,
                    Epetra_Vector *solnDot,
                    const double t,
                    const double delta_t);

  //! Base class for saving the solution on the domain in an xml node.
  /*!
   *
   * @param oNode                Reference to the XML_Node
   * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
   * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
   * @param t                    time
   *
   * @param duplicateOnAllProcs  If this is true, all processors will include
   *                             the same XML_Node information as proc 0. If
   *                             false, the xml_node info will only exist on proc 0.
   */
  virtual void
  saveDomain(ZZCantera::XML_Node& oNode,
             const Epetra_Vector *soln_GlAll_ptr,
             const Epetra_Vector *solnDot_GlAll_ptr,
             const double t,
             bool duplicateOnAllProcs = false);

  //! Base Class for reading the solution from the saved file
  /*!
   *  This class assumes that the XML_Node is <domain> in the example below.
   *
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
   *
   *
   * @param domainNode                Reference to the XML_Node to read the solution from
   * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
   * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
   *
   */
  virtual void 
  readDomain(const ZZCantera::XML_Node& domainNode,
             Epetra_Vector * const soln_GlAll_ptr,
             Epetra_Vector * const solnDot_GlAll_ptr, double globalTimeRead);


  //!  Fill the vector isAlgebraic with the values from the DomainDescription
  /*!
   * @param isAlgebraic  Epetra_IntVector to be filled with the IsAlgebraic values
   */
  virtual void fillIsAlgebraic(Epetra_IntVector  & isAlgebraic);

  //!  Fill the vector isArithmeticScaled with the values from the DomainDescription
  /*!
   * @param isArithmeticScaled  Epetra_IntVector to be filled with the IsArithmeticScaled values
   */
  virtual void fillIsArithmeticScaled(Epetra_IntVector & isArithmeticScaled);

  virtual void
  calcDeltaSolnVariables(const double t, const Epetra_Vector& soln,
			 const Epetra_Vector* solnDot_ptr, Epetra_Vector& deltaSoln,
                         const Epetra_Vector* const atolVector_ptr, 
                         const Solve_Type_Enum solveType = TimeDependentAccurate_Solve,
                         const  Epetra_Vector* solnWeights=0);

  //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
  /*!
   * @param atolDefault             Default atol value
   * @param soln                    Solution vector. This is a constant
   *                                the residual calculation.
   * @param atolVector              (OUTPUT) Reference for the atol vector to fill up
   */
  virtual void setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
			     Epetra_Vector_Ghosted & atolVector,
			     const Epetra_Vector_Ghosted * const atolV = 0);

  //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
  /*!
   * @param atolDefault             Default atol value
   * @param soln                    Solution vector. This is a constant
   *                                the residual calculation.
   * @param atolVector              (OUTPUT) Reference for the atol vector to fill up
   */
  virtual void setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln, 
				     const Epetra_Vector_Ghosted & solnDot, 
				     Epetra_Vector_Ghosted & atolVector_DAEInit,
				     const Epetra_Vector_Ghosted * const atolV = 0);

  //! Evaluates the atol vector used in the delta damping process.
  /*!
   *   @param relcoeff     Relative constant to multiply all terms by
   *   @param soln         current solution vector.
   *   @param atolDeltaDamping      If non-zero, this copies the vector into the object as input
   *                      The default is zero.
   */
  virtual void
  setAtolDeltaDamping(double atolDefault, double relcoeff, 
		      const Epetra_Vector_Ghosted & soln, 
		      Epetra_Vector_Ghosted & atolDeltaDamping,
		      const Epetra_Vector_Ghosted * const atolV = 0);
		    

 //! Evaluates the atol vector used in the delta damping process for the DAE problem
  /*!
   *   @param relcoeff     Relative constant to multiply all terms by
   *   @param soln         current solution vector.
   *   @param solnDot      Current solutionDot vector.
   *   @param atolDeltaDamping       If non-zero, this copies the vector into the object as input
   *                       The default is zero.
   */
  virtual void
  setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, 
			      const Epetra_Vector_Ghosted & soln,
			      const Epetra_Vector_Ghosted & solnDot,
			      Epetra_Vector_Ghosted & atolDeltaDamping,
			      const Epetra_Vector_Ghosted * const atolV = 0);

   
    
  // Method for writing the header for the surface domain to a tecplot file.
  /*
   * Only proc0 will write tecplot files.
   */
  virtual void writeSolutionTecplotHeader();

  //! Method for writing the solution on the surface domain to a tecplot file.
  /*!
   * Only proc0 will write tecplot files. Therefore, we must be sure to always
   * use a Epetra_Vector_GlAll solution type.
   *
   * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
   * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
   * @param t                    time
   *
   */
  virtual void writeSolutionTecplot(const Epetra_Vector_GlAll *soln_GlAll_ptr,
			    const Epetra_Vector_GlAll *solnDot_GlAll_ptr,
			    const double t );


  //! Base class for writing the solution on the domain to a logfile.
  /*!
   *
   * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
   * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
   * @param soln_ptr             Pointer to the solution vector
   * @param solnDot_ptr          Pointer to the time-derivative of the solution vector
   * @param solnOld_ptr          Pointer to the solution vector at the old time step
   * @param residInternal _ptr   Pointer to the current value of the residual just calculated
   *                             by a special call to the residEval()
   * @param t                    time
   * @param rdelta_t             The inverse of the value of delta_t
   * @param indentSpaces         Indentation that all output should have as a starter
   * @param duplicateOnAllProcs  If this is true, all processors will include
   *                             the same log information as proc 0. If
   *                             false, the loginfo will only exist on proc 0.
   */
  virtual void
  showSolution(const Epetra_Vector *soln_GlAll_ptr,
               const Epetra_Vector *solnDot_GlAll_ptr,
               const Epetra_Vector *soln_ptr,
               const Epetra_Vector *solnDot_ptr,
               const Epetra_Vector *solnOld_ptr,
               const Epetra_Vector_Owned *residInternal_ptr,
               const double t,
               const double rdelta_t,
               int indentSpaces,
               bool duplicateOnAllProcs = false);

  //! Base class for writing a solution vector, not the solution, on the domain to a logfile.
  /*!
   * @param solnVecName          Name of the Solution Vector
   * @param solnVector_GlALL_ptr       Pointer to the Global-All solution vector
   * @param solnVector_ptr             Pointer to the solution vector
   * @param t                    time
   * @param rdelta_t             The inverse of the value of delta_t
   * @param indentSpaces         Indentation that all output should have as a starter
   * @param duplicateOnAllProcs  If this is true, all processors will include
   *                             the same log information as proc 0. If
   *                             false, the loginfo will only exist on proc 0.
   */
  virtual void
  showSolutionVector(std::string& solnVecName,
		     const Epetra_Vector *solnVector_GlAll_ptr,
		     const Epetra_Vector *solnVector_ptr,
		     const double t,
		     const double rdelta_t,
		     int indentSpaces,
		     bool duplicateOnAllProcs = false,
		     FILE *of = stdout);

  //! Base class for writing a solution vector, not the solution, on the domain to a logfile.
  /*!
   * @param solnVecName          Name of the Solution Vector
   * @param solnVector_GlALL_ptr       Pointer to the Global-All solution vector
   * @param solnVector_ptr             Pointer to the solution vector
   * @param t                    time
   * @param rdelta_t             The inverse of the value of delta_t
   * @param indentSpaces         Indentation that all output should have as a starter
   * @param duplicateOnAllProcs  If this is true, all processors will include
   *                             the same log information as proc 0. If
   *                             false, the loginfo will only exist on proc 0.
   */
  virtual void
  showSolutionIntVector(std::string& solnVecName,
			const Epetra_IntVector *solnVector_GlAll_ptr,
			const Epetra_IntVector *solnVector_ptr,
			const double t,
			const double rdelta_t,
			int indentSpaces,
			bool duplicateOnAllProcs = false,
			FILE *of = stdout);

  //! Base class for writing the solution on the domain to a logfile.
  /*!
   *
   * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
   * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
   * @param soln_ptr             Pointer to the solution vector
   * @param solnDot_ptr          Pointer to the time-derivative of the solution vector
   * @param solnOld_ptr          Pointer to the solution vector at the old time step
   * @param residInternal _ptr   Pointer to the current value of the residual just calculated
   *                             by a special call to the residEval()
   * @param t                    time
   * @param rdelta_t             The inverse of the value of delta_t
   * @param indentSpaces         Indentation that all output should have as a starter
   * @param duplicateOnAllProcs  If this is true, all processors will include
   *                             the same log information as proc 0. If
   *                             false, the loginfo will only exist on proc 0.
   */
  void
  showSolution0All(const Epetra_Vector *soln_GlAll_ptr,
                   const Epetra_Vector *solnDot_GlAll_ptr,
                   const Epetra_Vector *soln_ptr,
                   const Epetra_Vector *solnDot_ptr,
                   const Epetra_Vector *solnOld_ptr,
                   const Epetra_Vector_Owned *residInternal_ptr,
                   const double t,
                   const double rdelta_t,
                   int indentSpaces,
                   bool duplicateOnAllProcs = false);

  //! Extract the double value out of a solution vector for a particular
  //! variable defined at a node corresponding to the surface domain
  /*!
   *  @param  soln_ptr   Pointer to the ghosted solution vector
   *  @param  v1         VarType of the variable to be returned
   *
   *  @return           Returns the value of the variable. If the variable doesn't exist
   *                    on the processor this routine returns the value of
   *                    -1.0E300.
   */
   double
   extractSolnValue(Epetra_Vector_Ghosted *soln_ptr, VarType v1); 

  //! Transfer the bulk flux vectors to the surface flux vectors
  /*!
   *  This routine update the following right flux vectors from the neighboring right bulk domain
   *      DiffFluxRightBulkDomain_LastResid_NE[i] 
   *	  TotalFluxRightBulkDomain_LastResid_NE[i] 
   *	  VarVectorRightBulkDomain_LastResid_NE[i]
   *
   *  It then updates the corresponding left flux vectors from the left bulk domain
   *
   *  This routine is typically called from the residEval routine. However, it's modular
   *  enough to be carved out as its own routine. We make it virtual because it may be
   *  overridden in other routines.
   */
  virtual void updateBulkFluxVectors();

  // *******************************************************************************
  //   End of functions - start of member data
  // *******************************************************************************

  //! Light description of what this domain is about
  /*!
   *   This is currently a 1 to 1 mapping between SurfDomainDescription structures
   *   and SurfDomain1D structures. This may be relaxed in the future.
   */
  SurfDomainDescription &SDD_;

  //! Number of owned nodes in this domain
  /*!
   * For surfaces, this is either 1 or 0
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
   * This is set even if the processor doesn't own the node.
   */
  NodalVars *NodalVarPtr;

  //! This is the index of this surface domain in the list of surface domains
  //! that are located at the current node.
  int Index_NodalSD_;

  LocalNodeIndices *LI_ptr_;

  //! Number of equations defined at the current node.
  /*!
   *  This will be equal to or greater than the number of unknowns assigned to the
   *  surface domain.
   *   This is set even if this processor doesn't own the node.
   */
  int NumNodeEqns;

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



  int  NumDomainEqnsLeft_ ;
  int NumDomainEqnsRight_ ;
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

  //! Variable Vector from left bulk domain
  /*!
   *
   */
  std::vector<double> VarVectorLeftBulkDomain_LastResid_NE;

  //! Variable Vector from left bulk domain
  /*!
   *
   */
  std::vector<double> VarVectorRightBulkDomain_LastResid_NE;

  //! Residual fed into the surface domain before application of BCs
  /*!
   *   This may be used to close balances for cases where Dirichlet conditions
   *   are applied at this surface.
   *
   *   Length = Number of equations defined at this surface node
   */
  std::vector<double> Resid_BeforeSurDomain_NE;


  std::vector<double> DomainResidVector_LastResid_NE;

private:

  //! Call an error exit condition
  /*!
   *
   * @param msg  Message to be printed out
   */
  void
  err(const char *msg);

};
//==================================================================================
//! Specification of a set of simple Dirichlet conditions on a surface domain
/*!
 *  This boundary condition specifies one or more dirichlet conditions on variables
 *  at an interface
 */
class SurBC_Dirichlet : public SurDomain1D
{
public:


  //! Constructor
  /*!
   *
   * @param sdd   Contains the surface domain description.
   */
  SurBC_Dirichlet(m1d::SurfDomainDescription &sdd);

  SurBC_Dirichlet(const SurBC_Dirichlet &r);

  //! Destructor
  virtual
  ~SurBC_Dirichlet();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  SurBC_Dirichlet &
  operator=(const SurBC_Dirichlet &r);

  //! Prepare all of the indices for fast calculation of the residual
  /*!
   *  Here we collect all of the information necessary to
   *  speedily implement SpecFlag_NE and Value_NE within the
   *  residual calculation.
   *  We transfer the information from SDT_Dirichlet structure to
   * this structure for quick processing.
   */
  virtual void
  domain_prep(LocalNodeIndices *li_ptr);

  //! Change the value of an existing dirichlet condition
  /*!
   *  @param vtDir    Variable type that is to be changed
   *  @param newVal   new value of dirichlet condition
   *
   *  @return Number of boundary conditions that are changed.
   */
  virtual int
  changeDirichletConditionValue(VarType vtDir,  double newVal);

  //! Change the boundary condition applied to a variable
  /*!
   *   @param vtDir    Variable type class. Note, general matches are allowed with this parameter
   *   @param BC_Type  Type of the boundary condition
   *   @param value    Value of the dirichlet condition or flux - default 0.0
   *   @param BC_TimeDep BoundaryCondition Pointers for time dependent BC for BC_Tppe = 3,4,5
   *                   (default 0)
   *   @param TimeDep  Function pointer to a function that returns a double given a single parameter (the time).
   *                   Defaults to a NULL pointer.
   *
   *   @return  Returns the number of boundary conditions matched.
   *            A negative number means that an error has been encountered
   */
  int 
  changeBoundaryCondition(VarType vtDir, int BC_Type, double value = 0.0, BoundaryCondition * BC_TimeDep = 0,
			  TimeDepFunctionPtr TimeDep = 0);

  //! Report on the boundary condition applied on the first match to VarType
  /*!
   *   Inputs:
   *   @param time     Current time for evaluating time dependent BC
   *   @param vtDir    Variable type class. Note, general matches are allowed with this parameter
   *   Outputs
   *   @param BC_Type  Type of the boundary condition
   *   @param value    Value of the Dirichlet condition or flux - default 0.0
   *   @param BC_TimeDep BoundaryCondition Pointers for time dependent BC for BC_Type = 3,4,5
   *                   (default 0)
   *   @param TimeDep  Function pointer to a function that returns a double given a single parameter (the time).
   *                   Defaults to a NULL pointer.
   *
   *   @return  Returns the number of boundary conditions matched.
   *            A negative number means that an error has been encountered
   */
  int
  reportBoundaryCondition(double time, const VarType vtDir, int &BC_Type, 
			  double &value,
			  BoundaryCondition * &BC_TimeDep,
			  TimeDepFunctionPtr &TimeDep) const;

  //! Basic function to calculate the residual for the domain.
  /*!
   *  We calculate the additions and/or replacement of the
   *  residual here for the equations that this dirichlet condition
   *  is responsible for.
   *
   * @param res           Output vector containing the residual
   * @param doTimeDependentResid  boolean indicating whether the time
   *                         dependent residual is requested
   * @param soln_ptr     solution vector at which the residual should be
   *                     evaluated
   * @param solnDot_ptr  solution dot vector at which the residual should
   *                     be evaluated.
   *  @param t           time
   *  @param rdelta_t    inverse of delta_t
   */
  virtual void
  residEval(Epetra_Vector &res,
            const bool doTimeDependentResid,
            const Epetra_Vector *soln_ptr,
            const Epetra_Vector *solnDot_ptr,
            const Epetra_Vector *solnOld_ptr,
            const double t,
            const double rdelta_t,
            const ResidEval_Type_Enum residType = Base_ResidEval,
	    const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);

  //! Base class for saving the solution on the domain in an xml node.
  /*!
   *
   * @param oNode                Reference to the XML_Node
   * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
   * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
   * @param t                    time
   *
   * @param duplicateOnAllProcs  If this is true, all processors will include
   *                             the same XML_Node information as proc 0. If
   *                             false, the xml_node info will only exist on proc 0.
   */
  virtual void
  saveDomain(ZZCantera::XML_Node& oNode,
             const Epetra_Vector *soln_GlAll_ptr,
             const Epetra_Vector *solnDot_GlAll_ptr,
             const double t,
             bool duplicateOnAllProcs = false);

  //! Base class for writing the solution on the domain to a logfile.
  /*!
   *
   * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
   * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
   * @param soln_ptr             Pointer to the solution vector
   * @param solnDot_ptr          Pointer to the time derivative of the solution vector
   * @param solnOld_ptr          Pointer to the solution vector at the old time step
   * @param residInternal _ptr   Pointer to the current value of the residual just calculated
   *                             by a special call to the residEval()
   * @param t                    time
   * @param rdelta_t             The inverse of the value of delta_t
   * @param indentSpaces         Indentation that all output should have as a starter
   * @param duplicateOnAllProcs  If this is true, all processors will include
   *                             the same log information as proc 0. If
   *                             false, the loginfo will only exist on proc 0.
   */
  virtual void
  showSolution(const Epetra_Vector *soln_GlAll_ptr,
               const Epetra_Vector *solnDot_GlAll_ptr,
               const Epetra_Vector *soln_ptr,
               const Epetra_Vector *solnDot_ptr,
               const Epetra_Vector *solnOld_ptr,
               const Epetra_Vector_Owned *residInternal_ptr,
               const double t,
               const double rdelta_t,
               int indentSpaces,
               bool duplicateOnAllProcs = false);

  //! Generate the initial conditions
  /*!
   *   For surface dirichlet conditions, we impose the t = 0- condition.
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
  virtual void
  initialConditions(const bool doTimeDependentResid,
                    Epetra_Vector *soln,
                    Epetra_Vector *solnDot,
                    const double t,
                    const double delta_t);

  // ******************************************************************************
  //  Member Data for this boundary condition
  // ******************************************************************************

  //! Number of equations defined at the current node.
  /*!
   *  This will be equal to or greater than the number of unknowns assigned to the
   *  surface domain.
   *   This is set even if this processor doesn't own the node.
   */
  //int NumNodeEqns;

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
  std::vector<BoundaryCondition *> BC_TimeDep_NE;



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

//==================================================================================
} /* End of namespace */
//==================================================================================
#endif /* M1D_SURDOMAIN1D_H_ */
