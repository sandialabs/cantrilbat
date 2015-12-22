/*
 * m1d_BulkDomain1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#ifndef M1D_BULKDOMAIN1D_H_
#define M1D_BULKDOMAIN1D_H_
//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_BulkDomainDescription.h"
#include "m1d_Domain1D.h"

#include "Epetra_Vector.h"

namespace m1d
{
class LocalNodeIndices;

//! Base class for solving residuals for bulk domains
/*!
 * There is a 1 to 1 mapping between the local control volume indexing and
 * the Local Consecutive Ordering indexing
 *
 * There is a 1 to 1 mapping between the global control volume indexing
 * and the Global node number indexing that is given by a single offset.
 *
 *
 *
 */
class BulkDomain1D : public Domain1D
{

public:

    //! Constructor
    /*!
     * @param bdd   Contains the bulk domain description.
     */
    BulkDomain1D(m1d::BulkDomainDescription &bdd);

    //! Copy Constructor
    /*!
     * @param r      Object to be copied into the current object
     */
    BulkDomain1D(const BulkDomain1D &r);

    //! Destructor
    virtual ~BulkDomain1D();

    //! Assignment operator
    /*!
     * @param r      Object to be copied into the current object
     * @return       Returns a changeable reference to the current object
     */
    BulkDomain1D &
    operator=(const BulkDomain1D &r);

    //! Returns the identifying tag for the domain
    virtual std::string
    id() const;

    //! Prepare all of the indices for fast calculation of the residual
    /*!
     *  Ok, at this point, we will have figured out the number of equations
     *  to be calculated at each node point. The object NodalVars will have
     *  been fully formed.
     *
     *  We use this to figure out what local node numbers/ cell numbers are
     *  needed and to set up indices for their efficient calling.
     *
     *  Child objects of this one will normally call this routine in a
     *  recursive fashion.
     */
    virtual void
    domain_prep(LocalNodeIndices *li_ptr);

    //! Basic function to calculate the residual for the current domain.
    /*!
     *  This base class is used just for volumetric domains.
     *
     *  All residual terms are written with the following sign convention
     *  based on keeping the time derivative term positive.
     *
     *       res = dcdt - dc2 /dx2 - src = 0
     *
     *  Boundary conditions for residuals are written in such a way that
     *  bulk domains may be split and then recombined with no loss in accuracy.
     *
     * @param res  Output vector containing the residual
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

    //! Generate the initial conditions
    /*!
     *   The basic algorithm is to loop over the volume domains.
     *   Then, we loop over the surface domains. Within the domains, we
     *   use the virtual function structure to go from general to the more
     *   specific direction (i.e., parent to child calling).
     *
     *   In this routine, we make sure that if there are displacement unknowns
     *   then the initial solution holds the X0NodePos values.
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
    saveDomain(Cantera::XML_Node& oNode,
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
     *    <domain id="BulkDomain1D_0" numVariables="6" points="10" type="bulk">
     *      <grid_data>
     *        <floatArray size="10" title="X0" type="length" units="m">
     *          0.000000000000000E+00,   8.748888888888889E-05,   1.749777777777778E-04,
     *          2.624666666666667E-04,   3.499555555555555E-04,   4.374444444444444E-04,
     *          5.249333333333334E-04,   6.124222222222222E-04,   6.999111111111111E-04,
     *          7.873999999999999E-04
     *        </floatArray>
     *     </domain>
     *  </simulation>
     *
     * @param domainNode           Reference to the XML_Node, named domain, to read the solution from
     * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
     * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
     */
    virtual void 
      readDomain(const Cantera::XML_Node& domainNode,
		 Epetra_Vector* const soln_GlAll_ptr,
		 Epetra_Vector* const solnDot_GlAll_ptr, double globalTimeRead);

    //! Fill up a vector indicating whether an unknown is an algebraic condition or not
    /*!
     *  The source of the information is lcated in the BulkDomainDescription for the domain.
     */
    virtual void fillIsAlgebraic(Epetra_IntVector & isAlgebraic);
 

    //!  Fill the vector isArithmeticScaled with the values from the DomainDescription
    /*!
     * @param isArithmeticScaled  Epetra_IntVector to be filled with the IsArithmeticScaled values
     */
    virtual void fillIsArithmeticScaled(Epetra_IntVector & isArithmeticScaled);



    // Method for writing the header for the surface domain to a tecplot file.
    /*
     * Only proc0 will write tecplot files.
     */
    virtual void writeSolutionTecplotHeader();

    // Method for writing the solution on the surface domain to a tecplot file.
    /*
     * Only proc0 will write tecplot files.
     *
     * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
     * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
     * @param t                    time
     *
     */
    virtual void writeSolutionTecplot(const Epetra_Vector *soln_GlAll_ptr,
				      const Epetra_Vector *solnDot_GlAll_ptr,
				      const double t );




    //! Base class for writing the solution on the domain to a logfile.
    /*!
     *
     * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
     * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
     * @param soln_ptr             Pointer to the solution vector
     * @param solnDot_ptr          Pointer to the time derivative of the solution vector
     * @param solnOld_ptr          Pointer to the solution vector at the old time step
     * @param residInternal_ptr    Pointer to the current value of the residual just calculated
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
		 const Epetra_Vector *residInternal_ptr,
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
    
    //! Base class for writing a solution int vector, not the solution, on the domain to a logfile.
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
			  const Epetra_IntVector *solnIntVector_GlAll_ptr,
			  const Epetra_IntVector *solnIntVector_ptr,
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
     * @param solnDot_ptr          Pointer to the time derivative of the solution vector
     * @param solnOld_ptr          Pointer to the solution vector at the old time step
     * @param residInternal_ptr    Pointer to the current value of the residual just calculated
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
		     const Epetra_Vector *residInternal_ptr,
		     const double t,
		     const double rdelta_t,
		     int indentSpaces,
		     bool duplicateOnAllProcs = false);
    
    //! Get parameters specified by text strings
    /*!
     *  (virtual from Domain1D)
     *   @param[in]  paramID   String name for the item to be requested
     *   @param[out] paramVal    Vector of information returned.
     *
     *   @return  Returns the number of items returned. A value of -1 signifies a failure.
     *
     */
    virtual int 
    reportSolutionParam(const std::string& paramID, double* const paramVal) const;

    //! Get vectors of solution quantities requested by text strings
    /*!
     *
     *   @param[in]  requestID   String name for the item to be requested
     *   @param[in]  requestType Type of the request
     *                      0    solution variable
     *                      1    other variable
     *   @param[out] vecInfo     Vector of information returned.
     *
     *   @return  Returns the number of items returned. A value of -1 signifies a failure.
     */
    virtual int
    reportSolutionVector(const std::string& requestID, const int requestType, const Epetra_Vector* soln_ptr,
                         std::vector<double>& vecInfo) const;

#ifdef MECH_MODEL
    //! Get the local value of the stress, from the solution vector, 
    //! or a reference value if not part of the solution. 

    double getPointStress(const NodalVars * const nv,
			  const doublereal* const solutionPoint) const;

#endif

    //! Get the local value of the temperature at a node or control volume interface
    //! given the local solution vector at that point
    /*!
     *   This function checks to see if the temperature is part of the solution 
     *   vector. If it is not, it returns the TemperatureReference_ value. If 
     *   it is, it looks up the index into the solution vector and then returns 
     *   the value.
     *   
     *   @return Returns the temperature in Kelvin
     */
    double getPointTemperature(const NodalVars* const nv, 
			       const doublereal* const solutionPoint) const;
    
    //! Get the local value of the total pressure at a node or control volume interface
    //! given the local solution vector at that point
    /*!
     *   This function checks to see if the pressure is part of the solution 
     *   vector. If it is not, it returns the PressureReference_ value. If 
     *   it is, it looks up the index into the solution vector and then returns 
     *   the value of the total pressure based on the local condition
     *
     *     @return Returns the total pressure in Pascals
     */
    double getPointPressure(const NodalVars* const nv, 
			    const doublereal* const solutionPoint) const;

    // ===========================================================================

    //! Light description of what this domain is about
    m1d::BulkDomainDescription &BDD_;

    //! Number of owned nodes in this domain
    int NumOwnedNodes;

    //! First global node that is owned by this processor on this bulk domain
    int FirstOwnedGbNode;

    //! Last Global node that is owned by this processor on this bulk domain
    int LastOwnedGbNode;

    //! Number of owned cells on this domain
    /*!
     *   This number is used to loop over the number of cells in the
     *   domain from left to right. This loop is called LCO for local
     *   cell order.
     */
    int NumLcCells;

    //! True if this processor owns the left-most node  of this domain
    bool IOwnLeft;

    //! True if this processor owns the right-most node of this domain
    bool IOwnRight;

    //! Index of the local node that corresponds to the current cell number
    /*!
     *  This is the index of the local node number at the center of the cell
     */
    std::vector<int> Index_DiagLcNode_LCO;

    //! Index of the left node that corresponds to the current cell number
    /*!
     * This is the index of the local node that is to the left of the current
     * cell. Note, if this cell is the first cell in the bulk domain, then this
     * is set to -1, even if there is a node to the left from another domain.
     */
    std::vector<int> Index_LeftLcNode_LCO;

    //! Index of the right node that corresponds to the current cell number
    /*!
     * This is the index of the local node that is to the right of the current
     * cell. Note, if this cell is the last cell in the bulk domain, then this
     * is set to -1, even if there is a node to the right from another domain.
     */
    std::vector<int> Index_RightLcNode_LCO;

    //!  boolean indicating whether displacements are part of the solution vector
    /*!
     *   If displacements are part of the solution vector the position of nodes
     *   is given by Xpos = X0pos + d.
     *   If displacements are not part of the solution vector, the position of nodes
     *   is given by Xpos = X0pos
     */
    bool MeshInSolnVector;

    //! Diffusive fluxes at the left boundary of the domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the diffusive flux calculated
     *  at the left boundary during the last residual calculation
     *
     *  This quantities are useful (I think) for the specification of boundary
     *  conditions for multi domains.
     *
     *  They are also useful for the specification of global balances, when Dirichlet conditions are set
     *  on some variables. Then, we seek the specification of fluxes which preserve global conservation laws.
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> DiffFluxLeftBound_LastResid_NE;

    //! Diffusive Fluxes at the right boundary of the domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the diffusive flux calculated
     *  at the right boundary during the last residual calculation
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> DiffFluxRightBound_LastResid_NE;


    //! Total Fluxes at the left boundary of the domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the diffusive flux calculated
     *  at the left boundary during the last residual calculation
     *
     *  This quantities are useful (I think) for the specification of boundary
     *  conditions for multi domains.
     *
     *  They are also useful for the specification of global balances, when Dirichlet conditions are set
     *  on some variables. Then, we seek the specification of fluxes which preserve global conservation laws.
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> TotalFluxLeftBound_LastResid_NE;

    //! Total Fluxes at the right boundary of the domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the total flux calculated
     *  at the right boundary during the last residual calculation
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> TotalFluxRightBound_LastResid_NE;

    //! Value of the variable as seen from the left side of domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the  left side variable value
     *  at the left boundary during the last residual calculation
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> VarVectorLeftBound_LastResid_NE;

    //! Value of the variable as seen from the right side of domain from the last residual calculation
    /*!
     *  This is a temporary variable that holds the right side variable value
     *  at the right boundary during the last residual calculation
     *
     *  Length = number of equations defined on this domain
     */
    std::vector<double> VarVectorRightBound_LastResid_NE;

    //! Value of Residual contributions at the left side of the domain from this domain only
    /*!
     *  this is the contributions to the residual from this domain.
     *  This is only calculated for show Solution cases.
     *
     *  Length = number of equations defined at the left boundary node
     */
    std::vector<double> DomainResidVectorLeftBound_LastResid_NE;

    //! Value of Residual contributions at the right side of the domain from this domain only
    /*!
     *  this is the contributions to the residual from this domain.
     *  This is only calculated for show Solution cases.
     *
     *  Length = number of equations defined at the left boundary node
     */
    std::vector<double> DomainResidVectorRightBound_LastResid_NE;


    //! Pointer to the local node indices for this processor
    LocalNodeIndices *LI_ptr_;

private:
    //! local error routine
    /*!
     *  @param msg error message
     */
    void err(const char *msg);

};

}

#endif /* M1D_BulkDomain1D_H_ */
