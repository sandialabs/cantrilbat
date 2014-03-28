/*
 * m1d_BulkDomain1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 * $Id: m1d_Domain1D.h 592 2013-05-13 16:57:58Z hkmoffa $
 */

#ifndef M1D_BSDOMAIN1D_H_
#define M1D_BSDOMAIN1D_H_
//! This is a heavyweight base class that provides the function
//! evaluation for a single domain whether a surface or a bulk.

#include "m1d_DomainDescription.h"
#include "m1d_ProblemResidEval.h"

#include "Epetra_Vector.h"

#include "cantera/base/global.h"

namespace m1d
{
class LocalNodeIndices;
//! Base class for solving residuals for bulk and surface domains
/*!
 *
 *
 * There is a 1 to 1 mapping between the local control volume indexing and
 * the Local Consecutive Ordering indexing
 *
 * There is a 1 to 1 mapping between the global control volume indexing
 * and the Global node number indexing that is given by a single offset.
 *
 *
 *
 */
class Domain1D
{

public:

  //! Constructor
  /*!
   *
   * @param bdd   Contains the bulk domain description.
   */
  Domain1D();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  Domain1D(const Domain1D &r);

  //!  Destructor
  virtual
  ~Domain1D();

  //! Assignment operator
  /*!
   *  @param r  object to be copied
   */
  Domain1D &
  operator=(const Domain1D &r);

  //! Specify an identifying tag for this domain.
  /*!
   *  The identifying tag for the domain is the name of the domain. It will appear on all
   *  output files.
   *
   * @param s Name of the domain
   */
  virtual void
  setID(const std::string& s);

  //! Returns the identifying tag for the domain
  virtual std::string
  id() const;

  //! Prepare all of the indices for fast calculation of the residual
  /*!
   *  @param li_ptr   Pointer to the LocalNodeIndices Structure that contains information
   *                  about how the mesh is layed out within this domain and other domains
   *                  in the problem.
   */
  virtual void
  domain_prep(m1d::LocalNodeIndices *li_ptr);

  //! Basic function to calculate the residual for the current domain.
  /*!
   *  This base class is used just for volumetric domains.
   *
   *  All residual terms are written with the following sign convention
   *  based on keeping the time derivative term positive.
   *
   *       res = dcdt - dc2 /dx2 - src = 0
   *
   * @param res  Output vector containing the residual
   * @param doTimeDependentResid  boolean indicating whether the time
   *                         dependent residual is requested
   * @param soln_ptr     solution vector at which the residual should be
   *                     evaluated
   * @param solnDot_ptr  solution dot vector at which the residual should
   *                     be evaluated.
   * @param solnOld_ptr  Pointer to the solution vector at the old time step
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


  //! Utility function to calculate quantities before the main residual routine.
  /*!
   *  This is used for a loop over nodes. All calculated quantities must be internally storred.
   *
   *  Currently this is called during the residual evalultion of the problem.
   * 
   * @param res  Output vector containing the residual
   * @param doTimeDependentResid  boolean indicating whether the time
   *                         dependent residual is requested
   * @param soln_ptr     solution vector at which the residual should be
   *                     evaluated
   * @param solnDot_ptr  solution dot vector at which the residual should
   *                     be evaluated.
   * @param solnOld_ptr  Pointer to the solution vector at the old time step
   *  @param t           time
   *  @param rdelta_t    inverse of delta_t
   */
  virtual void 
  residEval_PreCalc(const bool doTimeDependentResid,
		    const Epetra_Vector *soln_ptr,
		    const Epetra_Vector *solnDot_ptr,
		    const Epetra_Vector *solnOld_ptr,
		    const double t,
		    const double rdelta_t,
		    const ResidEval_Type_Enum residType,
		    const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);

  //! Auxiliary function to calculate the residual for the current domain.
  /*!
   *  By default this function does nothing.
   *
   * @param res  Output vector containing the residual
   * @param doTimeDependentResid  boolean indicating whether the time
   *                         dependent residual is requested
   * @param soln_ptr     solution vector at which the residual should be
   *                     evaluated
   * @param solnDot_ptr  solution dot vector at which the residual should
   *                     be evaluated.
   * @param solnOld_ptr  Pointer to the solution vector at the old time step
   *  @param t           time
   *  @param rdelta_t    inverse of delta_t
   */
  virtual void
  residEval_PostCalc(Epetra_Vector &res,
		     const bool doTimeDependentResid,
		     const Epetra_Vector *soln_ptr,
		     const Epetra_Vector *solnDot_ptr,
		     const Epetra_Vector *solnOld_ptr,
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
   *  After this call interrogation, of the previous time step's results will not be valid.
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
   *   @param  t                     Current time to be accepted, n
   *   @param  t_old                 Previous time step value, t_old may be equal to t, 
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
   *
   * @param simulationNode       Reference to the XML_Node named simulation
   * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
   * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
   */
  virtual void
  readSimulation(const Cantera::XML_Node& simulationNode,
		 Epetra_Vector * const soln_GlAll_ptr,
		 Epetra_Vector * const solnDot_GlAll_ptr);


  //! Base Class for reading the solution from the saved file
  /*!
   *
   * @param domainNode          Reference to the XML_Node to read the solution from
   * @param soln_GLALL_ptr      Pointer to the Global-All solution vector
   * @param solnDot_ptr         Pointer to the time derivative of the Global-All solution vector
   *
   */
  virtual void 
  readDomain(const Cantera::XML_Node& domainNode,
             Epetra_Vector * const soln_GlAll_ptr,
             Epetra_Vector * const solnDot_GlAll_ptr);


  //! Method for writing the header for the surface domain to a tecplot file.
  /**
   * Only proc0 will write tecplot files.
   */
  virtual void writeSolutionTecplotHeader();

  //! Method for writing the solution on the surface domain to a tecplot file.
  /**
   * Only proc0 will write tecplot files.
   *
   * @param soln_GlAll_ptr       Pointer to the Global-All solution vector
   * @param solnDot_GlAll_ptr    Pointer to the time derivative of the Global-All solution vector
   * @param t                    time
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
   * @param solnVecName           String name of the solution vector
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

  //! Base class for writing an int solution vector, not the solution, on the domain to a logfile.
  /*!
   * @param solnVecName           String name of the solution vector
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


  //! Generate the initial conditions
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


  void
  incrementCounters(const ResidEval_Type_Enum residType);

  //! Number of equations associated with this domain,
  //! whether it be a bulk or surface domain
  int NumDomainEqns;

  //! Identifying tag for the domain
  /*!
   *
   */
  std::string m_id;

  //! Coordinate system of the domain
  /*!
   *   0 = cartesian
   *   1 = cylindrical  - the axial dimension is r in cylindrical coordinates
   *   2 = spherical    -the axial dimension is r in spherical coordinates
   */
  int CoordinateSystem_;

  //! CrossSectional Area of the domain
  /*!
   *     The default is 1 m**2
   */
  double CellArea_;

  //! Reference Temperature (Kelvin)
  /*!
   *  For each domain, we have a reference temperature. This temperature will be used for property
   *  evaluation as the default temperature within the domain whenever there isn't another source for 
   *  the value of the temperature.
   *     The default is 298.15 Kelvin
   */
  double TemperatureReference_;

  //! Reference Pressure  (pascal)
  /*!
   *  For each domain, we have a reference thermodynamic pressure. This pressure will be used for property
   *  evaluation as the default pressure within the domain whenever there isn't another source for 
   *  the value of the thermodynamic pressure.
   *    The default is 1.01325E5 Pascal = OneAtm
   */
  double PressureReference_;

protected:

  //! Current value of the residual type
  ResidEval_Type_Enum residType_Curr_;

private:
  //! Local error routine
  /*!
   *
   * @param msg Meesage indicating where the error message is originating from
   */
  void
  err(const char *msg) const;

public:

  //! Counter for the total number of base residual calculations undertaken
  int counterResBaseCalcs_;

  //! Counter for the total number of base Jacobian residual calculations undertaken
  int counterJacBaseCalcs_;  

  //! Counter for the total number of Jacobian delta residual calculations undertaken
  int counterJacDeltaCalcs_;

  //! Counter for the total number of show residual calculations undertaken
  int counterResShowSolutionCalcs_;
};
}
#endif /* */
