/*
 * m1d_BulkDomain1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 * $Id: m1d_Domain1D.h 592 2013-05-13 16:57:58Z hkmoffa $
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

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"

#include "cantera/base/xml.h"

#include <string>

namespace m1d
{

class LocalNodeIndices;

class globalHeatBalVals
{
  public:
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

     virtual ~globalHeatBalVals()
     {

     }

     virtual void zero() 
     {
          totalHeatCapacity = 0.0;
          HeatFluxRight = 0.0;
          HeatFluxLeft = 0.0;
          oldNEnthalpy = 0.0;
          newNEnthalpy = 0.0;
     }

     double totalHeatCapacity;
     double HeatFluxRight;
     double HeatFluxLeft;
     double oldNEnthalpy;
     double newNEnthalpy;
     double moleFluxOut;
     std::vector<double>speciesMoleFluxOut;
     double enthFluxOut;
};


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


  virtual void
  eval_PostSoln(
            const bool doTimeDependentResid,
            const Epetra_Vector *soln_ptr,
            const Epetra_Vector *solnDot_ptr,
            const Epetra_Vector *solnOld_ptr,
            const double t,
            const double rdelta_t);

   virtual void eval_HeatBalance(const int ifunc,
				  const double t,
				  const double deltaT,
				  const Epetra_Vector *soln_ptr,
				  const Epetra_Vector *solnDot_ptr,
				  const Epetra_Vector *solnOld_ptr,
				  struct globalHeatBalVals& dVals);

   virtual void
   eval_SpeciesElemBalance(const int ifunc, const double t, const double deltaT,
	                   const Epetra_Vector *soln_ptr, const Epetra_Vector *solnDot_ptr, const Epetra_Vector *solnOld_ptr,
	                   class globalHeatBalVals& dVals);


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

    //! Get solution parameters specified by text strings
    /*!
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
     *   @param[in]  soln_ptr    Current solution vector (if null, not available)
     *   @param[out] vecInfo     Vector of information returned.
     *
     *   @return  Returns the number of items returned. A value of -1 signifies a failure.
     */
    virtual int
    reportSolutionVector(const std::string& requestID, const int requestType, const Epetra_Vector *soln_ptr,
                         std::vector<double>& vecInfo) const;


  //! Set the underlying state of the system from the solution vector
  /*!
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
   * @param doTimeDependentResid
   * @param soln
   * @param solnDot
   * @param t
   * @param delta_t delta t. If zero then delta_t equals 0.
   * @param t_old
   */
  virtual void
  setStateFromSolution(const bool doTimeDependentResid, const Epetra_Vector_Ghosted *soln, const Epetra_Vector_Ghosted *solnDot,
                       const double t, const double delta_t, const double t_old);

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

  // ------------------------------------------------  DATA ---------------------------------------------------------------

  //! Number of equations associated with this domain,
  //! whether it be a bulk or surface domain
  int NumDomainEqns;

  //! Identifying tag for the domain
  /*!
   *
   */
  std::string m_id;

  //! The type of coordinate system that is used
  /*!
   *  There are two that are envisioned: Rectinear_Coordinates and Cylindrical_Coordinates
   */
  CoordinateSystem_Type_Enum coordinateSystemType_;

  //! CrossSectional Area of the domain
  /*!
   *     This get's copied to the object from the problem statement.
   *     The default is 1 m**2
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

  // The initial state of the battery has some initial stress due to the mfg process. This produces no strain/displacement. 
  // While this is being imput as a scalar, it in principal could be an initial condition of each node, making it a function
  // of position in the battery layers. 

#ifdef MECH_MODEL
  double SolidStressAxialRef_;
#endif 

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
}
#endif /* */
