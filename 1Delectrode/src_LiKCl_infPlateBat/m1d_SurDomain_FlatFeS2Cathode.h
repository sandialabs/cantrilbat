/**
 * @file m1d_SurDomain_FlatFeS2Cathode.h
 *
 */
/*
 *    $Id: m1d_SurDomain_FlatFeS2Cathode.h 363 2012-08-22 03:37:42Z hkmoffa $
 */

#ifndef M1D_SURDOMAIN_FLATFE2SCATHODE_H_
#define M1D_SURDOMAIN_FLATFE2SCATHODE_H_
//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_ReactingSurDomain.h"

#include "Electrode.h"

#include "Epetra_Vector.h"

//======================================================================================================================
namespace m1d
{
// Forward declarations
class NodalVars;
class LocalNodeIndices;
class BulkDomainDescription;

//=====================================================================================================================
//! Specification of a set of boundary conditions on the top of the Cu2S surface
/*!
 *
 */
class SurDomain_FlatFeS2Cathode : public ReactingSurDomain
{
public:
  //! Constructor
  /*!
   *
   * @param sdd   Contains the surface domain description.
   */
  SurDomain_FlatFeS2Cathode(m1d::SurfDomainDescription &sdd, int problemType);

  //! Copy Constructor
  /*!
   *
   * @param r  Item to be copied
   */
  SurDomain_FlatFeS2Cathode(const SurDomain_FlatFeS2Cathode &r);

  //! Destructor
  virtual
  ~SurDomain_FlatFeS2Cathode();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  SurDomain_FlatFeS2Cathode &
  operator=(const SurDomain_FlatFeS2Cathode &r);

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

  //! Basic function to calculate the residual for the domain.
  /*!
   *  We calculate the additions and/or replacement of the
   *  residual here for the equations that this dirichlet condition
   *  is responsible for.
   *
   * @param res           Output vector containing the residual
   * @param doTimeDependentResid  boolean indicating whether the time
   *                         dependent residual is requested
   * @param soln_ptr     Solution vector at which the residual should be evaluated
   * @param solnDot_ptr  Pointer to the time-derivative of the solution vector
   * @param solnOld_ptr  Pointer to the solution vector at the old time step
   * @param t           time
   * @param rdelta_t    inverse of delta_t
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

  //! utility routine to update the objects used to calculate quantities at the surface
  /*!
   *
   * @param soln_ptr     solution vector at which the residual should be
   *                     evaluated
   * @param t            time
   * @param residType    Residual evaluation type
   */
  virtual void
  updateDependencies(const Epetra_Vector *soln_ptr, const double t, const ResidEval_Type_Enum residType =
      Base_ResidEval);
  void
  getMFElectrolyte_soln(const double * const solnBulk);

  /*
   *  Get the voltages for the metal and solution
   */
  void
  getVoltages(const double * const solnElectrolyte, const double * const solnSolid);
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

  //! Base class for writing the solution on the domain to a logfile.
  /*!
   *
   * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
   * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
   * @param soln_ptr             Pointer to the solution vector
   * @param solnDot_ptr          Pointer to the time-derivative of the solution vector
   * @param solndOld_ptr         Pointer to the solution vector at the old time step
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

  // Method for writing the header for the surface domain to a tecplot file.
  /*
   * Only proc0 will write tecplot files.
   */
  void writeSolutionTecplotHeader();

  // Method for writing the solution on the surface domain to a tecplot file.
  /*
   * Only proc0 will write tecplot files.
   *
   * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
   * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
   * @param t                    time
   *
   */
  void writeSolutionTecplot(const Epetra_Vector *soln_GlAll_ptr,
			    const Epetra_Vector *solnDot_GlAll_ptr,
			    const double t );

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

  // ******************************************************************************
  //  Member Data for this boundary condition
  // ******************************************************************************
protected:

  //! Pointer to the electrode object
  /*!
   *   We do not own the electrode object
   */
  Cantera::Electrode *ElectrodeC_;

  //! Thermodynamics object for the liquid electrolyte
  /*!
   *  This is a shallow pointer
   */
  Cantera::ThermoPhase *electrolyteThermo_;

  //! Pointer to the bulk domain description object
  //! for the electrolyte
  BulkDomainDescription *bedd_;

  //! mole fraction Solution
  std::vector<double> mfElectrolyte_Soln;

  //! mole fraction Thermo
  std::vector<double> mfElectrolyte_Thermo;

  //! voltage electrolyte
  double phiElectrolyte_;

  //! voltage electrode
  double phiCathode_;

  //CanteraLite::ReactionRate * RRTop_;
  //! Number of equations defined at the current node.
  /*!
   * This is set even if this processor doesn't own the node.
   */
  int NumNodeEqns;

  //! Boolean flag indicating which variables at the node are being specified
  //! with Dirichlet conditions
  /*!
   *   Vector has length equal to the number of equations defined at the node
   */
  std::vector<int> SpecFlag_NE;

  //! Value of the variable
  /*!
   *   Vector has length equal to the number of equations defined at the node
   */
  std::vector<double> Value_NE;

  std::vector<double> electrodeSpeciesProdRates_;

  std::vector<double> phaseMoleFlux_;

  double surfaceArea_;

  double concTot_Curr_;

  //! Type of the boundary condition specified on the cathode
  /*!
   *   0 specify the voltage
   *   1 specify the current
   */
  int voltageVarBCType_;

  //! Specified current in the cathode
  /*!
   *  This is actually the current from the cathode into the electrolyte.
   *  Therefore, during a normal discharge operation of the battery, this will be a
   *  negative quantity.
   *
   *   Note, this is only relevant when voltageVarBCType_ = 1
   */
  double icurrCathodeSpecified_;

  friend class FlatBatteryResidEval;

};

//==================================================================================
} /* End of namespace */
//==================================================================================
#endif /* M1D_SURDOMAIN1D_H_ */
