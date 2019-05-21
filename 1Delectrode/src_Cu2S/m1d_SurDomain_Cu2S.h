/**
 * @file m1d_SurDomain_Cu2S.h
 *  Basic object to calculate the surface residuals for surface domains.
 */
/*
 *    $Id: m1d_SurDomain_Cu2S.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_SURDOMAIN_CU2S_H_
#define M1D_SURDOMAIN_CU2S_H_
//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_ReactingSurDomain.h"
#include "ReactionRate.h"

#include "Epetra_Vector.h"

//======================================================================================================================
namespace m1d
{
// Forward declarations
class NodalVars;
class LocalNodeIndices;

//=====================================================================================================================
//! Specification of a set of boundary conditions on the top of the Cu2S surface
/*!
 *
 */
class Cu2S_TopSurface : public ReactingSurDomain
{
public:
  //! Constructor
  /*!
   *
   * @param sdd   Contains the surface domain description.
   */
  Cu2S_TopSurface(m1d::SurfDomainDescription &sdd, int problemType);

  //! Copy Constructor
  /*!
   *
   * @param r  Item to be copied
   */
  Cu2S_TopSurface(const Cu2S_TopSurface &r);

  //! Destructor
  virtual
  ~Cu2S_TopSurface();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  Cu2S_TopSurface &
  operator=(const Cu2S_TopSurface &r);

  void
  setParams(double x_h2s = 1.41E-7,
            double a1 = 2.71E3,
            double e1 = 6.30,
            double aneg1 = 1.88E3,
            double eneg1 = 6.3,
            double pres_atm = 1.0,
            double temperature = 300.);

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
   * @param soln_ptr     solution vector at which the residual should be
   *                     evaluated
   * @param solnDot_ptr  solution dot vector at which the residual should
   *                     be evaluated.
   *  @param t           time
   *  @param rdelta_t    inverse of delta_t
   */
  virtual void
  residEval(Epetra_Vector &res, const bool doTimeDependentResid,
            const Epetra_Vector* const soln_ptr,
            const Epetra_Vector* const solnDot_ptr,
            const Epetra_Vector* const solnOld_ptr,
            const double t, const double rdelta_t,
            const Zuzax::ResidEval_Type residType = Zuzax::ResidEval_Type::Base_ResidEval,
	    const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve) override;

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
  saveDomain(Zuzax::XML_Node& oNode,
             const Epetra_Vector* const soln_GlAll_ptr,
             const Epetra_Vector* const solnDot_GlAll_ptr,
             const double t,
             bool duplicateOnAllProcs = false) override;

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
               bool duplicateOnAllProcs = false) override;

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
                    Epetra_Vector* const soln,
                    Epetra_Vector* const solnDot,
                    const double t,
                    const double delta_t) override;

  // ******************************************************************************
  //  Member Data for this boundary condition
  // ******************************************************************************
protected:

  //! Mole fraction for H2S gas
  /*!
   *   default is X = 1.0E-8
   */
  double m_X_H2S_g;
  //double m_C_H2S_g;  -> don't think I need this

  //!  Problem 1 - CS Mechanism A1
  /*! Preexponential for the reaction on the gas Cu2S surface
   *     involving the creation of Cu2S lattice sites and Cv
   *     from H2s.
   *       units = cm/sec default = 2.71E5 cm/sec.
   *         but we are using 4.00 E7
   */
  double m_A1;

  //! Activation energy for the Cu2S surface reaction
  /*!
   *    We have fixed this at 6.3 kcal gmol-1
   */
  double m_E1;

  //! Problem 1 - CS Mechanism Aneg11
  /*!
   *    Prexponentialf for the reverse reaction that creates
   *        Cu2S lattice sites from H2S.
   *        units = cm4/mol/s
   *         default = 1.01E5 cm4/mol/s
   *
   *   we will be using  4.27E3
   */
  double m_Aneg1;

  //! Activation energy for the Cu2S reverse surface reaction
  /*!
   *    We have fixed this at 6.3 kcal gmol-1
   */
  double m_Eneg1;

  //! Molar Concentration of Lattice
  /*!
   *  This has units of mol/cm3 and is units of a full lattice
   *  constant (i.e., Cu2S counts as one mol)
   *     3.52E-2 gmol cm-3  - the value of chalcocite
   */
  double m_Conc_Solid;

  double m_k1;

  CanteraLite::ReactionRate * RRTop_;
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
};

//=====================================================================================================================
//! Specification of a set of boundary conditions on the bottom of the Cu2S surface
/*!
 *
 */
class Cu2S_BotSurface : public ReactingSurDomain
{
public:
  //! Constructor
  /*!
   *
   * @param sdd   Contains the surface domain description.
   */
  Cu2S_BotSurface(m1d::SurfDomainDescription &sdd, int problemType);

  //! Copy Constructor
  /*!
   *
   * @param r  Item to be copied
   */
  Cu2S_BotSurface(const Cu2S_BotSurface &r);

  //! Destructor
  virtual
  ~Cu2S_BotSurface();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  Cu2S_BotSurface &
  operator=(const Cu2S_BotSurface &r);

  void
  setParams();

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
            const Zuzax::ResidEval_Type residType = Zuzax::ResidEval_Type::Base_ResidEval,
	    const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve);

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
  saveDomain(Zuzax::XML_Node& oNode,
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
protected:

  CanteraLite::ReactionRate * RRBot_;

  double m_A2;

  double m_E2;

  double m_k1;
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
};

//==================================================================================
} /* End of namespace */
//==================================================================================
#endif /* M1D_SURDOMAIN1D_H_ */
