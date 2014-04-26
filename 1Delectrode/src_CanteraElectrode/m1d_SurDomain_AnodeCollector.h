/**
 * @file m1d_SurDomain_AnodeCollector.h
 *  Basic object to calculate the surface residuals for surface domains.
 */
/*
 *    $Id: m1d_SurDomain_AnodeCollector.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_SURDOMAIN_ANODECOLLECTOR_H_
#define M1D_SURDOMAIN_ANODECOLLECTOR_H_
//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_SurDomain1D.h"

#include "Electrode.h"

#include "Epetra_Vector.h"

using namespace Cantera;

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
class SurDomain_AnodeCollector : public SurBC_Dirichlet
{
public:
  //! Constructor
  /*!
   *
   * @param sdd   Contains the surface domain description.
   */
  SurDomain_AnodeCollector(m1d::SurfDomainDescription &sdd, int problemType);

  //! Copy Constructor
  /*!
   *
   * @param r  Item to be copied
   */
  SurDomain_AnodeCollector(const SurDomain_AnodeCollector &r);

  //! Destructor
  virtual
  ~SurDomain_AnodeCollector();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  SurDomain_AnodeCollector &
  operator=(const SurDomain_AnodeCollector&r);

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
   *  @param residType   Residual evaluation type
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

  //!  Get the voltages for the metal and solution
  /*!
   * @param solnElectrolyte  solution at the current node
   */
  void
  getVoltages(const double * const solnElectrolyte);

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

  //! Pointer to the bulk domain description object
  //! for the porous anode-electrolyte bulk region
  BulkDomainDescription *bedd_;

public:
  //! voltage electrolyte
  double phiElectrolyte_;

  //! voltage electrode
  double phiAnode_;

  //! Voltage at the edge of the anode current collector
  double phiAnodeCC_;

  //! current at the collector
  double icurrCollector_;

};

//==================================================================================
} /* End of namespace */
//==================================================================================
#endif /* M1D_SURDOMAIN1D_H_ */
