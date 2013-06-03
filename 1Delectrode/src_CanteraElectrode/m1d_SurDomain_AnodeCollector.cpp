/**
 * @file m1d_SurDomain_Cu2S.cpp
 *  object to calculate the  surface domains in the Cu2S problem
 */

/*
 *  $Id: m1d_SurDomain_AnodeCollector.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_SurDomain_AnodeCollector.h"

#include "m1d_NodalVars.h"
#include "m1d_SurfDomainTypes.h"

#include "m1d_exception.h"
#include "m1d_GlobalIndices.h"
#include "m1d_BulkDomainDescription.h"
#include "m1d_BulkDomain1D.h"

#include "m1d_SDT_AnodeCollector.h"

#include "Epetra_Comm.h"
#include "Epetra_Vector.h"

using namespace std;

//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
SurDomain_AnodeCollector::SurDomain_AnodeCollector(SurfDomainDescription &sdd, int problemType) :
  SurBC_Dirichlet(sdd), bedd_(0), phiElectrolyte_(0.0), phiAnode_(0.0), icurrCollector_(0.0)

{
  //! Determine bedd_
  //       use SurfDomainDescription &SDD_;
  //                    LeftBulk or RightBulk
  bedd_ = SDD_.RightBulk;
  if (!bedd_) {
    bedd_ = SDD_.LeftBulk;
  }
  if (!bedd_) {
    throw m1d_Error("SurDomain_FlatLiSiAnode::SurDomain_FlatLiSiAnode", "Can't find adjoining bulk electrolyte domain");
  }
}
//=====================================================================================================================
SurDomain_AnodeCollector::SurDomain_AnodeCollector(const SurDomain_AnodeCollector &r) :
  SurBC_Dirichlet(r.SDD_), bedd_(0), phiElectrolyte_(0.0), phiAnode_(0.0), icurrCollector_(0.0)
{
  operator=(r);
}
//=====================================================================================================================
// Destructor
SurDomain_AnodeCollector::~SurDomain_AnodeCollector()
{
}
//=====================================================================================================================
// Assignment operator
/*
 * @param r      Object to be copied into the current object
 * @return       Returns a changeable reference to the current object
 */
SurDomain_AnodeCollector &
SurDomain_AnodeCollector::operator=(const SurDomain_AnodeCollector &r)
{
  if (this == &r) {
    return *this;
  }
  SurBC_Dirichlet::operator=(r);

  SpecFlag_NE = r.SpecFlag_NE;
  Value_NE = r.Value_NE;
  bedd_ = r.bedd_;
  phiElectrolyte_ = r.phiElectrolyte_;
  phiAnode_ = r.phiAnode_;
  icurrCollector_ = r.icurrCollector_;

  return *this;
}
//=====================================================================================================================
// Prepare all of the indices for fast calculation of the residual
/*
 *  Here we collect all of the information necessary to
 *  speedily implement SpecFlag_NE and Value_NE within the
 *  residual calculation.
 *  We transfer the information from SDT_Dirichlet structure to
 * this structure for quick processing.
 */
void
SurDomain_AnodeCollector::domain_prep(LocalNodeIndices *li_ptr)
{
  /*
   * First call the parent domain prep to get the node information
   *   -  Index_LcNode
   *   -  NumEqns
   *   - NodalVarPtr
   */
  SurBC_Dirichlet::domain_prep(li_ptr);

}
//=====================================================================================================================
// Basic function to calculate the residual for the domain.
/*
 *  We calculate the additions and/or replacement of the
 *  residual here for the equations that this Dirichlet condition
 *  is responsible for.
 *
 * @param res           Output vector containing the residual
 * @param doTimeDependentResid  boolean indicating whether the time
 *                         dependent residual is requested
 * @param soln         Solution vector at which the residual should be
 *                     evaluated
 */
void
SurDomain_AnodeCollector::residEval(Epetra_Vector &res,
                                    const bool doTimeDependentResid,
                                    const Epetra_Vector *soln_ptr,
                                    const Epetra_Vector *solnDot_ptr,
                                    const Epetra_Vector *solnOld_ptr,
                                    const double t,
                                    const double rdelta_t,
                                    const ResidEval_Type_Enum residType,
				    const Solve_Type_Enum solveType)
{
 residType_Curr_ = residType;
  int ieqn;
  const Epetra_Vector &soln = *soln_ptr;
  incrementCounters(residType);
  /*
   *  Quick return if we don't own the node that the boundary condition
   *  is to be applied on.
   */
  if (NumOwnedNodes == 0) {
    return;
  }
  /*
   *  Figure out the equation start for this node
   *   We start at the start of the equations for this node
   *   because we will be applying dirichlet conditions on the bulk
   *   equations.
   */
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

  /*
   * get the offsets for the BulkDomain and the surface domain.
   */
  int offsetBD = NodalVarPtr->OffsetIndex_BulkDomainEqnStart_BDN[0];

  int EQ_Current_offset_BD = offsetBD + bedd_->EquationIndexStart_EqName[Current_Conservation];
  int EQ_Current_offset_ED = EQ_Current_offset_BD + 1;
  /*
   *  Loop over the equations that the boundary conditions are going to be applied to
   *    -> This takes care of the current surface domain equation
   */
  for (int i = 0; i < NumNodeEqns; i++) {
    if (SpecFlag_NE[i]) {

      ieqn = index_EqnStart + i;
      double solnVal = soln[ieqn];
      double val = Value_NE[i];
      switch (BC_Type_NE[i])
	{
	case 0:
	  /*
	   *  For Dirichlet equations, replace the equation
	   */
	  res[ieqn] = val - solnVal;
	  break;
	case 1:
	  /*
	   *  For flux boundary conditions, subtract from equation indicating a flux out
	   */
	  res[ieqn] -= val;
	  // CAL: WARNING m1d_SurDomain_CathodeCollector.cpp has += val
	  break;
	case 2:
	  /*
	   *  For Dirichlet boundary conditions with oscillation, replace the equation
	   */
	  res[ieqn] = val * TimeDep_NE[i](t) - solnVal;
	  break;
	case 3:
	  /*
	   *  For flux boundary conditions with oscillation, replace the equation
	   */
	  res[ieqn] -= val * TimeDep_NE[i](t);
	  break;
	case 4: // voltage BCconstant
	case 6: // voltage BCsteptable
	case 8: // voltage BClineartable
	  /*
	   *  For time dependent Dirichlet boundary condition using BoundaryCondition class 
	   */
	  res[ieqn] = BC_TimeDep_NE[i]->value(t) - solnVal;
	  break;
	case 5: // current BCconstant
	case 7: // current BCsteptable
	case 9: // current BClineartable
	  /*
	   *  For time dependent flux boundary condition using BoundaryCondition class 
	   */
	  res[ieqn] -= BC_TimeDep_NE[i]->value(t);
	  break;
	default:
	  throw m1d_Error("SurDomain_AnodeCollector::residEval", 
			  "BC_Type_NE[i] 0-9 for Dirichlet, Neumann, and Time Dependence"); 
	}
    }
  }
#ifdef DEBUG_HKM
  if (doTimeDependentResid) {
    if (residType == Base_ResidEval) {
    }
  }
#endif

  getVoltages(&(soln[index_EqnStart]));
  /*
   * get the consistent currents
   */
  BulkDomain1D *bd = bedd_->BulkDomainPtr_;
  icurrCollector_ = bd->DiffFluxLeftBound_LastResid_NE[EQ_Current_offset_ED];

}
//=====================================================================================================================
void
SurDomain_AnodeCollector::getVoltages(const double * const solnElectrolyte)
{
  int indexVS = bedd_->VariableIndexStart_VarName[Voltage];
  phiElectrolyte_ = solnElectrolyte[indexVS];
  phiAnode_ = solnElectrolyte[indexVS + 1];
}
//=====================================================================================================================
// Generate the initial conditions
/*
 *   For surface Dirichlet conditions, we impose the t = 0- condition.
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
void
SurDomain_AnodeCollector::initialConditions(const bool doTimeDependentResid,
                                            Epetra_Vector *soln_ptr,
                                            Epetra_Vector *solnDot,
                                            const double t,
                                            const double delta_t)
{
  int ieqn;
  Epetra_Vector &soln = *soln_ptr;
  /*
   *  Quick return if we don't own the node that the boundary condition
   *  is to be applied on.
   */
  if (NumOwnedNodes == 0) {
    return;
  }

  /*
   *  Figure out the equation start for this node
   *   We start at the start of the equations for this node
   *   because we will be applying Dirichlet conditions on the bulk
   *   equations.
   */
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

  /*
   *  Loop over the equations that the boundary conditions are going to be applied to
   */
  for (int i = 0; i < NumNodeEqns; i++) {
    if (SpecFlag_NE[i]) {
      /*
       *  For Dirichlet equations, replace the solution
       */
      ieqn = index_EqnStart + i;
      double val = Value_NE[i];
      if (BC_Type_NE[i] == 0) {
        soln[ieqn] = val;
      }
      else if (BC_Type_NE[i] == 2) {
	soln[ieqn] = val * TimeDep_NE[i](t);
      }
      else if (BC_Type_NE[i] == 4 || BC_Type_NE[i] == 6 || BC_Type_NE[i] == 8) {
	soln[ieqn] = BC_TimeDep_NE[i]->value(t);
      }
    }
  }
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
