/*
 * m1d_BulkDomain1D.cpp
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#include "m1d_Domain1D.h"

#include "m1d_DomainDescription.h"
#include "m1d_exception.h"

//#include "cantera/base/ctml.h"
#include "m1d_globals.h"
#include "m1d_ProblemStatement.h"

#include "mdp_stringUtils.h"


#include <cstdio>
#include <cstdlib>

using namespace std;
using namespace m1d;

namespace m1d
{
//=====================================================================================================================
// Constructor
Domain1D::Domain1D() :
    NumDomainEqns(0), 
    coordinateSystemType_(Rectilinear_Coordinates),
    crossSectionalArea_(1.0), 
    cylinderLength_(1.0),
    TemperatureReference_(298.15), 
    PressureReference_(1.01325e5),
    energyEquationProbType_(0),
    solidMechanicsProbType_(0),
    porosityEquationProbType_(  Porosity_EqnType_Status::None ),
    residType_Curr_(Base_ResidEval), 
    counterResBaseCalcs_(0), 
    counterJacBaseCalcs_(0),
    counterJacDeltaCalcs_(0), 
    counterResShowSolutionCalcs_(0)
{
    coordinateSystemType_ = PSinput_ptr->coordinateSystemType_;
    crossSectionalArea_ = PSinput_ptr->crossSectionalArea_;
    cylinderLength_ = PSinput_ptr->cylinderLength_;
    TemperatureReference_ = PSinput_ptr->TemperatureReference_;
    PressureReference_ = PSinput_ptr->PressureReference_;
    energyEquationProbType_ = PSinput_ptr->Energy_equation_prob_type_;
#ifdef MECH_MODEL
    solidMechanicsProbType_ = PSinput_ptr->Solid_Mechanics_prob_type_;
#else
    if (PSinput_ptr->Solid_Mechanics_prob_type_ != 0) {
        throw m1d_Error("Domain1D::Domain1D()", "MECH_MODEL define not set, but Solid_Mechanics_prob_type_ nonzero, " 
                        + mdpUtil::int2str(PSinput_ptr->Solid_Mechanics_prob_type_));
    }
#endif
}
//=====================================================================================================================
Domain1D::Domain1D(const Domain1D &r) :
    NumDomainEqns(0), 
    coordinateSystemType_(Rectilinear_Coordinates),
    crossSectionalArea_(1.0), 
    cylinderLength_(1.0),
    TemperatureReference_(298.15), 
    PressureReference_(1.01325e5),
    energyEquationProbType_(0),
    solidMechanicsProbType_(0),
    residType_Curr_(Base_ResidEval), 
    counterResBaseCalcs_(0), 
    counterJacBaseCalcs_(0),
    counterJacDeltaCalcs_(0),
    counterResShowSolutionCalcs_(0)
{
  *this = r;
}
//=====================================================================================================================
Domain1D::~Domain1D()
{
}
//=====================================================================================================================
Domain1D &
Domain1D::operator=(const Domain1D &r)
{
  if (this == &r) {
    return *this;
  }

  NumDomainEqns                   = r.NumDomainEqns;
  coordinateSystemType_           = r.coordinateSystemType_;
  crossSectionalArea_             = r.crossSectionalArea_;
  cylinderLength_                 = r.cylinderLength_;
  residType_Curr_ = r.residType_Curr_;
  counterResBaseCalcs_ = r.counterResBaseCalcs_;
  counterJacBaseCalcs_ = r.counterJacBaseCalcs_;
  counterJacDeltaCalcs_ = r.counterJacDeltaCalcs_;
  counterResShowSolutionCalcs_ = r.counterResShowSolutionCalcs_;
  TemperatureReference_ = r.TemperatureReference_;
  PressureReference_ = r.PressureReference_;
  energyEquationProbType_ = r.energyEquationProbType_;
  solidMechanicsProbType_ = r.solidMechanicsProbType_;
  return *this;
}
//=====================================================================================================================
/*
 * Specify an identifying tag for this domain.
 */
void
Domain1D::setID(const std::string& s)
{
  m_id = s;
}
//=====================================================================================================================
/*
 * Return the identifying tag for this domain.
 */
std::string
Domain1D::id() const
{
  err("id()");
  return string("");
}
//=====================================================================================================================
// Prepare all of the indices for fast calculation of the residual
void
Domain1D::domain_prep(m1d::LocalNodeIndices *li_ptr)
{
  err("domain_prep()");
}
//=====================================================================================================================
// Basic function to calculate the residual for the current domain.
/*
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
 *  @param solnOld_ptr  Pointer to the solution vector at the old time step
 *  @param t           time
 *  @param rdelta_t    inverse of delta_t
 */
void
Domain1D::residEval(Epetra_Vector &res,
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
  err("residEval()");
}
//=====================================================================================================================
void
Domain1D::eval_PostSoln(
                    const bool doTimeDependentResid,
                    const Epetra_Vector *soln_ptr,
                    const Epetra_Vector *solnDot_ptr,
                    const Epetra_Vector *solnOld_ptr,
                    const double t,
                    const double rdelta_t)
{
}
//=====================================================================================================================
void
Domain1D::eval_HeatBalance(const int ifunc,
			  const double t,
			  const double deltaT,
			  const Epetra_Vector *soln_ptr,
			  const Epetra_Vector *solnDot_ptr,
			  const Epetra_Vector *solnOld_ptr,
			  struct globalHeatBalVals& dVals)
{
}
//=====================================================================================================================
void
Domain1D::eval_SpeciesElemBalance(const int ifunc,
			  const double t,
			  const double deltaT,
			  const Epetra_Vector *soln_ptr,
			  const Epetra_Vector *solnDot_ptr,
			  const Epetra_Vector *solnOld_ptr,
			  struct globalHeatBalVals& dVals)
{
}
//=====================================================================================================================
// Utility function to calculate quantities before the main residual routine.
/*
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
void
Domain1D::residEval_PreCalc(const bool doTimeDependentResid,
			    const Epetra_Vector *soln_ptr,
			    const Epetra_Vector *solnDot_ptr,
			    const Epetra_Vector *solnOld_ptr,
			    const double t,
			    const double rdelta_t,
			    const ResidEval_Type_Enum residType,
			    const Solve_Type_Enum solveType)
{
}
  //=====================================================================================================================
  // Auxiliary function to calculate the residual for the current domain.
  /*
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
  void
  Domain1D::residEval_PostCalc(Epetra_Vector &res,
			       const bool doTimeDependentResid,
			       const Epetra_Vector *soln_ptr,
			       const Epetra_Vector *solnDot_ptr,
			       const Epetra_Vector *solnOld_ptr,
			       const double t,
			       const double rdelta_t,
			       const ResidEval_Type_Enum residType,
			       const Solve_Type_Enum solveType)
  {
  }
//===================================================================================================================
// Function that gets called at end or the start of every time step
/*
 *  This function provides a hook for a residual that gets called whenever a
 *  time step is accepted. The call is made with the current time as the time
 *  that is accepted. The old time may be obtained from t and rdelta_t_accepted.
 *
 *   @param  doTimeDependentResid  This is true if we are solving a time dependent
 *                                 problem.
 *   @param  soln_ptr              Solution value at the current time
 *   @param  solnDot_ptr           derivative of the solution at the current time.
 *   @param  solnOld_ptr           Solution value at the old time step, n-1
 *   @param  t                     current time to be accepted, n
 *   @param  t_old                 previous time step value
 */ 
void
Domain1D::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector *soln_ptr,
                              const Epetra_Vector *solnDot_ptr, const Epetra_Vector *solnOld_ptr,
                              const double t, const double t_old) 
{

}
//===================================================================================================================
// Revert the Residual object's conditions to the conditions at the start of the global time step
/*
 *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init_init.
 *  We get rid of the pendingIntegratedFlags_ flag here as well.
 */
void
Domain1D::revertToInitialGlobalTime()
{
    // There is nothing to do in the base case
}
//=====================================================================================================================
// Base class for saving the solution on the domain in an xml node.
/*
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
void
Domain1D::saveDomain(Cantera::XML_Node& oNode,
                     const Epetra_Vector *soln_GLALL_ptr,
                     const Epetra_Vector *solnDot_GLALL_ptr,
                     const double t,
                     bool duplicateOnAllProcs)
{
  err("saveDomain()");
}
//====================================================================================================================
void
Domain1D::readSimulation(const Cantera::XML_Node& simulationNode,
			 Epetra_Vector * const soln_GLALL_ptr,
			 Epetra_Vector * const solnDot_GLALL_ptr)
{
    string ida = id();
    Cantera::XML_Node* domainNode = simulationNode.findNameID("domain", ida);
    if (!domainNode) {
	throw m1d_Error("Domain1D::readSimulation()", "cant find domain node" + ida);
    }
    readDomain(*domainNode, soln_GLALL_ptr, solnDot_GLALL_ptr, -1.0);
}
//====================================================================================================================
void
Domain1D::readDomain(const Cantera::XML_Node& oNode,
                     Epetra_Vector * const soln_GLALL_ptr,
                     Epetra_Vector * const solnDot_GLALL_ptr, double globalTimeRead)
{
  err("readDomain()");
}
//====================================================================================================================
  //! Method for writing the header for the surface domain to a tecplot file.
  /*
   * Only proc0 will write tecplot files.
   */
void Domain1D::writeSolutionTecplotHeader() {
  err("writeSolutionTecplotHeader()");
}
//====================================================================================================================
//! Method for writing the solution on the surface domain to a tecplot file.
/*
 * Only proc0 will write tecplot files.
 *
 * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
 * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
 * @param t                    time
 *
 */
void Domain1D::writeSolutionTecplot(const Epetra_Vector *soln_GlAll_ptr,
				    const Epetra_Vector *solnDot_GlAll_ptr,
				    const double t ){
  err("writeSolutionTecplot()");
}
//====================================================================================================================
// Base class for writing the solution on the domain to a logfile.
/*
 *
 * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
 * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
 * @param soln_ptr             Pointer to the solution vector
 * @param solnDot_ptr          Pointer to the time-derivative of the solution vector
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
Domain1D::showSolution(const Epetra_Vector *soln_GlAll_ptr,
                       const Epetra_Vector *solnDot_GlAll_ptr,
                       const Epetra_Vector *soln_ptr,
                       const Epetra_Vector *solnDot_ptr,
                       const Epetra_Vector *solnOld_ptr,
                       const Epetra_Vector_Owned *residInternal_ptr,
                       const double t,
                       const double rdelta_t,
                       const int indentSpaces,
                       bool duplicateOnAllProcs)
{
  err("showSolution()");
}
//====================================================================================================================
  // Base class for writing a solution vector, not the solution, on the domain to a logfile.
  /*
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
  void
  Domain1D::showSolutionVector(std::string& solnVecName,
			       const Epetra_Vector *solnVector_GlAll_ptr,
			       const Epetra_Vector *solnVector_ptr,
			       const double t,
			       const double rdelta_t,
			       int indentSpaces,
			       bool duplicateOnAllProcs,
			       FILE *of) 
  {
    err("showSolutionVector()");
  }
//====================================================================================================================
  // Base class for writing an int solution vector, not the solution, on the domain to a logfile.
  /*
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
  void
  Domain1D::showSolutionIntVector(std::string& solnVecName,
			       const Epetra_IntVector *solnIntVector_GlAll_ptr,
			       const Epetra_IntVector *solnIntVector_ptr,
			       const double t,
			       const double rdelta_t,
			       int indentSpaces,
			       bool duplicateOnAllProcs,
			       FILE *of) 
  {
    err("showSolutionIntVector()");
  }
//==================================================================================================================================
int 
Domain1D::reportSolutionParam(const std::string& paramID, double* const paramVal) const
{
    paramVal[0] = 0.0;
    return -1;
}
//==================================================================================================================================
int 
Domain1D::reportSolutionVector(const std::string& requestID, const int requestType, const Epetra_Vector *soln_ptr,
                               std::vector<double>& vecInfo) const
{
    vecInfo.clear();
    return -1;
}
//==================================================================================================================================
// Set the underlying state of the system from the solution vector
/*
 *   Note this is an important routine for the speed of the solution.
 *   It would be great if we could supply just exactly what is changing here.
 *   This routine is always called at the beginning of the residual evaluation process.
 *
 *   This is a natural place to put any precalculations of nodal quantities that
 *   may be needed by the residual before its calculation.
 *
 *   Also, this routine is called with rdelta_t = 0. This implies that a step isn't being taken. However, the
 *   the initial conditions must be propagated.
 *
 * @param doTimeDependentResid
 * @param soln
 * @param solnDot
 * @param t
 * @param delta_t 
 @ @param t_old
 */
void
Domain1D::setStateFromSolution(const bool doTimeDependentResid, const Epetra_Vector_Ghosted *soln, 
			       const Epetra_Vector_Ghosted *solnDot, const double t, const double delta_t, const double t_old)
{
}
//====================================================================================================================
// Generate the initial conditions
/*
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
void
Domain1D::initialConditions(const bool doTimeDependentResid,
                            Epetra_Vector *soln,
                            Epetra_Vector *solnDot,
                            const double t,
                            const double delta_t)
{

}
//=====================================================================================================================
void Domain1D:: calcDeltaSolnVariables(const double t, const Epetra_Vector& soln,
				       const Epetra_Vector* solnDot_ptr, Epetra_Vector& deltaSoln,
                                       const Epetra_Vector* const atolVector_ptr,
				       const Solve_Type_Enum solveType,
				       const  Epetra_Vector* solnWeights)
{
  printf("Domain1D: function not implemented\n");
  exit(-1);
}
//=====================================================================================================================
void Domain1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
			     Epetra_Vector_Ghosted & atolVector,
			     const Epetra_Vector_Ghosted * const atolV)
{
  printf("Domain1D: function not implemented\n");
  exit(-1);
}
//================================================================================================================
void
Domain1D::setAtolDeltaDamping(double atolDefault, double relcoeff,  const Epetra_Vector_Ghosted & soln, 
			      Epetra_Vector_Ghosted & atolDeltaDamping,
			      const Epetra_Vector_Ghosted * const atolV)
{
  printf("Domain1D: function not implemented\n");
  exit(-1);
}
//================================================================================================================
void 
Domain1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln,
				const Epetra_Vector_Ghosted & solnDot,
				Epetra_Vector_Ghosted & atolVector_DAEInit,
				const Epetra_Vector_Ghosted * const atolV)
{
  printf("Domain1D: function not implemented\n");
  exit(-1);
}
//================================================================================================================
void
Domain1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, 
				      const Epetra_Vector_Ghosted & soln,
				      const Epetra_Vector_Ghosted & solnDot,
				      Epetra_Vector_Ghosted & atolDeltaDamping,
				      const Epetra_Vector_Ghosted * const atolV)
{
  printf("Domain1D: function not implemented\n");
  exit(-1);
}
//=====================================================================================================================
void
Domain1D::incrementCounters(const ResidEval_Type_Enum residType)
{
  if (residType == Base_ResidEval) {
    counterResBaseCalcs_++;
  } else if (residType == JacBase_ResidEval) {
    counterJacBaseCalcs_++;
  } else if (residType == Base_ShowSolution) {
    counterResShowSolutionCalcs_++;
  } else {
    counterJacDeltaCalcs_++;
  }
}
//=====================================================================================================================
void
Domain1D::err(const char *msg) const
{
  printf("BulkDomain1D: function not implemented: %s\n", msg);
  exit(-1);
}
//=====================================================================================================================
}
//=====================================================================================================================
