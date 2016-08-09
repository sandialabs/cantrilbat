/**
 * @file m1d_SurDomain_Cu2S.cpp
 *  object to calculate the  surface domains in the Cu2S problem
 */

/*
 *  $Id: m1d_SurDomain_Cu2S.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_SurDomain_Cu2S.h"

#include "m1d_NodalVars.h"
#include "m1d_SurfDomainTypes.h"

#include "m1d_exception.h"
#include "m1d_GlobalIndices.h"

#include "Cu2S_models.h"

#include "Epetra_Comm.h"
#include "Epetra_Vector.h"

#include "m1d_Comm.h"

using namespace std;
using namespace CanteraLite;
//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
//
//       m_A1 = 4.0E7   m_Aneg1 = 4.27E3  orig ->
//       m_A1 = 4.0E5   m_Aneg1 = 4.27E1  next ->
Cu2S_TopSurface::Cu2S_TopSurface(SurfDomainDescription &sdd, int problemType) :
  ReactingSurDomain(sdd), m_X_H2S_g(1.0E-8), m_A1(4.0E5), m_E1(6.3), m_Aneg1(4.27E1), m_Eneg1(6.3),
      m_Conc_Solid(3.52E-2), m_k1(1.0)
{
  if (problemType == 0) {
    RRTop_ = new CanteraLite::ConstantLinearReactionRate(m_k1);
  } else if (problemType == 1) {
    RRTop_ = new CanteraLite::Cu2SReactionRate(m_X_H2S_g, m_A1, m_E1, m_Aneg1, m_Eneg1, 620. / 760., 300.);
  } else if (problemType == 2) {
    RRTop_ = new CanteraLite::ConstantReactionRate(m_k1);
  } else if (problemType == 3) {
    RRTop_ = new CanteraLite::Cu2SReactionRate(m_X_H2S_g, m_A1, m_E1, m_Aneg1, m_Eneg1, 620. / 760., 300.);
  }
}
//=====================================================================================================================
Cu2S_TopSurface::Cu2S_TopSurface(const Cu2S_TopSurface &r) :
  ReactingSurDomain(r.SDD_), m_X_H2S_g(1.0E-8), m_A1(4.0E7), m_E1(6.3), m_Aneg1(4.27E3), m_Eneg1(6.3),
      m_Conc_Solid(3.52E-2), m_k1(1.0)
{
  operator=(r);
}
//=====================================================================================================================
// Destructor
Cu2S_TopSurface::~Cu2S_TopSurface()
{
}
//=====================================================================================================================
// Assignment operator
/*
 * @param r      Object to be copied into the current object
 * @return       Returns a changeable reference to the current object
 */
Cu2S_TopSurface &
Cu2S_TopSurface::operator=(const Cu2S_TopSurface &r)
{
  if (this == &r) {
    return *this;
  }
  ReactingSurDomain::operator=(r);

  SpecFlag_NE = r.SpecFlag_NE;
  Value_NE = r.Value_NE;

  m_X_H2S_g = r.m_X_H2S_g;
  m_A1 = r.m_A1;
  m_E1 = r.m_E1;
  m_Aneg1 = r.m_Aneg1;
  m_Eneg1 = r.m_Eneg1;
  m_Conc_Solid = r.m_Conc_Solid;
  m_k1 = r.m_k1;

  // Do a shallow copy because we don't have a duplicator function
  RRTop_ = r.RRTop_;
  if (RRTop_) {
    throw m1d_Error("Cu2S_TopSurface::operator=()", "incomplete");
  }
  NumNodeEqns = r.NumNodeEqns;
  return *this;
}
//=====================================================================================================================
void
Cu2S_TopSurface::setParams(double x_h2s,
                           double a1,
                           double e1,
                           double aneg1,
                           double eneg1,
                           double pres_atm,
                           double temperature)
{
  Cu2SReactionRate *cRR = dynamic_cast<Cu2SReactionRate *> (RRTop_);
  if (!cRR) {
    throw m1d_Error("Cu2S_TopSurface::setParams", "");
  }
  cRR->m_X_H2S_g = x_h2s;
  cRR->m_A1 = a1;
  cRR->m_E1 = e1;
  cRR->m_Aneg1 = aneg1;
  cRR->m_Eneg1 = eneg1;
  cRR->m_Temperature = temperature;

  double ctot = pres_atm / (82.05 * temperature);
  cRR->m_C_H2S_g = m_X_H2S_g * ctot;
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
Cu2S_TopSurface::domain_prep(LocalNodeIndices *li_ptr)
{
  /*
   * First call the parent domain prep to get the node information
   *   -  Index_LcNode
   *   -  NumEqns
   *   - NodalVarPtr
   */
  SurDomain1D::domain_prep(li_ptr);
  /*
   *  Now figure out:
   *    - SpecFlag_NE[i]
   *    - Value_NE[i];
   */

  /*
   *  Make sure that the SurfDomainType associated with this domain
   *  is a straight surface Dirichlet condition
   */
  SDT_Dirichlet *sdt = dynamic_cast<SDT_Dirichlet *> (&SDD_);
  AssertThrow(sdt, "bad cast");
  /*
   * Zero out the current domain's setup section.
   */
  NumNodeEqns = NodalVarPtr->NumEquations;

  SpecFlag_NE.resize(NumNodeEqns);
  Value_NE.resize(NumNodeEqns);
  for (int j = 0; j < NumNodeEqns; j++) {
    SpecFlag_NE[j] = 0;
    Value_NE[j] = 0.0;
  }
  NumBCs = 0;

  /*
   * Loop over the Dirichlet conditions defined in the SurfDomainType
   * Structure. Transfer the information to this structure for quick
   * processing.
   */
  for (int i = 0; i < sdt->NumConditions; i++) {
    VarType vtDir = sdt->VariableID[i];
    VAR_TYPE vtmDir = vtDir.VariableType;
    VAR_TYPE_SUBNUM stDir = vtDir.VariableSubType;
    // If the subtype is -1, then we apply to all equations of the
    // variable main type
    /*
     *
     */
    if (stDir == -1) {
      for (int j = 0; j < NumNodeEqns; j++) {
        VarType eqt = NodalVarPtr->VariableNameList_EqnNum[j];
        if ((vtmDir == Variable_Type_Any) || (eqt.VariableType == vtmDir)) {
          NumBCs++;
          SpecFlag_NE[j] = 1;
          Value_NE[j] = sdt->Value[i];
        } else {
          SpecFlag_NE[j] = 0;
        }
      }
    } else {
      for (int j = 0; j < NumNodeEqns; j++) {
        VarType vtNode = NodalVarPtr->VariableNameList_EqnNum[j];
        if (vtDir == vtNode) {
          NumBCs++;
          SpecFlag_NE[j] = 1;
          Value_NE[j] = sdt->Value[i];
        } else {
          SpecFlag_NE[j] = 0;
        }
      }
    }
  }
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
Cu2S_TopSurface::residEval(Epetra_Vector &res,
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
   *   Figure out the node to apply the equation on
   *   It is one node -> This can be done using base class level
   */
  // NodalVars* nodeCent = LI_ptr->NodalVars_LcNode[Index_LcNode];

  /*
   *  Figure out the equation start for this node
   *   We start at the start of the equations for this node
   *   because we will be applying dirichlet conditions on the bulk
   *   equations.
   */
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

  /*
   *  Figure out the surface domain on which these boundary conditions
   *  will be applied. There is one controlling one
   */
  //NodalVars *nv = NodalVarPtr;
  /*
   *   Figure out the node to apply the equation on
   *   It is one node -> This can be done using base class level
   */
  // NodalVars* nodeCent = LI_ptr->NodalVars_LcNode[Index_LcNode];


  /*
   *  Get the equation offsets at the node
   */
  int x_offset = NodalVarPtr->Offset_VarType[Displacement_Axial];
  int c_offset = NodalVarPtr->Offset_VarType[Concentration_Species];
  /*
   *   Calculate the state of the system on that controlling surface domain
   */

  /*
   *  Loop over the equations that the boundary conditions are going to be applied to
   */
  for (int i = 0; i < NumNodeEqns; i++) {
    if (SpecFlag_NE[i]) {
      /*
       *  For Dirichlet equations, replace the equation
       */
      ieqn = index_EqnStart + i;
      double solnVal = soln[ieqn];
      double val = Value_NE[i];
      res[ieqn] = val - solnVal;
    }
  }
  /*
   * Boundary conditions on top
   */
  double area_cvb = 1.0;
  if (coordinateSystemType_ == Cylindrical_Coordinates) {
    //area_cvb = Pi * m_cellBound[NumLcCells - 1] * cylinderLength_;
    area_cvb = cylinderLength_;
    exit(-1);
  } else if (coordinateSystemType_ == Spherical_Coordinates) {
    //area_cvb = 2 * Pi * m_cellBound[m_nodes - 1] * m_cellBound[m_nodes - 1];
    area_cvb = 2 * ZZCantera::Pi;
    exit(-1);
  }
  double stoicCoeffC0;
  double stoicCoeffCS;
  double molarReactionRate;

  if (RRTop_) {
    const double *concVec = &(soln[index_EqnStart + c_offset]);
    molarReactionRate = RRTop_->calculateRate(concVec);
    stoicCoeffC0 = RRTop_->getStoichSpec(0);
    stoicCoeffCS = RRTop_->getStoichSpec(1);
  } else {
    stoicCoeffC0 = -1.0;
    stoicCoeffCS = 1.0;
    molarReactionRate = 0.0;
  }
  /*
   * The solid is the product of the reaction
   */
  double solidCreationRate = stoicCoeffCS * molarReactionRate * area_cvb;
  /*
   * The reactant gets used up during the reaction
   */
  res[index_EqnStart + c_offset] -= stoicCoeffC0 * molarReactionRate * area_cvb;
  /*
   *  The node velocity is calculated on the fly
   */
  double nodeVeloc = 0.0;
  if (solnDot_ptr) {
    nodeVeloc = (*solnDot_ptr)[index_EqnStart + x_offset];
    /*
     *  The mesh movement is determined by the total continuity equation
     */
    res[index_EqnStart + x_offset] = m_Conc_Solid * nodeVeloc * area_cvb - solidCreationRate;
  }
  /*
   *  For the steady state equations case, we are counting on the fact that the mesh movement residual
   *  within the volume domain handles the case.
   */

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
Cu2S_TopSurface::saveDomain(Zuzax::XML_Node& oNode,
                            const Epetra_Vector *soln_GLALL_ptr,
                            const Epetra_Vector *solnDot_GLALL_ptr,
                            const double t,
                            bool duplicateOnAllProcs)
{
  // const doublereal* s = soln_GLALL_ptr + loc();
  // Find the number of global equations on this domain, whether it's local or not
  //int numEquationsGb = SDD_.NumEquationsPerNode;
  // Find the global node number of the node where this domain resides
  int locGbNode = SDD_.LocGbNode;

  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  NodalVars *nv = gi->NodalVars_GbNode[locGbNode];
  int eqnStart = nv->EqnStart_GbEqnIndex;
  //XML_Node& inlt = o.addChild("inlet");
  Zuzax::XML_Node& inlt = oNode.addChild("domain");
  int numVar = nv->NumEquations;
  inlt.addAttribute("id", id());
  inlt.addAttribute("points", 1);
  inlt.addAttribute("type", "surface");
  inlt.addAttribute("numVariables", numVar);
  double x0pos = nv->x0NodePos();
  double xpos = nv->xNodePos();
  ZZctml::addFloat(inlt, "X0", x0pos, "", "", Zuzax::Undef, Zuzax::Undef);
  ZZctml::addFloat(inlt, "X", xpos, "", "", Zuzax::Undef, Zuzax::Undef);

  for (int k = 0; k < numVar; k++) {
    double sval = (*soln_GLALL_ptr)[eqnStart + k];
    string nm = nv->VariableName(k);
    VarType vv = nv->VariableNameList_EqnNum[k];
    string type = VarType::VarMainName(vv.VariableType);
    ZZctml::addFloat(inlt, nm, sval, "", "", Zuzax::Undef, Zuzax::Undef);
  }
}
//=====================================================================================================================
// Base class for writing the solution on the domain to a logfile.
/*
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
void
Cu2S_TopSurface::showSolution(const Epetra_Vector *soln_GlAll_ptr,
                              const Epetra_Vector *solnDot_GlAll_ptr,
                              const Epetra_Vector *soln_ptr,
                              const Epetra_Vector *solnDot_ptr,
                              const Epetra_Vector *solnOld_ptr,
                              const Epetra_Vector_Owned *residInternal_ptr,
                              const double t,
                              const double rdelta_t,
                              int indentSpaces,
                              bool duplicateOnAllProcs)
{
  //char buf[132];
  int locGbNode = SDD_.LocGbNode;
 // int mypid = LI_ptr_->Comm_ptr_->MyPID();
 // bool doWrite = !mypid || duplicateOnAllProcs;
  bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
  string indent = "";
  for (int i = 0; i < indentSpaces; i++) {
    indent += " ";
  }
  stream0 ss;
  print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
  const char *ind = indent.c_str();
  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  NodalVars *nv = gi->NodalVars_GbNode[locGbNode];
  int eqnStart = nv->EqnStart_GbEqnIndex;
  //std::vector<VarType> &variableNameList = SDD_.VariableNameList;
  std::vector<VarType> &variableNameListNode = nv->VariableNameList_EqnNum;
  int numVar = nv->NumEquations;
  string sss = id();
  if (doWrite) {
    ss.drawline0(indentSpaces, 80);
    ss.print0("%s  Solution on Surface Domain %10s : Number of variables = %d\n", ind, sss.c_str(), numVar);
    //Zuzax::writelog(buf);
    ss.print0("%s                                           : Number of boundary conditions = %d\n", ind, NumBCs);
    //Zuzax::writelog(buf);
    doublereal x0 = nv->x0NodePos();
    ss.print0("%s                                           : Node %d at pos %g\n", ind, locGbNode, x0);
    //Zuzax::writelog(buf);
    ss.drawline0(indentSpaces, 80);
    ss.print0("%s     VariableName         Value        DirichletCondition\n", ind);
    //Zuzax::writelog(buf);
    ss.drawline0(indentSpaces + 2, 60);
    int jDir = 0;
    for (int k = 0; k < numVar; k++) {
      VarType &vt = variableNameListNode[k];
      string name = vt.VariableName(15);
      double sval = (*soln_GlAll_ptr)[eqnStart + k];
      ss.print0("%s   %-15s   %-10.4E ", ind, name.c_str(), sval);
      //Zuzax::writelog(buf);
      if (SpecFlag_NE[k] != 0) {
        ss.print0(" (Dir %d val = %-10.4E)", jDir, Value_NE[jDir]);
        //Zuzax::writelog(buf);
        jDir++;
      }
      ss.print0("\n");
      //Zuzax::writelog(buf);
    }
    ss.drawline0(indentSpaces + 2, 60);
    ss.drawline0(indentSpaces, 80);
  }
  print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
}
//=====================================================================================================================
// Generate the initial conditions
/*
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
void
Cu2S_TopSurface::initialConditions(const bool doTimeDependentResid,
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
   *   Figure out the node to apply the equation on
   *   It is one node -> This can be done using base class level
   */
  // NodalVars* nodeCent = LI_ptr->NodalVars_LcNode[Index_LcNode];

  /*
   *  Figure out the equation start for this node
   *   We start at the start of the equations for this node
   *   because we will be applying dirichlet conditions on the bulk
   *   equations.
   */
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

  /*
   *  Figure out the surface domain on which these boundary conditions
   *  will be applied. There is one controlling one
   */

  /*
   *   Calculate the state of the system on that controlling surface domain
   */

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
      soln[ieqn] = val;
    }
  }
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================

Cu2S_BotSurface::Cu2S_BotSurface(SurfDomainDescription &sdd, int problemType) :
  ReactingSurDomain(sdd), m_A2(9.6), m_E2(0.0), m_k1(100000.)
{
  if (problemType == 0) {
    // Set the bottom to a zero dirichlet condition
    RRBot_ = new CanteraLite::ConstantLinearReactionRate(m_k1);
  } else if (problemType == 1) {
    RRBot_ = new CanteraLite::Cu2SCuRR(m_A2, m_E2, m_Temperature);
  } else if (problemType == 2) {
    RRBot_ = new CanteraLite::Cu2SCuRR(m_A2, m_E2, m_Temperature);
  } else if (problemType == 3) {
    RRBot_ = new CanteraLite::Cu2SCuRR(m_A2, m_E2, m_Temperature);
  }
}
//=====================================================================================================================
Cu2S_BotSurface::Cu2S_BotSurface(const Cu2S_BotSurface &r) :
  ReactingSurDomain(r.SDD_), m_A2(9.6), m_E2(0.0), m_k1(100000.)
{
  Cu2S_BotSurface::operator=(r);
}
//=====================================================================================================================
// Destructor
Cu2S_BotSurface::~Cu2S_BotSurface()
{
}
//=====================================================================================================================
// Assignment operator
/*
 * @param r      Object to be copied into the current object
 * @return       Returns a changeable reference to the current object
 */
Cu2S_BotSurface &
Cu2S_BotSurface::operator=(const Cu2S_BotSurface &r)
{
  if (this == &r) {
    return *this;
  }
  ReactingSurDomain::operator=(r);

  SpecFlag_NE = r.SpecFlag_NE;
  Value_NE = r.Value_NE;

  m_A2 = r.m_A2;
  m_E2 = r.m_E2;

  // Do a shallow copy because we don't have a duplicator function
  RRBot_ = r.RRBot_;
  if (RRBot_) {
    throw m1d_Error("Cu2S_BotSurface::operator=()", "incomplete");
  }
  NumNodeEqns = r.NumNodeEqns;
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
Cu2S_BotSurface::domain_prep(LocalNodeIndices *li_ptr)
{
  /*
   * First call the parent domain prep to get the node information
   *   -  Index_LcNode
   *   -  NumEqns
   *   - NodalVarPtr
   */
  SurDomain1D::domain_prep(li_ptr);
  /*
   *  Now figure out:
   *    - SpecFlag_NE[i]
   *    - Value_NE[i];
   */

  /*
   *  Make sure that the SurfDomainType associated with this domain
   *  is a straight surface Dirichlet condition
   */
  SDT_Dirichlet *sdt = dynamic_cast<SDT_Dirichlet *> (&SDD_);
  AssertThrow(sdt, "bad cast");
  /*
   * Zero out the current domain's setup section.
   */
  NumNodeEqns = NodalVarPtr->NumEquations;

  SpecFlag_NE.resize(NumNodeEqns);
  Value_NE.resize(NumNodeEqns);
  for (int j = 0; j < NumNodeEqns; j++) {
    SpecFlag_NE[j] = 0;
    Value_NE[j] = 0.0;
  }
  NumBCs = 0;

  /*
   * Loop over the Dirichlet conditions defined in the SurfDomainType
   * Structure. Transfer the information to this structure for quick
   * processing.
   */
  for (int i = 0; i < sdt->NumConditions; i++) {
    VarType vtDir = sdt->VariableID[i];
    VAR_TYPE vtmDir = vtDir.VariableType;
    VAR_TYPE_SUBNUM stDir = vtDir.VariableSubType;
    // If the subtype is -1, then we apply to all equations of the
    // variable main type
    /*
     *
     */
    if (stDir == -1) {
      for (int j = 0; j < NumNodeEqns; j++) {
        VarType eqt = NodalVarPtr->VariableNameList_EqnNum[j];
        if ((vtmDir == Variable_Type_Any) || (eqt.VariableType == vtmDir)) {
          NumBCs++;
          SpecFlag_NE[j] = 1;
          Value_NE[j] = sdt->Value[i];
        } else {
          SpecFlag_NE[j] = 0;
        }
      }
    } else {
      for (int j = 0; j < NumNodeEqns; j++) {
        VarType vtNode = NodalVarPtr->VariableNameList_EqnNum[j];
        if (vtDir == vtNode) {
          NumBCs++;
          SpecFlag_NE[j] = 1;
          Value_NE[j] = sdt->Value[i];
        } else {
          SpecFlag_NE[j] = 0;
        }
      }
    }
  }
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
Cu2S_BotSurface::residEval(Epetra_Vector &res,
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
   *   Figure out the node to apply the equation on
   *   It is one node -> This can be done using base class level
   */
  // NodalVars* nodeCent = LI_ptr->NodalVars_LcNode[Index_LcNode];

  /*
   *  Figure out the equation start for this node
   *   We start at the start of the equations for this node
   *   because we will be applying dirichlet conditions on the bulk
   *   equations.
   */
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

  /*
   *  Figure out the surface domain on which these boundary conditions
   *  will be applied. There is one controlling one
   */
  //NodalVars *nv = NodalVarPtr;
  /*
   *   Figure out the node to apply the equation on
   *   It is one node -> This can be done using base class level
   */
  // NodalVars* nodeCent = LI_ptr->NodalVars_LcNode[Index_LcNode];


  /*
   *  Get the equation offsets at the node
   */
  int x_offset = NodalVarPtr->Offset_VarType[Displacement_Axial];
  int c_offset = NodalVarPtr->Offset_VarType[Concentration_Species];
  /*
   *   Calculate the state of the system on that controlling surface domain
   */

  /*
   *  Loop over the equations that the boundary conditions are going to be applied to
   */
  for (int i = 0; i < NumNodeEqns; i++) {
    if (SpecFlag_NE[i]) {
      /*
       *  For Dirichlet equations, replace the equation
       */
      ieqn = index_EqnStart + i;
      double solnVal = soln[ieqn];
      double val = Value_NE[i];
      res[ieqn] = val - solnVal;
    }
  }
  /*
   * Boundary conditions on top
   */
  double area_cvb = 1.0;
  if (coordinateSystemType_ == Cylindrical_Coordinates) {
    //area_cvb = Pi * m_cellBound[m_nodes - 1] * cylinderLength_;
    area_cvb = cylinderLength_;
    // the whole program needs to get fixed for different coordinate systems
    exit(-1);
  } else if (coordinateSystemType_ == Spherical_Coordinates) {
    //area_cvb = 2 * Pi * m_cellBound[m_nodes - 1] * m_cellBound[m_nodes - 1];
    area_cvb = 2 * Zuzax::Pi;
    // the whole program needs to get fixed for different coordinate systems
    exit(-1);
  }
  double stoicCoeffC0;
  double molarReactionRate;

  if (RRBot_) {
    const double *concVec = &(soln[index_EqnStart + c_offset]);
    molarReactionRate = RRBot_->calculateRate(concVec);
    stoicCoeffC0 = RRBot_->getStoichSpec(0);
  } else {
    stoicCoeffC0 = -1.0;
    molarReactionRate = 0.0;
  }
  /*
   * The reactant gets used up during the reaction
   */
  res[index_EqnStart + c_offset] -= stoicCoeffC0 * molarReactionRate * area_cvb;
  /*
   *  The mesh movement is determined by the total continuity equation
   */
  res[index_EqnStart + x_offset] = (*soln_ptr)[index_EqnStart + x_offset];
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
Cu2S_BotSurface::saveDomain(Zuzax::XML_Node& oNode,
                            const Epetra_Vector *soln_GLALL_ptr,
                            const Epetra_Vector *solnDot_GLALL_ptr,
                            const double t,
                            bool duplicateOnAllProcs)
{
  // const doublereal* s = soln_GLALL_ptr + loc();
  // Find the number of global equations on this domain, whether it's local or not
  //int numEquationsGb = SDD_.NumEquationsPerNode;
  // Find the global node number of the node where this domain resides
  int locGbNode = SDD_.LocGbNode;

  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  NodalVars *nv = gi->NodalVars_GbNode[locGbNode];
  int eqnStart = nv->EqnStart_GbEqnIndex;
  //XML_Node& inlt = o.addChild("inlet");
  Zuzax::XML_Node& inlt = oNode.addChild("domain");
  int numVar = nv->NumEquations;
  inlt.addAttribute("id", id());
  inlt.addAttribute("points", 1);
  inlt.addAttribute("type", "surface");
  inlt.addAttribute("numVariables", numVar);
  double x0pos = nv->x0NodePos();
  double xpos = nv->xNodePos();
  ZZctml::addFloat(inlt, "X0", x0pos, "", "", Zuzax::Undef, Zuzax::Undef);
  ZZctml::addFloat(inlt, "X", xpos, "", "", Zuzax::Undef, Zuzax::Undef);

  for (int k = 0; k < numVar; k++) {
    double sval = (*soln_GLALL_ptr)[eqnStart + k];
    string nm = nv->VariableName(k);
    VarType vv = nv->VariableNameList_EqnNum[k];
    string type = VarType::VarMainName(vv.VariableType);
    ZZctml::addFloat(inlt, nm, sval, "", "", Zuzax::Undef, Zuzax::Undef);
  }
}
//=====================================================================================================================
// Base class for writing the solution on the domain to a logfile.
/*
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
void
Cu2S_BotSurface::showSolution(const Epetra_Vector *soln_GlAll_ptr,
                              const Epetra_Vector *solnDot_GlAll_ptr,
                              const Epetra_Vector *soln_ptr,
                              const Epetra_Vector *solnDot_ptr,
                              const Epetra_Vector *solnOld_ptr,
	                      const Epetra_Vector_Owned *residInternal_ptr,
                              const double t,
                              const double rdelta_t,
                              int indentSpaces,
                              bool duplicateOnAllProcs)
{
  //char buf[132];
  //buf[0] = '\0';
  int locGbNode = SDD_.LocGbNode;
  //int mypid = LI_ptr_->Comm_ptr_->MyPID();
 // bool doWrite = !mypid || duplicateOnAllProcs;
  bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
  string indent = "";
  for (int i = 0; i < indentSpaces; i++) {
    indent += " ";
  }
  const char *ind = indent.c_str();
  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  NodalVars *nv = gi->NodalVars_GbNode[locGbNode];
  int eqnStart = nv->EqnStart_GbEqnIndex;
  //std::vector<VarType> &variableNameList = SDD_.VariableNameList;
  std::vector<VarType> &variableNameListNode = nv->VariableNameList_EqnNum;
  int numVar = nv->NumEquations;
  string sss = id();
  stream0 ss;
  print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
  if (doWrite) {
    ss.drawline0(indentSpaces, 80);
    ss.print0("%s  Solution on Surface Domain %10s : Number of variables = %d\n", ind, sss.c_str(), numVar);
 
    //Zuzax::writelog(buf);
    ss.print0("%s                                           : Number of boundary conditions = %d\n", ind, NumBCs);
    //Zuzax::writelog(buf);
    doublereal x0 = nv->x0NodePos();
    ss.print0("%s                                           : Node %d at pos %g\n", ind, locGbNode, x0);
    //Zuzax::writelog(buf);
    ss.drawline0(indentSpaces, 80);
    ss.print0("%s     VariableName         Value        DirichletCondition\n", ind);
    //Zuzax::writelog(buf);
    ss.drawline0(indentSpaces + 2, 60);
    int jDir = 0;
    for (int k = 0; k < numVar; k++) {
      VarType &vt = variableNameListNode[k];
      std::string name = vt.VariableName(20);
      double sval = (*soln_GlAll_ptr)[eqnStart + k];
      ss.print0("%s   %-20s   %-10.4E ", ind, name.c_str(), sval);
      //Zuzax::writelog(buf);
      if (SpecFlag_NE[k] != 0) {
        ss.print0(" (Dir %d val = %-10.4E)", jDir, Value_NE[jDir]);
        //Zuzax::writelog(buf);
        jDir++;
      }
      ss.print0("\n");
      //Zuzax::writelog(buf);
    }
    ss.drawline0(indentSpaces + 2, 60);
    ss.drawline0(indentSpaces, 80);
  }
  print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
}
//=====================================================================================================================
// Generate the initial conditions
/*
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
void
Cu2S_BotSurface::initialConditions(const bool doTimeDependentResid,
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
   *   Figure out the node to apply the equation on
   *   It is one node -> This can be done using base class level
   */
  // NodalVars* nodeCent = LI_ptr->NodalVars_LcNode[Index_LcNode];

  /*
   *  Figure out the equation start for this node
   *   We start at the start of the equations for this node
   *   because we will be applying dirichlet conditions on the bulk
   *   equations.
   */
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

  /*
   *  Figure out the surface domain on which these boundary conditions
   *  will be applied. There is one controlling one
   */

  /*
   *   Calculate the state of the system on that controlling surface domain
   */

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
      soln[ieqn] = val;
    }
  }
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
