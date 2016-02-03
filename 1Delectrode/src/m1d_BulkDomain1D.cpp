/**
 * m1d_BulkDomain1D.cpp
 */
/*
 *  $Id: m1d_BulkDomain1D.cpp 504 2013-01-07 22:32:48Z hkmoffa $
 */

//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.
#include "m1d_BulkDomain1D.h"

#include "m1d_LocalNodeIndices.h"
#include "m1d_DomainDescription.h"
#include "m1d_GlobalIndices.h"
#include "m1d_NodalVars.h"

#include "m1d_exception.h"
#include "m1d_Comm.h"

#include "Epetra_Comm.h"

#include "mdp_stringUtils.h"

#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"

#include "stdio.h"
#include "stdlib.h"

using namespace std;
using namespace m1d;

namespace m1d
{

//======================================================================================================================
BulkDomain1D::BulkDomain1D(m1d::BulkDomainDescription &bdd) :
  Domain1D(), BDD_(bdd), NumOwnedNodes(0), FirstOwnedGbNode(-1), LastOwnedGbNode(-1), IOwnLeft(false),
      IOwnRight(false), MeshInSolnVector(false), LI_ptr_(0)
{

}
//=====================================================================================
BulkDomain1D::BulkDomain1D(const BulkDomain1D &r) :
  Domain1D(), BDD_(r.BDD_), NumOwnedNodes(0), FirstOwnedGbNode(-1), LastOwnedGbNode(-1), NumLcCells(-1),
      IOwnLeft(false), IOwnRight(false), MeshInSolnVector(false), LI_ptr_(0)
{
  BulkDomain1D::operator=(r);
}
//=====================================================================================
BulkDomain1D::~BulkDomain1D()
{
}
//=====================================================================================
BulkDomain1D &
BulkDomain1D::operator=(const BulkDomain1D &r)
{
  if (this == &r) {
    return *this;
  }
  Domain1D::operator=(r);

  /*
   *   We do a shallow copy here until we know better how to do this
   */
  BDD_ = r.BDD_;
  NumOwnedNodes = r.NumOwnedNodes;
  FirstOwnedGbNode = r.FirstOwnedGbNode;
  LastOwnedGbNode = r.LastOwnedGbNode;
  NumLcCells = r.NumLcCells;
  IOwnLeft = r.IOwnLeft;
  IOwnRight = r.IOwnRight;
  Index_DiagLcNode_LCO = r.Index_DiagLcNode_LCO;
  Index_LeftLcNode_LCO = r.Index_LeftLcNode_LCO;
  Index_RightLcNode_LCO = r.Index_RightLcNode_LCO;
  MeshInSolnVector = r.MeshInSolnVector;
  DiffFluxLeftBound_LastResid_NE = r.DiffFluxLeftBound_LastResid_NE;
  DiffFluxRightBound_LastResid_NE = r.DiffFluxRightBound_LastResid_NE;
  TotalFluxLeftBound_LastResid_NE = r.TotalFluxLeftBound_LastResid_NE;
  TotalFluxRightBound_LastResid_NE = r.TotalFluxRightBound_LastResid_NE;
  VarVectorLeftBound_LastResid_NE = r.VarVectorLeftBound_LastResid_NE;
  VarVectorRightBound_LastResid_NE = r.VarVectorRightBound_LastResid_NE;

  LI_ptr_ = r.LI_ptr_;

  return *this;
}
//=====================================================================================================================
/*
 * Return the identifying tag for this domain.
 */
std::string
BulkDomain1D::id() const
{
  int id = BDD_.ID();
  if (m_id != "") {
    return m_id;
  } else {
    return std::string("BulkDomain1D_") + Cantera::int2str(id);
  }
}
//=====================================================================================================================
// Basic function to calculate the residual for the domain.
void
BulkDomain1D::residEval(Epetra_Vector &res,
                        const bool doTimeDependentResid,
                        const Epetra_Vector *soln_ptr,
                        const Epetra_Vector *solnDot_ptr,
                        const Epetra_Vector *solnOld_ptr,
                        const double t,
                        const double rdelta_t,
                        ResidEval_Type_Enum residType,
			const Solve_Type_Enum solveType)
{
  residType_Curr_ = residType;
  err("residEval()");
}
//=====================================================================================================================
void
BulkDomain1D::domain_prep(LocalNodeIndices *li_ptr)
{
  LI_ptr_ = li_ptr;
  /*
   * The first thing to do is to find the intersection of this domain with the
   * owned nodes on this processor.
   */

  // Find the global node value of the left boundary of this domain
  int firstGbNodeDomain = BDD_.FirstGbNode;

  //  find the Global node of the right boundary of this domain
  int lastGbNodeDomain = BDD_.LastGbNode;

  int leftOwnedLcNode = 0;
  int leftOwnedGbNode = LI_ptr_->IndexGbNode_LcNode[leftOwnedLcNode];

  // Find the node on the processor
  int totalOwned_nodes = LI_ptr_->NumOwnedLcNodes;
  int rightOwnedLcNode = totalOwned_nodes - 1;
  int rightOwnedGbNode = LI_ptr_->IndexGbNode_LcNode[rightOwnedLcNode];

  // Do the left intersection
  FirstOwnedGbNode = leftOwnedGbNode;
  int firstOwnedLcNode = 0;
  if (firstGbNodeDomain > FirstOwnedGbNode) {

    if (rightOwnedLcNode < firstGbNodeDomain) {
      // We are here when there is no intersection of this processor
      // with this bulk domain
      NumOwnedNodes = 0;
      FirstOwnedGbNode = -1;
      LastOwnedGbNode = -1;
      NumLcCells = 0;
      IOwnLeft = false;
      IOwnRight = false;
      return;
    }
    FirstOwnedGbNode = firstGbNodeDomain;
    firstOwnedLcNode = LI_ptr_->GbNodeToLcNode(FirstOwnedGbNode);
    AssertThrow(firstOwnedLcNode >= 0, "");
  }

  // Do the right intersection - assume first that it's the last node on the processor
  LastOwnedGbNode = rightOwnedGbNode;
  int lastOwnedLcNode = rightOwnedLcNode;
  if (lastGbNodeDomain < LastOwnedGbNode) {
    if (lastGbNodeDomain < FirstOwnedGbNode) {
      // We are here when there is no intersection of this processor
      // with this bulk domain
      NumOwnedNodes = 0;
      FirstOwnedGbNode = -1;
      LastOwnedGbNode = -1;
      NumLcCells = 0;
      IOwnLeft = false;
      IOwnRight = false;
      return;
    }

    LastOwnedGbNode = lastGbNodeDomain;
    lastOwnedLcNode = LI_ptr_->GbNodeToLcNode(LastOwnedGbNode);
    AssertThrow(lastOwnedLcNode >= 0, "");
  }

  // Determine how many owned nodes of the bulk domain are on the processor
  NumOwnedNodes = lastOwnedLcNode - firstOwnedLcNode + 1;
  AssertThrow(NumOwnedNodes > 0, "NumOwnedNodes > 0");
  AssertThrow(NumOwnedNodes == (LastOwnedGbNode - FirstOwnedGbNode + 1), "");

  // Determine the number of owned cells in this domain on this processor
  NumLcCells = NumOwnedNodes;

  // True if this processor owns the left-most node of this domain
  if (firstGbNodeDomain == FirstOwnedGbNode) {
    IOwnLeft = true;
  } else {
    IOwnLeft = false;
  }

  // True if this processor owns the right-most node of this domain
  if (lastGbNodeDomain == LastOwnedGbNode) {
    IOwnRight = true;
  } else {
    IOwnRight = false;
  }

  Index_DiagLcNode_LCO.resize(NumLcCells);
  Index_LeftLcNode_LCO.resize(NumLcCells);
  Index_RightLcNode_LCO.resize(NumLcCells);
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int inode = firstOwnedLcNode + iCell;
    Index_DiagLcNode_LCO[iCell] = inode;
    if (iCell == 0) {
      if (IOwnLeft) {
        Index_LeftLcNode_LCO[iCell] = -1;
      } else {
        Index_LeftLcNode_LCO[iCell] = LI_ptr_->IDLeftLcNode_LcNode[inode];
      }
    } else {
      Index_LeftLcNode_LCO[iCell] = LI_ptr_->IDLeftLcNode_LcNode[inode];
    }

    if (iCell == (NumLcCells - 1)) {
      if (IOwnRight) {
        Index_RightLcNode_LCO[iCell] = -1;
      } else {
        Index_RightLcNode_LCO[iCell] = LI_ptr_->IDRightLcNode_LcNode[inode];
      }
    } else {
      Index_RightLcNode_LCO[iCell] = LI_ptr_->IDRightLcNode_LcNode[inode];
    }
  }

  // Download the number of equations in this domain
  NumDomainEqns = BDD_.NumEquationsPerNode;

  
  // Find the number of equations at the first node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;

  int iGbNode = BDD_.FirstGbNode;
  NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
  size_t numEqFirst = nv->NumEquations;
  iGbNode = BDD_.LastGbNode;
  nv = gi->NodalVars_GbNode[iGbNode];
  size_t numEqLast = nv->NumEquations;

  DomainResidVectorLeftBound_LastResid_NE.resize(numEqFirst, 0.0);
  DomainResidVectorRightBound_LastResid_NE.resize(numEqLast, 0.0);
 
  //  Size the flux vectors at the ends of the domains
  DiffFluxLeftBound_LastResid_NE.resize( numEqFirst, 0.0);
  TotalFluxLeftBound_LastResid_NE.resize(numEqFirst, 0.0);
  VarVectorLeftBound_LastResid_NE.resize(numEqFirst, 0.0);
  DiffFluxRightBound_LastResid_NE.resize(numEqLast, 0.0);
  TotalFluxRightBound_LastResid_NE.resize(numEqLast, 0.0);
  VarVectorRightBound_LastResid_NE.resize(numEqLast, 0.0);



}
//=====================================================================================================================
// Generate the initial conditions
/*
 *   The basic algorithm is to loop over the volume domains.
 *   Then, we loop over the surface domains. Within the domains, we
 *   use the virtual function structure to go from general to the more
 *   specific direction (i.e., parent to child calling).
 *
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
//
//   HKM This function should be filled in with a capability to initialize constant values for 
//       field variables, at least!
void
BulkDomain1D::initialConditions(const bool doTimeDependentResid,
                                Epetra_Vector *soln,
                                Epetra_Vector *solnDot,
                                const double t,
                                const double delta_t)
{
  Domain1D::initialConditions(doTimeDependentResid, soln, solnDot, t, delta_t);
}
//=====================================================================================================================
static void
drawline(int sp, int ll)
{
  for (int i = 0; i < sp; i++) {
    Cantera::writelog(" ");
  }
  for (int i = 0; i < ll; i++) {
    Cantera::writelog("-");
  }
  Cantera::writelog("\n");
}

//=====================================================================================================================
static void
drawline0(stream0 &ss, int sp, int ll)
{
  for (int i = 0; i < sp; i++) {
    ss.fprintf0only(" ");
  }
  for (int i = 0; i < ll; i++) {
    ss.fprintf0only("-");
  }  
  ss.fprintf0only("\n");   
}


//=====================================================================================================================
// saving the solution on the domain in an xml node.
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
BulkDomain1D::saveDomain(Cantera::XML_Node& oNode,
                         const Epetra_Vector *soln_GLALL_ptr,
                         const Epetra_Vector *solnDot_GLALL_ptr,
                         const double t,
                         bool duplicateOnAllProcs)
{
  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;

  // Add a child for this domain
  Cantera::XML_Node& bdom = oNode.addChild("domain");

  // Number of equations per node
  int numEquationsPerNode = BDD_.NumEquationsPerNode;

  // Vector containing the variable names as they appear in the solution vector
  std::vector<VarType> &variableNameList = BDD_.VariableNameList;

  //! First global node of this bulk domain
  int firstGbNode = BDD_.FirstGbNode;

  //! Last Global node of this bulk domain
  int lastGbNode = BDD_.LastGbNode;
  int numNodes = lastGbNode - firstGbNode + 1;

  bdom.addAttribute("id", id());
  bdom.addAttribute("points", (size_t)numNodes);
  bdom.addAttribute("type", "bulk");
  bdom.addAttribute("numVariables", (size_t) numEquationsPerNode);

  // Dump out the coordinates
  Cantera::XML_Node& gv = bdom.addChild("grid_data");

  std::vector<double> varContig(numNodes);

  int i = 0;
  for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
    NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
    varContig[i] = nv->x0NodePos();
  }
  ctml::addNamedFloatArray(gv, "X0", varContig.size(), &(varContig[0]), "m", "length");

  for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
    VarType vt = variableNameList[iVar];
    i = 0;
    std::string nmm = vt.VariableName(200);
    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
      int istart = nv->EqnStart_GbEqnIndex;
      varContig[i] = (*soln_GLALL_ptr)[istart + ibulk + iVar];
    }
    ctml::addNamedFloatArray(gv, nmm, varContig.size(), &(varContig[0]), "kmol/m3", "concentration");

  }
}

//====================================================================================================================
//
//  This treatment assumes that the problem size stays constant. If this is not the case, the routine will
//  error exit. If we need to change the problem size, then we will need to reinitialize a lot more that just
//  the solution vector. This is best done after we have put in grid refinement.
//
//  Also, we assume that the number of variables stays the same. This may be fiddled with sooner than the number
//  of grid points. There are probably some interesting possibilities here with just initializing a subset of 
//  variables. However, right now, if the number of variables aren't equal and in the same order, the
//  routine will error exit.
//
//  Also we don't consider any interpolation in between time steps. We enter with an XML_Node specific
//  to a particular time step. And then populate the solution vector with the saved solution.
//
//  MP Implementation
//     We've punted on this for now. The basic strategy will be to identity which of the three situations
//     we are currently doing:
//     
//          1) global data into local-processor data structure
//          2) global data into global data structure
//          3) local-processor data into local-processor data structure
// 
//     We are currently set up for #1. However, that may change. Even #1 will fail
//     
void
BulkDomain1D::readDomain(const Cantera::XML_Node& SimulationNode,
                         Epetra_Vector * const soln_GLALL_ptr, Epetra_Vector * const solnDot_GLALL_ptr, double globalTimeRead)
{
    // get the NodeVars object pertaining to this global node
    GlobalIndices *gi = LI_ptr_->GI_ptr_;

    string ids = id();
    Cantera::XML_Node *domainNode_ptr = SimulationNode.findNameID("domain", ids);

    // Number of equations per node
    int numEquationsPerNode = BDD_.NumEquationsPerNode;

    // Vector containing the variable names as they appear in the solution vector
    std::vector<VarType> &variableNameList = BDD_.VariableNameList;

    //! First global node of this bulk domain
    int firstGbNode = BDD_.FirstGbNode;

    //! Last Global node of this bulk domain
    int lastGbNode = BDD_.LastGbNode;
    int numNodes = lastGbNode - firstGbNode + 1;

    string iidd      = (*domainNode_ptr)["id"]; 
    string s_points  = (*domainNode_ptr)["points"]; 
    int points = atoi(s_points.c_str());
    if (points != numNodes) {
        printf("we have an unequal number of points\n");
        exit(-1);
    }
    string ttype    = (*domainNode_ptr)["type"]; 
    string snumVar  = (*domainNode_ptr)["numVariables"]; 
    int numVar = atoi(snumVar.c_str());
    if (numVar != numEquationsPerNode) {
       printf("we have an unequal number of equations\n");
       exit(-1);
    }

    //
    //  Go get the grid data XML node and read it in
    //
    const Cantera::XML_Node* gd_ptr = (*domainNode_ptr).findByName("grid_data");

    std::vector<double> varContig(numNodes);
    ctml::getFloatArray(*gd_ptr, varContig, true, "", "X0");
    int i = 0;
    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      varContig[i] = nv->x0NodePos();
      nv->setupInitialNodePosition(varContig[i], 0.0);
    }

    for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
       VarType vt = variableNameList[iVar];
       i = 0;
       std::string nmm = vt.VariableName(200);
       ctml::getFloatArray(*gd_ptr, varContig, true, "", nmm);
       for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
          NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	  size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
          int istart = nv->EqnStart_GbEqnIndex;
          (*soln_GLALL_ptr)[istart + offset] =  varContig[i];
       }
    }

/*

    for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
       VarType vt = variableNameList[iVar];
       i = 0;
       std::string nmm = vt.VariableName(200);
       ctml::getFloatArray(*gd_ptr, varContig, true, "", nmm);
       for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
          NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
          int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
          int istart = nv->EqnStart_GbEqnIndex;
          (*soln_GLALL_ptr)[istart + ibulk + iVar] =  varContig[i];
       }
    }
*/

}
//====================================================================================================================
//  Fill the vector isAlgebraic with the values from the DomainDescription
void
BulkDomain1D::fillIsAlgebraic(Epetra_IntVector  & isAlgebraic)     
{
  int myBDD_ID = BDD_.ID();
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars *nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    
    int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];

    int bmatch = -1;
    for (int i = 0; i < nodeCent->NumBulkDomains; i++) {
      int id = nodeCent->BulkDomainIndex_BDN[i];
      if (id == myBDD_ID) {
	bmatch = i;
	break;
      }
    }
    if (bmatch != nodeCent->BulkDomainIndex_fromID[myBDD_ID]) {
      printf("we have a prob\n");
      exit(-1);
    }
    int indexBulkDomainOffset = nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[bmatch];
    

    int numVar = BDD_.NumEquationsPerNode;

    for (int k = 0; k < numVar; k++) {
      int isA = BDD_.IsAlgebraic_NE[k];
      isAlgebraic[indexCent_EqnStart + indexBulkDomainOffset + k] = isA;
    }
    
  }
}
//====================================================================================================================
//  Fill the vector isArithmeticScaled with the values from the DomainDescription
void
BulkDomain1D::fillIsArithmeticScaled(Epetra_IntVector  & isArithmeticScaled)     
{
  int myBDD_ID = BDD_.ID();
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars *nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    
    int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];

    int bmatch = -1;
    for (int i = 0; i < nodeCent->NumBulkDomains; i++) {
      int id = nodeCent->BulkDomainIndex_BDN[i];
      if (id == myBDD_ID) {
	bmatch = i;
	break;
      }
    }
    if (bmatch != nodeCent->BulkDomainIndex_fromID[myBDD_ID]) {
      printf("we have a prob\n");
      exit(-1);
    }
    size_t indexBulkDomainOffset = nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[bmatch];
    size_t numVar = BDD_.NumEquationsPerNode;
    for (size_t k = 0; k < numVar; k++) {
      int isA = BDD_.IsArithmeticScaled_NE[k];
      isArithmeticScaled[indexCent_EqnStart + indexBulkDomainOffset + k] = isA;
    }
   
  }
}
//====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 * @param atolVector Reference for the atol vector to fill up
 */
void BulkDomain1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
				 Epetra_Vector_Ghosted & atolVector,
				 const Epetra_Vector_Ghosted * const atolV)
{
  int myBDD_ID = BDD_.ID();
  for (size_t iCell = 0; iCell < (size_t) NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars *nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
    int bmatch = -1;
    for (int i = 0; i < nodeCent->NumBulkDomains; i++) {
      int id = nodeCent->BulkDomainIndex_BDN[i];
      if (id == myBDD_ID) {
	bmatch = i;
	break;
      }
    }
    if (bmatch != nodeCent->BulkDomainIndex_fromID[myBDD_ID]) {
      printf("we have a prob\n");
      exit(-1);
    }
    size_t indexBulkDomainOffset = nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[bmatch];
    size_t numVar = BDD_.NumEquationsPerNode;
    for (size_t k = 0; k < numVar; k++) {
      atolVector[indexCent_EqnStart + indexBulkDomainOffset + k] = atolDefault;
    } 
  } 
}
//================================================================================================================
void
BulkDomain1D::setAtolDeltaDamping(double atolDefault, double relcoeff,  const Epetra_Vector_Ghosted & soln, 
				  Epetra_Vector_Ghosted & atolDeltaDamping,
				  const Epetra_Vector_Ghosted * const atolV)
{
  int myBDD_ID = BDD_.ID();
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars *nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
    int bmatch = -1;
    for (int i = 0; i < nodeCent->NumBulkDomains; i++) {
      int id = nodeCent->BulkDomainIndex_BDN[i];
      if (id == myBDD_ID) {
	bmatch = i;
	break;
      }
    }
    if (bmatch != nodeCent->BulkDomainIndex_fromID[myBDD_ID]) {
      printf("we have a prob\n");
      exit(-1);
    }
    size_t indexBulkDomainOffset = nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[bmatch];
    size_t numVar = BDD_.NumEquationsPerNode;
    for (size_t k = 0; k < numVar; k++) {
      atolDeltaDamping[indexCent_EqnStart + indexBulkDomainOffset + k] = atolDefault * relcoeff;
    }  
  }
}
//================================================================================================================
void 
BulkDomain1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln, 
				    const Epetra_Vector_Ghosted & solnDot, 
				    Epetra_Vector_Ghosted & atolVector_DAEInit,
				    const Epetra_Vector_Ghosted * const atolV)
{
  int myBDD_ID = BDD_.ID();
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars *nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
    int bmatch = -1;
    for (int i = 0; i < nodeCent->NumBulkDomains; i++) {
      int id = nodeCent->BulkDomainIndex_BDN[i];
      if (id == myBDD_ID) {
	bmatch = i;
	break;
      }
    }
    if (bmatch != nodeCent->BulkDomainIndex_fromID[myBDD_ID]) {
      printf("we have a prob\n");
      exit(-1);
    }
    int indexBulkDomainOffset = nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[bmatch];
    int numVar = BDD_.NumEquationsPerNode;
    for (int k = 0; k < numVar; k++) {
      atolVector_DAEInit[indexCent_EqnStart + indexBulkDomainOffset + k] = atolDefault;
    } 
  } 
}
//================================================================================================================
void
BulkDomain1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, 
					  const Epetra_Vector_Ghosted & soln,
					  const Epetra_Vector_Ghosted & solnDot,
					  Epetra_Vector_Ghosted & atolDeltaDamping,
					  const Epetra_Vector_Ghosted * const atolV)
{
  int myBDD_ID = BDD_.ID();
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars *nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
    int bmatch = -1;
    for (int i = 0; i < nodeCent->NumBulkDomains; i++) {
      int id = nodeCent->BulkDomainIndex_BDN[i];
      if (id == myBDD_ID) {
	bmatch = i;
	break;
      }
    }
    if (bmatch != nodeCent->BulkDomainIndex_fromID[myBDD_ID]) {
      printf("we have a prob\n");
      exit(-1);
    }
    int indexBulkDomainOffset = nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[bmatch];
    int numVar = BDD_.NumEquationsPerNode;
    for (int k = 0; k < numVar; k++) {
      atolDeltaDamping[indexCent_EqnStart + indexBulkDomainOffset + k] = atolDefault * relcoeff;
    }  
  }
}
//=====================================================================================================================
// Method for writing the header for the bulk domain to a tecplot file.
void
BulkDomain1D::writeSolutionTecplotHeader()
{
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = !mypid ; //only proc 0 should write

  if (doWrite) {
    
    //open tecplot file
    FILE* ofp;
    string sss = id();
    char filename[20];
    sprintf(filename,"%s%s",sss.c_str(),".dat");
    ofp = fopen( filename, "w");

    //write title and variable list
    fprintf( ofp, "TITLE = \"Solution on Domain %s\"\n",sss.c_str() );

    // Number of equations per node
    int numVar = BDD_.NumEquationsPerNode;
    // Vector containing the variable names as they appear in the solution vector
    std::vector<VarType> &variableNameList = BDD_.VariableNameList;
    //! First global node of this bulk domain

    fprintf( ofp, "VARIABLES = ");
    fprintf( ofp, "\"x [m]\"  \n" );

    for (int k = 0; k < numVar; k++) {
      VarType &vt = variableNameList[k];
      string name = vt.VariableName(15);
      fprintf( ofp, "\"%s\" \t", name.c_str() );
    }
    fprintf(ofp, "\n" );
    fclose(ofp);
  }

}
//=====================================================================================================================
// Method for writing the solution on the surface domain to a tecplot file.
/*
 *
 * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
 * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
 * @param t                    time
 *
 * @param duplicateOnAllProcs  If this is true, all processors will include
 *                             the same log information as proc 0. If
 *                             false, the loginfo will only exist on proc 0.
 */
void
BulkDomain1D::writeSolutionTecplot(const Epetra_Vector *soln_GlAll_ptr, const Epetra_Vector *solnDot_GlAll_ptr,
			           const double t)
{
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = !mypid ; //only proc 0 should write
  if (doWrite) {
    // get the NodeVars object pertaining to this global node
    GlobalIndices *gi = LI_ptr_->GI_ptr_;
    // Number of equations per node
    int numVar = BDD_.NumEquationsPerNode;
    //! First global node of this bulk domain
    int firstGbNode = BDD_.FirstGbNode;
    //! Last Global node of this bulk domain
    int lastGbNode = BDD_.LastGbNode;
    int numNodes = lastGbNode - firstGbNode + 1;
    
    //open tecplot file
    FILE* ofp;
    string sss = id();
    char filename[20];
    sprintf(filename,"%s%s",sss.c_str(),".dat");
    ofp = fopen(filename, "a");
    
    /*
     *  Write out the Heading for the solution at the current time. It's put in a ZONE structure
     *  with T being the heading and SOLUTIONTIME being the value of the time
     */
    fprintf(ofp, "ZONE T = \"t = %g [s]\" I = %d SOLUTIONTIME = %19.11E\n", t, numNodes, t);

    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++) {
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      //x-position
      fprintf(ofp, "%g \t", nv->xNodePos());
      for (int iVar = 0; iVar < numVar; iVar++) {
	//other variables
	int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
	int istart = nv->EqnStart_GbEqnIndex;
	fprintf(ofp, "%g \t", (*soln_GlAll_ptr)[istart + ibulk + iVar]);
      }
      fprintf(ofp, "\n");
    }
    fclose(ofp);

#undef DAKOTAOUT
#ifdef  DAKOTAOUT
    double firstOutputTime = 0.1;
    static double outputTime = 0.1;
    if ( t > outputTime ) {
      std::ofstream dfp;
      if ( outputTime > firstOutputTime )
	dfp.open( "results.out", std::ios_base::app );
      else
	dfp.open( "results.out", std::ios_base::out );
      int indexVS = BDD_.VariableIndexStart_VarName[Voltage];
      double potentialE_left = (*soln_GlAll_ptr)[indexVS];
      double lithium_left = (*soln_GlAll_ptr)[1];
      
      NodalVars *nv = gi->NodalVars_GbNode[lastGbNode-1];
      int istart = nv->EqnStart_GbEqnIndex;
      double potentialE_right = (*soln_GlAll_ptr)[istart+indexVS];
      double lithium_right = (*soln_GlAll_ptr)[istart+1];
      double ohmicLoss =  potentialE_left - potentialE_right;
      double moleFracDelta = lithium_left - lithium_right;

      // uncomment next line to get ohmic loss in electrolyte
      //dfp << ohmicLoss << "  t" << outputTime << std::endl;

      // uncomment next line to get change in litium mole fraction
      //dfp << moleFracDelta << "  t" << outputTime << std::endl;
      
      double Li_eutectic = 0.585/2.0;
      double Rgas_local = 8.3144;
      double Faraday_local = 9.6485e4;
      double concentrationPotential = Rgas_local * TemperatureReference_ / Faraday_local
	* log( lithium_left / lithium_right ) ;

      //uncomment following line to get concentration potential loss
      dfp << concentrationPotential << "  t" << outputTime << std::endl;

      //uncomment following line to get combined ohmic and concentration losses
      outputTime *= 10.0;
    }
      
#endif
#undef DAKOTAOUT
  }
}



//=====================================================================================================================
// Base class for writing the solution on the domain to a logfile.
/*
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
BulkDomain1D::showSolution(const Epetra_Vector *soln_GlAll_ptr,
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
#ifdef DO_OLD_WAY
  showSolution0All(soln_GlAll_ptr, solnDot_GlAll_ptr, soln_ptr, solnDot_ptr, solnOld_ptr,
      residInternal_ptr,t, rdelta_t, indentSpaces, duplicateOnAllProcs);
  return;
#else
  // nn is the number of block rows in the printout
  int nn = NumDomainEqns / 5;
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
  std::vector<VarType> &variableNameList = BDD_.VariableNameList;
  int iBlock;
  int iGbNode;
  int n;
  int nPoints = BDD_.LastGbNode - BDD_.FirstGbNode + 1;
  stream0 ss;

  string indent = "";
  for (int i = 0; i < indentSpaces; i++) {
    indent += " ";
  }
  const char *ind = indent.c_str();
  doublereal v;
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  // Number of points in each vector
  string sss = id();
  bool do0Write = (!mypid || duplicateOnAllProcs);
  print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
  if (do0Write) {
    drawline0(ss, indentSpaces, 80);
    ss.print0("%s  Solution on Bulk Domain %12s : Number of variables = %d\n", ind, sss.c_str(), NumDomainEqns);
    ss.print0("%s                                         : Number of Nodes = %d\n", ind, nPoints);
#ifdef MECH_MODEL
    ss.print0("%s                                         : Beginning pos %g\n", ind, BDD_.Xpos_start);
    ss.print0("%s                                         : Ending    pos %g\n", ind, BDD_.Xpos_end);
#endif
  }
  print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));

  for (iBlock = 0; iBlock < nn; iBlock++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (do0Write) {
      drawline0(ss, indentSpaces, 80);
      ss.print0("%s        z   ", ind);
      for (n = 0; n < 5; n++) {
        int ivar = iBlock * 5 + n;
        VarType vt = variableNameList[ivar];
        string name = vt.VariableName(15);
        ss.print0(" %15s", name.c_str());
      }
      ss.print0("\n");
      drawline0(ss, indentSpaces, 80);
    }
    for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
      doWrite = ((iGbNode >= FirstOwnedGbNode) && (iGbNode <= LastOwnedGbNode));
      if (doWrite) {
        NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
        doublereal x = nv->xNodePos();
        ss.print0("\n%s    %-10.4E ", ind, x);
        int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
        int istart = nv->EqnStart_GbEqnIndex;
        for (n = 0; n < 5; n++) {
          int ivar = iBlock * 5 + n;
          VarType vt = variableNameList[ivar];
          v = (*soln_GlAll_ptr)[istart + ibulk + iBlock * 5 + n];
          ss.print0(" %-10.4E ", v);
        }
      }
    }

    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    if (do0Write) {
      printf("\n");
    }
  }

  int nrem = NumDomainEqns - 5 * nn;
  if (nrem > 0) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (do0Write) {
      drawline0(ss, indentSpaces, 80);
      ss.print0("%s        z   ", ind);
      for (n = 0; n < nrem; n++) {
        int ivar = nn * 5 + n;
        VarType vt = variableNameList[ivar];
        string name = vt.VariableName(15);
        ss.print0(" %15s", name.c_str());
      }
      ss.print0("\n");

      drawline0(ss, indentSpaces, 80);
    }
    int iCell = 0;
    for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
      doWrite = ((iGbNode >= FirstOwnedGbNode) && (iGbNode <= LastOwnedGbNode));
      if (doWrite) {
        int ilc = Index_DiagLcNode_LCO[iCell];
        if (ilc < 0) {
          throw m1d_Error("BulkDomain1D::showSolution pid= " + Cantera::int2str(mypid), "confused");
        }
        int igb = LI_ptr_->IndexGbNode_LcNode[ilc];
        if (iGbNode != igb) {
          throw m1d_Error("BulkDomain1D::showSolution pid= " + Cantera::int2str(mypid), "confused"
              + Cantera::int2str(iGbNode) + " vs " + Cantera::int2str(igb));
        }
        NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
        doublereal x = nv->xNodePos();
        ss.print0("%s %-10.4E ", ind, x);
        int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
        int istart = nv->EqnStart_GbEqnIndex;
        for (n = 0; n < nrem; n++) {
          int ivar = iBlock * 5 + n;
          VarType vt = variableNameList[ivar];
          v = (*soln_GlAll_ptr)[istart + ibulk + nn * 5 + n];
          ss.print0(" %-10.4E ", v);
        }
        ss.print0("\n");
        iCell++;
      }
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    if (do0Write) {
      drawline0(ss, indentSpaces, 80);
    }
  }
#endif
}
//=====================================================================================================================
  
  // Base class for writing a solution vector, not the solution, on the domain to a logfile.
  /*
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
  void
  BulkDomain1D::showSolutionVector(std::string& solnVecName,
				  const Epetra_Vector *solnVector_GlAll_ptr,
				  const Epetra_Vector *solnVector_ptr,
				  const double t,
				  const double rdelta_t,
				  int indentSpaces,
		 	          bool duplicateOnAllProcs, FILE *of)
  {

    // nn is the number of block rows in the printout
    int nn = NumDomainEqns / 5;
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
    std::vector<VarType> &variableNameList = BDD_.VariableNameList;
    int iBlock;
    int iGbNode;
    int n;
    int nPoints = BDD_.LastGbNode - BDD_.FirstGbNode + 1;
    stream0 ss(of);

    string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
      indent += " ";
    }
    const char *ind = indent.c_str();
    doublereal v;
    GlobalIndices *gi = LI_ptr_->GI_ptr_;
    // Number of points in each vector
    string sss = id();
    bool do0Write = (!mypid || duplicateOnAllProcs);
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (do0Write) {
      drawline0(ss, indentSpaces, 100);
      ss.print0("%s  %s Vector on Bulk Domain %12s : Number of variables = %d\n", ind, solnVecName.c_str(),
		sss.c_str(), NumDomainEqns);
      ss.print0("%s                                         : Number of Nodes = %d\n", ind, nPoints);

      ss.print0("%s                                         : Beginning pos %g\n", ind, BDD_.Xpos_start);
      ss.print0("%s                                         : Ending    pos %g\n", ind, BDD_.Xpos_end);
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));

    for (iBlock = 0; iBlock < nn; iBlock++) {
      print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
      if (do0Write) {
	drawline0(ss, indentSpaces, 100);
	ss.print0("%s       z    GblEqnIndx ", ind);
	for (n = 0; n < 5; n++) {
	  int ivar = iBlock * 5 + n;
	  VarType vt = variableNameList[ivar];
	  string name = vt.VariableName(15);
	  ss.print0("%-15.15s", name.c_str());
	}
	ss.print0("\n");
	drawline0(ss, indentSpaces, 100);
      }
      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
	doWrite = ((iGbNode >= FirstOwnedGbNode) && (iGbNode <= LastOwnedGbNode));
	if (doWrite) {
	  NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	  doublereal x = nv->xNodePos();
	  ss.print0("\n%s % -11.4E ", ind, x);
	  int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
	  int istart = nv->EqnStart_GbEqnIndex;
	  ss.print0("%7d", istart + ibulk + iBlock * 5);
	  for (n = 0; n < 5; n++) {
	    int ivar = iBlock * 5 + n;
	    VarType vt = variableNameList[ivar];
	    v = (*solnVector_GlAll_ptr)[istart + ibulk + iBlock * 5 + n];
	    ss.print0("  % -11.4E  ", v);
	  }
	}
      }

      print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
      if (do0Write) {
	ss.fprintf0only("\n");
      }
    }

    int nrem = NumDomainEqns - 5 * nn;
    if (nrem > 0) {
      print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
      if (do0Write) {
	drawline0(ss, indentSpaces, 100);
	ss.print0("%s       z    GblEqnIndx ", ind);
	for (n = 0; n < nrem; n++) {
	  int ivar = nn * 5 + n;
	  VarType vt = variableNameList[ivar];
	  string name = vt.VariableName(15);
	  ss.print0("%-15.15s", name.c_str());
	}
	ss.print0("\n");

	drawline0(ss, indentSpaces, 100);
      }
      int iCell = 0;
      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
	doWrite = ((iGbNode >= FirstOwnedGbNode) && (iGbNode <= LastOwnedGbNode));
	if (doWrite) {
	  int ilc = Index_DiagLcNode_LCO[iCell];
	  if (ilc < 0) {
	    throw m1d_Error("BulkDomain1D::showSolutionVector pid= " + Cantera::int2str(mypid), "confused");
	  }
	  int igb = LI_ptr_->IndexGbNode_LcNode[ilc];
	  if (iGbNode != igb) {
	    throw m1d_Error("BulkDomain1D::showSolutionVector pid= " + Cantera::int2str(mypid), "confused"
			    + Cantera::int2str(iGbNode) + " vs " + Cantera::int2str(igb));
	  }
	  NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	  doublereal x = nv->xNodePos();
	  ss.print0("%s % -11.4E ", ind, x);
	  int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
	  int istart = nv->EqnStart_GbEqnIndex;
	  ss.print0("%7d", istart + ibulk + nn * 5);
	  for (n = 0; n < nrem; n++) {
	    int ivar = iBlock * 5 + n;
	    VarType vt = variableNameList[ivar];
	    v = (*solnVector_GlAll_ptr)[istart + ibulk + nn * 5 + n];
	    ss.print0("  % -11.4E  ", v);
	  }
	  ss.print0("\n");
	  iCell++;
	}
      }
      print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
      if (do0Write) {
	drawline0(ss, indentSpaces, 100);
      }
    }
  }
//=====================================================================================================================
  
  // Base class for writing a solution int vector, not the solution, on the domain to a logfile.
  /*
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
  void
  BulkDomain1D::showSolutionIntVector(std::string& solnVecName,
				      const Epetra_IntVector *solnIntVector_GlAll_ptr,
				      const Epetra_IntVector *solnIntVector_ptr,
				      const double t,
				      const double rdelta_t,
				      int indentSpaces,
				      bool duplicateOnAllProcs, FILE *of)
  {

    // nn is the number of block rows in the printout
    int nn = NumDomainEqns / 5;
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
    std::vector<VarType> &variableNameList = BDD_.VariableNameList;
    int iBlock;
    int iGbNode;
    int n;
    int nPoints = BDD_.LastGbNode - BDD_.FirstGbNode + 1;
    stream0 ss(of);

    string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
      indent += " ";
    }
    const char *ind = indent.c_str();
    int v;
    GlobalIndices *gi = LI_ptr_->GI_ptr_;
    // Number of points in each vector
    string sss = id();
    bool do0Write = (!mypid || duplicateOnAllProcs);
    if (do0Write) {
      drawline0(ss, indentSpaces, 100);
      ss.print0("%s  %s Vector on Bulk Domain %12s : Number of variables = %d\n", ind, solnVecName.c_str(),
		sss.c_str(), NumDomainEqns);
      ss.print0("%s                                         : Number of Nodes = %d\n", ind, nPoints);

      ss.print0("%s                                         : Beginning pos %g\n", ind, BDD_.Xpos_start);
      ss.print0("%s                                         : Ending    pos %g\n", ind, BDD_.Xpos_end);
    }

    for (iBlock = 0; iBlock < nn; iBlock++) {
      print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
      if (do0Write) {
	drawline0(ss, indentSpaces, 100);
	ss.print0("%s       z    GblEqnIndx ", ind);
	for (n = 0; n < 5; n++) {
	  int ivar = iBlock * 5 + n;
	  VarType vt = variableNameList[ivar];
	  string name = vt.VariableName(15);
	  ss.print0("%-15.15s", name.c_str());
	}
	ss.print0("\n");
	drawline0(ss, indentSpaces, 100);
      }
      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
	doWrite = ((iGbNode >= FirstOwnedGbNode) && (iGbNode <= LastOwnedGbNode));
	if (doWrite) {
	  NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	  doublereal x = nv->xNodePos();
	  ss.print0("\n%s % -11.4E ", ind, x);
	  int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
	  int istart = nv->EqnStart_GbEqnIndex;
	  ss.print0("%7d", istart + ibulk + iBlock * 5);
	  for (n = 0; n < 5; n++) {
	    int ivar = iBlock * 5 + n;
	    VarType vt = variableNameList[ivar];
	    v = (*solnIntVector_GlAll_ptr)[istart + ibulk + iBlock * 5 + n];
	    ss.print0("  % -11d  ", v);
	  }
	}
      }

      print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
      if (do0Write) {
	ss.print0("\n");
      }
    }

    int nrem = NumDomainEqns - 5 * nn;
    if (nrem > 0) {
      print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
      if (do0Write) {
	drawline0(ss, indentSpaces, 100);
	ss.print0("%s       z    GblEqnIndx ", ind);
	for (n = 0; n < nrem; n++) {
	  int ivar = nn * 5 + n;
	  VarType vt = variableNameList[ivar];
	  string name = vt.VariableName(15);
	  ss.print0("%-15.15s", name.c_str());
	}
	ss.print0("\n");

	drawline0(ss, indentSpaces, 100);
      }
      int iCell = 0;
      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
	doWrite = ((iGbNode >= FirstOwnedGbNode) && (iGbNode <= LastOwnedGbNode));
	if (doWrite) {
	  int ilc = Index_DiagLcNode_LCO[iCell];
	  if (ilc < 0) {
	    throw m1d_Error("BulkDomain1D::showSolutionVector pid= " + Cantera::int2str(mypid), "confused");
	  }
	  int igb = LI_ptr_->IndexGbNode_LcNode[ilc];
	  if (iGbNode != igb) {
	    throw m1d_Error("BulkDomain1D::showSolutionVector pid= " + Cantera::int2str(mypid), "confused"
			    + Cantera::int2str(iGbNode) + " vs " + Cantera::int2str(igb));
	  }
	  NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	  doublereal x = nv->xNodePos();
	  ss.print0("%s % -11.4E ", ind, x);
	  int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
	  int istart = nv->EqnStart_GbEqnIndex;
	  ss.print0("%7d", istart + ibulk + nn * 5);
	  for (n = 0; n < nrem; n++) {
	    int ivar = iBlock * 5 + n;
	    VarType vt = variableNameList[ivar];
	    v = (*solnIntVector_GlAll_ptr)[istart + ibulk + nn * 5 + n];
	    ss.print0("  % -11d  ", v);
	  }
	  ss.print0("\n");
	  iCell++;
	}
      }
      print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
      if (do0Write) {
	drawline0(ss, indentSpaces, 100);
      }
    }
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
BulkDomain1D::showSolution0All(const Epetra_Vector *soln_GlAll_ptr,
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
 m1d::m1d_Error(" BulkDomain1D showSolution0All","drop dead");

  // nn is the number of block rows in the printout
  int nn = NumDomainEqns / 5;
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = !mypid || duplicateOnAllProcs;
  std::vector<VarType> &variableNameList = BDD_.VariableNameList;
  int iBlock;
  int iGbNode;
  int n;
  int nPoints = BDD_.LastGbNode - BDD_.FirstGbNode + 1;
  //stream0 ss;
  // print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
  char buf[100];
  string indent = "";
  for (int i = 0; i < indentSpaces; i++) {
    indent += " ";
  }
  const char *ind = indent.c_str();
  doublereal v;
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  // Number of points in each vector
  string sss = id();
  if (doWrite) {
    drawline(indentSpaces, 80);
    sprintf(buf, "%s  Solution on Bulk Domain %12s : Number of variables = %d\n", ind, sss.c_str(), NumDomainEqns);
    Cantera::writelog(buf);
    sprintf(buf, "%s                                         : Number of Nodes = %d\n", ind, nPoints);
    Cantera::writelog(buf);
    sprintf(buf, "%s                                         : Beginning pos %g\n", ind, BDD_.Xpos_start);
    Cantera::writelog(buf);
    sprintf(buf, "%s                                         : Ending    pos %g\n", ind, BDD_.Xpos_end);
    Cantera::writelog(buf);
  }
  if (doWrite) {
    for (iBlock = 0; iBlock < nn; iBlock++) {
      drawline(indentSpaces, 80);
      sprintf(buf, "%s        z   ", ind);
      Cantera::writelog(buf);
      for (n = 0; n < 5; n++) {
        int ivar = iBlock * 5 + n;
        VarType vt = variableNameList[ivar];
        string name = vt.VariableName(15);
        sprintf(buf, " %15s", name.c_str());
        Cantera::writelog(buf);
      }
      sprintf(buf, "\n");
      Cantera::writelog(buf);
      drawline(indentSpaces, 80);

      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
        NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
        doublereal x = nv->xNodePos();
        sprintf(buf, "\n%s    %-10.4E ", ind, x);
        Cantera::writelog(buf);
        int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
        int istart = nv->EqnStart_GbEqnIndex;
        for (n = 0; n < 5; n++) {
          int ivar = iBlock * 5 + n;
          VarType vt = variableNameList[ivar];
          v = (*soln_GlAll_ptr)[istart + ibulk + iBlock * 5 + n];
          sprintf(buf, " %-10.4E ", v);
          Cantera::writelog(buf);
        }
      }
      Cantera::writelog("\n");
    }

    int nrem = NumDomainEqns - 5 * nn;
    if (nrem > 0) {
      drawline(indentSpaces, 80);
      sprintf(buf, "%s        z   ", ind);
      Cantera::writelog(buf);
      for (n = 0; n < nrem; n++) {
        int ivar = nn * 5 + n;
        VarType vt = variableNameList[ivar];
        string name = vt.VariableName(15);
        sprintf(buf, " %15s", name.c_str());
        Cantera::writelog(buf);
      }
      sprintf(buf, "\n");
      Cantera::writelog(buf);
      drawline(indentSpaces, 80);

      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
        NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
        doublereal x = nv->xNodePos();
        sprintf(buf, "%s    %-10.4E ", ind, x);
        Cantera::writelog(buf);
        int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
        int istart = nv->EqnStart_GbEqnIndex;
        for (n = 0; n < nrem; n++) {
          int ivar = iBlock * 5 + n;
          VarType vt = variableNameList[ivar];
          v = (*soln_GlAll_ptr)[istart + ibulk + nn * 5 + n];
          sprintf(buf, " %-10.4E ", v);
          Cantera::writelog(buf);
        }
        sprintf(buf, "\n");
        Cantera::writelog(buf);
      }
    }
    drawline(indentSpaces, 80);
  }
}
//=====================================================================================================================
double
BulkDomain1D::getPointTemperature(const NodalVars* const nv,
                                  const doublereal* const solutionPoint) const
{
    /*
     * Get the temperature: Check to see if the temperature is in the solution vector.
     *   If it is not, then use the reference temperature
     */
    std::map<m1d::Var_Name_Enum, size_t>::const_iterator it = nv->Offset_VarType.find(Temperature);
    size_t ieqnTemp = it->second;
    if (ieqnTemp != npos) {
        return solutionPoint[ieqnTemp];
    }
    return TemperatureReference_;
}
//=====================================================================================================================
double
BulkDomain1D::getPointPressure(const NodalVars* const nv,
                               const doublereal* const solutionPoint) const
{
    // 
    //  Note there are two pressure variables listed, Pressure_Axial and Pressure_Radial
    //  Pressure_Radial is really an eigenvalue, so it doesn't factor in here
    //  Pressure_Axial is not live, yet. When it is, will have to reevaluate this.
    //
    const static size_t presS = (size_t) Pressure_Axial;
    int iPres = nv->indexBulkDomainVar0(presS);
    if (iPres >= 0) {
        return solutionPoint[iPres];
    }
    return PressureReference_;
}
//=====================================================================================================================
// Get parameters specified by text strings
int 
BulkDomain1D::reportSolutionParam(const std::string& paramName, double* const paramVal) const
{
  paramVal[0] = 0.0;
  return -1;
}
//=====================================================================================================================
int
BulkDomain1D::reportSolutionVector(const std::string& requestID, const int requestType, const Epetra_Vector *soln_ptr,
				   std::vector<double>& vecInfo) const
{
    //GlobalIndices *gi = LI_ptr_->GI_ptr_;
    size_t findV = npos;
    VarType vtF;
    string name;
    vecInfo.clear();
    if (requestType == 1) {
	if (requestID != "x [m]") {
	    throw m1d_Error("reportSolutionVector", "don't know how to respond to request: " + requestID);
	}
    } else if (requestType == 0) {
    for (size_t i = 0; i < (size_t) NumDomainEqns; ++i) {
	const VarType& vt = BDD_.VariableNameList[i];
        const string& varName = vt.VariableName();
	if (mdpUtil::LowerCaseStringEquals(varName, requestID)) {
            findV = i;
	    vtF = vt;
	    break;
	}
    }
    if (findV == npos) {
	return -1;
    }
    }
    const Epetra_Vector& soln = *soln_ptr;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
	int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
	int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];


	NodalVars *nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];




        const double *solnCentStart = &(soln[indexCent_EqnStart]);
	//
        // Find the start of the solution at the current node
        //
	if (requestType == 0) {
	    size_t iOffSet = nodeCent->indexBulkDomainVar(vtF);
	    if (iOffSet == npos) {
		return -1;
	    }
	    vecInfo.push_back(solnCentStart[iOffSet]);
	} else {
	    if (requestID == "x [m]") {
		double x = nodeCent->xNodePos();
		vecInfo.push_back(x);
	    }
	}
    }

    return NumLcCells;
}
//=====================================================================================================================
void
BulkDomain1D::err(const char *msg)
{
  printf("BulkDomain1D: function not implemented: %s\n", msg);
  exit(-1);
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================
