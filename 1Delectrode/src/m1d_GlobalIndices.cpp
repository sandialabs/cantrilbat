/**
 * @file m1d_GlobalIndices.cpp
 *
 */

/*
 *  $Id: m1d_GlobalIndices.cpp 504 2013-01-07 22:32:48Z hkmoffa $
 */
#include "m1d_GlobalIndices.h"
#include "m1d_DomainLayout.h"
#include "m1d_NodalVars.h"
#include "m1d_exception.h"
#include "m1d_EpetraExtras.h"


#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "cantera/base/stringUtils.h" 

#define safeDelete(ptr)  if (ptr) { delete ptr; ptr = 0; }

namespace m1d
{

//GlobalIndices *GI_ptr = 0;

//===========================================================================
GlobalIndices::GlobalIndices(Epetra_Comm *comm_ptr) :
  Comm_ptr_(comm_ptr), NumProc(1), MyProcID(0), NumGbNodes(0), NumGbEqns(0), NumOwnedLcNodes(0), NumOwnedLcEqns(0),
  GbNodetoOwnedLcNodeMap(0),
  GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap(0),
      XNodePos_GbNode(0), GbEqnstoAllMap(0), SolnAll(0), SolnIntAll(0), SolnDotAll(0), DL_ptr_(0)
{
}
//===========================================================================
GlobalIndices::~GlobalIndices() {
  safeDelete(GbNodetoOwnedLcNodeMap);
  safeDelete(GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap);
  for (int iLcNodes = 0; iLcNodes < NumGbNodes; iLcNodes++) {
    delete NodalVars_GbNode[iLcNodes];
    NodalVars_GbNode[iLcNodes] = 0;
  }
  safeDelete(XNodePos_GbNode);
  safeDelete(GbEqnstoAllMap);
  safeDelete(SolnAll);
  safeDelete(SolnIntAll);
  safeDelete(SolnDotAll);
  DL_ptr_ = 0;
}
//===========================================================================
GlobalIndices::GlobalIndices(const GlobalIndices &r) :
  Comm_ptr_(0), NumProc(1), MyProcID(0), NumGbNodes(0), NumGbEqns(0), NumOwnedLcNodes(0), NumOwnedLcEqns(0),
  GbNodetoOwnedLcNodeMap(0), 
  GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap(0),
      XNodePos_GbNode(0), GbEqnstoAllMap(0), SolnAll(0), SolnIntAll(0), SolnDotAll(0), DL_ptr_(0)
{
  *this = r;
}
//===========================================================================
GlobalIndices &
GlobalIndices::operator=(const GlobalIndices &r)
{
  if (this == &r)
    return *this;

  for (int iLcNodes = 0; iLcNodes < NumGbNodes; iLcNodes++) {
    safeDelete(NodalVars_GbNode[iLcNodes]);
  }

  Comm_ptr_ = r.Comm_ptr_;
  NumGbNodes = r.NumGbNodes;
  NumGbEqns = r.NumGbEqns;
  IndexStartGbNode_Proc = r.IndexStartGbNode_Proc;
  IndexStartGbEqns_Proc = r.IndexStartGbEqns_Proc;
  NumOwnedLcNodes_Proc = r.NumOwnedLcNodes_Proc;
  NumOwnedLcEqns_Proc = r.NumOwnedLcEqns_Proc;
  NumEqns_GbNode = r.NumEqns_GbNode;
  IndexStartGbEqns_GbNode = r.IndexStartGbEqns_GbNode;

  NumOwnedLcNodes = r.NumOwnedLcNodes;
  NumOwnedLcEqns = r.NumOwnedLcEqns;

  safeDelete(GbNodetoOwnedLcNodeMap);
  GbNodetoOwnedLcNodeMap = new Epetra_Map(*r.GbNodetoOwnedLcNodeMap);

  safeDelete(GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap);
  GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap = new Epetra_BlockMap(*r.GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap);

  /*
   * Create a deep copy of the NodalVars vector using the copy constructor
   * for NodalVars
   */
  NodalVars_GbNode.resize(NumGbNodes);
  for (int iLcNodes = 0; iLcNodes < NumGbNodes; iLcNodes++) {
    NodalVars_GbNode[iLcNodes] = new NodalVars(*(r.NodalVars_GbNode[iLcNodes]));
  }

  safeDelete(XNodePos_GbNode);
  XNodePos_GbNode = new Epetra_Vector(*(r.XNodePos_GbNode));

  safeDelete(GbEqnstoAllMap);
  GbEqnstoAllMap = new Epetra_BlockMap(*r.GbEqnstoAllMap);

  safeDelete(SolnAll);
  SolnAll = new Epetra_Vector(*r.SolnAll);

  safeDelete(SolnIntAll);
  SolnIntAll = new Epetra_IntVector(*r.SolnIntAll);


  safeDelete(SolnDotAll);
  SolnDotAll = new Epetra_Vector(*r.SolnDotAll);

  DL_ptr_ = r.DL_ptr_;

  return *this;
}
//=====================================================================================================================
// Initialize the number of global nodes
/*
 * Size all arrays based on the total number of global nodes
 */
void
GlobalIndices::init(DomainLayout *dl_ptr)
{
  DomainLayout &DL = *dl_ptr;
  DL_ptr_ = dl_ptr;
  // We get the total number of global nodes from the DomainLayout object
  NumGbNodes = DL.NumGbNodes;
  NodalVars_GbNode.resize(NumGbNodes);
  NumEqns_GbNode.resize(NumGbNodes);
  IndexStartGbEqns_GbNode.resize(NumGbNodes);
  int *tmpIndex = new int[NumGbNodes];
  for (int iGbNode = 0; iGbNode < NumGbNodes; iGbNode++) {
    AssertTrace(NodalVars_GbNode[iGbNode] == 0);
    NodalVars_GbNode[iGbNode] = new NodalVars(iGbNode, dl_ptr);
    tmpIndex[iGbNode] = iGbNode;
  }
  /*
   * Create a map with all of the unknowns on all processors
   */
  Epetra_Map *eAll_map = new Epetra_Map(NumProc * NumGbNodes, NumGbNodes, tmpIndex, 0, *Comm_ptr_);
  XNodePos_GbNode = new Epetra_Vector(*eAll_map, true);
  safeDelete(eAll_map);
  delete [] tmpIndex;
}
//=====================================================================================================================
int
GlobalIndices::discoverNumEqnsPerNode()
{
  int numEqn = 0;
  for (int iGbNode = 0; iGbNode < NumGbNodes; iGbNode++) {
#ifdef DEBUG_MODE
  if (iGbNode == 18) {
       //printf("we are here\n");
  }
#endif
    NodalVars *nv = NodalVars_GbNode[iGbNode];
    nv->DiscoverDomainsAtThisNode();
    nv->GenerateEqnOrder();

    NumEqns_GbNode[iGbNode] = nv->NumEquations;
    IndexStartGbEqns_GbNode[iGbNode] = numEqn;
    nv->EqnStart_GbEqnIndex = numEqn;
    numEqn += nv->NumEquations;
  }
  NumGbEqns = numEqn;
  return numEqn;
}
//=====================================================================================================================
void
GlobalIndices::procDivide()
{
  int iGbNode, iProc;

  // Resize arrays based on the Number of processors
  IndexStartGbNode_Proc.resize(NumProc);
  IndexStartGbEqns_Proc.resize(NumProc);
  NumOwnedLcNodes_Proc.resize(NumProc);
  NumOwnedLcEqns_Proc.resize(NumProc);

  //  Break up the problem amongst the processors, based on
  //  distributing the nodes evenly. Ok, this perhaps needs work,
  //  but let's use it to start, because it probably doesn't matter
  //  much anyway.
  //
  //  NumLCNodes_Proc[] will contain the breakup, in terms of distributing
  //  the nodes. Numbering of global nodes will be contiguous from proc 0
  //  to proc NumProc-1 based on this vector.
  //
  int base = NumGbNodes / NumProc;
  int remainder = NumGbNodes - NumProc * base;

  for (iProc = 0; iProc < NumProc; iProc++) {
    NumOwnedLcNodes_Proc[iProc] = base;
    if (iProc < remainder) {
      NumOwnedLcNodes_Proc[iProc]++;
    }
  }

  IndexStartGbNode_Proc[0] = 0;
  for (iProc = 1; iProc < NumProc; iProc++) {
    IndexStartGbNode_Proc[iProc] = IndexStartGbNode_Proc[iProc - 1] + NumOwnedLcNodes_Proc[iProc - 1];
  }

  IndexStartGbEqns_Proc[0] = 0;
  for (iProc = 0; iProc < NumProc; iProc++) {
    int iGbNode_start = IndexStartGbNode_Proc[iProc];
    //int iGbNode_end  = iGbNode_start  + numLnNodes_Proc[iProc] - 1;
    if (iProc > 0) {
      IndexStartGbEqns_Proc[iProc] = IndexStartGbEqns_GbNode[iGbNode];
    }

    int iEqnSum = 0;
    iGbNode = iGbNode_start;
    for (int iLnNode = 0; iLnNode < NumOwnedLcNodes_Proc[iProc]; iLnNode++, iGbNode++) {
      iEqnSum += NumEqns_GbNode[iGbNode];
    }
    NumOwnedLcEqns_Proc[iProc] = iEqnSum;
  }

  NumOwnedLcNodes = NumOwnedLcNodes_Proc[MyProcID];
  NumOwnedLcEqns = NumOwnedLcEqns_Proc[MyProcID];

}
//=====================================================================================================================
//   Form a Map to describe the distribution of our owned equations
//  and our owned nodes on the processors
void
GlobalIndices::initNodeMaps()
{
  AssertTrace(Comm_ptr_);
  GbNodetoOwnedLcNodeMap = new Epetra_Map(NumGbNodes, NumOwnedLcNodes, 0, *Comm_ptr_);
}
//=====================================================================================================================
void
GlobalIndices::initBlockNodeMaps(int *numEqns_LcNode)
{
  int * myGlobalNodes = GbNodetoOwnedLcNodeMap->MyGlobalElements();

  GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap = new Epetra_BlockMap(NumGbNodes, NumOwnedLcNodes, myGlobalNodes,
      numEqns_LcNode, 0, *Comm_ptr_);

  std::vector<int> listI(NumGbNodes);
  for (int i = 0; i < NumGbNodes; i++) {
    listI[i] = i;
  }

  GbEqnstoAllMap = new Epetra_BlockMap(NumProc * NumGbNodes, NumGbNodes, &(listI[0]), &(NumEqns_GbNode[0]), 0, *Comm_ptr_);

  SolnAll = new Epetra_Vector(*GbEqnstoAllMap, true);
  SolnDotAll = new Epetra_Vector(*GbEqnstoAllMap, true);
  SolnIntAll = new Epetra_IntVector(*GbEqnstoAllMap, true);
}
//=====================================================================================================================

//=====================================================================================================================
void
GlobalIndices::InitMesh()
{
  DL_ptr_->InitializeXposNodes(this);
}
//=====================================================================================================================
// This utility function will return the global node number
// given the global equation number
/*
 * @param   rowEqnNum returns the node equation number
 * @return  Returns the global node number. Will return -1 if there
 *          is a problem.
 */
int
GlobalIndices::GbEqnToGbNode(const int GbEqnNum, int & rowEqnNum) const
{
#ifdef DEBUG_MODE
  if ((GbEqnNum < 0) || GbEqnNum >= NumGbEqns) {
    throw m1d_Error("GlobalIndices::GbEqnToGbNode", "GbEqnNum out of bounds: " + Cantera::int2str(GbEqnNum));
  }
#endif
  int nfound = -1;
  for (int i = 1; i < NumGbNodes; i++) {
    if (IndexStartGbEqns_GbNode[i] > GbEqnNum) {
      nfound = i - 1;
      break;
    }
  }
  nfound = NumGbNodes - 1;
  if (nfound >= 0) {
    rowEqnNum = GbEqnNum - IndexStartGbEqns_GbNode[nfound];
  } else {
    rowEqnNum = -1;
  }
  return nfound;
}
//=====================================================================================================================
void
GlobalIndices::updateGlobalPositions(Epetra_Vector *Xpos_LcOwnedNode_p)
{
  gather_nodeV_OnAll(*XNodePos_GbNode, *Xpos_LcOwnedNode_p, Comm_ptr_);
  for (int iGbNode = 0; iGbNode < NumGbNodes; iGbNode++) {
    NodalVars *nv = NodalVars_GbNode[iGbNode];
    nv->changeNodePosition((*XNodePos_GbNode)[iGbNode]);
  }
}
//=====================================================================================================================
}
//=====================================================================================================================
