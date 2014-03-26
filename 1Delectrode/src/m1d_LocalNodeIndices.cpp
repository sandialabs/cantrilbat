/**
 * @file m1d_LocalRowNodeIndices.cpp
 *
 */

/*
 *  $Id: m1d_LocalNodeIndices.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"

#include "m1d_defs.h"
#include "m1d_GlobalIndices.h"
#include "m1d_LocalNodeIndices.h"
#include "m1d_exception.h"
#include "m1d_NodalVars.h"
#include "m1d_Comm.h"


namespace m1d
{

#define safeDelete(ptr)  if (ptr) { delete ptr; ptr = 0; }

//===========================================================================
// Establish global storage for the global pointer
//LocalNodeIndices *LI_ptr = 0;

//===========================================================================
LocalNodeIndices::LocalNodeIndices(Epetra_Comm *comm_ptr, GlobalIndices *gi_ptr) :
  Comm_ptr_(comm_ptr), MyProcID(0), NumLcNodes(0), NumOwnedLcNodes(0), NumExtNodes(0), NumLcRowNodes(0), NumLcEqns(0),
      NumLcOwnedEqns(0), GCNIndexGbNode(-1), GCNIndexLcNode(-1), RightLcNode(-1), LeftLcNode(-1), 
      GbNodetoLcNodeColMap(0), 
      GbNodetoOwnedLcNodeMap(0), 
      GbBlockNodeEqnstoLcBlockNodeEqnsColMap(0), GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap(0),
      Importer_GhostEqns(0), Importer_NodalValues(0), NumNodeColors(0), NodeColorMap(0), NumEqnColors(0),
      EqnColorMap(0), EqnColors(0), Xpos_LcNode_p(0), Xpos_LcOwnedNode_p(0), GI_ptr_(gi_ptr)
{
}
//===========================================================================
LocalNodeIndices::~LocalNodeIndices() {

  //  for (int iLcNode = 0; iLcNode < NumOwnedLcNodes; iLcNode++) {
  // safeDelete(NodalVars_LcNode[iLcNode]);
  //}
  safeDelete(NodeColorMap);
  safeDelete(GbNodetoLcNodeColMap);
  safeDelete(GbNodetoOwnedLcNodeMap);
  safeDelete(GbBlockNodeEqnstoLcBlockNodeEqnsColMap);
  safeDelete(GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap);
  safeDelete(EqnColorMap);
  safeDelete(EqnColors);
  safeDelete(Importer_GhostEqns);
  safeDelete(Importer_NodalValues);
  safeDelete(Xpos_LcNode_p);
  safeDelete(Xpos_LcOwnedNode_p);
}
//===========================================================================
LocalNodeIndices::LocalNodeIndices(const LocalNodeIndices &r) :
  MyProcID(0), NumLcNodes(0), NumOwnedLcNodes(0), NumExtNodes(0), NumLcRowNodes(0), NumLcEqns(0), NumLcOwnedEqns(0),
      GCNIndexGbNode(-1), GCNIndexLcNode(-1), RightLcNode(-1), LeftLcNode(-1), 
      GbNodetoLcNodeColMap(0),
      GbNodetoOwnedLcNodeMap(0),
      GbBlockNodeEqnstoLcBlockNodeEqnsColMap(0), GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap(0), Importer_GhostEqns(0),
      Importer_NodalValues(0), NumNodeColors(0), NodeColorMap(0), EqnColorMap(0), EqnColors(0), Xpos_LcNode_p(0),
      Xpos_LcOwnedNode_p(0), GI_ptr_(r.GI_ptr_)
{
  *this = r;
}
//===========================================================================
// assignment operator
/*
 *  Do a deep copy of all Maps within the operator
 */
LocalNodeIndices &
LocalNodeIndices::operator=(const LocalNodeIndices &r)
{
  if (this == &r)
    return *this;

  Comm_ptr_ = r.Comm_ptr_;
  MyProcID = r.MyProcID;
  NumLcNodes = r.NumLcNodes;
  NumOwnedLcNodes = r.NumOwnedLcNodes;
  NumExtNodes = r.NumExtNodes;
  NumLcRowNodes = r.NumLcRowNodes;
  NumLcOwnedEqns = r.NumLcOwnedEqns;
  IndexGbNode_LcNode = r.IndexGbNode_LcNode;
  NumEqns_LcNode = r.NumEqns_LcNode;
  IndexLcEqns_LcNode = r.IndexLcEqns_LcNode;
  IndexGbEqns_LcEqns = r.IndexGbEqns_LcEqns;
  IsExternal_LcNode = r.IsExternal_LcNode;
  GCNIndexGbNode = r.GCNIndexGbNode;
  GCNIndexLcNode = r.GCNIndexLcNode;
  RightLcNode = r.RightLcNode;
  LeftLcNode = r.LeftLcNode;

  safeDelete(GbNodetoLcNodeColMap);
  GbNodetoLcNodeColMap = new Epetra_Map(*r.GbNodetoLcNodeColMap);

  safeDelete(GbNodetoOwnedLcNodeMap);
  GbNodetoOwnedLcNodeMap = new Epetra_Map(*r.GbNodetoOwnedLcNodeMap);

  safeDelete(GbBlockNodeEqnstoLcBlockNodeEqnsColMap);
  GbBlockNodeEqnstoLcBlockNodeEqnsColMap = new Epetra_BlockMap(*r.GbBlockNodeEqnstoLcBlockNodeEqnsColMap);

  safeDelete(GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap);
  GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap = new Epetra_BlockMap(*r.GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap);

  GCNIndexGbNode = r.GCNIndexGbNode;
  GCNIndexLcNode = r.GCNIndexLcNode;

  // Make a deep copy of the importer
  safeDelete(Importer_GhostEqns);
  if (r.Importer_GhostEqns) {
    Importer_GhostEqns = new Epetra_Import(*(r.Importer_GhostEqns));
  }

  // Make a deep copy of the importer
  safeDelete(Importer_NodalValues);
  if (r.Importer_NodalValues) {
    Importer_NodalValues = new Epetra_Import(*(r.Importer_NodalValues));
  }

  NumNodeColors = r.NumNodeColors;
  safeDelete(NodeColorMap);
  if (r.NodeColorMap) {
    NodeColorMap = new Epetra_MapColoring(*(r.NodeColorMap));
  }

  NumEqnColors = r.NumEqnColors;
  safeDelete(EqnColorMap);
  if (r.EqnColorMap) {
    EqnColorMap = new Epetra_MapColoring(*(r.EqnColorMap));
  }

  safeDelete(EqnColors);
  if (r.EqnColors) {
    EqnColors = new Epetra_IntVector(*(r.EqnColors));
  }

  NodalVars_LcNode = r.NodalVars_LcNode;

  safeDelete(Xpos_LcNode_p);
  Xpos_LcNode_p = new Epetra_Vector(*(r.Xpos_LcNode_p));

  safeDelete(Xpos_LcOwnedNode_p);
  double *V;
  Epetra_DataAccess eee = View;
  Xpos_LcNode_p->ExtractView(&V);
  Xpos_LcOwnedNode_p = new Epetra_Vector(eee, *(GI_ptr_->GbNodetoOwnedLcNodeMap), V);

  GI_ptr_ = r.GI_ptr_;

  return *this;
}
//===========================================================================
void
LocalNodeIndices::determineLcNodeMaps(DomainLayout *dl_ptr)
{
  MyProcID = GI_ptr_->MyProcID;
  NumLcRowNodes = GI_ptr_->NumOwnedLcNodes_Proc[MyProcID];
  NumOwnedLcNodes = GI_ptr_->NumOwnedLcNodes_Proc[MyProcID];

  NumExtNodes = 2;
  if (GI_ptr_->MyProcID == 0) {
    NumExtNodes--;
  }
  if (GI_ptr_->MyProcID == (GI_ptr_->NumProc - 1)) {
    NumExtNodes--;
  }
  NumLcNodes = GI_ptr_->NumOwnedLcNodes + NumExtNodes;
  NumLcRowNodes = GI_ptr_->NumOwnedLcNodes;

  /*
   * Ok, we now know the number of local nodes. So, let's
   * resize all vectors based on the number of local nodes.
   */
  IndexGbNode_LcNode.resize(NumLcNodes);
  NumEqns_LcNode.resize(NumLcNodes);
  IndexLcEqns_LcNode.resize(NumLcNodes);
  IsExternal_LcNode.resize(NumLcNodes);
  IDLeftLcNode_LcNode.resize(NumLcNodes);
  IDLeftLcNode_LcNode[0] = -1;
  IDRightLcNode_LcNode.resize(NumLcNodes);
  IDRightLcNode_LcNode[NumOwnedLcNodes - 1] = -1;
  IDLeftGbNode_LcNode.resize(NumLcNodes, -1);
  IDRightGbNode_LcNode.resize(NumLcNodes, -1);
  X0pos_LcNode.resize(NumLcNodes, -1.0);

  NodalVars_LcNode.resize(NumLcNodes);
  for (int iLcNode = 0; iLcNode < NumOwnedLcNodes; iLcNode++) {
    // get the global node number
    int gn = IndexGbNode_LcNode[iLcNode];
    NodalVars_LcNode[iLcNode] = GI_ptr_->NodalVars_GbNode[gn];
  }

  int *MyGlobalElements = (GI_ptr_->GbNodetoOwnedLcNodeMap)->MyGlobalElements();
  for (int iLcNode = 0; iLcNode < NumOwnedLcNodes; iLcNode++) {
    IndexGbNode_LcNode[iLcNode] = MyGlobalElements[iLcNode];
    int numEqn = (GI_ptr_->NumEqns_GbNode)[IndexGbNode_LcNode[iLcNode]];
    NumEqns_LcNode[iLcNode] = numEqn;
    IsExternal_LcNode[iLcNode] = false;
    if (iLcNode > 0) {
      IDLeftLcNode_LcNode[iLcNode] = iLcNode - 1;
    }
    if (iLcNode < NumOwnedLcNodes - 1) {
      IDRightLcNode_LcNode[iLcNode] = iLcNode + 1;
    }
  }

  int index = NumOwnedLcNodes;

  // Identify the Right-most node on the processor (higher global node number)
  /*
   *  Assign the local index to the index of the highest owned node first.
   */
  RightLcNode = NumOwnedLcNodes - 1;
  /*
   *  If there is a right ghost node, reassign it to that
   */
  if (GI_ptr_->MyProcID < (GI_ptr_->NumProc - 1)) {
    /*
     * Get the global ID of the node to the right of the current nodes
     * on the processor
     */
    int rightGlob = GI_ptr_->IndexStartGbNode_Proc[GI_ptr_->MyProcID + 1];
    RightLcNode = index;
    IndexGbNode_LcNode[index] = rightGlob;
    NumEqns_LcNode[index] = GI_ptr_->NumEqns_GbNode[rightGlob];
    /*
     * Assign the right node of the existing end to the current node
     */
    IDRightLcNode_LcNode[NumOwnedLcNodes - 1] = index;
    /*
     * Assign the right node to the new end to -1
     *
     */
    IDRightLcNode_LcNode[index] = -1;
    IsExternal_LcNode[index] = true;
    /*
     * Increment the index for the next block
     */
    index++;
  }

  // Identify the Left-most node on the processor (lower global node number)
  /*
   *  Assign the local index to the index of the lowest owned node first.
   */
  LeftLcNode = 0;
  if (GI_ptr_->MyProcID > 0) {
    /*
     * Get the global ID of the node to the left of the current nodes
     * on the processor
     */
    int leftGlob = GI_ptr_->IndexStartGbNode_Proc[GI_ptr_->MyProcID] - 1;
    LeftLcNode = index;
    IndexGbNode_LcNode[index] = leftGlob;
    NumEqns_LcNode[index] = GI_ptr_->NumEqns_GbNode[IndexGbNode_LcNode[index]];
    /*
     * Assign the left node of the first local node to the current node
     */
    IDLeftLcNode_LcNode[0] = index;
    /*
     * Assign the left node of the new left end to -1
     */
    IDLeftLcNode_LcNode[index] = -1;
    IsExternal_LcNode[index] = true;
    /*
     * Increment the index for the next block
     */
    index++;
  }
  /*
   * Use an assert for sanity checking
   */
  AssertThrow(index == NumLcNodes, "LocalNodeIndices::initSize()");

  int igb;
  for (int iLcNode = 0; iLcNode < NumOwnedLcNodes; iLcNode++) {
    int index = IDLeftLcNode_LcNode[iLcNode];
    if (index >= 0) {
      igb = IndexGbNode_LcNode[index];
      IDLeftGbNode_LcNode[iLcNode] = igb;
    }
    index = IDRightLcNode_LcNode[iLcNode];
    if (index >= 0) {
      igb = IndexGbNode_LcNode[index];
      IDRightGbNode_LcNode[iLcNode] = igb;
    }
  }

  initLcNodeMaps();

  initLcBlockNodeMaps();

  makeNodeColors();

  Xpos_LcNode_p = new Epetra_Vector(*GbNodetoLcNodeColMap, false);
  double *V;
  Epetra_DataAccess eee = View;
  Xpos_LcNode_p->ExtractView(&V);
  safeDelete(Xpos_LcOwnedNode_p);
  Xpos_LcOwnedNode_p = new Epetra_Vector(eee, *(GI_ptr_->GbNodetoOwnedLcNodeMap), V);
}
//===========================================================================
void
LocalNodeIndices::determineLcEqnMaps()
{
  makeEqnColors();

  makeImporters();
}
//===========================================================================
void
LocalNodeIndices::initLcNodeMaps()
{
  /*
   *  Create a map of the local nodes and ghost nodes. The global ids are the global node numbers.
   *  We rely on Epetra to tally up the total number of global elements produced. This map will included ghosted nodes.
   */
  GbNodetoLcNodeColMap = new Epetra_Map(-1, NumLcNodes, DATA_PTR(IndexGbNode_LcNode), 0, *Comm_ptr_);
#ifdef DEBUG_MATRIX_STRUCTURE
  print0_epBlockMap(*GbNodetoLcNodeColMap);
#endif

  /*
   *  Create a map of the local nodes. The global ids are the global node numbers.
   *  We rely on Epetra to tally up the total number of global elements produced. This map will not include ghosted nodes.
   */
   GbNodetoOwnedLcNodeMap = new Epetra_Map(-1, NumOwnedLcNodes, DATA_PTR(IndexGbNode_LcNode), 0, *Comm_ptr_);

}
//===========================================================================
void
LocalNodeIndices::initLcBlockNodeMaps()
{
  /*
   *  Create a block map of the local eqns and ghost eqns. The global element ids are the global node numbers.
   *  The number of rows in the block are the number of equations defined at each node.
   *  We rely on Epetra to tally up the total number of global elements produced. This map will included ghosted nodes
   *  and equations
   */
  GbBlockNodeEqnstoLcBlockNodeEqnsColMap = new Epetra_BlockMap(-1, NumLcNodes, DATA_PTR(IndexGbNode_LcNode), DATA_PTR(
                                                               NumEqns_LcNode), 0, *Comm_ptr_);

  /*
   *  Create a map of the local equations only. The global ids are the global node numbers.
   *  We relie on Epetra to tally up the total number of global elements produced.  this map will
   *  not include ghosted nodes, and therefore will be 1 to 1.
   */
  GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap = new Epetra_BlockMap(-1, NumOwnedLcNodes, DATA_PTR(IndexGbNode_LcNode),
                                                                    DATA_PTR(NumEqns_LcNode), 0, *Comm_ptr_);
}
//===========================================================================
// Construct a coloring map for the LcNodes on this processor.
/*
 *   This map includes the external nodes that are defined on the processor
 *   as well.
 *
 *   This coloring map will be used to calculate the jacobian. It's based on
 *   an assumption about the matrix stencil. Every node that is adjacent to
 *   another node or is involved with its residual evaluation receives
 *   a different color. Right now we assume the matrix stencil is strictly
 *   an adjacent node operation.   If this assumption is not
 *   correct, this is the location to change the algorithm.
 *
 *   This routine may be called after the LcNodes map is created.
 */
void
LocalNodeIndices::makeNodeColors()
{
  /*
   *   Making node colors is pretty easy for one d problems. We will assign
   *   colors based on the global node number's remainder wrt division
   *   of 3. Additionally, any global node gets an additional color for a
   *   total of 4 nodal colors.
   */

  std::vector<int> nodeColor(NumLcNodes);
  NumNodeColors = 3;
  for (int iNode = 0; iNode < NumLcNodes; iNode++) {
    /*
     *  Get the global node number
     */
    int gbNode = IndexGbNode_LcNode[iNode];

    int rmd = gbNode % NumNodeColors;
    nodeColor[iNode] = rmd;
  }
  if (GCNIndexLcNode != -1) {
    nodeColor[GCNIndexLcNode] = 3;
    NumNodeColors = 4;
  }

  /*
   *  Now, call the constructor that creates the coloring map
   */
  NodeColorMap = new Epetra_MapColoring(*(GbNodetoLcNodeColMap), DATA_PTR(nodeColor));
}
//===========================================================================
// Construct a coloring map for the LcEqns on this processor.
/*
 *   This map includes the external nodes that are defined on the processor
 *   as well.
 *
 *   This coloring map will be used to calculate the jacobian. It's based on
 *   an assumption about the matrix stencil. Every node that is adjacent to
 *   another node or is involved with its residual evaluation receives
 *   a different color. Right now we assume the matrix stencil is strictly
 *   an adjacent node operation.   If this assumption is not
 *   correct, this is the location to change the algorithm.
 *
 *   This routine may be called after the LcNodes map is created.
 */
void
LocalNodeIndices::makeEqnColors()
{

  /*
   * Get the maximum number of colors in the node map
   */
  int maxNodeColor = NodeColorMap->MaxNumColors();
  AssertTrace(NumNodeColors == maxNodeColor);

  std::vector<int> colorMaxEqnsLocal(maxNodeColor, 0);
  std::vector<int> colorMaxEqns(maxNodeColor, 0);
  std::vector<int> colorStart(maxNodeColor, 0);

  for (int iNode = 0; iNode < NumLcNodes; iNode++) {
    int ncolor = (*NodeColorMap)[iNode];
    NodalVars * nv = NodalVars_LcNode[iNode];
    int neqns = nv->NumEquations;
    colorMaxEqnsLocal[ncolor] = MAX(colorMaxEqnsLocal[ncolor],neqns);
  }

  const Epetra_Comm &cc = NodeColorMap->Comm();
  cc.MaxAll(DATA_PTR(colorMaxEqnsLocal), DATA_PTR(colorMaxEqns), maxNodeColor);
  NumEqnColors = 0;
  for (int i = 0; i < maxNodeColor; i++) {
    colorStart[i] = NumEqnColors;
    NumEqnColors += colorMaxEqns[i];
  }

  std::vector<int> eqnColor(NumLcEqns, -1);

  int iLcEqn = 0;
  for (int iNode = 0; iNode < NumLcNodes; iNode++) {
    NodalVars * nv = NodalVars_LcNode[iNode];
    int numEqns = nv->NumEquations;
    int iNodeColor = (*NodeColorMap)[iNode];
    int cs = colorStart[iNodeColor];
    for (int iEqn = 0; iEqn < numEqns; iEqn++) {
      eqnColor[iLcEqn] = cs + iEqn;
      iLcEqn++;
    }
  }
  /*
   *  Now, call the constructor that creates the coloring map
   */
  EqnColorMap = new Epetra_MapColoring(*(GbBlockNodeEqnstoLcBlockNodeEqnsColMap), DATA_PTR(eqnColor));

  EqnColors = new Epetra_IntVector(*(GbBlockNodeEqnstoLcBlockNodeEqnsColMap));

  int n = EqnColors->MyLength();
  AssertTrace(n == NumLcEqns);
  for (int iEqn = 0; iEqn < NumLcEqns; iEqn++) {
    (*EqnColors)[iEqn] = eqnColor[iEqn];
  }
}
//===========================================================================
// global node to local node mapping
/*
 * Given a global node, this function will return the local node value.
 * If the global node is not on this processor, then this function returns -1.
 *
 * @param gbNode  global node
 * @return returns the local node value
 */
int
LocalNodeIndices::GbNodeToLcNode(const int gbNode) const
{
  // we can do better, but let's get on with it.
  for (int iNode = 0; iNode < NumLcNodes; iNode++) {
    /*
     *  Get the global node number
     */
    int gbNode1 = IndexGbNode_LcNode[iNode];
    if (gbNode == gbNode1) {
      return iNode;
    }
  }
  return -1;
}
//=====================================================================================================================
// Global eqn to local eqn mapping
/*
 * Given a global eqn, this function will return the local eqn value.
 * If the global eqn is not on this processor, then this function returns -1.
 * The local eqn may or may not be owned by this processor.
 *
 * @param gbEqn  global node
 * @return returns the local eqn value
 */
int
LocalNodeIndices::GbEqnToLcEqn(const int gbEqn) const
{
  int rowEqnNum;
  int gbNode = GI_ptr_->GbEqnToGbNode(gbEqn, rowEqnNum);
  int lcnode = GbNodeToLcNode(gbNode);
  if (lcnode >= 0) {
    int lceqn = IndexLcEqns_LcNode[lcnode];
    return (rowEqnNum + lceqn);
  }
  return -1;
}

//=====================================================================================================================
void
LocalNodeIndices::InitializeLocalNodePositions()
{
  /*
   * Here we assume that the positions are correct in the nodal vars position.
   * We call the next routine to make sure that they are propagated to Xpos_LcNode[];
   */
  UpdateNodalVarsPositions();

  // Check to make sure we have set up the vectors correctly
  AssertTrace((*Xpos_LcNode_p)[0] == (*Xpos_LcOwnedNode_p)[0]);
  AssertTrace((*Xpos_LcNode_p)[1] == (*Xpos_LcOwnedNode_p)[1]);

  GI_ptr_->updateGlobalPositions(Xpos_LcOwnedNode_p);
}
//=====================================================================================================================
void
LocalNodeIndices::UpdateNodalVarsPositions()
{
  for (int iNode = 0; iNode < NumLcNodes; iNode++) {
    NodalVars *nv = NodalVars_LcNode[iNode];
    // Initial Spatial position of the node
    X0pos_LcNode[iNode] = nv->x0NodePos();
    (*Xpos_LcNode_p)[iNode] = nv->xNodePos();
  }
}
//=====================================================================================================================
//  Extract the positions from the solution vector and propagate them into all other structures
/*
 *  The following member data are updated:
 *      this->Xpos_LcNode_p[];
 *      nv->XNodePos
 * @param soln
 */
void
LocalNodeIndices::ExtractPositionsFromSolution(const Epetra_Vector * const soln_p)
{
  for (int iNode = 0; iNode < NumLcNodes; iNode++) {
    NodalVars *nv = NodalVars_LcNode[iNode];
    int startN = IndexLcEqns_LcNode[iNode];
    int offset = (nv->Offset_VarType)[Displacement_Axial];
    if (offset >= 0) {
      (*Xpos_LcNode_p)[iNode] = (*soln_p)[startN + offset] + (nv->x0NodePos());
      nv->changeNodePosition((*Xpos_LcNode_p)[iNode]);
    }
  }
}
//=====================================================================================================================
//  We set initial conditions here that make sense from a global perspective.
//  This should be done as a starting point. If there are better answers, it should be overridden.
//
//    Set Displacement_Axial unknowns to 0.
void
LocalNodeIndices::setInitialConditions(const bool doTimeDependentResid,
                                       Epetra_Vector *soln,
                                       Epetra_Vector *solnDot,
                                       const double t,
                                       const double delta_t)
{

  for (int iNode = 0; iNode < NumLcNodes; iNode++) {
    NodalVars *nv = NodalVars_LcNode[iNode];
    double xpos = nv->x0NodePos();
    int startN = IndexLcEqns_LcNode[iNode];
    int offset = (nv->Offset_VarType)[Displacement_Axial];
    if (offset >= 0) {
      (*soln)[startN + offset] = 0.0;
    }
    AssertTrace(xpos == (*Xpos_LcNode_p)[iNode]);
  }
}
//=====================================================================================================================
// Generate the nodal variables structure
/*
 *   This routine will update the pointer to the NodalVars structure
 */
void
LocalNodeIndices::GenerateNodalVars()
{
  for (int iNode = 0; iNode < NumLcNodes; iNode++) {
    int gbnode = IndexGbNode_LcNode[iNode];
    NodalVars_LcNode[iNode] = GI_ptr_->NodalVars_GbNode[gbnode];
  }
  /*
   * Update the positions information in the NodalVars structure
   */
  UpdateNodalVarsPositions();
}
//=====================================================================================================================
int
LocalNodeIndices::UpdateEqnCount()
{
  int numOwnedEqns = 0;
  int numLcEqns = 0;
  int ieqn = 0;
  for (int iLcNode = 0; iLcNode < NumLcNodes; iLcNode++) {
    NodalVars *nv = NodalVars_LcNode[iLcNode];
    NumEqns_LcNode[iLcNode] = nv->NumEquations;
    IndexLcEqns_LcNode[iLcNode] = ieqn;
    if (!IsExternal_LcNode[iLcNode]) {
      numOwnedEqns += NumEqns_LcNode[iLcNode];
    }
    numLcEqns += NumEqns_LcNode[iLcNode];
    ieqn += NumEqns_LcNode[iLcNode];
  }
  NumLcOwnedEqns = numOwnedEqns;
  NumLcEqns = numLcEqns;
  return numOwnedEqns;
}
//=====================================================================================================================
// Generate Equation mapping vectors
/*
 *  This must be called after the total number of local equations on a processor, NumLcEqns, is found out
 */
void
LocalNodeIndices::generateEqnMapping()
{
  IndexGbEqns_LcEqns.resize(NumLcEqns);
  for (int iLcNode = 0; iLcNode < NumLcNodes; iLcNode++) {
    NodalVars *nv = NodalVars_LcNode[iLcNode];
    int nEqns = nv->NumEquations;
    int ilstart = IndexLcEqns_LcNode[iLcNode];
    int gStart = nv->EqnStart_GbEqnIndex;
    for (int i = 0; i < nEqns; i++) {
      IndexGbEqns_LcEqns[ilstart + i] = gStart + i;
    }
  }
}
//=====================================================================================================================
void
LocalNodeIndices::updateGhostEqns(Epetra_Vector * const targSolnV, const Epetra_Vector * const srcSolnV)
{
  // Do compatibility checks to make sure that the Epetra_Vectors are overlap
  // and owned vectors specifically
  bool checkStruct = true;
  if (checkStruct) {
    const Epetra_BlockMap &src_map = *GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap;
    const Epetra_BlockMap &srcS_map = srcSolnV->Map();

    if (src_map.NumMyElements() != srcS_map.NumMyElements()) {
      throw m1d_Error("LocalNodeIndices::updateGhostEqns Error", "srcSoln not compatible");
    }
    if (src_map.NumMyPoints() != srcS_map.NumMyPoints()) {
      throw m1d_Error("LocalNodeIndices::updateGhostEqns Error", "srcSoln not compatible");
    }
    if (src_map.MaxMyGID() != srcS_map.MaxMyGID()) {
      throw m1d_Error("LocalNodeIndices::updateGhostEqns Error", "srcSoln not compatible");
    }
    if (src_map.MinMyGID() != srcS_map.MinMyGID()) {
      throw m1d_Error("LocalNodeIndices::updateGhostEqns Error", "srcSoln not compatible");
    }

    const Epetra_BlockMap &tar_map = *(GbBlockNodeEqnstoLcBlockNodeEqnsColMap);
    const Epetra_BlockMap &tarS_map = targSolnV->Map();
    if (tar_map.NumMyElements() != tarS_map.NumMyElements()) {
      throw m1d_Error("LocalNodeIndices::updateGhostEqns Error", "tarSoln not compatible");
    }
    if (tar_map.NumMyPoints() != tarS_map.NumMyPoints()) {
      throw m1d_Error("LocalNodeIndices::updateGhostEqns Error", "tarSoln not compatible");
    }
    if (tar_map.MaxMyGID() != tarS_map.MaxMyGID()) {
      throw m1d_Error("LocalNodeIndices::updateGhostEqns Error", "tarSoln not compatible");
    }
    if (tar_map.MinMyGID() != tarS_map.MinMyGID()) {
      throw m1d_Error("LocalNodeIndices::updateGhostEqns Error", "tarSoln not compatible");
    }
  }

  /*
   * Use the Ghost Eqn importer created earlier 
   */
  if (Importer_GhostEqns) {
    int code = targSolnV->Import(*srcSolnV, *Importer_GhostEqns, Insert, 0);
    if (code) {
      throw m1d_Error("LocalNodeIndices::updateGhostEqns Error", "import operation failed");
    }
  } else {
    throw m1d_Error("LocalNodeIndices::updateGhostEqns", "ghost importer not made");
  }

  return;
}
//=====================================================================================================================
void
LocalNodeIndices::makeImporters()
{
  /*
   * Create a map to do the ghost node exchanges
   */
  AssertTrace(GbBlockNodeEqnstoLcBlockNodeEqnsColMap != 0);
  const Epetra_BlockMap &targ_map = *GbBlockNodeEqnstoLcBlockNodeEqnsColMap;

  /*
   *  Extract the map for the distributed vector from the object
   */
  AssertTrace(GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap!= 0);
  const Epetra_BlockMap &src_map = *GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap;

  /*
   *  Create an import object that describes the communication
   *  to update ghost equations
   */
  safeDelete(Importer_GhostEqns);
  Importer_GhostEqns = new Epetra_Import(targ_map, src_map);

  /*
   * Example of how to use the importer to updates ghost nodes
   */
  //
  // int code = solnV->Import(*src,
  //                          *Importer_GhostEqns, Insert, 0);
  // if (code) {
  //    throw m1d_Error("Error", "Error");
  // }

  /*
   *  Create an import object that describes the communication
   *  to update ghost node vectors
   */
  AssertTrace(GbNodetoLcNodeColMap != 0);
  const Epetra_BlockMap &targN_map = *GbNodetoLcNodeColMap;

  /*
   *  Extract the map for the distributed vector from the object
   */
  AssertTrace(GI_ptr_->GbNodetoOwnedLcNodeMap != 0);
  const Epetra_BlockMap &srcN_map = *(GI_ptr_->GbNodetoOwnedLcNodeMap);

  /*
   *  Create an import object that describes the communication to bring
   *  all of the unknowns from the distributed object onto processor 0
   */
  safeDelete(Importer_NodalValues);
  Importer_NodalValues = new Epetra_Import(targN_map, srcN_map);

  return;
}
//=====================================================================================================================
//=====================================================================================================================
}
/* End namespace m1d() */
//=====================================================================================================================
