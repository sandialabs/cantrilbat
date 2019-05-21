/**
 * @file m1d_GlobalIndices.cpp
 *  Definitions for the structure that fills up and maintains the global indexing for the problem
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_GlobalIndices.h"
#include "m1d_DomainLayout.h"
#include "m1d_NodalVars.h"
#include "m1d_exception.h"
#include "m1d_EpetraExtras.h"


#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "zuzax/base/stringUtils.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
GlobalIndices::GlobalIndices(Epetra_Comm* comm_ptr) :
    Comm_ptr_(comm_ptr), 
    NumProc(1),
    MyProcID(0), 
    NumGbNodes(0), 
    NumGbNodes_s(0), 
    NumGbEqns(0), 
    NumGbEqns_s(0), 
    NumOwnedLcNodes(0), 
    NumOwnedLcEqns(0),
    GbNodetoOwnedLcNodeMap(0),
    GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap(0),
    XNodePos_GbNode(0),
    GbEqnstoAllMap(0),
    SolnAll(0),
    SolnIntAll(0),
    SolnDotAll(0),
    DL_ptr_(nullptr)
{
}
//==================================================================================================================================
GlobalIndices::~GlobalIndices()
{
    safeDelete(GbNodetoOwnedLcNodeMap);
    safeDelete(GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap);
    for (size_t iLcNodes = 0; iLcNodes < NumGbNodes_s; iLcNodes++) {
        delete NodalVars_GbNode[iLcNodes];
        NodalVars_GbNode[iLcNodes] = nullptr;
    }
    safeDelete(XNodePos_GbNode);
    safeDelete(GbEqnstoAllMap);
    safeDelete(SolnAll);
    safeDelete(SolnIntAll);
    safeDelete(SolnDotAll);
    DL_ptr_ = nullptr;
}
//==================================================================================================================================
GlobalIndices::GlobalIndices(const GlobalIndices& r) :
    GlobalIndices(r.Comm_ptr_)
{
    *this = r;
}
//==================================================================================================================================
GlobalIndices& GlobalIndices::operator=(const GlobalIndices& r)
{
    if (this == &r) {
        return *this;
    }
   
    for (size_t iLcNodes = 0; iLcNodes < NumGbNodes_s; iLcNodes++) {
        safeDelete(NodalVars_GbNode[iLcNodes]);
    }

    Comm_ptr_ = r.Comm_ptr_;
    NumProc = r.NumProc;
    MyProcID = r.MyProcID;
    NumGbNodes = r.NumGbNodes;
    NumGbNodes_s = r.NumGbNodes_s;
    NumGbEqns = r.NumGbEqns;
    NumGbEqns_s = r.NumGbEqns_s;
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
    for (size_t iLcNodes = 0; iLcNodes < NumGbNodes_s; iLcNodes++) {
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
//==================================================================================================================================
void
GlobalIndices::init(DomainLayout* dl_ptr)
{
    DomainLayout& DL = *dl_ptr;
    DL_ptr_ = dl_ptr;
    // We get the total number of global nodes from the DomainLayout object
    NumGbNodes = DL.NumGbNodes;
    NumGbNodes_s = DL.NumGbNodes;
    NodalVars_GbNode.resize(NumGbNodes_s);
    NumEqns_GbNode.resize(NumGbNodes_s);
    IndexStartGbEqns_GbNode.resize(NumGbNodes_s);
    int* tmpIndex = new int[NumGbNodes];
    for (size_t iGbNode = 0; iGbNode < NumGbNodes_s; iGbNode++) {
        AssertTrace(NodalVars_GbNode[iGbNode] == 0);
        NodalVars_GbNode[iGbNode] = new NodalVars(iGbNode, dl_ptr);
        tmpIndex[iGbNode] = (int) iGbNode;
    }
    /*
     * Create a map with all of the doubles on all processors
     */
    Epetra_Map* eAll_map = new Epetra_Map(NumProc * NumGbNodes, NumGbNodes, tmpIndex, 0, *Comm_ptr_);
    XNodePos_GbNode = new Epetra_Vector(*eAll_map, true);
    safeDelete(eAll_map);
    delete [] tmpIndex;
}
//==================================================================================================================================
size_t GlobalIndices::discoverNumEqnsPerNode()
{
    size_t numTotEqns = 0;
    for (size_t iGbNode = 0; iGbNode < NumGbNodes_s; iGbNode++) {
        NodalVars* nv = NodalVars_GbNode[iGbNode];
        nv->DiscoverDomainsAtThisNode();
        nv->GenerateEqnOrder();

        NumEqns_GbNode[iGbNode] = nv->NumEquations;
        IndexStartGbEqns_GbNode[iGbNode] = numTotEqns;
        nv->EqnStart_GbEqnIndex = numTotEqns;
        numTotEqns += nv->NumEquations;
    }
    NumGbEqns = numTotEqns;
    NumGbEqns_s = numTotEqns;
    return numTotEqns;
}
//==================================================================================================================================
void GlobalIndices::procDivide()
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
//==================================================================================================================================
//   Form a Map to describe the distribution of our owned nodes on the processors
void
GlobalIndices::initNodeMaps()
{
    AssertTrace(Comm_ptr_);
    GbNodetoOwnedLcNodeMap = new Epetra_Map(NumGbNodes, NumOwnedLcNodes, 0, *Comm_ptr_);
}
//==================================================================================================================================
void GlobalIndices::initBlockNodeMaps()
{ 
    int igNodeStart = IndexStartGbNode_Proc[MyProcID];
    int *numEqns_LcNode = &NumEqns_GbNode[igNodeStart];
    int* myGlobalNodes = GbNodetoOwnedLcNodeMap->MyGlobalElements();

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
//==================================================================================================================================
void
GlobalIndices::InitMesh()
{
    DL_ptr_->InitializeXposNodes(this);
}
//==================================================================================================================================
// This utility function will return the global node number given the global equation number
int
GlobalIndices::GbEqnToGbNode(const int GbEqnNum, int& rowEqnNum) const
{
#ifdef DEBUG_MODE
    if ((GbEqnNum < 0) || GbEqnNum >= NumGbEqns) {
        throw m1d_Error("GlobalIndices::GbEqnToGbNode", "GbEqnNum out of bounds: " + Zuzax::int2str(GbEqnNum));
    }
#endif
    int nfound = -1;
    for (int i = 1; i < NumGbNodes; i++) {
        if (IndexStartGbEqns_GbNode[i] > GbEqnNum) {
            nfound = i - 1;
            break;
        }
    }
    if (nfound >= 0) {
        rowEqnNum = GbEqnNum - IndexStartGbEqns_GbNode[nfound];
    } else {
        rowEqnNum = -1;
    }
    return nfound;
#ifdef BINSEARCH
    int jBot, jTop, jTrial;
    int nfound = -1;
    jBot = 0;
    jTop = NumGbNodes - 1;
    while ((jTop - jBot) > 1) {
        jTrial = (jTop + jBot) / 2;
        if (GbEqnNum < IndexStartGbEqns_GbNode[jTrial]) {
            jTop = jTrial;
        } else {
            jBot = jTrial;
        }
    }
    nfound = jBot;
    rowEqnNum = GbEqnNum - IndexStartGbEqns_GbNode[nfound];
#endif
}
//==================================================================================================================================
void GlobalIndices::updateGlobalPositions(const Epetra_Vector* const Xpos_LcOwnedNode_p)
{
    gather_nodeV_OnAll(*XNodePos_GbNode, *Xpos_LcOwnedNode_p, Comm_ptr_);
    for (int iGbNode = 0; iGbNode < NumGbNodes; iGbNode++) {
        NodalVars* nv = NodalVars_GbNode[iGbNode];
        nv->changeNodePosition((*XNodePos_GbNode)[iGbNode]);
    }
}
//==================================================================================================================================
int GlobalIndices::GbNodeToLocalNodeNum( const int GbNodeNum, int& procNum) const
{
    int iLocalNodeNum;
    for (int iproc = 0; iproc <  NumProc - 1; ++iproc) {
       if (GbNodeNum < IndexStartGbNode_Proc[iproc+1]) {
           procNum = iproc;
           iLocalNodeNum = GbNodeNum - IndexStartGbNode_Proc[procNum];
           return iLocalNodeNum;
       }
    }
    procNum = NumProc - 1;
    iLocalNodeNum = GbNodeNum - IndexStartGbNode_Proc[procNum];
    return iLocalNodeNum; 
}
//==================================================================================================================================
int GlobalIndices::GbEqnNumToLocEqnNum( const int GbEqnNum, int& procNum) const
{
    int iLocalEqnNum;
    for (int iproc = 0; iproc <  NumProc - 1; ++iproc) {
       if (GbEqnNum < IndexStartGbEqns_Proc[iproc+1]) {
           procNum = iproc;
           iLocalEqnNum = GbEqnNum - IndexStartGbEqns_Proc[procNum];
           return iLocalEqnNum;
       }
    }
    procNum = NumProc - 1;
    iLocalEqnNum = GbEqnNum - IndexStartGbEqns_Proc[procNum];
    return iLocalEqnNum; 
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
