/**
 * @file m1d_NodalVars.cpp
 *
 */

/*
 *  $Id: m1d_NodalVars.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_defs.h"
#include "m1d_NodalVars.h"
#include "m1d_LocalNodeIndices.h"
#include "m1d_DomainLayout.h"
#include "m1d_BulkDomainDescription.h"
#include "m1d_SurfDomainDescription.h"

namespace m1d
{

//===========================================================================
NodalVars::NodalVars(DomainLayout *dl_ptr) :
  GbNode(-1), NumEquations(0), NumBulkDomains(0), NumSurfDomains(0), EqnStart_GbEqnIndex(-1),
      XNodePos(M1D_DOUBLE_NOTSET), X0NodePos(M1D_DOUBLE_NOTSET), XFracNodePos(M1D_DOUBLE_NOTSET), DL_ptr_(dl_ptr)
{
}
//===========================================================================
NodalVars::NodalVars(const int gbnode, DomainLayout *dl_ptr) :
  GbNode(gbnode), NumEquations(0), NumBulkDomains(0), NumSurfDomains(0), EqnStart_GbEqnIndex(-1),
      XNodePos(M1D_DOUBLE_NOTSET), X0NodePos(M1D_DOUBLE_NOTSET), XFracNodePos(M1D_DOUBLE_NOTSET), DL_ptr_(dl_ptr)

{
}
//===========================================================================
NodalVars::~NodalVars()
{
}
//===========================================================================
NodalVars::NodalVars(const NodalVars &r) :
  NumEquations(0), NumBulkDomains(0), NumSurfDomains(0), EqnStart_GbEqnIndex(-1), XNodePos(M1D_DOUBLE_NOTSET),
      X0NodePos(M1D_DOUBLE_NOTSET), XFracNodePos(M1D_DOUBLE_NOTSET)
{
  *this = r;
}
//===========================================================================
NodalVars &
NodalVars::operator=(const NodalVars &r)
{
  if (this == &r)
    return *this;

  GbNode = r.GbNode;
  NumEquations = r.NumEquations;
  NumBulkDomains = r.NumBulkDomains;
  NumSurfDomains = r.NumSurfDomains;
  EqnStart_GbEqnIndex = r.EqnStart_GbEqnIndex;
  BulkDomainIndex_BDN = r.BulkDomainIndex_BDN;
  OffsetIndex_BulkDomainEqnStart_BDN = r.OffsetIndex_BulkDomainEqnStart_BDN;
  SurfDomainIndex_SDN = r.SurfDomainIndex_SDN;
  OffsetIndex_SurfDomainEqnStart_SDN = r.OffsetIndex_SurfDomainEqnStart_SDN;
  VariableNameList_EqnNum = r.VariableNameList_EqnNum;
  VariableSubType_EqnNum = r.VariableSubType_EqnNum;
  EquationNameList_EqnNum = r.EquationNameList_EqnNum;
  EquationSubType_EqnNum = r.EquationSubType_EqnNum;
  Offset_VarType = r.Offset_VarType;
  XNodePos = r.XNodePos;
  X0NodePos = r.X0NodePos;
  XFracNodePos = r.XFracNodePos;
  //XcellBoundRight = r.XcellBoundRight;
  //XcellBoundLeft = r.XcellBoundLeft;

  return *this;
}
//===========================================================================

void
NodalVars::DiscoverDomainsAtThisNode()
{
  DomainLayout &DL = *DL_ptr_;
  /*
   * Fill in the vector, BulkDomainIndex_BDN;
   */
  NumBulkDomains = 0;
  for (int ibdd = 0; ibdd < DL.NumBulkDomains; ibdd++) {
    BulkDomainDescription *bdd = DL.BulkDomainDesc_global[ibdd];
    int bid = bdd->ID();
    int firstGbNode = bdd->FirstGbNode;
    int lastGbNode = bdd->LastGbNode;
    if (firstGbNode <= GbNode && lastGbNode >= GbNode) {
      NumBulkDomains++;
      BulkDomainIndex_BDN.resize(NumBulkDomains);
      BulkDomainIndex_BDN[NumBulkDomains - 1] = bid;
      BulkDomainIndex_fromID[bid] = NumBulkDomains - 1;
    } else {
      BulkDomainIndex_fromID[bid] = -1;
    }
  }

  /*
   * Fill in the vector, SurfDomainIndex_SDN;
   */
  NumSurfDomains = 0;
  for (int isdd = 0; isdd < DL.NumSurfDomains; isdd++) {
    SurfDomainDescription *sdd = DL.SurfDomainDesc_global[isdd];
    int sid = sdd->ID();
    int firstGbNode = sdd->LocGbNode;
    if (firstGbNode == GbNode) {
      NumSurfDomains++;
      SurfDomainIndex_SDN.resize(NumSurfDomains);
      SurfDomainIndex_SDN[NumSurfDomains - 1] = sid;
      SurfDomainIndex_fromID[sid] = NumSurfDomains - 1;
    } else {
      SurfDomainIndex_fromID[sid] = -1;
    }
  }

}
//===========================================================================
// This routine discovers and order the equation unknowns at a node
/*
 * This is the ONE PLACE IN THE CODE that discovers and orders the
 * equations at a node.
 */
void
NodalVars::GenerateEqnOrder()
{
  //LocalNodeIndices &LI = *LI_ptr;
  DomainLayout &DL = *DL_ptr_;
  /*
   * First we decide whether simple order is true. Simple ordering
   * would be to order the equations first by bulk domain and then
   * by surface domain unknowns.
   */
  bool simpleOrdering = true;
  if (NumSurfDomains > 0 && NumBulkDomains > 1) {
    int sdi = SurfDomainIndex_SDN[0];
    SurfDomainDescription * sdd = DL.SurfDomainDesc_global[sdi];
    int bdi0 = BulkDomainIndex_BDN[0];
    int bdi1 = BulkDomainIndex_BDN[1];
    BulkDomainDescription * bdd0 = DL.BulkDomainDesc_global[bdi0];
    BulkDomainDescription * bdd1 = DL.BulkDomainDesc_global[bdi1];
    DomainDescription * ddL = sdd->LeftDomain;
    BulkDomainDescription *bddL = dynamic_cast<BulkDomainDescription *> (ddL);
    if (!bddL) {
      throw m1d_Error("NodalVars::GenerateEqnOrder", "Not handled yet");
    }
    BulkDomainDescription * bddR = sdd->RightBulk;
    if (bdd0 != bddL) {
      printf("confused\n");
      throw m1d_Error("ReconcileEqnOrder", "confused");
    }
    if (bdd1 != bddR) {
      printf("confused\n");
      throw m1d_Error("ReconcileEqnOrder", "confused");
    }
    int numEqnsRightBulk = bdd1->NumEquationsPerNode;
    for (int iEqn = 0; iEqn < numEqnsRightBulk; iEqn++) {
      int j = sdd->RightToLeftBulkEqnMapping[iEqn];
      if (j != -1) {
        simpleOrdering = false;
      }
    }
  }

  /*
   * Resize the offset indices
   */
  OffsetIndex_BulkDomainEqnStart_BDN.resize(NumBulkDomains, -1);
  OffsetIndex_SurfDomainEqnStart_SDN.resize(NumSurfDomains, -1);

  int iOffset = 0;
  VariableNameList_EqnNum.clear();
  EquationNameList_EqnNum.clear();
  VariableSubType_EqnNum.clear();
  EquationSubType_EqnNum.clear();

  if (simpleOrdering) {

    //! Number of bulk domains located at the node

    for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
      int ibdd = BulkDomainIndex_BDN[ibd];
      BulkDomainDescription *bdd = DL.BulkDomainDesc_global[ibdd];
      int numEqnsBulk = bdd->NumEquationsPerNode;
      OffsetIndex_BulkDomainEqnStart_BDN[ibd] = iOffset;
      for (int iEqn = 0; iEqn < numEqnsBulk; iEqn++) {
        VariableNameList_EqnNum.push_back(bdd->VariableNameList[iEqn]);
        EquationNameList_EqnNum.push_back(bdd->EquationNameList[iEqn]);
        VariableSubType_EqnNum.push_back(VariableNameList_EqnNum[iOffset].VariableSubType);
        EquationSubType_EqnNum.push_back(EquationNameList_EqnNum[iOffset].EquationSubType);
        iOffset++;
      }
    }

    for (int isd = 0; isd < NumSurfDomains; isd++) {
      int isdd = SurfDomainIndex_SDN[isd];
      SurfDomainDescription *sdd = DL.SurfDomainDesc_global[isdd];
      int numEqnsSurf = sdd->NumEquationsPerNode;
      OffsetIndex_SurfDomainEqnStart_SDN[isd] = iOffset;
      for (int iEqn = 0; iEqn < numEqnsSurf; iEqn++) {
        VariableNameList_EqnNum.push_back(sdd->VariableNameList[iEqn]);
        EquationNameList_EqnNum.push_back(sdd->EquationNameList[iEqn]);
        VariableSubType_EqnNum.push_back(VariableNameList_EqnNum[iOffset].VariableSubType);
        EquationSubType_EqnNum.push_back(EquationNameList_EqnNum[iOffset].EquationSubType);
        iOffset++;
      }
    }

  } else {
    SurfDomainDescription * sdd = 0;
    for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
      int ibdd = BulkDomainIndex_BDN[ibd];
      BulkDomainDescription *bdd = DL.BulkDomainDesc_global[ibdd];
      int numEqnsBulk = bdd->NumEquationsPerNode;
      OffsetIndex_BulkDomainEqnStart_BDN[ibd] = iOffset;
      bool domainEqnContract = false;
      if (ibd == 1) {
        domainEqnContract = true;
        int sdi = SurfDomainIndex_SDN[0];
        sdd = DL.SurfDomainDesc_global[sdi];
      }
      for (int iEqn = 0; iEqn < numEqnsBulk; iEqn++) {
        /*
         *  If the right bulk equation is to be mapped into a
         *  left bulk equation at the domain, then we do not add a separate
         *  equation for it at the interface. Instead, we skip to the 
         *  next equation
         */
        if (domainEqnContract) {
          int j = sdd->RightToLeftBulkEqnMapping[iEqn];
          if (j != -1)
            continue;
        }
        VariableNameList_EqnNum.push_back(bdd->VariableNameList[iEqn]);
        EquationNameList_EqnNum.push_back(bdd->EquationNameList[iEqn]);
        VariableSubType_EqnNum.push_back(VariableNameList_EqnNum[iOffset].VariableSubType);
        EquationSubType_EqnNum.push_back(EquationNameList_EqnNum[iOffset].EquationSubType);
        iOffset++;
      }
    }

    for (int isd = 0; isd < NumSurfDomains; isd++) {
      int isdd = SurfDomainIndex_SDN[isd];
      SurfDomainDescription *sdd = DL.SurfDomainDesc_global[isdd];
      int numEqnsSurf = sdd->NumEquationsPerNode;
      OffsetIndex_SurfDomainEqnStart_SDN[isd] = iOffset;
      for (int iEqn = 0; iEqn < numEqnsSurf; iEqn++) {
        VariableNameList_EqnNum.push_back(sdd->VariableNameList[iEqn]);
        EquationNameList_EqnNum.push_back(sdd->EquationNameList[iEqn]);
        VariableSubType_EqnNum.push_back(VariableNameList_EqnNum[iOffset].VariableSubType);
        EquationSubType_EqnNum.push_back(EquationNameList_EqnNum[iOffset].EquationSubType);

        iOffset++;
      }
    }
  }
  /*
   *  Ok, now we can finally calculate how many equations are at this node
   *  We still don't know what EqnStart_LcEqnIndex is until we sum up the
   *  equations on a local processor and global unknown basis.
   */
  NumEquations = iOffset;

  Offset_VarType.clear();
  for (VAR_TYPE iv = Variable_Type_Unspecified; iv < Max_Var_Name; ++iv) {
    Offset_VarType[iv] = -1;
    if (iv == Variable_Type_Unspecified) {
      Offset_VarType[iv] = 0;
    } else {
      for (int iOffset = 0; iOffset < NumEquations; iOffset++) {
        VAR_TYPE cv = VariableNameList_EqnNum[iOffset].VariableType;
        if (cv == iv) {
          if (Offset_VarType[iv] == -1) {
            Offset_VarType[iv] = iOffset;
          }
        }
      }
    }
  }
}
//=====================================================================================================================
// Generate a name for the kth variable at the node
/*
 *
 * @param k  input the index of the variable on the local node
 * @return  Returns a string representing the variable.
 */
std::string
NodalVars::VariableName(int k)
{
  VarType kv = VariableNameList_EqnNum[k];
  return kv.VariableName(200);
}
//=====================================================================================================================
//! Change the node position
/*!
 * @param xNodePos  new value of the node position
 */
void
NodalVars::changeNodePosition(double xNodePos)
{
  XNodePos = xNodePos;
}
//=====================================================================================================================
//! Set up the initial node positions
/*!
 * @param xNodePos
 * @param x0NodePos
 * @param xFracNodePos
 */
void
NodalVars::setupInitialNodePosition(double x0NodePos, double xFracNodePos)
{
  XNodePos = x0NodePos;
  X0NodePos = x0NodePos;
  XFracNodePos = xFracNodePos;
}
//=====================================================================================================================
//! returns the node position
double
NodalVars::xNodePos() const
{
  return XNodePos;
}
//=====================================================================================================================
//! returns the initial node position
double
NodalVars::x0NodePos() const
{
  return X0NodePos;
}
//=====================================================================================================================
//! Return the fraction node position from left boundary
double
NodalVars::xFracNodePos() const
{
  return XFracNodePos;
}
//=====================================================================================================================
}
//=====================================================================================================================
