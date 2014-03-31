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
#include "m1d_globals.h"

#include <list> 

namespace m1d
{

//===========================================================================
NodalVars::NodalVars(DomainLayout *dl_ptr) :
      GbNode(-1), NumEquations(0), NumBulkDomains(0), NumSurfDomains(0), EqnStart_GbEqnIndex(-1),
      XNodePos(M1D_DOUBLE_NOTSET), 
      X0NodePos(M1D_DOUBLE_NOTSET), 
      XFracNodePos(M1D_DOUBLE_NOTSET), DL_ptr_(dl_ptr)
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

  GbNode                                  = r.GbNode;
  NumEquations                            = r.NumEquations;
  NumBulkDomains                          = r.NumBulkDomains;
  NumSurfDomains                          = r.NumSurfDomains;
  EqnStart_GbEqnIndex                     = r.EqnStart_GbEqnIndex;
  BulkDomainIndex_BDN                     = r.BulkDomainIndex_BDN;
  OffsetIndex_BulkDomainEqnStart_BDN      = r.OffsetIndex_BulkDomainEqnStart_BDN;
  BulkDomainIndex_fromID                  = r.BulkDomainIndex_fromID;
  SurfDomainIndex_SDN                     = r.SurfDomainIndex_SDN;
  OffsetIndex_SurfDomainEqnStart_SDN      = r.OffsetIndex_SurfDomainEqnStart_SDN;
  SurfDomainIndex_fromID                  = r.SurfDomainIndex_fromID;
  VariableNameList_EqnNum                 = r.VariableNameList_EqnNum;
  VariableSubType_EqnNum                  = r.VariableSubType_EqnNum;
  EquationNameList_EqnNum                 = r.EquationNameList_EqnNum;
  EquationSubType_EqnNum                  = r.EquationSubType_EqnNum;
  Offset_VarType                          = r.Offset_VarType;
  Offset_VarTypeVector                    = r.Offset_VarTypeVector;
  Number_VarType                          = r.Number_VarType;
  Number_VarTypeVector                    = r.Number_VarTypeVector;
  XNodePos                                = r.XNodePos;
  X0NodePos                               = r.X0NodePos;
  XFracNodePos                            = r.XFracNodePos;
  //
  //  Do a shallow pointer copy here. This is appropriate in most cases
  //
  DL_ptr_                                 = r.DL_ptr_;

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
//============================================================================================================
void  insert_into_list(std::list<VarType>& varList, VarType& v1)
{
    //   int ss = varList.size();
    std::list<VarType>::iterator it; 
    for (it = varList.begin(); it != varList.end(); ++it) {
	VarType vtL = *it; 
      if (v1 == vtL) {
	    printf("shouldn't be here\n");
    return;
        }
        if (vtL > v1) {
	    varList.insert(it,v1);
            return;
        }
    }
    varList.push_back(v1);
}
//============================================================================================================
static int lookupPosition(int val, std::vector<int>& yy) {
    for (size_t i = 0; i < yy.size(); i++) {
	if (val == yy[i]) {
	    return (int) i;
	}
    }
    return -1;
}
static void print_line(const char *str, int n)
{
    for (int i = 0; i < n; i++) {
        printf("%s", str);
    }
    printf("\n");
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
    NumEquations = iOffset;

  } else {

      if (NumBulkDomains < 2) {
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
      } else if ( NumBulkDomains == 2) {
 
	  // Two bulk domains.
	  std::list<VarType> varList;
	  std::list<EqnType> eqnList;
	  std::vector<int> EqntoLeftMapping;
	  std::vector<int> EqntoRightMapping;
	  std::vector<int> RightToEqnMapping;
	  std::vector<int> LeftToEqnMapping;
	  int leftBulkDomEqn = 0;
	  int rightBulkDomEqn = 0;
	  bool endOfVars = false;
	  int ibddL = BulkDomainIndex_BDN[0];
	  int ibddR = BulkDomainIndex_BDN[1];
	  int iEqn = 0;
	  BulkDomainDescription *bddL = DL.BulkDomainDesc_global[ibddL];
	  BulkDomainDescription *bddR = DL.BulkDomainDesc_global[ibddR];
	  OffsetIndex_BulkDomainEqnStart_BDN[0] = 0;
	  OffsetIndex_BulkDomainEqnStart_BDN[1] = 0;
	  int numEqnsBulkL = bddL->NumEquationsPerNode;
	  int numEqnsBulkR = bddR->NumEquationsPerNode;
	  int sdi = SurfDomainIndex_SDN[0];
	  SurfDomainDescription *sdd = DL.SurfDomainDesc_global[sdi];
	  do {
	      VarType varL = VarType(Max_Var_Name);
	      EqnType eqnL = EqnType(Max_Eqn_Name);
	      if (leftBulkDomEqn < numEqnsBulkL) {
		  varL = bddL->VariableNameList[leftBulkDomEqn];
		  eqnL = bddL->EquationNameList[leftBulkDomEqn];
	      }
	      VarType varR = VarType(Max_Var_Name);
	      EqnType eqnR = EqnType(Max_Eqn_Name);
	      if (rightBulkDomEqn < numEqnsBulkR) {
		  varR = bddR->VariableNameList[rightBulkDomEqn];
		  eqnR = bddR->EquationNameList[rightBulkDomEqn];
	      }
	      if (varL.VariableType <= varR.VariableType) {
		  {
		     
		      //! Add the left equation and variable
		      varList.push_back(varL);
		      eqnList.push_back(eqnL);
		      LeftToEqnMapping.push_back(iEqn);
		      EqntoLeftMapping.push_back(leftBulkDomEqn);
		      EqntoRightMapping.push_back(-1);
		      //m1d::VAR_TYPE_SUBNUM varSubL = varL.VariableSubType;
		      // increment counters
		      iEqn++;
		      leftBulkDomEqn++;

		  }
	
	      } else {
		  int jLEqn = sdd->RightToLeftBulkEqnMapping[rightBulkDomEqn];
		  if (jLEqn < 0) {
		      varList.push_back(varR);
		      eqnList.push_back(eqnR);
		      RightToEqnMapping.push_back(iEqn);
		      EqntoLeftMapping.push_back(-1);
		      EqntoRightMapping.push_back(rightBulkDomEqn);
		      // increment counters
		      iEqn++;
		      rightBulkDomEqn++;
		  } else {
		      // Will have to fix up cases later
		      int pos = lookupPosition(jLEqn, LeftToEqnMapping);
		      if (pos >= 0) {
			  RightToEqnMapping.push_back(pos);
			  EqntoRightMapping[pos] = rightBulkDomEqn;
		      } else {  
			  RightToEqnMapping.push_back(-1);
		      }
		      rightBulkDomEqn++;
		  }
	      }
	      if ((leftBulkDomEqn >= numEqnsBulkL) &&  (rightBulkDomEqn >= numEqnsBulkR)) {
		  endOfVars = true;
	      }
	  } while (!endOfVars);
	  NumEquations = iEqn;
	  for (int rightBulkDomEqn = 0; rightBulkDomEqn < numEqnsBulkR;  rightBulkDomEqn++) {
	      if (RightToEqnMapping[rightBulkDomEqn] == -1) {
		  int leftBulkDomEqn = sdd->RightToLeftBulkEqnMapping[rightBulkDomEqn];
		  if (leftBulkDomEqn < 0) {
		      printf("logic error\n");
		      exit(-1);
		  }
		  iEqn = LeftToEqnMapping[leftBulkDomEqn];
		  RightToEqnMapping[rightBulkDomEqn] = iEqn;
		  EqntoRightMapping[iEqn] = rightBulkDomEqn;
	      }
	  }
  
	  std::list<VarType>::iterator itVList;
	  std::list<EqnType>::iterator itEList;
	  for (itVList = varList.begin(), itEList = eqnList.begin(); itVList != varList.end(); ++itVList, ++itEList) {
	      VariableNameList_EqnNum.push_back(*itVList);
	      EquationNameList_EqnNum.push_back(*itEList);
	      VariableSubType_EqnNum.push_back(itVList->VariableSubType);
	      EquationSubType_EqnNum.push_back(itEList->EquationSubType);
	  }

	  for (int isd = 0; isd < NumSurfDomains; isd++) {
	      int isdd = SurfDomainIndex_SDN[isd];
	      SurfDomainDescription *sdd = DL.SurfDomainDesc_global[isdd];
	      int numEqnsSurf = sdd->NumEquationsPerNode;
	      OffsetIndex_SurfDomainEqnStart_SDN[isd] = NumEquations;
	      for (int iEqnSurf = 0; iEqnSurf < numEqnsSurf; iEqnSurf++) {
		  VariableNameList_EqnNum.push_back(sdd->VariableNameList[iEqnSurf]);
		  EquationNameList_EqnNum.push_back(sdd->EquationNameList[iEqnSurf]);
		  VariableSubType_EqnNum.push_back(VariableNameList_EqnNum[NumEquations].VariableSubType);
		  EquationSubType_EqnNum.push_back(EquationNameList_EqnNum[NumEquations].EquationSubType);
		  EqntoRightMapping.push_back(-1);
		  EqntoLeftMapping.push_back(-1);
		  NumEquations++;
	      }
	  }
	  bool printResult = false;
#ifdef DEBUG_MODE
          if (m1d::s_printLvl_DebugTables) {
              printResult = true;
          }
#endif
	  //
	  //  With all complicated logic we best have to have a print out section to see the results
	  //
	  if (printResult) {
	      print_line("-", 152);
	      printf("    Equation Order for Global Node %d", GbNode);
	      printf("          numbulkdomains = %d", NumBulkDomains);
	      printf("          numsurfdomains = %d\n", NumSurfDomains);
	      printf("                 EQUATION                         |               "
		     "  LEFT BULK DOMAIN                |              RIGHT BULK DOMAIN                   |\n");
	      printf("     eqnNum             Var                 Eqn   |    eqnNum          Var  "
		     "                  Eqn   |     eqnNum           Var               Eqn       |\n");
	      for (int iEqn = 0; iEqn <  NumEquations; iEqn++) {
		  VarType var =  VariableNameList_EqnNum[iEqn];
		  EqnType eqn =  EquationNameList_EqnNum[iEqn];
		  string varS = var.VariableName(20);
		  string eqnS = eqn.EquationName(20);
		  int rightEqn = EqntoRightMapping[iEqn];
		  int leftEqn =  EqntoLeftMapping[iEqn];

		  printf(" %5d %20.20s  %20.20s |", iEqn, varS.c_str(), eqnS.c_str());
		  if (leftEqn >= 0) {
		      VarType varL =  bddL->VariableNameList[leftEqn];
		      EqnType eqnL =  bddL->EquationNameList[leftEqn];
		      string varLS = varL.VariableName(20);
		      string eqnLS = eqnL.EquationName(20);
		      printf("%5d %20.20s  %20.20s |", leftEqn, varLS.c_str(), eqnLS.c_str());
		  } else {
		      printf("                                                 |");
		  }
		  if (rightEqn >= 0) {
		      VarType varR =  bddR->VariableNameList[rightEqn];
		      EqnType eqnR =  bddR->EquationNameList[rightEqn];
		      string varRS = varR.VariableName(20);
		      string eqnRS = eqnR.EquationName(20);
		      printf(" %5d %20.20s  %20.20s |", rightEqn, varRS.c_str(), eqnRS.c_str());
		  } else {
		      printf("                                                  |");
		  }
		  printf("\n");
	      }
	      print_line("-", 152);
	  }

      } else {
	  printf("NUmbers of bulk domains greater than 2 hasn't been coded yet\n");
	  exit(-1);
      }

  } // End of Complex method
  /*
   *  Ok, now we have finally calculated how many equations are at this node.
   *  We still don't know what EqnStart_LcEqnIndex is until we sum up the
   *  equations on a local processor and global unknown basis.
   */
 
  //
  //  Create Offset_VarType array
  //
  Offset_VarType.clear();
  for (VAR_TYPE iv = Variable_Type_Unspecified; iv < Max_Var_Name; ++iv) {
    Offset_VarType[iv] = npos;
    if (iv == Variable_Type_Unspecified) {
      Offset_VarType[iv] = 0;
    } else {
      for (int iOffset = 0; iOffset < NumEquations; iOffset++) {
        VAR_TYPE cv = VariableNameList_EqnNum[iOffset].VariableType;
        if (cv == iv) {
          if (Offset_VarType[iv] == npos) {
            Offset_VarType[iv] = iOffset;
          }
        }
      }
    }
  }

  //
  //  Create Offset_EqnType array
  //
  Offset_EqnType.clear();
  for (EQ_TYPE iv = Equation_Type_Unspecified; iv < Max_Eqn_Name; ++iv) {
    Offset_EqnType[iv] = npos;
    if (iv == Equation_Type_Unspecified) {
      Offset_EqnType[iv] = 0;
    } else {
      for (int iOffset = 0; iOffset < NumEquations; iOffset++) {
        EQ_TYPE cv = EquationNameList_EqnNum[iOffset].EquationType;
        if (cv == iv) {
            if (Offset_EqnType[iv] == npos) {
              Offset_EqnType[iv] = iOffset;
            }
        }
        if (iv == Species_Eqn_Offset) {
          if ((cv == Species_Conservation) || 
	      (cv == MoleFraction_Summation) ||
              (cv == ChargeNeutrality_Summation)) {
              if (Offset_EqnType[iv] == npos) {
                  Offset_EqnType[iv] = iOffset;
              }
          }
        }
      }
    }
  }
  // 
  //  Create Number_VarType
  //
  for (VAR_TYPE iv = Variable_Type_Unspecified; iv < Max_Var_Name; ++iv) {
    if (iv == Variable_Type_Unspecified) {
        Number_VarType[iv] = 0;
    } else {
        int voff = Offset_VarType[iv];
        if (voff < 0) {
            Number_VarType[iv] = 0;
        } else {
            int count = 0;
            for (size_t index = voff; index < VariableNameList_EqnNum.size(); index++) {
                VarType vt = VariableNameList_EqnNum[index];
                if (vt.VariableType != iv) {
                    break;
                } else {
                    count++;
                }
            }
            Number_VarType[iv] = count;
        }
    }
  }
  // 
  //  Create Number_EqnType
  //
  for (EQ_TYPE iv = Equation_Type_Unspecified; iv < Max_Eqn_Name; ++iv) {
    if (iv == Equation_Type_Unspecified) {
        Number_EqnType[iv] = 0;
    } else {
        int voff = Offset_EqnType[iv];
        if (voff < 0) {
            Number_EqnType[iv] = 0;
        } else {
            int count = 0;
            for (size_t index = voff; index < EquationNameList_EqnNum.size(); index++) {
                EqnType vt = EquationNameList_EqnNum[index];
                if (vt.EquationType != iv) {
                    break;
                } else {
                    count++;
                }
            }
            Number_EqnType[iv] = count;
        }
    }
  }
  //
  //  Create Offset_VarTypeVector and Number_VarTypeVector
  //
  size_t ss = Max_Var_Name;
  Offset_VarTypeVector.resize(ss, npos);
  Number_VarTypeVector.resize(ss, npos);
  for (VAR_TYPE iv = Displacement_Axial; iv < Max_Var_Name; ++iv) {
    size_t siv = (size_t) iv;
    if (Offset_VarType[iv] >= 0) {
         Offset_VarTypeVector[siv] = (size_t) Offset_VarType[iv];
    } else {
         Offset_VarTypeVector[siv] = npos;
    }
    if (Number_VarType[iv] >= 0) {
         Number_VarTypeVector[siv] = (size_t) Number_VarType[iv];
    } else {
         Number_VarTypeVector[siv] = npos;
    }
  }

  //
  //  Create Offset_EqnTypeVector and Number_EqnTypeVector
  //
  ss = Max_Eqn_Name;
  Offset_EqnTypeVector.resize(ss, npos);
  Number_EqnTypeVector.resize(ss, npos);
  for (EQ_TYPE iv = Momentum_Axial; iv < Max_Eqn_Name; ++iv) {
    size_t siv = (size_t) iv;
    if (Offset_EqnType[iv] >= 0) {
         Offset_EqnTypeVector[siv] = (size_t) Offset_EqnType[iv];
    } else {
         Offset_EqnTypeVector[siv] = npos;
    }
    if (Number_EqnType[iv] >= 0) {
         Number_EqnTypeVector[siv] = (size_t) Number_EqnType[iv];
    } else {
         Number_EqnTypeVector[siv] = npos;
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
//
//  variableType and variableSubType completely determine a variable within the code.
//  This means that mole fractions in different domains that are different should have different
//  VAR_TYPE_SUBNUM values. This has to be detemined by the application.
//
size_t NodalVars::indexBulkDomainVar(VAR_TYPE variableType, VAR_TYPE_SUBNUM variableSubType) const
{
    std::map<VAR_TYPE, size_t>::const_iterator it;
    if ((it = Offset_VarType.find(variableType)) == Offset_VarType.end()) {
      return npos;
    } 
    size_t start = it->second;
#ifdef DEBUG_MODE
    if ((it = Number_VarType.find(variableType)) == Number_VarType.end()) {
        return npos;
    } 
    size_t numV = it->second;
#else
    size_t numV = Number_VarType.find(variableType)->second;
#endif
    //
    //  We do a quick trial for a common case.
    if (variableSubType < (int) numV) {
	VarType vt = VariableNameList_EqnNum[start + variableSubType];
	if (variableSubType == vt.VariableSubType) {
	    return start + variableSubType;
	}
    }
    for (size_t i = 0; i < numV; i++) {
	VarType vt = VariableNameList_EqnNum[start + i];
	if (variableSubType == vt.VariableSubType) {
	    return start + i;
	}
    }
    return npos;
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
