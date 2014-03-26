/**
 * @file m1d_BulkDomainDescription.cpp
 *
 */

/*
 *  $Id: m1d_BulkDomainDescription.cpp 567 2013-03-21 23:03:11Z hkmoffa $
 */
#include "m1d_BulkDomainDescription.h"
#include "m1d_NodalVars.h"
#include "m1d_GlobalIndices.h"
#include "m1d_exception.h"
#include "m1d_DomainLayout.h"

#ifndef MIN
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif
//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
  BulkDomainDescription::BulkDomainDescription(DomainLayout *dl_ptr, std::string domainName) :
    DomainDescription(dl_ptr, domainName), IDBulkDomain(-1), FirstGbNode(-1), LastGbNode(-1), Xpos_end(0.0), BulkDomainPtr_(0)
{

}
//=====================================================================================================================
BulkDomainDescription::~BulkDomainDescription()
{
}
//=====================================================================================================================
BulkDomainDescription::BulkDomainDescription(const BulkDomainDescription &r) :
  DomainDescription(DL_ptr_), IDBulkDomain(r.IDBulkDomain)
{
  *this = r;
}
//=====================================================================================================================
BulkDomainDescription &
BulkDomainDescription::operator=(const BulkDomainDescription &r)
{
  if (this == &r)
    return *this;

  DomainDescription::operator=(r);

  IDBulkDomain = r.IDBulkDomain;
  FirstGbNode = r.FirstGbNode;
  LastGbNode = r.LastGbNode;
  Xpos_start = r.Xpos_start;
  Xpos_end = r.Xpos_end;
  // These are shallow copies and may be wrong. It depends on the
  // context
  LeftSurf = r.LeftSurf;
  RightSurf = r.RightSurf;
  DL_ptr_ = r.DL_ptr_;
  BulkDomainPtr_ = r.BulkDomainPtr_;

  return *this;
}
//=====================================================================================================================
// Sets the id of the surface domain
void
BulkDomainDescription::setID(int id)
{
  IDBulkDomain = id;
}
//=====================================================================================================================
// Reports the id of the surface domain
int
BulkDomainDescription::ID() const
{
  return IDBulkDomain;
}
//=====================================================================================================================
// Specify the left and right surface domains for this
// bulk domain
void
BulkDomainDescription::setAdjSurfDomains(SurfDomainDescription *leftSurf, SurfDomainDescription *rightSurf)
{
  LeftSurf = leftSurf;
  RightSurf = rightSurf;
}
//=====================================================================================================================
// Specify the left and right global nodes for this bulk domain
/*
 *
 * @param leftGbNode   Left value of the node
 * @param rightGbNode  Right value of the node
 */
void
BulkDomainDescription::setGbNodeBounds(const int leftGbNode, const int rightGbNode)
{
  FirstGbNode = leftGbNode;
  LastGbNode = rightGbNode;
}
//=====================================================================================================================
// Specify the left and right X Positions of the domain
/*
 * @param xleft   Left value of the domain boundary
 * @param xright  right value of the domain boundary
 */
void
BulkDomainDescription::setXposBounds(const double xleft, const double xright)
{
  Xpos_start = xleft;
  Xpos_end = xright;
}
//=====================================================================================================================
//! Initialize the positions of the nodes in the domain
void
BulkDomainDescription::InitializeXposNodes(GlobalIndices *gi_ptr)
{
  int numnodes = LastGbNode - FirstGbNode + 1;
  double L = Xpos_end - Xpos_start;
  double delta = L / (numnodes - 1.0);
  for (int inode = 0; inode < numnodes; inode++) {
    double xrel = delta * inode;
    double xpos = Xpos_start + xrel;
    int iGbNode = FirstGbNode + inode;
    NodalVars *nv = gi_ptr->NodalVars_GbNode[iGbNode];
    double xfrac = (xpos - DL_ptr_->XLoc_LeftBoundary) / (DL_ptr_->XLoc_RightBoundary - DL_ptr_->XLoc_LeftBoundary);
    nv->setupInitialNodePosition(xpos, xfrac);
    (*gi_ptr->XNodePos_GbNode)[iGbNode] = xpos;
    if (inode == (numnodes - 1.0)) {
      if (fabs(xpos - Xpos_end) > 1.0E-7) {
	throw m1d_Error("BulkDomainDescription::InitializeXposNodes()",
			"Inconsistency in mesh positions");
      }
    }
  }
}
//=====================================================================================================================
// Set the equation description
/*
 *  This routine is responsible for setting the variables:
 *    - NumEquationsPerNode
 *    - VariableNameList
 *    - EquationNameList
 *    - EquationIndexStart_EqName
 */
void
BulkDomainDescription::SetEquationDescription()
{
  /*
   * We will assume that EquationNameList is fully populated.
   * Then populate everything else using that List here 
   */
  if (EquationNameList.size() == 0) {
    throw m1d_Error("BulkDomainDescription::SetEquationDescription()", "No Equations are defined");
  }
  NumEquationsPerNode = EquationNameList.size();
  // Assertion to make to make it easier
  AssertTrace((int) VariableNameList.size() == NumEquationsPerNode);
  int nvar = 0;
  int iVarTPos;
  int iEqnTPos;
  int iSpeciesConservationTPos = (int) Species_Conservation;
  IsAlgebraic_NE.resize(NumEquationsPerNode, 0);
  IsArithmeticScaled_NE.resize(NumEquationsPerNode, 0);

  for (int iEqn = 0; iEqn < NumEquationsPerNode; iEqn++) {
    const EqnType &eqnT = EquationNameList[iEqn];
    const VarType &varT = VariableNameList[iEqn];
    iEqnTPos = (int) eqnT.EquationType;
    iVarTPos = (int) varT.VariableType;
    if (eqnT.EquationType == Momentum_Axial) {
      VarType varN(Velocity_Axial);
      AssertTrace(varT == varN);
      if (EquationIndexStart_EqName[iEqnTPos] == -1) {
        EquationIndexStart_EqName[iEqnTPos] = iEqn;
      }
      if (VariableIndexStart_VarName[iVarTPos] == -1) {
        VariableIndexStart_VarName[iVarTPos] = iEqn;
      }
    } else if (eqnT.EquationType == Momentum_Radial) {
      VarType varN(Velocity_Axial);
      AssertTrace(varT == varN);
      EquationIndexStart_EqName[iEqnTPos] = iEqn;
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Momentum_Swirl) {
      VarType varN(Velocity_Swirl);
      AssertTrace(varT == varN);
      EquationIndexStart_EqName[iEqnTPos] = iEqn;
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == MeshDisplacement_Axial) {
      VarType varN(Displacement_Axial);
      AssertTrace(varT == varN);
      EquationIndexStart_EqName[iEqnTPos] = iEqn;
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Enthalpy_Conservation) {
      VarType varN(Temperature);
      AssertTrace(varT == varN);
      EquationIndexStart_EqName[iEqnTPos] = iEqn;
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Continuity) {
      VarType varN(Velocity_Axial);
      AssertTrace(varT == varN);
      EquationIndexStart_EqName[iEqnTPos] = iEqn;
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Continuity_Global) {
      VarType varN(PressureGradient_Radial);
      AssertTrace(varT == varN);
      EquationIndexStart_EqName[iEqnTPos] = iEqn;
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Species_Conservation) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Concentration_Species, subT, eqnT.EquationSubTypeName);
      if ((varT.VariableType != Concentration_Species) && (varT.VariableType != MoleFraction_Species)) {
        throw m1d_Error("error", "error");
      }
      EquationIndexStart_EqName[iEqnTPos] = MIN(iEqn, EquationIndexStart_EqName[iEqnTPos]);
    } else if (eqnT.EquationType == Current_Conservation) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Voltage, subT);
      AssertTrace(varT == varN);
      EquationIndexStart_EqName[iEqnTPos] = MIN(iEqn, EquationIndexStart_EqName[iEqnTPos]);
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == MoleFraction_Summation) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Concentration_Species, subT, eqnT.EquationSubTypeName);
      if ((varN.VariableType != varT.VariableType) && (MoleFraction_Species != varT.VariableType)) {
        throw m1d_Error("error", "error");
      }
      EquationIndexStart_EqName[iEqnTPos] = MIN(iEqn, EquationIndexStart_EqName[iEqnTPos]);
      /*
       *  Special line -> want the species to point to to MoleFraction_Summation equation if it is the first species!
       */
      EquationIndexStart_EqName[iSpeciesConservationTPos] = MIN(iEqn, EquationIndexStart_EqName[iSpeciesConservationTPos]);
      if (EquationIndexStart_EqName[iSpeciesConservationTPos] == -1) {
	EquationIndexStart_EqName[iSpeciesConservationTPos] = iEqn;
      }

    } else if (eqnT.EquationType == ChargeNeutrality_Summation) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Concentration_Species, subT, eqnT.EquationSubTypeName);
      if ((varN.VariableType != varT.VariableType) && (MoleFraction_Species != varT.VariableType)) {
        throw m1d_Error("error", "error");
      }
      EquationIndexStart_EqName[iEqnTPos] = MIN(iEqn, EquationIndexStart_EqName[iEqnTPos]);
      /*
       *  Special line -> want the species to point to to ChargeNeutrality_Summation equation if it is the first species!
       */
      EquationIndexStart_EqName[iSpeciesConservationTPos] = MIN(iEqn, EquationIndexStart_EqName[iSpeciesConservationTPos]);
      if (EquationIndexStart_EqName[iSpeciesConservationTPos] == -1) {
	EquationIndexStart_EqName[iSpeciesConservationTPos] = iEqn;
      }
    } else if (eqnT.EquationType == Current_Specification) {
	throw m1d_Error("BulkDomainDescription::SetEquationDescription()", "Current Specification for bulk domain");
    } else if (eqnT.EquationType == Voltage_Specification) {
	throw m1d_Error("BulkDomainDescription::SetEquationDescription()", "Voltage Specification for bulk domain");
    } else if (eqnT.EquationType == Dirichlet_Specification) {
	//   We can have different variables with this specification. The first one is temperature
	// VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
	EquationIndexStart_EqName[iEqnTPos] = MIN(iEqn, EquationIndexStart_EqName[iEqnTPos]);
    } else {
      throw m1d_Error("BulkDomainDescription::SetEquationDescription()", "Unknown Conversion");
    }
    nvar++;

    if (EquationIndexStart_EqName[iEqnTPos] == -1) {
      EquationIndexStart_EqName[iEqnTPos] = iEqn;
    }
    if (VariableIndexStart_VarName[iVarTPos] == -1) {
      VariableIndexStart_VarName[iVarTPos] = iEqn;
    }
    /*
     *  PROBABLY need some logic here to figure out if species related offsets are valid
     */
  }
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D *
BulkDomainDescription::mallocDomain1D()
{
  throw m1d_Error("BulkDomainDescription::mallocDomain1D()", "base class called");
  BulkDomainPtr_ = 0;
  return BulkDomainPtr_;
}

//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
