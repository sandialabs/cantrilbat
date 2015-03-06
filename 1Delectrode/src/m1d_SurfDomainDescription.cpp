/**
 * @file m1d_SurfDomainDescription.cpp
 *
 */

/*
 *  $Id: m1d_SurfDomainDescription.cpp 567 2013-03-21 23:03:11Z hkmoffa $
 */
#include "m1d_SurfDomainDescription.h"

#include "m1d_NodalVars.h"
#include "m1d_GlobalIndices.h"
#include "m1d_exception.h"
#include "m1d_DomainLayout.h"
#ifndef MIN
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif

namespace m1d
{
//=====================================================================================================================
  SurfDomainDescription::SurfDomainDescription(DomainLayout *dl_ptr, std::string domainName) :
    DomainDescription(dl_ptr,domainName), IDSurfDomain(-1), LocGbNode(-1), LeftDomain(0), LeftBulk(0), RightDomain(0), RightBulk(0)
{


}
//=====================================================================================================================
SurfDomainDescription::~SurfDomainDescription()
{
}
//=====================================================================================================================
SurfDomainDescription::SurfDomainDescription(const SurfDomainDescription &r) :
  DomainDescription(r.DL_ptr_), IDSurfDomain(-1), LocGbNode(-1), LeftDomain(0), LeftBulk(0), RightDomain(0),
      RightBulk(0)
{
  *this = r;
}
//=====================================================================================================================
SurfDomainDescription &
SurfDomainDescription::operator=(const SurfDomainDescription &r)
{
  if (this == &r) {
    return *this;
  }

  DomainDescription::operator=(r);

  IDSurfDomain = r.IDSurfDomain;
  LocGbNode = r.LocGbNode;
  LeftDomain = r.LeftDomain;
  LeftBulk = r.LeftBulk;
  RightDomain = r.RightDomain;
  RightBulk = r.RightBulk;

  return *this;
}
//=====================================================================================================================
// Sets the id of the surface domain
void
SurfDomainDescription::setID(int id)
{
  IDSurfDomain = id;
}
//=====================================================================================================================
// Reports the id of the surface domain
int
SurfDomainDescription::ID() const
{
  return IDSurfDomain;
}
//=====================================================================================================================
// set the adjacent domains
void
SurfDomainDescription::setAdjBulkDomains(BulkDomainDescription *leftBulk, BulkDomainDescription *rightBulk)
{
  LeftDomain = (DomainDescription *) leftBulk;
  LeftBulk = leftBulk;
  RightDomain = (DomainDescription *) rightBulk;
  RightBulk = rightBulk;

}
//=====================================================================================================================
void
SurfDomainDescription::setGbNode(const int locGbNode)
{
  LocGbNode = locGbNode;
}
//====================================================================================================================
// Determine the list of Equations and Variables
/*
 *  This routine is responsible for setting the variables:
 *    - VariableNameList
 *    - EquationNameList
 *  We do nothing here, so that the base class is overriden, but nothing is done
 */
void
SurfDomainDescription::SetEquationsVariablesList()
{
}
//=====================================================================================================================
// Set the equation description
/*
 *  This routine is responsible for setting the variables:
 *    - NumEquationsPerNode
 *    - EquationIndexStart_EqName 
 *    - VariableIndexStart_VarName
 *
 *    For surfaces, there very well may be no extra equations assigned on
 *    the surface domain.
 */
void
SurfDomainDescription::SetEquationDescription()
{
  /*
   * We will assume that EquationNameList is fully populated.
   * Then populate everything else from there.
   *
   *  It is not an error for a surface domain to have no extra
   *  equations assigned to it.
   */
  if (EquationNameList.size() == 0) {
    return;
  }
  NumEquationsPerNode = EquationNameList.size();
  // Assertion to make to make it easier
  AssertTrace((int) VariableNameList.size() == NumEquationsPerNode);
  int nvar = 0;
  int iVarTPos;
  int iEqnTPos;
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
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Momentum_Radial) {
      VarType varN(Velocity_Axial);
      AssertTrace(varT == varN);
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Momentum_Swirl) {
      VarType varN(Velocity_Swirl);
      AssertTrace(varT == varN);
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == MeshDisplacement_Axial) {
      VarType varN(Displacement_Axial);
      AssertTrace(varT == varN);
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Enthalpy_Conservation) {
      VarType varN(Temperature);
      AssertTrace(varT == varN);
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Continuity) {
      VarType varN(Velocity_Axial);
      AssertTrace(varT == varN);
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Continuity_Global) {
      VarType varN(PressureGradient_Radial);
      AssertTrace(varT == varN);
      IsArithmeticScaled_NE[iEqn] = 1;
    } else if (eqnT.EquationType == Species_Conservation) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Concentration_Species, subT, eqnT.EquationSubTypeName);
      if ((varT.VariableType != Concentration_Species) && (varT.VariableType != MoleFraction_Species)) {
        throw m1d_Error("error", "error");
      }
    } else if (eqnT.EquationType == Current_Conservation) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Voltage, subT);
      AssertTrace(varT == varN);

    } else if (eqnT.EquationType == MoleFraction_Summation) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Concentration_Species, subT, eqnT.EquationSubTypeName);
      if ((varN.VariableType != varT.VariableType) && (MoleFraction_Species != varT.VariableType)) {
        throw m1d_Error("error", "error");
      }

    } else if (eqnT.EquationType == ChargeNeutrality_Summation) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Concentration_Species, subT, eqnT.EquationSubTypeName);
      if ((varN.VariableType != varT.VariableType) && (MoleFraction_Species != varT.VariableType)) {
        throw m1d_Error("error", "error");
      }

    } else if (eqnT.EquationType == Voltage_Specification) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Voltage, subT);
      AssertTrace(varT == varN);
      IsArithmeticScaled_NE[iEqn] = 1;

    } else if (eqnT.EquationType == Current_Specification) {
      VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
      VarType varN(Voltage, subT);
      AssertTrace(varT == varN);
      IsArithmeticScaled_NE[iEqn] = 1;

    } else {
      throw m1d_Error("SurfDomainDescription::SetEquationDescription()", "Unknown Conversion");
    }
    nvar++;
    if (EquationIndexStart_EqName[iEqnTPos] == -1) {
      EquationIndexStart_EqName[iEqnTPos] = iEqn;
    }
    if (VariableIndexStart_VarName[iVarTPos] == -1) {
      VariableIndexStart_VarName[iVarTPos] = iEqn;
    }
  }

}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
SurDomain1D *
SurfDomainDescription::mallocDomain1D()
{
  throw m1d_Error("SurfDomainDescription::mallocDomain1D()", "base class called");
  return (SurDomain1D *) 0;
}
//=====================================================================================================================
// Sets the mapping between right and left domains to either
// of two extremes
/*
 * This sets two common mapping extremes. A different function
 * sets the inbetween cases.
 *
 * @param mapType     0 If the map type is 0, this means that
 *                      the right equations are mapped into the
 *                      left equations to the full extent possible
 *                    1 If the map type is 1, his means that
 *                      the right equations are never mapped into the
 *                      left equations at all. They are separate equations
 */
void
SurfDomainDescription::setRLMapping(int mapType)
{
  int i;
  if (RightBulk == 0) {
    return;
  }
  /*
   * get the number of equations in the right domain.
   */
  int neRight = RightBulk->NumEquationsPerNode;
  RightToLeftBulkEqnMapping.resize(neRight);
  for (i = 0; i < neRight; i++) {
    RightToLeftBulkEqnMapping[i] = -1;
  }

  if (LeftDomain == 0) {
    return;
  }
  int neLeft = LeftDomain->NumEquationsPerNode;
  /*
   * We are done if the mapType is one.
   */
  if (mapType == 1) {
    return;
  }
  if (mapType != 0) {
    throw m1d_Error("SurfDomainDescription::setRLMapping", "unknown map type");
  }
  /*
   * Go find the mapping. Mappings for the right which are not available
   * remain set to -1
   */
  for (i = 0; i < neRight; i++) {
    /*
     * get the equation type and subtype of the right side.
     */
    EQ_Name_Enum etr = (RightBulk->EquationNameList[i]).EquationType;
    int subtyper = (RightBulk->EquationNameList[i]).EquationSubType;
    /*
     * Go find it on the left side
     */
    int ifound = 0;
    for (int j = 0; j < neLeft; j++) {
      EQ_Name_Enum etl = (LeftBulk->EquationNameList[j]).EquationType;
      int subtypel = (LeftBulk->EquationNameList[j]).EquationSubType;
      if (etr == etl && subtyper == subtypel) {
	if (ifound == 1) {
	  m1d_Error("SurfDomainDescription::setRLMapping", "confused");
	}
	RightToLeftBulkEqnMapping[i] = j;
        ifound = 1;
      }
    }

  }
}
//=====================================================================================================================
// Sets the mapping between right and left domains to the
// arbitrary case
/*
 *   @param rightConnectivity map of the right into the left
 *  Length is equal to the number of equations in the right domain
 */
void
SurfDomainDescription::setRLMapping(const int * const rightConnectivity)
{
  int ne = RightBulk->NumEquationsPerNode;
  int neLeft = LeftBulk->NumEquationsPerNode;
  RightToLeftBulkEqnMapping.resize(ne);
  for (int i = 0; i < ne; i++) {
    int j = rightConnectivity[i];
    if (j < -1 || j >= neLeft) {
      throw m1d_Error("SurfDomainDescription::setRLMapping", "unknown map type");
    }
    RightToLeftBulkEqnMapping[i] = rightConnectivity[i];
  }
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
