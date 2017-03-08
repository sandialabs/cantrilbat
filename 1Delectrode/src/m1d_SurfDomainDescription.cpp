/**
 * @file m1d_SurfDomainDescription.cpp
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_SurfDomainDescription.h"

#include "m1d_NodalVars.h"
#include "m1d_GlobalIndices.h"
#include "m1d_exception.h"
#include "m1d_DomainLayout.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
SurfDomainDescription::SurfDomainDescription(DomainLayout* dl_ptr, std::string domainFunctionName,
                                             std::string domainName) :
    DomainDescription(dl_ptr,domainName),
    IDSurfDomain(-1), LocGbNode(-1),
    LeftDomain(0),
    LeftBulk(0),
    RightDomain(nullptr),
    RightBulk(nullptr),
    SurDomain1DPtr_(nullptr)
{
}
//==================================================================================================================================
SurfDomainDescription::~SurfDomainDescription()
{
}
//==================================================================================================================================
SurfDomainDescription::SurfDomainDescription(const SurfDomainDescription& r) :
    DomainDescription(r.DL_ptr_), IDSurfDomain(-1), LocGbNode(-1), LeftDomain(0), LeftBulk(0), RightDomain(0),
    RightBulk(0),
    SurDomain1DPtr_(0)
{
    *this = r;
}
//==================================================================================================================================
SurfDomainDescription&
SurfDomainDescription::operator=(const SurfDomainDescription& r)
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
    SurDomain1DPtr_ = r.SurDomain1DPtr_;

    return *this;
}
//==================================================================================================================================
void
SurfDomainDescription::setID(int id)
{
    IDSurfDomain = id;
}
//==================================================================================================================================
int
SurfDomainDescription::ID() const
{
    return IDSurfDomain;
}
//==================================================================================================================================
void
SurfDomainDescription::setAdjBulkDomains(BulkDomainDescription* leftBulk, BulkDomainDescription* rightBulk)
{
    LeftDomain = (DomainDescription*) leftBulk;
    LeftBulk = leftBulk;
    RightDomain = (DomainDescription*) rightBulk;
    RightBulk = rightBulk;
}
//==================================================================================================================================
void
SurfDomainDescription::setGbNode(const int locGbNode)
{
    LocGbNode = locGbNode;
}
//==================================================================================================================================
void
SurfDomainDescription::SetEquationsVariablesList()
{
}
//==================================================================================================================================
/*
 *  This routine is responsible for setting the variables:
 *    - NumEquationsPerNode
 *    - EquationIndexStart_EqName
 *    - VariableIndexStart_VarName
 *
 *    For surfaces, there very well may be no extra equations assigned on the surface domain.
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
        const EqnType& eqnT = EquationNameList[iEqn];
        const VarType& varT = VariableNameList[iEqn];
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
#ifdef MECH_MODEL
        } else if (eqnT.EquationType ==  Mechanical_Model_Axial) {
            VarType varN(Displacement_Axial);
            AssertTrace(varT == varN);
            IsArithmeticScaled_NE[iEqn] = 1;
        } else if (eqnT.EquationType == Mechanical_Stress_Axial) {
            VarType varN(Solid_Stress_Axial);
            AssertTrace(varT == varN);
            IsArithmeticScaled_NE[iEqn] = 1;
#endif
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
//==================================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
SurDomain1D*
SurfDomainDescription::mallocDomain1D()
{
    throw m1d_Error("SurfDomainDescription::mallocDomain1D()", "base class called");
    SurDomain1DPtr_ = 0;
    return SurDomain1DPtr_;
}
//==================================================================================================================================
void SurfDomainDescription::setRLMapping(int mapType)
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
//==================================================================================================================================
void
SurfDomainDescription::setRLMappingVector(const int* const rightConnectivity)
{
    int ne = RightBulk->NumEquationsPerNode;
    int neLeft = LeftBulk->NumEquationsPerNode;
    RightToLeftBulkEqnMapping.resize(ne);
    for (int i = 0; i < ne; i++) {
        int j = rightConnectivity[i];
        if (j < -1 || j >= neLeft) {
            throw m1d_Error("SurfDomainDescription::setRLMappingVector", "unknown map type");
        }
        RightToLeftBulkEqnMapping[i] = rightConnectivity[i];
    }
}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
