/**
 * @file m1d_SDD_Mixed.cpp
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#include "m1d_SDD_Mixed.h"
#include "m1d_SurDomain1D.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
SDD_Mixed::SDD_Mixed(DomainLayout* dl_ptr, std::string domainFunctionName, std::string domainName) :
    SurfDomainDescription(dl_ptr, domainFunctionName, domainName),
    NumConditions(0),
    SBC_Type_(0)
{
}
//==================================================================================================================================
SDD_Mixed::SDD_Mixed(DomainLayout* dl_ptr, double value, std::string domainFunctionName, std::string domainName) :
    SurfDomainDescription(dl_ptr, domainFunctionName, domainName),
    NumConditions(1),
    SBC_Type_(0)
{
    EquationID.resize(1);
    VariableID.resize(1);
    Value.resize(1);
    EquationID[0].setID(Equation_Type_Any, -1, "");
    VariableID[0].setID(Variable_Type_Any, -1, "");
    Value[0] = value;
    TimeDep.push_back(nullptr);
    BC_TimeDep_.push_back(nullptr);
    BC_Type_.push_back(0);
}
//==================================================================================================================================
SDD_Mixed::SDD_Mixed(const SDD_Mixed& r) :
    SurfDomainDescription(r.DL_ptr_),
    NumConditions(0), 
    SBC_Type_(0)
{
    *this = r;
}
//==================================================================================================================================
SDD_Mixed::~SDD_Mixed()
{
}
//==================================================================================================================================
SDD_Mixed& SDD_Mixed::operator=(const SDD_Mixed& r)
{
    if (this == &r) {
        return *this;
    }
    SurfDomainDescription::operator=(r);
    NumConditions = r.NumConditions;
    EquationID = r.EquationID;
    Value = r.Value;
    TimeDep = r.TimeDep;  // Caution -> shallow pointer copy
    BC_TimeDep_ = r.BC_TimeDep_;   // Caution -> shallow pointer copy
    BC_Type_ = r.BC_Type_;
    SBC_Type_ = r.SBC_Type_;

    return *this;
}
//==================================================================================================================================
void
SDD_Mixed::addDirichletCondition(EqnType equationID, VarType variableID, double value)
{
    NumConditions++;
    EquationID.push_back(equationID);
    VariableID.push_back(variableID);
    Value.push_back(value);
    TimeDep.push_back(nullptr);
    BC_TimeDep_.push_back(nullptr);
    BC_Type_.push_back(0);
}
//==================================================================================================================================
void
SDD_Mixed::addDirichletCondition(VarType variableID, double value)
{
    NumConditions++;
    VariableID.push_back(variableID);
    int eqInt = int(variableID.VariableType);
    EQ_TYPE eq = (EQ_TYPE) eqInt;
    EqnType equationID(eq, variableID.VariableType);
    EquationID.push_back(equationID);
    Value.push_back(value);
    TimeDep.push_back(nullptr);
    BC_TimeDep_.push_back(nullptr);
    BC_Type_.push_back(0);
}
//==================================================================================================================================
void
SDD_Mixed::addDirichletCondition(EqnType equationID, VarType variableID, double value, double (*timeDep)(double))
{
    NumConditions++;
    EquationID.push_back(equationID);
    VariableID.push_back(variableID);
    Value.push_back(value);
    TimeDep.push_back(timeDep);
    BC_TimeDep_.push_back(nullptr);
    BC_Type_.push_back(2);
}
//==================================================================================================================================
void
SDD_Mixed::addDirichletCondition(EqnType equationID, VarType variableID, int BC_Type, BoundaryCondition* BC_timeDep)
{
    NumConditions++;
    EquationID.push_back(equationID);
    VariableID.push_back(variableID);
    Value.push_back(0);
    TimeDep.push_back(nullptr);
    BC_TimeDep_.push_back(BC_timeDep);
    BC_Type_.push_back(BC_Type);
}
//==================================================================================================================================
void
SDD_Mixed::addFluxCondition(const EqnType& equationID, const VarType& variableID, double value)
{
    NumConditions++;
    EquationID.push_back(equationID);
    VariableID.push_back(variableID);
    Value.push_back(value);
    TimeDep.push_back(nullptr);
    BC_TimeDep_.push_back(nullptr);
    BC_Type_.push_back(1);
}
//==================================================================================================================================
void
SDD_Mixed::addFluxCondition(const EqnType& equationID, const VarType& variableID, double value, double (*timeDep)(double))
{
    NumConditions++;
    EquationID.push_back(equationID);
    VariableID.push_back(variableID);
    Value.push_back(value);
    TimeDep.push_back(timeDep);
    BC_TimeDep_.push_back(nullptr);
    BC_Type_.push_back(3);
}
//==================================================================================================================================
void
SDD_Mixed::addFluxCondition(const EqnType& equationID, const VarType& variableID, int BC_Type, BoundaryCondition* BC_timeDep)
{
    NumConditions++;
    EquationID.push_back(equationID);
    VariableID.push_back(variableID);
    Value.push_back(0.0);
    TimeDep.push_back(nullptr);
    BC_TimeDep_.push_back(BC_timeDep);
    BC_Type_.push_back(BC_Type);
}
//==================================================================================================================================
void
SDD_Mixed::addRobinCondition(EqnType equationID, VarType variableID, BoundaryCondition* BC_timeDep, int bc_type)
{
    NumConditions++;
    EquationID.push_back(equationID);
    VariableID.push_back(variableID);
    Value.push_back(0.0);
    TimeDep.push_back(nullptr);
    BC_TimeDep_.push_back(BC_timeDep);
    BC_Type_.push_back(bc_type);
}
//==================================================================================================================================
void
SDD_Mixed::SetEquationDescription()
{
    /*
     * Dirichlet equations don't have any extra equations
     * that are solved on surface domains.
     */
    EquationNameList.clear();

    /*
     * Set the policy for connecting bulk domains
     * This really isn't set yet.
     */
    setRLMapping(0);
    /*
     * Fill in the rest of the information
     */
    SurfDomainDescription::SetEquationDescription();
}
//==================================================================================================================================
SurDomain1D*
SDD_Mixed::mallocDomain1D()
{
    SurDomain1DPtr_ = new SurBC_Dirichlet(*this);
    return SurDomain1DPtr_;
}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------

