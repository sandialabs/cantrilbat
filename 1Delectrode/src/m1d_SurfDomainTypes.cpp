/**
 * @file m1d_SurfDomainTypes.cpp
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#include "m1d_SurfDomainTypes.h"
#include "m1d_SurDomain1D.h"

namespace m1d
{

//=====================================================================================================================
  SDT_Dirichlet::SDT_Dirichlet(DomainLayout *dl_ptr, std::string domainFunctionName, std::string domainName) :
    SurfDomainDescription(dl_ptr, domainFunctionName, domainName), 
    NumConditions(0), 
    EquationID(0), 
    Value(0), 
    TimeDep(0), 
    BC_TimeDep_(0), 
    BC_Type_(0)
{
}
//=====================================================================================================================
SDT_Dirichlet::SDT_Dirichlet(DomainLayout *dl_ptr, double value, std::string domainFunctionName, std::string domainName) :
    SurfDomainDescription(dl_ptr, domainFunctionName, domainName), 
    NumConditions(1), 
    TimeDep(0), 
    BC_TimeDep_(0), 
    BC_Type_(0)
{
  EquationID.resize(1);
  VariableID.resize(1);
  Value.resize(1);
  EquationID[0].setID(Equation_Type_Any, -1, "");
  VariableID[0].setID(Variable_Type_Any, -1, "");
  Value[0] = value;
  TimeDep.push_back(0);
  BC_TimeDep_.push_back(0);
  BC_Type_.push_back(0);
}
//=====================================================================================================================
SDT_Dirichlet::SDT_Dirichlet(const SDT_Dirichlet &r) :
  SurfDomainDescription(r.DL_ptr_), 
  NumConditions(0), EquationID(0), Value(0), TimeDep(0), BC_TimeDep_(0), BC_Type_(0)
{
  *this = r;
}
//=====================================================================================================================
SDT_Dirichlet::~SDT_Dirichlet()
{
}
//=====================================================================================================================
SDT_Dirichlet &
SDT_Dirichlet::operator=(const SDT_Dirichlet &r)
{
  if (this == &r) {
    return *this;
  }
  SurfDomainDescription::operator=(r);
  NumConditions = r.NumConditions;
  EquationID = r.EquationID;
  Value = r.Value;
  TimeDep = r.TimeDep;
  BC_TimeDep_ = r.BC_TimeDep_;
  BC_Type_ = r.BC_Type_;

  return *this;
}
//=====================================================================================================================
// Add a Dirichlet Condition
/*
 *
 * @param  equationID  Equation ID to apply the Dirichlet condition to
 * @param  variableID  VariableID to apply the Dirichlet condition to
 * @param  value  Value to apply
 */
void
SDT_Dirichlet::addDirichletCondition(EqnType equationID, VarType variableID, double value) {
  NumConditions++;
  EquationID.push_back(equationID);
  VariableID.push_back(variableID);
  Value.push_back(value);
  TimeDep.push_back(0);
  BC_TimeDep_.push_back(0);
  BC_Type_.push_back(0);
}
//=====================================================================================================================
// Add a Dirichlet Condition assuming the default mapping between variable and equation ID
/*
 * @param  variableID  VariableID to apply the Dirichlet condition to
 * @param  value  Value to apply
 */
void
SDT_Dirichlet::addDirichletCondition(VarType variableID, double value)
{
  NumConditions++;
  VariableID.push_back(variableID);
  Value.push_back(value);
  int eqInt = int(variableID.VariableType);
  EQ_TYPE eq = (EQ_TYPE) eqInt;
  EqnType equationID(eq, variableID.VariableType);
  EquationID.push_back(equationID);
  TimeDep.push_back(0);
  BC_TimeDep_.push_back(0);
  BC_Type_.push_back(0);
}

//=====================================================================================================================
// Add a Dirichlet Condition with time dependence
/*
 *
 * @param  equationID  Equation ID to apply the Dirichlet condition to
 * @param  variableID  VariableID to apply the Dirichlet condition to
 * @param  value  Value to apply
 * @param  timeDep  Function pointer for time dependence: normalized by value
 */
void
SDT_Dirichlet::addDirichletCondition(EqnType equationID, VarType variableID, double value, double (*timeDep)(double))
{
  NumConditions++;
  EquationID.push_back(equationID);
  VariableID.push_back(variableID);
  Value.push_back(value);
  TimeDep.push_back(timeDep);
  BC_Type_.push_back(2);
}
//=====================================================================================================================
// Add a Dirichlet Condition with time dependence of BoundaryCondition class
/*
 *
 * @param  equationID  Equation ID to apply the Dirichlet condition to
 * @param  variableID  VariableID to apply the Dirichlet condition to
 * @param  value  Value to apply
 * @param  timeDep  Function pointer for time dependence: normalized by value
 */
void
SDT_Dirichlet::addDirichletCondition(EqnType equationID, VarType variableID, int BC_Type, BoundaryCondition *BC_timeDep)
{
  NumConditions++;
  EquationID.push_back(equationID);
  VariableID.push_back(variableID);
  Value.push_back(0);
  TimeDep.push_back(0);
  BC_TimeDep_.push_back(BC_timeDep);
  BC_Type_.push_back(BC_Type);
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
SDT_Dirichlet::SetEquationDescription()
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
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
SurDomain1D *
SDT_Dirichlet::mallocDomain1D()
{
  SurDomain1DPtr_ = new SurBC_Dirichlet(*this);
  return SurDomain1DPtr_;
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
SDT_Mixed::SDT_Mixed(DomainLayout *dl_ptr, std::string domainFunctionName, std::string domainName) :
    SDT_Dirichlet(dl_ptr, domainFunctionName, domainName), 
    SBC_Type_(0)
{
}
//=====================================================================================================================
SDT_Mixed::SDT_Mixed(const SDT_Mixed &r) :
  SDT_Dirichlet(r.DL_ptr_), SBC_Type_(0)
{
  *this = r;
}
//=====================================================================================================================
SDT_Mixed::~SDT_Mixed()
{
}
//=====================================================================================================================
SDT_Mixed &
SDT_Mixed::operator=(const SDT_Mixed &r)
{
  if (this == &r) {
    return *this;
  }

  SDT_Dirichlet::operator=(r);
  SBC_Type_ = r.SBC_Type_;
  return *this;
}
//=====================================================================================================================
// Add a flux Condition
/*
 *
 * @param  equationID  Equation ID to apply the flux condition to
 * @param  variableID  VariableID to apply the flux condition to
 * @param  value  Value to apply
 */
void
SDT_Mixed::addFluxCondition(const EqnType& equationID, const VarType& variableID, double value)
{
  NumConditions++;
  EquationID.push_back(equationID);
  VariableID.push_back(variableID);
  Value.push_back(value);
  TimeDep.push_back(0);
  BC_TimeDep_.push_back(0);
  BC_Type_.push_back(1);
}
//=====================================================================================================================
// Add a flux Condition using time dependent continuous function
/*
 *
 * @param  equationID  Equation ID to apply the flux condition to
 * @param  variableID  VariableID to apply the flux condition to
 * @param  value  Value to apply
 */
void
SDT_Mixed::addFluxCondition(const EqnType& equationID, const VarType& variableID, double value, double (*timeDep)(double))
{
  NumConditions++;
  EquationID.push_back(equationID);
  VariableID.push_back(variableID);
  Value.push_back(value);
  TimeDep.push_back(timeDep);
  BC_TimeDep_.push_back(0);
  BC_Type_.push_back(3);
}
//=====================================================================================================================
// Add a flux Condition using Boundary Condition class
/*
 *
 * @param  equationID  Equation ID to apply the flux condition to
 * @param  variableID  VariableID to apply the flux condition to
 * @param  value  Value to apply
 */
void
SDT_Mixed::addFluxCondition(const EqnType & equationID, const VarType& variableID, int BC_Type, BoundaryCondition *BC_timeDep)
{
  NumConditions++;
  EquationID.push_back(equationID);
  VariableID.push_back(variableID);
  Value.push_back(0.0);
  TimeDep.push_back(0);
  BC_TimeDep_.push_back(BC_timeDep);
  BC_Type_.push_back(BC_Type);
}
//=====================================================================================================================
// Add a Robin Mixed boundary Condition
/*
 *
 * @param  equationID  Equation ID to apply the flux condition to
 * @param  variableID  VariableID to apply the flux condition to
 * @param  value  Value to apply
 */
void
SDT_Mixed::addRobinCondition(EqnType equationID, VarType variableID, BoundaryCondition *BC_timeDep, int bc_type)
{
  NumConditions++;
  EquationID.push_back(equationID);
  VariableID.push_back(variableID);
  Value.push_back(0.0);
  TimeDep.push_back(0);
  BC_TimeDep_.push_back(BC_timeDep);
  BC_Type_.push_back(bc_type);
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
SurDomain1D *
SDT_Mixed::mallocDomain1D()
{
  SurDomain1DPtr_ = new SurBC_Dirichlet(*this);
  return SurDomain1DPtr_;
}
//=====================================================================================================================

} /* End of Namespace */
//=====================================================================================================================

