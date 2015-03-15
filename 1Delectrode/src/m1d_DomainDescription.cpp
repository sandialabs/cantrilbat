/**
 * @file m1d_DomainDescription.cpp
 *
 */

/*
 *  $Id: m1d_DomainDescription.cpp 567 2013-03-21 23:03:11Z hkmoffa $
 */
#include "m1d_DomainDescription.h"
#include "m1d_NodalVars.h"
#include "m1d_GlobalIndices.h"
#include "m1d_exception.h"
#include "m1d_DomainLayout.h"
#include "m1d_ProblemStatement.h"
#ifndef MIN
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif

namespace m1d
{

//===========================================================================
  DomainDescription::DomainDescription(DomainLayout *dl_ptr, std::string domainName) :
  NumEquationsPerNode(0), IsArithmeticScaled_NE(0), DomainName(domainName),
  DL_ptr_(dl_ptr),
  SolutionBehavior_printLvl_(0),
  Residual_printLvl_(0)
{
    EquationIndexStart_EqName.resize((int) Max_Eqn_Name, -1);
    VariableIndexStart_VarName.resize((int) Max_Var_Name, -1);

    /*
     * Get problem
     */
    // pointer to the Problem residual -> not assigned yet
    // ProblemResidEval *problemResid = DL_ptr_->problemResid_;

    // Pointer to the input file
    ProblemStatement *psInput_ptr = DL_ptr_->psInput_ptr_;

    SolutionBehavior_printLvl_ = psInput_ptr->SolutionBehavior_printLvl_;
    Residual_printLvl_ = psInput_ptr->Residual_printLvl_;
}
//===========================================================================
DomainDescription::~DomainDescription()
{
}
//===========================================================================
DomainDescription::DomainDescription(const DomainDescription &r) :
  NumEquationsPerNode(0), DomainName(""), DL_ptr_(0)
{
  *this = r;
}
//===========================================================================
DomainDescription &
DomainDescription::operator=(const DomainDescription &r)
{
  if (this == &r) {
    return *this;
  }

  NumEquationsPerNode = r.NumEquationsPerNode;
  VariableNameList = r.VariableNameList;
  EquationNameList = r.EquationNameList;
  EquationIndexStart_EqName = r.EquationIndexStart_EqName;
  VariableIndexStart_VarName = r.VariableIndexStart_VarName;
  IsAlgebraic_NE = r.IsAlgebraic_NE;
  IsArithmeticScaled_NE = r.IsArithmeticScaled_NE;
  DomainName = r.DomainName;
  DL_ptr_ = r.DL_ptr_;
  SolutionBehavior_printLvl_ = r.SolutionBehavior_printLvl_;
  Residual_printLvl_ = r.Residual_printLvl_;

  return *this;
}
//===================================================================================================================================
void
DomainDescription::ReadModelDescriptions()
{
}
//===================================================================================================================================
// Determine the list of Equations and Variables
/*
 *  This routine is responsible for setting the variables:
 *    - VariableNameList
 *    - EquationNameList
 */
void
DomainDescription::SetEquationsVariablesList() {
    throw m1d_Error("DomainDescription::SetEquationVariablesList()", "Base class implementation called");
}
//===================================================================================================================================
// Set the equation description
/*
 *  This routine is responsible for setting the variables:
 *    - NumEquationsPerNode
 *    - EquationIndexStart_EqName
 */
void
DomainDescription::SetEquationDescription()
{
  throw m1d_Error("DomainDescription::SetEquationDescription()", "Base class implementation called");
}
//===================================================================================================================================
void
DomainDescription::DetermineConstitutiveModels()
{
}
//===========================================================================
} /* End of Namespace */
//===========================================================================
