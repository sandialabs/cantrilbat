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
  DomainDescription::DomainDescription(DomainLayout *dl_ptr, std::string domainFunctionName, std::string domainName) :
  NumEquationsPerNode(0), 
  IsArithmeticScaled_NE(0), 
  DomainName(domainName),
  DomainFunctionName_(domainFunctionName),
  DL_ptr_(dl_ptr),
  SolutionBehavior_printLvl_(0),
  Residual_printLvl_(0),
  porosityEquationProbType_( Porosity_EqnType_Status::None)
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
  NumEquationsPerNode(0), 
  DomainName(""), 
  DomainFunctionName_(""), 
  DL_ptr_(0)
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
  DomainFunctionName_ = r.DomainFunctionName_;
  DL_ptr_ = r.DL_ptr_;
  SolutionBehavior_printLvl_ = r.SolutionBehavior_printLvl_;
  Residual_printLvl_ = r.Residual_printLvl_;
  porosityEquationProbType_ = r.porosityEquationProbType_;

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
//! Reorder the variables and equations on this domain
/*!
 *  this is a necessary step
 */
void DomainDescription::ReorderVariablesEquations()
{
    VarType varJ, varBest;
    size_t nvars = VariableNameList.size();
    if (nvars > 1) {
	std::vector<size_t> order(nvars, npos);
	std::vector<int> taken(nvars, 0);
	for (size_t i = 0; i < nvars; ++i) {
	    varBest = VarType(Max_Var_Name);
	    VarType varI = VarType(Max_Var_Name);
	    size_t iBest = npos;
	    for (size_t j = 0; j < nvars; ++j) {
		if (! taken[j]) {
		    varJ = VariableNameList[j];
		    if (varJ.VariableType < varBest.VariableType) {
			iBest = j;
			varBest = varJ;
		    }
		}
	    }
	    if (iBest != npos) {
		taken[iBest] = 1;
		order[i] = iBest;
	    } else {
		throw m1d_Error("DomainDescription::ReorderVariablesEquations()", "unknown index problem");
	    }
	}
	std::vector<int> IsArithmeticScaled_NE_copy = IsArithmeticScaled_NE;
	std::vector<int> IsAlgebraic_NE_copy = IsAlgebraic_NE;
	std::vector<VarType>  VariableNameList_copy = VariableNameList;
	std::vector<EqnType> EquationNameList_copy = EquationNameList;
	
	for (size_t i = 0; i < nvars; ++i) {
	    size_t iC = order[i];
	    VariableNameList[i] = VariableNameList_copy[iC];
	    EquationNameList[i] = EquationNameList_copy[iC];
	    IsAlgebraic_NE[i] =  IsAlgebraic_NE_copy[iC];
	    IsArithmeticScaled_NE[i] = IsArithmeticScaled_NE_copy[iC];
	}
    }
}
//===================================================================================================================================
void
DomainDescription::DetermineConstitutiveModels()
{
}
//===========================================================================
} /* End of Namespace */
//===========================================================================
