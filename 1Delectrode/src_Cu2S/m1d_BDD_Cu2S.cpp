/**
 * @file m1d_BDD_Cu2S.cpp
 */

/*
 *  $Id: m1d_BDD_Cu2S.cpp 361 2012-08-21 00:39:02Z hkmoffa $
 */

#include "m1d_BDD_Cu2S.h"


#include "m1d_TDGrowingFilm_dom1D.h"

namespace m1d
{

//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
/*
 * In the constructor, we have typically been laying out what the unknowns are
 * and what the equations are, that are solved within the domain.
 *
 */
BDD_Cu2S::BDD_Cu2S(DomainLayout *dl_ptr, int id) :
  BulkDomainDescription(dl_ptr)
{
  IsAlgebraic_NE.resize(7,0);
  IsArithmeticScaled_NE.resize(7,0);
  int eqnIndex = 0;
  EquationNameList.push_back(EqnType(MeshDisplacement_Axial, 0, "MeshDisplacement"));
  IsArithmeticScaled_NE[eqnIndex] = 1;

  EquationNameList.push_back(EqnType(Species_Conservation, 0, "Species0"));

  VariableNameList.push_back(VarType(Displacement_Axial));
  VariableNameList.push_back(VarType(Concentration_Species, 0, "Species0"));
}
//===========================================================================
BDD_Cu2S::BDD_Cu2S(const BDD_Cu2S &r) :
  BulkDomainDescription(r.DL_ptr_)
{
  *this = r;
}
//===========================================================================
BDD_Cu2S::~BDD_Cu2S()
{
}

//===========================================================================
BDD_Cu2S &
BDD_Cu2S::operator=(const BDD_Cu2S &r)
{
  if (this == &r) {
    return *this;
  }

  BulkDomainDescription::operator=(r);

  EquationID = r.EquationID;

  return *this;
}
//=====================================================================================================================
void
BDD_Cu2S::SetEquationsVariablesList()
{
    int eqnIndex = 0;
    EquationNameList.clear();
    VariableNameList.clear();
  
    IsAlgebraic_NE.resize(7,0);
    IsArithmeticScaled_NE.resize(7,0);
  
    EquationNameList.push_back(EqnType(MeshDisplacement_Axial, 0, "MeshDisplacement"));
    VariableNameList.push_back(VarType(Displacement_Axial));
    IsArithmeticScaled_NE[eqnIndex] = 1;
    
    EquationNameList.push_back(EqnType(Species_Conservation, 0, "Species0")); 
    VariableNameList.push_back(VarType(Concentration_Species, 0, "Species0"));
}
//===========================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D *
BDD_Cu2S::mallocDomain1D()
{
   BulkDomainPtr_ = new TDGrowingFilm_dom1D(this);
   return  BulkDomainPtr_;
}
//===========================================================================
} /* End of Namespace */
//===========================================================================

