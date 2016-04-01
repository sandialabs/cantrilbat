/**
 * @file m1d_DomainDescription.cpp
 */

/*
 *  $Id: m1d_BulkDomainTypes.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_BulkDomainTypes.h"

#include "m1d_SimpleDiff_dom1D.h"
#include "m1d_SimpleTDDiff_dom1D.h"

namespace m1d
{
//===========================================================================
BDT_SimpleDiff::BDT_SimpleDiff(DomainLayout *dl_ptr) :
  BulkDomainDescription(dl_ptr)
{
}
//===========================================================================
BDT_SimpleDiff::BDT_SimpleDiff(DomainLayout *dl_ptr, int id) :
  BulkDomainDescription(dl_ptr)
{
  const EqnType eqnT(Species_Conservation, 0, "Species0");
  EquationNameList.push_back(eqnT);
  VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
  VariableNameList.push_back(VarType(Concentration_Species, subT, eqnT.EquationSubTypeName));
}
//===========================================================================
BDT_SimpleDiff::BDT_SimpleDiff(const BDT_SimpleDiff &r) :
  BulkDomainDescription(r.DL_ptr_)
{
  *this = r;
}
//===========================================================================
BDT_SimpleDiff::~BDT_SimpleDiff()
{
}
//===========================================================================
BDT_SimpleDiff &
BDT_SimpleDiff::operator=(const BDT_SimpleDiff &r)
{
  if (this == &r) {
    return *this;
  }

  BulkDomainDescription::operator=(r);

  EquationID = r.EquationID;

  return *this;
}
//=============================================================================================================
void BDT_SimpleDiff::SetEquationsVariablesList()
{
  EquationNameList.clear();
  VariableNameList.clear();
  const EqnType eqnT(Species_Conservation, 0, "Species0");
  EquationNameList.push_back(eqnT);
  VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
  VariableNameList.push_back(VarType(Concentration_Species, subT, eqnT.EquationSubTypeName));
}
//===========================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D *
BDT_SimpleDiff::mallocDomain1D()
{
   BulkDomainPtr_ = new SimpleDiff_dom1D(this);
   return BulkDomainPtr_;
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
BDT_SimpleTDDiff::BDT_SimpleTDDiff(DomainLayout *dl_ptr) :
  BulkDomainDescription(dl_ptr)
{
}
//===========================================================================
BDT_SimpleTDDiff::BDT_SimpleTDDiff(DomainLayout *dl_ptr, int id) :
  BulkDomainDescription(dl_ptr)
{
  const EqnType eqnT(Species_Conservation, 0, "Species0");
  EquationNameList.push_back(eqnT);
  VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
  VariableNameList.push_back(VarType(Concentration_Species, subT, eqnT.EquationSubTypeName));
}
//===========================================================================
BDT_SimpleTDDiff::BDT_SimpleTDDiff(const BDT_SimpleTDDiff &r) :
  BulkDomainDescription(r.DL_ptr_)
{
  *this = r;
}
//===========================================================================
BDT_SimpleTDDiff::~BDT_SimpleTDDiff()
{
}
//===========================================================================
BDT_SimpleTDDiff &
BDT_SimpleTDDiff::operator=(const BDT_SimpleTDDiff &r)
{
  if (this == &r) {
    return *this;
  }
  BulkDomainDescription::operator=(r);
  EquationID = r.EquationID;
  return *this;
}
//==============================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D *
BDT_SimpleTDDiff::mallocDomain1D()
{
   BulkDomainPtr_ = new SimpleTDDiff_dom1D(this);
   return BulkDomainPtr_;
}
//==============================================================================================================================
void BDT_SimpleTDDiff::SetEquationsVariablesList()
{
  EquationNameList.clear();
  VariableNameList.clear();
  const EqnType eqnT(Species_Conservation, 0, "Species0");
  EquationNameList.push_back(eqnT);
  VAR_TYPE_SUBNUM subT = eqnT.EquationSubType;
  VariableNameList.push_back(VarType(Concentration_Species, subT, eqnT.EquationSubTypeName));
}
//=============================================================================================================================
} /* End of Namespace */
//===========================================================================

