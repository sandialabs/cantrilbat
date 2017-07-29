/**
 * @file m1d_SurfDomainTypes.cpp
 */


#include "m1d_SDD_Mixed_Cu2S.h"
#include "m1d_SurDomain_Cu2S.h"

//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
SDD_Mixed_Cu2S::SDD_Mixed_Cu2S(DomainLayout *dl_ptr, int pos) :
    SDD_Mixed(dl_ptr),
    m_position(pos)
{
}
//=====================================================================================================================
SDD_Mixed_Cu2S::SDD_Mixed_Cu2S(const SDD_Mixed_Cu2S &r) :
  SDD_Mixed(r.DL_ptr_),
   m_position(0)
{
  *this = r;
}
//=====================================================================================================================
SDD_Mixed_Cu2S::~SDD_Mixed_Cu2S()
{
}
//=====================================================================================================================
SDD_Mixed_Cu2S &
SDD_Mixed_Cu2S::operator=(const SDD_Mixed_Cu2S &r)
{
  if (this == &r) {
    return *this;
  }

  SDD_Mixed::operator=(r);
  m_position = r.m_position;

  return *this;
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
SDD_Mixed_Cu2S::SetEquationDescription()
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
  setRLMapping(1);
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
SDD_Mixed_Cu2S::mallocDomain1D()
{
  if (m_position == 0) {
    return new Cu2S_TopSurface(*this, 1);
  }
  return new Cu2S_BotSurface(*this, 1);
}

//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================

