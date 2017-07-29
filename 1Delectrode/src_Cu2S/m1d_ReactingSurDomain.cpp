/*
 * @file m1d_SurReactionRate.h
 * Virtual base class 
 */

#include "m1d_ReactingSurDomain.h"
//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
ReactingSurDomain::ReactingSurDomain(SurfDomainDescription &sdd) :
  SurDomain1D(sdd), m_NumSpecies(0), m_Temperature(300.)
{
}
//=====================================================================================================================
ReactingSurDomain::ReactingSurDomain(const ReactingSurDomain &r) :
  SurDomain1D(r.SDD_)
{
  ReactingSurDomain::operator=(r);
}
//=====================================================================================================================
ReactingSurDomain::~ReactingSurDomain()
{
}
//=====================================================================================================================
ReactingSurDomain &
ReactingSurDomain::operator=(const ReactingSurDomain &r)
{
  if (this == &r) {
    return *this;
  }
  SurDomain1D::operator=(r);

  m_NumSpecies = r.m_NumSpecies;
  m_Temperature = r.m_Temperature;

  return *this;
}
//=====================================================================================================================
double
ReactingSurDomain::calculateRate(const double * const y)
{
  throw m1d_Error("SurReactionRate::calculateRate", "base class called -> error");
}
//=====================================================================================================================
}
//=====================================================================================================================
