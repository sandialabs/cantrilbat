/*
 * @file ReactionRate.h
 * Virtual base class 
 */
/*
 *$Id: ReactionRate.cpp 506 2013-01-07 22:43:59Z hkmoffa $
 */

#include "ReactionRate.h"
#include "cantera/base/ctexceptions.h"

namespace CanteraLite
{
/**************************************************************************
 *
 *
 */
ReactionRate::ReactionRate(double temperature) :
  m_NumSpecies(0), m_Temperature(temperature)
{
}

/*************************************************************************
 *
 *
 */
ReactionRate::~ReactionRate()
{
}

/*************************************************************************
 *
 *
 */
void
ReactionRate::setStoichSpec(int nspec, const double * const stoichCoeff)
{
  m_NumSpecies = nspec;
  m_stoichCoeff.resize(nspec, 0.0);
  for (int i = 0; i < nspec; i++) {
    m_stoichCoeff[i] = stoichCoeff[i];
  }
}

/*************************************************************************
 *
 *
 */
double
ReactionRate::getStoichSpec(int ispec)
{
  return (m_stoichCoeff[ispec]);
}

/**************************************************************************
 *
 *
 */
double
ReactionRate::calculateRate(const double * const y)
{
  throw ZZCantera::CanteraError("ReactionRate::calculateRate",
                              "base class called -> error");
}

}
