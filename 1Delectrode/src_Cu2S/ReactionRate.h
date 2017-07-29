/*
 * @file ReactionRate.h
 * Virtual base class 
 */
/*
 *$Id: ReactionRate.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef REACTIONRATE_H
#define REACTIONRATE_H

#include "m1d_defs.h"

//=====================================================================================================================
namespace CanteraLite
{
//=====================================================================================================================
//=====================================================================================================================
//! Small utility class to return the rate of progress of a reaction given
//! a vector of concentrations and the temperature
class ReactionRate
{
public:
  ReactionRate(double temperature);
  virtual
  ~ReactionRate();

  //! Return the rate of progress of a single reaction
  /*!
   *
   * @param y     Vector of concentrations
   * @return      returns the rate of progress (for surface reactions the units are kmol m-2 s-1
   */
  virtual double
  calculateRate(const double * const y);

  //! Get the stoichiometric coefficient for the particular species
  virtual double
  getStoichSpec(int ispec);

  //! set the stoichiometric coefficient for a species
  virtual void
  setStoichSpec(int nspec, const double * const stoichCoeff);

protected:
  int m_NumSpecies;
public:
  double m_Temperature;
protected:
  std::vector<double> m_stoichCoeff;
};
//=====================================================================================================================
}
//=====================================================================================================================
#endif
//=====================================================================================================================
