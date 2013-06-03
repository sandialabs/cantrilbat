/*
 * @file ReactionRate.h
 * Virtual base class 
 */
/*
 *$Id: m1d_ReactingSurDomain.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef M1D_REACTINGSURDOMAIN_H
#define M1D_REACTINGSURDOMAIN_H

#include "m1d_SurDomain1D.h"

#include <vector>

namespace m1d
{
//=====================================================================================================================
class ReactingSurDomain : public SurDomain1D
{
public:
  //! Constructor for the surface reaction base class
  /*!
   *
   * @param sdd
   * @param temperature
   */
  ReactingSurDomain(SurfDomainDescription &sdd);

  //! Copy Constructor
  /*!
   *
   * @param sdd
   * @param temperature
   */
  ReactingSurDomain(const ReactingSurDomain &r);

  //! Destructor
  virtual
  ~ReactingSurDomain();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  ReactingSurDomain &
  operator=(const ReactingSurDomain &r);

  virtual double
  calculateRate(const double * const y);

protected:
  //! Number of species
  int m_NumSpecies;

  //! Temperature
  double m_Temperature;

};
//=====================================================================================================================
}

#endif
