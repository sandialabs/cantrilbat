/**
 * @file SDT_ElectrodeSepInterface.h
 * Definition of an interface between anode-separator or cathode-separator
 * (see class \link m1d::SDT_ElectrodeSepInterface SDT_ElectrodeSepInterface\endlink).
 */

/*
 * Copywrite 2014 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_SDT_ELECTRODESEPINTERFACE_H_
#define M1D_SDT_ELECTRODESEPINTERFACE_H_

#include "m1d_SurfDomainTypes.h"

namespace m1d
{

//! This class specifies that all equations are handled
//! by a simple Dirichlet condition or simple flux conditions of the third type
/*!
 *
 */
class SDT_ElectrodeSepInterface : public SDT_Mixed
{
public:

  //! Constructor
  /*!
   *   We construct the object but don't actually specify any Dirichlet conditions.
   *   Later we can add dirichlet conditions into the object.
   *
   * In the constructor, we have typically been laying out what the unknowns are
   * and what the equations are, that are solved within the domain.
   *
   *
   * @param dl_ptr  Domain Layout object that owns this description.
   */
  SDT_ElectrodeSepInterface(DomainLayout *dl_ptr, const char *domainName = "");

  //! Destructor
  virtual
  ~SDT_ElectrodeSepInterface();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  SDT_ElectrodeSepInterface(const SDT_ElectrodeSepInterface &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  SDT_ElectrodeSepInterface &
  operator=(const SDT_ElectrodeSepInterface &r);

  // --------------------------------------------------------------------------------------------------------------
  //                                   DATA
  // --------------------------------------------------------------------------------------------------------------

};

}
#endif
