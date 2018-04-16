/**
 * @file m1d_SDD_ElectrodeSepInterface.h
 * Definition of an interface between anode-separator or cathode-separator
 * (see class \link m1d::SDD_ElectrodeSepInterface SDD_ElectrodeSepInterface\endlink).
 */

/*
 * Copywrite 2014 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_SDD_ELECTRODESEPINTERFACE_H_
#define M1D_SDD_ELECTRODESEPINTERFACE_H_

#include "m1d_SDD_Mixed.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! This class specifies that all equations are handled
//! by a simple Dirichlet condition or simple flux conditions of the third type
/*!
 *
 */
class SDD_ElectrodeSepInterface : public SDD_Mixed
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
   *  @param[in]             dl_ptr              Domain Layout object that owns this description.
   *  @param[in]             domainName          String containing the name of the domain 
   */
  SDD_ElectrodeSepInterface(DomainLayout *dl_ptr, const std::string& domainName = "");

  //! Destructor
  virtual
  ~SDD_ElectrodeSepInterface();

  //! Copy Constructor
  /*!
   *  @param[in]             r                   Object to be copied
   */
  SDD_ElectrodeSepInterface(const SDD_ElectrodeSepInterface &r);

  //! Assignment operator
  /*!
   *  @param[in]             r                   Object to be copied
   *
   *  @return                                    Returns a changeable reference to the current object
   */
  SDD_ElectrodeSepInterface& operator=(const SDD_ElectrodeSepInterface &r);

  // --------------------------------------------------------------------------------------------------------------
  //                                   DATA
  // --------------------------------------------------------------------------------------------------------------

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
