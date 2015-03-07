/**
 * @file m1d_BDT_porAnode_LiKCl.h
 */
/*
 *   $Id: m1d_BDT_porAnode_LiKCl.h 504 2013-01-07 22:32:48Z hkmoffa $
 */

#ifndef M1D_BDT_PORANODE_LIKCL_H_
#define M1D_BDT_PORANODE_LIKCL_H_

#include "m1d_BulkDomainTypes.h"

#include <cantera/transport.h>      // transport properties
#include <cantera/thermo.h>      // transport properties
#include <cantera/thermo/IonsFromNeutralVPSSTP.h>  // ion properties

#include "Electrode.h"

namespace m1d
{

//=====================================================================================================================

//! This class consists of multiple species diffusing in a time
//! dependent manner.  There is a net flow and a net electric current.
/*!
 * 
 */
class BDT_porAnode_LiKCl : public BulkDomainDescription
{
public:

  //! Constructor
  /*!
   * This constructor constructs the bulk domain from a MultiPhase object.
   *
   * In the constructor, we have typically been laying out what the unknowns are
   * and what the equations are, that are solved within the domain.
   *
   */
  BDT_porAnode_LiKCl (DomainLayout *dl_ptr);

  //! Destructor
  virtual
  ~BDT_porAnode_LiKCl(); 

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  BDT_porAnode_LiKCl (const BDT_porAnode_LiKCl &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  BDT_porAnode_LiKCl &
  operator=(const BDT_porAnode_LiKCl  &r);

  //! Determine the list of Equations and Variables
  /*!
   *  This routine is responsible for setting the variables:
   *    - VariableNameList
   *    - EquationNameList
   */
  virtual void
  SetEquationsVariablesList();

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual BulkDomain1D *
  mallocDomain1D();

  // --------------------------------------------------------------------------------------------

  //! Pointer to the thermo object for the molten salt
  /*!
   *   We own this object
   */
  Cantera::IonsFromNeutralVPSSTP *ionicLiquid_;

  //! Pointer to the transport object for the molten salt
  /*!
   * We own this object
   */
  Cantera::Transport* trans_;



  //! top or bottom of the domain
  /*!
   *   0 - top, right
   *   1 - bottom, left
   */
  int m_position;


  //! Pointer to the electrode object
  /*!
   * We own the electrode object.
   */
  Cantera::Electrode *Electrode_;


};
//=====================================================================================================================
}
//=====================================================================================================================
#endif /*   */
