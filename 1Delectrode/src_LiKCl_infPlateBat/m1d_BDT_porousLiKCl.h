/**
 * @file m1d_BDT_porousLiKCl.h
 */
/*
 *   $Id: m1d_BDT_porousLiKCl.h 506 2013-01-07 22:43:59Z hkmoffa $
 */

#ifndef M1D_BDT_POROUSLIKCL_H_
#define M1D_BDT_POROUSLIKCL_H_

#include "m1d_BulkDomainTypes.h"

#include <cantera/transport.h>      // transport properties
#include <cantera/thermo.h>      // transport properties
#include <cantera/thermo/IonsFromNeutralVPSSTP.h>  // ion properties
#include <cantera/thermo/StoichSubstance.h>  // separator
namespace m1d
{

//==================================================================
//==================================================================
//==================================================================

//! This class consists of multiple species diffusing in a time
//! dependent manner.  There is a net flow and a net electric current.
/*!
 *  This class is used to test the implementation
 */
class BDT_porousLiKCl : public BulkDomainDescription
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

  BDT_porousLiKCl(DomainLayout *dl_ptr);

  //! Destructor
  virtual
  ~BDT_porousLiKCl();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  BDT_porousLiKCl(const BDT_porousLiKCl &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  BDT_porousLiKCl &
  operator=(const BDT_porousLiKCl &r);

  virtual void SetEquationsVariablesList();

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual BulkDomain1D *
  mallocDomain1D();

  // --------------------------------------------------------------------------------------------

  //! Equation type and var type to apply them
  std::vector<VarType> EquationID;

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

};

}

#endif /* M1D_BULKDOMAINTYPES_H_ */
