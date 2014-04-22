/**
 * @file m1d_BDD_porousElectrode.h
 */

/*
 *   $Id: m1d_BDD_porousElectrode.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_BDD_POROUSELECTRODE_H_
#define M1D_BDD_POROUSELECTRODE_H_

#include "m1d_BulkDomainDescription.h"

namespace Cantera {
  class Electrode;
  class Transport;
  class ThermoPhase;
}

namespace m1d
{
//! This class consists of multiple species diffusing in a time
//! dependent manner.  There is a net flow and a net electric current.
/*!
 * 
 */
class BDD_porousElectrode : public BulkDomainDescription
{
public:

  //! Constructor
  /*!
   * This constructor constructs the bulk domain from a MultiPhase object.
   *
   * In the constructor, we have typically been laying out what the unknowns are
   * and what the equations are, that are solved within the domain.
   *
   * @param dl_ptr   Pointer to the domain layout object
   */
  BDD_porousElectrode(DomainLayout *dl_ptr, std::string domainName = "", ELECTRODE_KEY_INPUT *input_ptr );

  //! Destructor
  virtual
  ~BDD_porousElectrode(); 

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  BDD_porousElectrode(const BDD_porousElectrode &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  BDD_porousElectrode &
  operator=(const BDD_porousElectrode &r);

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   * @return  Returns a pointer to the object that will calculate the residual efficiently
   */
  virtual BulkDomain1D *
  mallocDomain1D();

  // --------------------------------------------------------------------------------------------
  //            DATA
  // --------------------------------------------------------------------------------------------

  //! Pointer to the thermo object for the electrolyte
  /*!
   *   We own this object
   */
  Cantera::ThermoPhase *ionicLiquid_;

  //! Pointer to the transport object for the electrolyte
  /*!
   * We own this object
   */
  Cantera::Transport* trans_;

  //! Pointer to the electrode object
  /*!
   * We own the electrode object.
   */
  Cantera::Electrode *Electrode_;
};
//=====================================================================================================================
}
//=====================================================================================================================
#endif 
