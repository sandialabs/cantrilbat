/**
 * @file m1d_BDT_infPorCathode_LiKCl.h
 */


#ifndef M1D_BDT_INFPORCATHODE_LIKCL_H_
#define M1D_BDT_INFPORCATHODE_LIKCL_H_

#include "m1d_BulkDomainTypes.h"
#include "m1d_BDD_porousElectrode.h"

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
class BDT_infPorCathode_LiKCl : public BDD_porousElectrode
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
  BDT_infPorCathode_LiKCl(DomainLayout *dl_ptr);

  //! Destructor
  virtual
  ~BDT_infPorCathode_LiKCl();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  BDT_infPorCathode_LiKCl(const BDT_infPorCathode_LiKCl &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  BDT_infPorCathode_LiKCl &
  operator=(const BDT_infPorCathode_LiKCl &r);

  //! Read in the possible models for each domain
  /*!
   *  This procedure is done before the Equations anv variable list are set up.
   *  Needed information about what is possible is input here.
   *  We read the Cantera ThermoPhase and transport object into DomainDescriptions here.
   *
   *   We loop over volume and then surface domains.
   */
  virtual void
  ReadModelDescriptions();

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
  virtual BulkDomain1D *mallocDomain1D();

  //! This is done after the equations are set up
  /*!
   *  We loop over volume and then surface domains here.
   */
  virtual void
  DetermineConstitutiveModels();

  // --------------------------------------------------------------------------------------------

  //! Pointer to the thermo object for the molten salt
  /*!
   *   We own this object
   */
  Cantera::IonsFromNeutralVPSSTP *ionicLiquidIFN_;

  //! Pointer to the transport object for the molten salt
  /*!
   * We own this object
   */
  //Cantera::Transport* trans_;

  //! top or bottom of the domain
  /*!
   *   0 - top, left
   *   1 - bottom, right
   */
  int m_position;

  //! Pointer to the electrode object
  /*!
   * We own the electrode object.
   */
  //Cantera::Electrode *Electrode_;

};

}

#endif /*   */