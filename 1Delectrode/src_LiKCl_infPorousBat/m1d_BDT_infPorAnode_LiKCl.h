/**
 * @file m1d_BDT_infPorAnode_LiKCl.h
 */


#ifndef M1D_BDT_INFPORANODE_LIKCL_H_
#define M1D_BDT_INFPORANODE_LIKCL_H_

#include "m1d_BulkDomainTypes.h"
#include "m1d_BDD_porousElectrode.h"

#include <zuzax/transport.h>      // transport properties
#include <zuzax/thermo.h>      // transport properties
#include <zuzax/thermo/IonsFromNeutralVPSSTP.h>  // ion properties

#include "Electrode.h"

namespace m1d
{

//=====================================================================================================================

//! This class consists of multiple species diffusing in a time
//! dependent manner.  There is a net flow and a net electric current.
/*!
 * 
 */
class BDT_infPorAnode_LiKCl : public BDD_porousElectrode
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
  BDT_infPorAnode_LiKCl(DomainLayout *dl_ptr);

  //! Destructor
  virtual
  ~BDT_infPorAnode_LiKCl(); 

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  BDT_infPorAnode_LiKCl (const BDT_infPorAnode_LiKCl &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  BDT_infPorAnode_LiKCl &
  operator=(const BDT_infPorAnode_LiKCl  &r);

  //! Read in the possible models for each domain
  /*!
   *  This procedure is done before the Equations anv variable list are set up.
   *  Needed information about what is possible is input here.
   *  We read the Zuzax ThermoPhase and transport object into DomainDescriptions here.
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
  virtual BulkDomain1D *
  mallocDomain1D();

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
  Zuzax::IonsFromNeutralVPSSTP *ionicLiquidIFN_;

  //! top or bottom of the domain
  /*!
   *   0 - top, right
   *   1 - bottom, left
   */
  int m_position;
};
//=====================================================================================================================
}
//=====================================================================================================================
#endif /*   */
