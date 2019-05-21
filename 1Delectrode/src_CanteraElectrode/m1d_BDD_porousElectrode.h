/**
 * @file m1d_BDD_porousElectrode.h
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_BDD_POROUSELECTRODE_H_
#define M1D_BDD_POROUSELECTRODE_H_

#include "m1d_BulkDomainDescription.h"
#include "m1d_porousElectrode_dom1D.h"
#include "m1d_CanteraElectrodeGlobals.h"
#include "m1d_BDD_porousFlow.h"
#include "Electrode.h"
#include "Electrode_input.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace Zuzax
{
  class Electrode;
  class Transport;
  class ThermoPhase;
}
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! This class consists of multiple species diffusing in a time
//! dependent manner.  There is a net flow and a net electric current.
/*!
 * 
 */
class BDD_porousElectrode : public BDD_porousFlow
{
public:

  //! Constructor
  /*!
   *  This constructor constructs the bulk domain from a MultiPhase object.
   *
   *  In the constructor, we have typically been laying out what the unknowns are
   *  and what the equations are, that are solved within the domain.
   *
   *  @param[in]             dl_ptr              Pointer to the domain layout object
   *  @param[in]             electrodeType       Type of electrode
   *  @param[in]             domainFunctionName  Name of the domain function
   *                                               Default = ""
   *  @param[in]             domainName          Name of the domain 
   *                                               Default = ""
   */
  BDD_porousElectrode(DomainLayout *dl_ptr, int electrodeType, std::string domainFunctionName = "", std::string domainName = "");

  //! Destructor
  virtual ~BDD_porousElectrode(); 

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
  BDD_porousElectrode& operator=(const BDD_porousElectrode &r);

  //! Read in the possible models for each domain
  /*!
   *  This procedure is done before the Equations anv variable list are set up.
   *  Needed information about what is possible is input here.
   *  We read the Zuzax ThermoPhase and transport object into DomainDescriptions here.
   *
   *   We loop over volume and then surface domains.
   */
  virtual void ReadModelDescriptions();

  //! Determine the list of Equations and Variables
  /*!
   *  This routine is responsible for setting the variables:
   *    - VariableNameList
   *    - EquationNameList
   */
  virtual void SetEquationsVariablesList();

  //! This is done after the equations are set up
  /*!
   *  We loop over volume and then surface domains here.
   */
  virtual void DetermineConstitutiveModels();


  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   * @return  Returns a pointer to the object that will calculate the residual efficiently
   */
  virtual BulkDomain1D *
  mallocDomain1D();

  // --------------------------------------------------------------------------------------------
  //            DATA
  // --------------------------------------------------------------------------------------------

  //! Pointer to the electrode object
  /*!
   *  We own the electrode object.
   */
  Zuzax::Electrode* Electrode_;

  //! Pointer to the metal phase that does electrical conduction within the solid
  /*!
   *  We own this object
   */
  Zuzax::ThermoPhase* metalPhase_;

  //! Type of the electrode
  /*!
   *     0 anode
   *     1 cathode
   *     2 reference electrode 
   */
  int electrodeType_;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif 
